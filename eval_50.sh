#!/usr/bin/env bash
# eval_50.sh -- THEOBROMA broad-spectrum stress-test / plausibility / boundary suite.
# 50 probes across: dedup/InChIKey identity, classification coherence, ADMET domain,
# similarity metrics, taxonomy reconciliation, licensing, stereochemistry, cross-refs,
# integrity invariants, and API robustness. Read-only except where noted (all SELECT).
# Run:  bash eval_50.sh 2>&1 | tee eval_50_out.txt
# Each block prints a labeled result; compare against the "expect" comment.
set +e
PSQL="psql -h localhost -U theobroma -d theobroma -tAc"
BASE="http://localhost:5000"
hr(){ echo "----- $1"; }
q(){ $PSQL "$1"; }        # scalar/table query
code(){ curl -s -o /dev/null -w "%{http_code}" "$1"; }  # http status only

echo "=================== THEOBROMA eval_50 ==================="
echo "date: $(date -u +%FT%TZ)"

# ============ A. CORPUS / COUNT RECONCILIATION INVARIANTS ============
hr "1. total compounds (anchor; expect 1132805)"
q "SELECT count(*) FROM compounds;"

hr "2. kingdom partition sums to total? (expect diff 0; secondary kingdoms may inflate slices, so test PRIMARY only)"
q "SELECT (SELECT count(*) FROM compounds) - (SELECT count(*) FROM compounds WHERE kingdom IN ('plant','animal','fungi','bacteria','unresolved')) AS should_be_0;"

hr "3. any compound with NULL/empty kingdom? (expect 0; every row should have a primary kingdom)"
q "SELECT count(*) FROM compounds WHERE kingdom IS NULL OR kingdom='';"

hr "4. distinct 14-char skeletons vs compounds (skeletons<=compounds; report both)"
q "SELECT count(DISTINCT substring(inchikey,1,14)) AS skeletons, count(*) AS compounds FROM compounds;"

hr "5. multi-member stereo families: every family has >=2 members sharing 14-char skeleton (expect 0 singletons among 'families')"
q "SELECT count(*) FROM (SELECT substring(inchikey,1,14) k, count(*) c FROM compounds GROUP BY 1 HAVING count(*)>1) t WHERE c<2;"

hr "6. license tier partition sums to total (expect diff 0)"
q "SELECT (SELECT count(*) FROM compounds) - (SELECT count(*) FROM compounds WHERE license_tier IN ('CC0','CC BY 4.0','CC BY-NC 4.0','CC BY-NC-SA 4.0','CC BY-NC-ND 4.0','Unspecified')) AS should_be_0;"

# ============ B. INChIKEY / IDENTITY / DEDUP INTEGRITY ============
hr "7. InChIKey format validity: 27 chars, 14-BLOCK-10-1-1 pattern (expect 0 malformed)"
q "SELECT count(*) FROM compounds WHERE inchikey !~ '^[A-Z]{14}-[A-Z]{10}-[A-Z]$';"

hr "8. duplicate full InChIKeys (comp_id should be 1:1 with inchikey; expect 0 dups)"
q "SELECT count(*) FROM (SELECT inchikey FROM compounds GROUP BY inchikey HAVING count(*)>1) t;"

hr "9. non-standard InChIKey second-block flag: keys ending in 'N' vs others (standard=N; count non-N = isotopic/nonstandard)"
q "SELECT right(inchikey,1) AS lastchar, count(*) FROM compounds GROUP BY 1 ORDER BY 2 DESC;"

hr "10. stereo-layer block distribution: how many share skeleton but differ only in 2nd block (true stereoisomers)"
q "SELECT count(*) AS multi_stereo_skeletons FROM (SELECT substring(inchikey,1,14) k FROM compounds GROUP BY 1 HAVING count(DISTINCT substring(inchikey,16,10))>1) t;"

hr "11. RED FLAG: stereo 'families' whose members differ in MOLECULAR FORMULA (skeleton-hash grouping error; expect small/0)"
q "SELECT count(*) FROM (SELECT substring(inchikey,1,14) k FROM compounds WHERE molecular_formula IS NOT NULL GROUP BY 1 HAVING count(DISTINCT molecular_formula)>1) t;"

hr "12. metal-containing entries (InChI disconnects metals -> identity risk). Count formulas with common metals"
q "SELECT count(*) FROM compounds WHERE molecular_formula ~ '(Fe|Cu|Zn|Mg|Mn|Co|Ni|Mo|Pt|Pd|Ru|Na|K|Ca)[0-9A-Z]' ;"

# ============ C. NPCLASSIFIER / CLASSIFICATION COHERENCE ============
hr "13. np_class populated coverage (report count + pct)"
q "SELECT count(*) FILTER (WHERE np_class IS NOT NULL AND np_class<>'') AS have, round(100.0*count(*) FILTER (WHERE np_class IS NOT NULL AND np_class<>'')/count(*),2) AS pct FROM compounds;"

hr "14. multi-label check: np_class values containing the ' | ' separator (NPClassifier is multi-label; expect >0)"
q "SELECT count(*) FROM compounds WHERE np_class LIKE '% | %';"

hr "15. pathway coherence probe (Quassinoids should be Terpenoids/Triterpenoids only; list pathway/superclass combos)"
q "SELECT np_pathway, np_superclass, count(*) FROM compounds WHERE np_class LIKE '%Quassinoids%' GROUP BY 1,2 ORDER BY 3 DESC;"

hr "16. class present but pathway empty (coherence gap; report count)"
q "SELECT count(*) FROM compounds WHERE np_class IS NOT NULL AND np_class<>'' AND (np_pathway IS NULL OR np_pathway='');"

hr "17. inferred_class (XGBoost backfill) coverage and overlap: rows with inferred_class but empty np_class"
q "SELECT count(*) FILTER (WHERE inferred_class IS NOT NULL AND inferred_class<>'') AS inferred_total, count(*) FILTER (WHERE (inferred_class IS NOT NULL AND inferred_class<>'') AND (np_class IS NULL OR np_class='')) AS inferred_fills_gap FROM compounds;"

hr "18. organometallic/inorganic residue in classification: metal-formula rows that still carry np_class (NPClassifier shouldn't classify these well)"
q "SELECT count(*) FROM compounds WHERE molecular_formula ~ '(Fe|Cu|Zn|Pt|Pd|Ru|Co|Ni|Mo)[0-9A-Z]' AND np_class IS NOT NULL AND np_class<>'';"

hr "19. classyfire_superclass coverage (carried-forward; report count + distinct values)"
q "SELECT count(*) FILTER (WHERE classyfire_superclass IS NOT NULL AND classyfire_superclass<>'') AS have, count(DISTINCT classyfire_superclass) AS distinct_vals FROM compounds;"

hr "20. distinct np_pathway values (expect ~7 canonical NPClassifier pathways, plus multi-label combos)"
q "SELECT count(DISTINCT np_pathway) FROM compounds WHERE np_pathway IS NOT NULL AND np_pathway<>'';"

# ============ D. ADMET APPLICABILITY / SANITY ============
hr "21. ADMET coverage: compounds with any persisted ADMET value (report count)"
q "SELECT count(*) FROM admet;"

hr "22. ADMET probability endpoints within [0,1] (DILI as probe; expect min>=0 max<=1)"
q "SELECT min(dili), max(dili) FROM admet WHERE dili IS NOT NULL;" 2>/dev/null || echo "(dili column name differs; adjust)"

hr "23. LogP physical sanity (RDKit MolLogP; flag extreme outliers |logp|>30)"
q "SELECT count(*) FROM compounds WHERE logp IS NOT NULL AND abs(logp)>30;"

hr "24. MW sanity: zero/negative or absurd (>5000 Da) molecular weights"
q "SELECT count(*) FILTER (WHERE mw<=0) AS nonpos, count(*) FILTER (WHERE mw>5000) AS over5000 FROM compounds WHERE mw IS NOT NULL;"

hr "25. TPSA/HBA/HBD non-negative (expect 0 negatives)"
q "SELECT count(*) FROM compounds WHERE (tpsa<0) OR (hba<0) OR (hbd<0);"

hr "26. ADMET applicability caveat probe: large NPs (MW>1000) still carry ADMET (out-of-domain; report count)"
q "SELECT count(*) FROM compounds c JOIN admet a ON a.comp_id=c.comp_id WHERE c.mw>1000;"

# ============ E. SIMILARITY METRICS / BOUNDARIES ============
hr "27. NaFM coverage vs corpus (fp16 index rows; expect ~1128910 = corpus minus OOV atoms)"
wc -l ~/theobroma/data/vectors/nafm_comp_ids.txt 2>/dev/null || echo "(nafm_comp_ids.txt not found)"

hr "28. NaFM OOV gap: corpus minus NaFM-covered (expect ~3895 skipped, out-of-vocab atoms)"
echo "corpus=$(q "SELECT count(*) FROM compounds;")  nafm=$(wc -l < ~/theobroma/data/vectors/nafm_comp_ids.txt 2>/dev/null)"

hr "29. NaFM self-similarity: curcumin top hit is itself at ~1.0 (metric sanity)"
curl -s "$BASE/api/similarity?smiles=curcumin&metric=nafm&top_n=3" | head -c 300; echo

hr "30. Morgan metric returns for a common query (expect 200 + results)"
code "$BASE/api/similarity?q=quercetin&metric=morgan&limit=5"; echo " <- morgan quercetin"

hr "31. MACCS metric returns (expect 200)"
code "$BASE/api/similarity?q=quercetin&metric=maccs&limit=5"; echo " <- maccs quercetin"

hr "32. ChemBERTa metric returns (expect 200)"
code "$BASE/api/similarity?q=quercetin&metric=chemberta&limit=5"; echo " <- chemberta quercetin"

hr "33. metric disagreement probe: do morgan vs nafm top-neighbor differ for curcumin? (heterogeneous metrics SHOULD differ)"
echo "morgan:"; curl -s "$BASE/api/similarity?q=curcumin&metric=morgan&limit=2" | grep -oE 'THEO_[0-9]+' | head -2 | tr '\n' ' '; echo
echo "nafm:  "; curl -s "$BASE/api/similarity?q=curcumin&metric=nafm&top_n=2" | grep -oE 'THEO_[0-9]+' | head -2 | tr '\n' ' '; echo

hr "34. NaFM off-corpus novel SMILES (lookup-only design: should return empty/message, NOT 500)"
code "$BASE/api/similarity?smiles=C1CC1CC1CC1CCCCCCNNNN&metric=nafm&top_n=5"; echo " <- novel SMILES nafm (expect 200, empty)"

# ============ F. TAXONOMY RECONCILIATION / HOMONYMS ============
hr "35. taxonomy resolution coverage: rows with >=1 resolved rank (report pct)"
q "SELECT round(100.0*count(*) FILTER (WHERE (genus IS NOT NULL AND genus<>'') OR (family IS NOT NULL AND family<>'') OR (taxorder IS NOT NULL AND taxorder<>'') OR (taxclass IS NOT NULL AND taxclass<>'') OR (phylum IS NOT NULL AND phylum<>''))/count(*),2) AS pct_resolved FROM resolved_taxonomy;"

hr "36. APG clade values in taxclass (post-APG-fix; expect eudicots/monocots/magnoliids/etc.)"
q "SELECT taxclass, count(*) FROM resolved_taxonomy WHERE taxclass IN ('eudicots','monocots','magnoliids','basal angiosperms','gymnosperms') GROUP BY 1 ORDER BY 2 DESC;"

hr "37. deprecated plant classes should be GONE (expect 0 Magnoliopsida/Liliopsida)"
q "SELECT count(*) FROM resolved_taxonomy WHERE taxclass IN ('Magnoliopsida','Liliopsida','magnoliopsida','liliopsida');"

hr "38. homonym risk: genus names appearing across multiple kingdoms (report top offenders)"
q "SELECT ct.genus, count(DISTINCT c.kingdom) AS kingdoms FROM compound_taxonomy ct JOIN compounds c ON c.comp_id=ct.comp_id WHERE ct.genus IS NOT NULL AND ct.genus<>'' GROUP BY 1 HAVING count(DISTINCT c.kingdom)>1 ORDER BY 2 DESC LIMIT 10;"

hr "39. order-unfilled-for-plants signature (NCBI gap): plant rows with class but no order (report count)"
q "SELECT count(*) FROM resolved_taxonomy rt JOIN compounds c ON c.comp_id=rt.comp_id WHERE c.kingdom='plant' AND (rt.taxclass IS NOT NULL AND rt.taxclass<>'') AND (rt.taxorder IS NULL OR rt.taxorder='');"

hr "40. taxclass casing consistency (pinopsida vs Pinopsida etc.; list any mixed-case survivors)"
q "SELECT taxclass, count(*) FROM resolved_taxonomy WHERE taxclass ~ '^[A-Z]' AND taxclass NOT IN ('') GROUP BY 1 ORDER BY 2 DESC LIMIT 10;"

# ============ G. LICENSING LOGIC / TWO-LEVEL RULE ============
hr "41. attestation-vs-ref consistency (post-reconcile; expect 0 disagreements)"
q "SELECT count(*) FROM per_source_license_attestation a JOIN source_license_ref r ON LOWER(r.src)=LOWER(a.source) WHERE a.license_tier<>r.license_tier;"

hr "42. no CC0 compound without a CC0 source attestation (structure-as-fact invariant; expect 0 violations)"
q "SELECT count(*) FROM compounds c WHERE c.license_tier='CC0' AND NOT EXISTS (SELECT 1 FROM per_source_license_attestation a JOIN source_license_ref r ON LOWER(r.src)=LOWER(a.source) WHERE a.comp_id=c.comp_id AND r.license_tier='CC0');"

hr "43. most-permissive-wins invariant: compound tier == min(tier_rank) across its sources (expect 0 violations)"
q "SELECT count(*) FROM (SELECT c.comp_id, c.tier_rank AS ct, min(r.tier_rank) FILTER (WHERE r.tier_rank<=4) AS best FROM compounds c JOIN per_source_license_attestation a ON a.comp_id=c.comp_id JOIN source_license_ref r ON LOWER(r.src)=LOWER(a.source) GROUP BY 1,2) t WHERE ct <> COALESCE(best,5);"

hr "44. Unspecified only if ALL sources unspecified (expect 0 violations)"
q "SELECT count(*) FROM compounds c WHERE c.license_tier='Unspecified' AND EXISTS (SELECT 1 FROM per_source_license_attestation a JOIN source_license_ref r ON LOWER(r.src)=LOWER(a.source) WHERE a.comp_id=c.comp_id AND r.tier_rank<=4);"

hr "45. every compound has >=1 attestation (expect 0 orphans)"
q "SELECT count(*) FROM compounds c WHERE NOT EXISTS (SELECT 1 FROM per_source_license_attestation a WHERE a.comp_id=c.comp_id);"

# ============ H. CROSS-REFS / API ROBUSTNESS / UX ============
hr "46. ChEBI linkage: compounds with chebi_id (expect ~19508) and all valid-looking"
q "SELECT count(*) FILTER (WHERE chebi_id IS NOT NULL AND chebi_id<>'') AS have, count(*) FILTER (WHERE chebi_id IS NOT NULL AND chebi_id !~ '^(CHEBI:)?[0-9]+$') AS malformed FROM compounds;"

hr "47. API robustness: empty query, SQL-ish injection, unicode, overlong (expect no 500s)"
for u in "$BASE/api/search?q=&type=name" "$BASE/api/search?q=%27%3B%20DROP%20TABLE%20compounds%3B--&type=name" "$BASE/api/search?q=%CE%B1-pinene&type=name" "$BASE/api/search?q=$(python3 -c 'print("A"*600)')&type=name"; do
  printf "  %s -> %s\n" "$(echo "$u" | cut -c1-55)" "$(code "$u")"
done

hr "48. filter_options + clade search regression (both must be 200 post apg_clade fix)"
echo "  filter_options -> $(code "$BASE/api/filter_options")"
echo "  clade=eudicots -> $(code "$BASE/search?type=clade&q=eudicots")"

hr "49. license-provenance metadata now most-permissive (expect the corrected rule string)"
curl -s "$BASE/api/compound/THEO_0854403/license-provenance" | grep -oE '"resolution_rule":"[^"]*"'; echo

hr "50. stereoisomer API + substructure + bulk (all 200)"
echo "  stereoisomers -> $(code "$BASE/api/stereoisomers/THEO_0854403")"
echo "  substructure  -> $(code "$BASE/search?type=smiles&smiles_mode=substructure&q=c1ccccc1")"
echo "  bulk(csv)     -> $(code "$BASE/api/bulk?cols=comp_id,name&format=csv")"

echo "=================== eval_50 DONE ==================="
