#!/usr/bin/env bash
# eval_50b.sh -- corrected re-runs of the eval_50 probes that errored or mis-fired,
# now using the real schema (no molecular_formula column; ADMET cols are quoted
# CamelCase; similarity API uses smiles=/top_n=). Also resolves the taxonomy
# coverage definitional question and the substructure 302.
set +e
PSQL="psql -h localhost -U theobroma -d theobroma -tAc"
BASE="http://localhost:5000"
hr(){ echo "----- $1"; }
q(){ $PSQL "$1"; }
code(){ curl -s -o /dev/null -w "%{http_code}" "$1"; }

echo "=================== eval_50b (corrections) ==================="

# --- 11'. stereo-family formula mismatch, via InChI formula layer (2nd '/'-field) ---
# InChI = "InChI=1S/<FORMULA>/c...": extract formula = split_part(inchi,'/',2).
hr "11'. stereo families whose members differ in InChI FORMULA layer (skeleton-hash error; expect small/0)"
q "SELECT count(*) FROM (
     SELECT substring(inchikey,1,14) k
     FROM compounds WHERE inchi LIKE 'InChI=%'
     GROUP BY 1 HAVING count(DISTINCT split_part(inchi,'/',2))>1
   ) t;"

hr "11b. show a few example collision skeletons (if any) with their differing formulas"
q "SELECT substring(inchikey,1,14) AS skel, string_agg(DISTINCT split_part(inchi,'/',2), ' vs ') AS formulas, count(*) AS members
   FROM compounds WHERE inchi LIKE 'InChI=%'
   GROUP BY 1 HAVING count(DISTINCT split_part(inchi,'/',2))>1
   ORDER BY members DESC LIMIT 5;"

# --- 12'. metal-containing entries via InChI formula layer ---
hr "12'. metal-containing entries (InChI formula layer contains a metal symbol)"
q "SELECT count(*) FROM compounds
   WHERE split_part(inchi,'/',2) ~ '(Fe|Cu|Zn|Mg|Mn|Co|Ni|Mo|Pt|Pd|Ru|Ag|Au|Hg|Pb|Cd|Cr|Sn|Sb|As)([0-9]|$|[A-Z])';"

hr "12b. of those metal entries, how many still carry np_class (NPClassifier weak on organometallics)"
q "SELECT count(*) FROM compounds
   WHERE split_part(inchi,'/',2) ~ '(Fe|Cu|Zn|Mg|Mn|Co|Ni|Mo|Pt|Pd|Ru|Ag|Au|Hg)([0-9]|$|[A-Z])'
     AND np_class IS NOT NULL AND np_class<>'';"

# --- 22'. ADMET probability endpoints within [0,1] (DILI, quoted CamelCase) ---
hr "22'. DILI in [0,1] (expect min>=0 max<=1)"
q "SELECT min(\"DILI\"), max(\"DILI\") FROM admet WHERE \"DILI\" IS NOT NULL;"
hr "22b. hERG + AMES + BBB ranges (probability endpoints, expect [0,1])"
q "SELECT 'hERG' m, min(\"hERG\"), max(\"hERG\") FROM admet WHERE \"hERG\" IS NOT NULL
   UNION ALL SELECT 'AMES', min(\"AMES\"), max(\"AMES\") FROM admet WHERE \"AMES\" IS NOT NULL
   UNION ALL SELECT 'BBB', min(\"BBB_Martins\"), max(\"BBB_Martins\") FROM admet WHERE \"BBB_Martins\" IS NOT NULL;"

# --- 30-33'. similarity API with the correct param convention (smiles=, top_n=) ---
hr "30'. Morgan via smiles= + top_n= (expect 200)"
code "$BASE/api/similarity?smiles=quercetin&metric=morgan&top_n=5"; echo " <- morgan"
hr "31'. MACCS (expect 200)"
code "$BASE/api/similarity?smiles=quercetin&metric=maccs&top_n=5"; echo " <- maccs"
hr "32'. ChemBERTa (expect 200)"
code "$BASE/api/similarity?smiles=quercetin&metric=chemberta&top_n=5"; echo " <- chemberta"
hr "33'. metric disagreement: morgan vs nafm top-2 neighbours for curcumin"
echo "morgan:"; curl -s "$BASE/api/similarity?smiles=curcumin&metric=morgan&top_n=3" | grep -oE 'THEO_[0-9]+' | head -3 | tr '\n' ' '; echo
echo "nafm:  "; curl -s "$BASE/api/similarity?smiles=curcumin&metric=nafm&top_n=3"  | grep -oE 'THEO_[0-9]+' | head -3 | tr '\n' ' '; echo
echo "maccs: "; curl -s "$BASE/api/similarity?smiles=curcumin&metric=maccs&top_n=3" | grep -oE 'THEO_[0-9]+' | head -3 | tr '\n' ' '; echo

hr "30c. does q= also work, or is smiles= required? (UX check: try q=)"
code "$BASE/api/similarity?q=quercetin&metric=morgan&top_n=5"; echo " <- q= param (400 => smiles= required)"

# --- 35'. taxonomy coverage definitional split: kingdom vs deep-rank ---
hr "35'. resolution INCLUDING kingdom as a rank (should be ~100, matching Help's 99.8%)"
q "SELECT round(100.0*count(*) FILTER (WHERE (kingdom IS NOT NULL AND kingdom<>'') OR (genus IS NOT NULL AND genus<>'') OR (family IS NOT NULL AND family<>'') OR (taxorder IS NOT NULL AND taxorder<>'') OR (taxclass IS NOT NULL AND taxclass<>'') OR (phylum IS NOT NULL AND phylum<>''))/count(*),2) AS pct_incl_kingdom FROM resolved_taxonomy;"
hr "35b. resolution at genus level only (the honest 'deep' number)"
q "SELECT round(100.0*count(*) FILTER (WHERE genus IS NOT NULL AND genus<>'')/count(*),2) AS pct_genus,
          round(100.0*count(*) FILTER (WHERE family IS NOT NULL AND family<>'')/count(*),2) AS pct_family,
          round(100.0*count(*) FILTER (WHERE phylum IS NOT NULL AND phylum<>'')/count(*),2) AS pct_phylum
   FROM resolved_taxonomy;"

# --- 38'. homonym 'Homo' across 5 kingdoms: legitimate or artifact? inspect ---
hr "38'. what kingdoms is genus 'Homo' attached to, and how many compounds each?"
q "SELECT c.kingdom, count(*) FROM compound_taxonomy ct JOIN compounds c ON c.comp_id=ct.comp_id WHERE ct.genus='Homo' GROUP BY 1 ORDER BY 2 DESC;"

# --- 50'. substructure 302: where does it redirect? ---
hr "50'. substructure endpoint headers (is 302 a sane redirect to results?)"
curl -s -o /dev/null -w "status=%{http_code} redirect=%{redirect_url}\n" "$BASE/search?type=smiles&smiles_mode=substructure&q=c1ccccc1"

echo "=================== eval_50b DONE ==================="
