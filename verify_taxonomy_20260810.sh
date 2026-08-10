#!/usr/bin/env bash
# Verify the taxonomy repairs of 2026-08-10 before commit. Read-only.
P="psql -h localhost -U theobroma -d theobroma -tAc"
ok(){ echo "  PASS  $1"; }
no(){ echo "  FAIL  $1"; FAILED=1; }
chk(){ [ "$2" = "$3" ] && ok "$1 ($2)" || no "$1: got $2, expected $3"; }
FAILED=0
echo "===== taxonomy repair verification $(date -u +%FT%TZ)"

chk "T01 corpus unchanged" "$($P "SELECT count(*) FROM compounds;")" "1132805"

chk "T02 no plant tracheophyta" \
    "$($P "SELECT count(*) FROM resolved_taxonomy WHERE kingdom='plant' AND phylum='tracheophyta';")" "0"

chk "T03 no monocot order under magnoliopsida" \
    "$($P "SELECT count(*) FROM resolved_taxonomy WHERE kingdom='plant' AND taxclass='magnoliopsida'
           AND taxorder IN ('acorales','alismatales','arecales','asparagales','commelinales',
                            'dioscoreales','liliales','pandanales','poales','zingiberales','petrosaviales');")" "0"

echo "  INFO  T04 liliopsida compounds = $($P "SELECT count(DISTINCT comp_id) FROM resolved_taxonomy WHERE taxclass='liliopsida';")"

chk "T05 every phylum mapped" \
    "$($P "SELECT count(*) FROM (SELECT DISTINCT rt.phylum FROM resolved_taxonomy rt
           LEFT JOIN phylum_kingdom_map m ON m.phylum=rt.phylum
           WHERE rt.phylum IS NOT NULL AND rt.phylum<>'' AND m.phylum IS NULL) x;")" "0"

chk "T06 no plant family split across phyla" \
    "$($P "SELECT count(*) FROM (SELECT family FROM resolved_taxonomy
           WHERE kingdom='plant' AND family IS NOT NULL AND family<>'' AND phylum IS NOT NULL AND phylum<>''
           GROUP BY family HAVING count(DISTINCT phylum)>1) x;")" "0"

chk "T07 no hyphen-truncated binomials" \
    "$($P "SELECT count(*) FROM (SELECT 1 FROM compounds c
           CROSS JOIN LATERAL unnest(string_to_array(c.source_organism,'; ')) AS r(tok)
           CROSS JOIN LATERAL unnest(string_to_array(c.source_organism_curated,'; ')) AS k(tok)
           WHERE c.source_organism ~ '[a-z]+-[a-z]+' AND r.tok LIKE k.tok||'-%' AND r.tok<>k.tok) x;")" "0"

echo "  INFO  T08 fenugreek=$($P "SELECT count(*) FROM compounds WHERE source_organism_curated LIKE '%foenum-graecum%';") nux-vomica=$($P "SELECT count(*) FROM compounds WHERE source_organism_curated LIKE '%nux-vomica%';")"

echo "  --- T09 Table S1 kingdoms (voted, must be unchanged)"
$P "SELECT '        '||kingdom||' '||count(*) FROM resolved_taxonomy GROUP BY kingdom ORDER BY count(*) DESC;"

echo "  --- T10 tree grouping: lineage vs voted disagreement"
$P "SELECT '        '||rt.kingdom||' -> '||m.lineage_kingdom||' : '||count(*)
    FROM resolved_taxonomy rt JOIN phylum_kingdom_map m ON m.phylum=rt.phylum
    WHERE rt.kingdom<>m.lineage_kingdom GROUP BY rt.kingdom, m.lineage_kingdom ORDER BY count(*) DESC LIMIT 5;"

chk "T11 one kingdom branch per Fusarium lineage" \
    "$($P "SELECT count(DISTINCT COALESCE(m.lineage_kingdom, rt.kingdom)) FROM resolved_taxonomy rt
           LEFT JOIN phylum_kingdom_map m ON m.phylum=rt.phylum WHERE rt.genus='Fusarium';")" "1"

echo "  --- T12 Germacrane exemplar under liliopsida"
$P "SELECT '        n='||count(DISTINCT c.comp_id)||' ord='||count(DISTINCT rt.taxorder)||
           ' fam='||count(DISTINCT rt.family)||' gen='||count(DISTINCT rt.genus)
    FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id=c.comp_id
    JOIN compound_region_map g ON g.comp_id=c.comp_id
    WHERE g.macro_region='Middle East' AND rt.taxclass='liliopsida'
      AND (c.np_class ILIKE '%Germacrane sesquiterpenoids%' OR c.inferred_class ILIKE '%Germacrane sesquiterpenoids%');"

echo "  --- live site"
curl -s https://theobroma.l3s.uni-hannover.de/ | grep -q "maint-banner" && ok "T13 maintenance banner on" || no "T13 banner missing"
curl -s https://theobroma.l3s.uni-hannover.de/compound/THEO_0858442 | grep -c "showmore" \
  | xargs -I{} sh -c '[ {} -ge 4 ] && echo "  PASS  T14 list expanders ({})" || echo "  FAIL  T14 expanders={}"'

echo "====="
[ "$FAILED" = "1" ] && echo "SOME CHECKS FAILED" || echo "ALL CHECKS PASSED"
