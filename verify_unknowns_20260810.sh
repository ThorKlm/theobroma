#!/usr/bin/env bash
# Probes for defects of the same class as those found today, in places not yet checked.
P="psql -h localhost -U theobroma -d theobroma -tAc"
B="https://theobroma.l3s.uni-hannover.de"
api(){ curl -s "$B/api/search?$1&limit=1" | python3 -c 'import sys,json
try: print(json.load(sys.stdin).get("total","ERR"))
except: print("ERR")'; }
echo "===== unknowns probe $(date -u +%FT%TZ)"

# U01 every filter type the search page offers actually filters
echo "  --- U01 filter types: does each narrow the corpus?"
for t in name organism kingdom source region genus family order tax_class clade phylum \
         npclassifier_class npclassifier_superclass npclassifier_pathway classyfire_superclass; do
  case $t in
    name) q=curcumin;; organism) q=Curcuma;; kingdom) q=fungi;; source) q=LOTUS;;
    region) q=Middle+East;; genus) q=Salvia;; family) q=Lamiaceae;; order) q=lamiales;;
    tax_class|clade) q=liliopsida;; phylum) q=streptophyta;;
    npclassifier_class) q=Flavonols;; npclassifier_superclass) q=Flavonoids;;
    npclassifier_pathway) q=Alkaloids;; classyfire_superclass) q=Flavonoids;;
  esac
  n=$(api "type=$t&q=$q")
  [ "$n" = "1132805" ] && echo "        FAIL $t=$q -> $n (no filtering)" || echo "        ok   $t=$q -> $n"
done

# U02 the same filter as primary and as extra must agree
echo "  --- U02 primary vs extra agreement"
for t in npclassifier_class family region source; do
  case $t in npclassifier_class) q=Flavonols;; family) q=Lamiaceae;; region) q=Middle+East;; source) q=LOTUS;; esac
  a=$(api "type=$t&q=$q"); b=$(api "type=name&q=&extra_type_1=$t&extra_q_1=$q")
  echo "        $t: primary=$a extra=$b"
done

# U03 AND semantics: composed result must be a subset of each component
echo "  --- U03 composition is intersective"
a=$(api "type=tax_class&q=liliopsida")
b=$(api "type=region&q=Middle+East")
c=$(api "type=tax_class&q=liliopsida&extra_type_1=region&extra_q_1=Middle+East")
echo "        liliopsida=$a  middleeast=$b  both=$c  (both must be <= min)"

# U04 model-tier compounds reachable by every classification rank
echo "  --- U04 tier-3 visibility"
$P "SELECT '        tier3 with effective_class: '||count(*) FROM compounds
    WHERE classification_tier=3 AND effective_class<>'';"
$P "SELECT '        tier3 with np_class:        '||count(*) FROM compounds
    WHERE classification_tier=3 AND COALESCE(np_class,'')<>'';"

# U05 columns the app filters on that are empty for a whole tier
echo "  --- U05 per-tier column population"
$P "SELECT '        tier '||classification_tier||': np_class='||count(*) FILTER (WHERE COALESCE(np_class,'')<>'')
        ||' effective='||count(*) FILTER (WHERE COALESCE(effective_class,'')<>'')
        ||' superclass='||count(*) FILTER (WHERE COALESCE(effective_superclass,'')<>'')
    FROM compounds WHERE classification_tier IS NOT NULL GROUP BY classification_tier ORDER BY 1;"

# U06 route smoke test: nothing 500s after today's schema change
echo "  --- U06 route smoke test"
cid=$($P "SELECT comp_id FROM compounds LIMIT 1;")
for r in "/" "/search" "/browse" "/tree" "/annotate" "/help" "/statistics" "/scaffolds" "/similarity" \
         "/compound/$cid" "/api/stats" "/api/compound/$cid" "/api/compound/$cid/license-provenance" \
         "/api/stereoisomers/$cid" "/api/depict?smiles=CCO" "/api/taxonomy_tree?type=genus&q=Salvia"; do
  code=$(curl -s -o /dev/null -w "%{http_code}" "$B$r")
  [ "$code" = "200" ] && echo "        ok   $code $r" || echo "        FAIL $code $r"
done

# U07 dropdown vocabularies match the data they filter
echo "  --- U07 distinct values the UI offers vs the data"
$P "SELECT '        effective_class distinct: '||count(DISTINCT effective_class) FROM compounds WHERE effective_class<>'';"
$P "SELECT '        np_class distinct:        '||count(DISTINCT np_class) FROM compounds WHERE COALESCE(np_class,'')<>'';"
$P "SELECT '        taxclass distinct:        '||count(DISTINCT taxclass) FROM resolved_taxonomy WHERE taxclass<>'';"

# U08 whitespace and separator artifacts in every text column the UI filters on
echo "  --- U08 stray whitespace in filterable columns"
for col in effective_class effective_superclass effective_pathway np_class np_superclass np_pathway classyfire_superclass; do
  n=$($P "SELECT count(*) FROM compounds WHERE $col IS NOT NULL AND $col <> btrim($col);")
  echo "        $col untrimmed: $n"
done
$P "SELECT '        resolved_taxonomy untrimmed: '||
    (SELECT count(*) FROM resolved_taxonomy WHERE genus<>btrim(genus) OR family<>btrim(family)
      OR taxorder<>btrim(taxorder) OR taxclass<>btrim(taxclass) OR phylum<>btrim(phylum));"

# U09 case consistency: are rank values stored in one case?
echo "  --- U09 case consistency in taxonomy"
$P "SELECT '        taxclass with uppercase: '||count(*) FROM resolved_taxonomy WHERE taxclass <> lower(taxclass);"
$P "SELECT '        taxorder with uppercase: '||count(*) FROM resolved_taxonomy WHERE taxorder <> lower(taxorder);"
$P "SELECT '        family all-lowercase:    '||count(*) FROM resolved_taxonomy WHERE family <> '' AND family = lower(family);"

# U10 orphan rows in every child table after today's updates
echo "  --- U10 referential integrity"
for t in compound_taxonomy compound_region_map per_source_license_attestation resolved_taxonomy admet scaffolds; do
  n=$($P "SELECT count(*) FROM $t t LEFT JOIN compounds c ON c.comp_id=t.comp_id WHERE c.comp_id IS NULL;" 2>/dev/null)
  echo "        $t orphans: ${n:-n/a}"
done
echo "====="
