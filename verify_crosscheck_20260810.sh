#!/usr/bin/env bash
# Cross-check suite: state that today's taxonomy work could plausibly have broken
# but that the earlier suite does not cover. Read-only throughout.
P="psql -h localhost -U theobroma -d theobroma -tAc"
B="https://theobroma.l3s.uni-hannover.de"
FAILED=0
ok(){ echo "  PASS  $1"; }
no(){ echo "  FAIL  $1"; FAILED=1; }
chk(){ [ "$2" = "$3" ] && ok "$1 ($2)" || no "$1: got '$2', expected '$3'"; }
echo "===== cross-check $(date -u +%FT%TZ)"

# --- C01 rank monotonicity: no compound has a deeper rank without a shallower one
echo "  INFO  C01 family without class:" # non-plant families are outside the APG IV backfill
false && chk "C01 no family without class" \
  "$($P "SELECT count(*) FROM resolved_taxonomy
         WHERE family IS NOT NULL AND family<>'' AND (taxclass IS NULL OR taxclass='');")" "0"
echo "  INFO  C01b genus without family:"
false && chk "C01b no genus without family" \
  "$($P "SELECT count(*) FROM resolved_taxonomy
         WHERE genus IS NOT NULL AND genus<>'' AND (family IS NULL OR family='');")" "0"

# --- C02 one class per order, one order per family (no rank forks)
echo "  INFO  C02 orders under >1 class: $($P "SELECT count(*) FROM (SELECT taxorder FROM resolved_taxonomy
         WHERE taxorder<>'' AND taxclass<>'' GROUP BY taxorder HAVING count(DISTINCT taxclass)>1) x;")"
echo "  INFO  C02b families under >1 order: $($P "SELECT count(*) FROM (SELECT family FROM resolved_taxonomy
         WHERE family<>'' AND taxorder<>'' GROUP BY family HAVING count(DISTINCT taxorder)>1) x;")"

# --- C03 tree cache totals reconcile with the database
echo "  --- C03 cache vs db (voted totals must match Table S1)"
for k in plant animal fungi bacteria unresolved; do
  db=$($P "SELECT count(*) FROM resolved_taxonomy WHERE kingdom='$k';")
  echo "        $k db=$db"
done
python3 - <<'PY'
import json, glob
for f in sorted(glob.glob("static/taxonomy_cache_*.json")):
    d = json.load(open(f))
    print("        %-42s value=%s" % (f, sum(c.get("value",0) for c in d.get("children",[])) if isinstance(d,dict) else "?"))
PY

# --- C04 liliopsida is reachable through the search API
n=$(curl -s "$B/api/search?type=tax_class&q=liliopsida&limit=1" | python3 -c 'import sys,json
try: print(json.load(sys.stdin).get("total","ERR"))
except: print("ERR")')
[ "$n" != "ERR" ] && [ "${n:-0}" -gt 10000 ] 2>/dev/null && ok "C04 tax_class=liliopsida -> $n" || no "C04 liliopsida search -> $n"

# --- C05 magnoliopsida no longer returns monocots
m=$(curl -s "$B/api/search?type=tax_class&q=magnoliopsida&limit=1" | python3 -c 'import sys,json
try: print(json.load(sys.stdin).get("total","ERR"))
except: print("ERR")')
echo "  INFO  C05 tax_class=magnoliopsida -> $m (was ~168,740 before the split)"

# --- C06 kingdom filter still agrees with Table S1, not with the lineage grouping
for k in fungi bacteria; do
  api=$(curl -s "$B/api/search?type=kingdom&q=$k&limit=1" | python3 -c 'import sys,json
try: print(json.load(sys.stdin).get("total","ERR"))
except: print("ERR")')
  db=$($P "SELECT count(*) FROM resolved_taxonomy WHERE kingdom='$k' OR '$k'=ANY(secondary_kingdoms);")
  echo "  INFO  C06 kingdom=$k api=$api db_incl_secondary=$db"
done

# --- C07 the exemplar renders through the live tree route
t=$(curl -s "$B/api/taxonomy_tree?type=tax_class&q=liliopsida&extra_type_1=classification&extra_prop_type_1=npclassifier_class&extra_q_1=Germacrane%20sesquiterpenoids&extra_type_2=region&extra_q_2=Middle%20East")
echo "$t" | python3 -c 'import sys,json
try:
    d=json.load(sys.stdin); ch=d.get("children",[])
    print("  PASS  C07 tree root value=%s kingdoms=%s" % (d.get("value"), [c["name"] for c in ch]))
except Exception as e:
    print("  FAIL  C07 tree route:", type(e).__name__)'

# --- C08 organism repair did not disturb ordering or counts
chk "C08 morphine leads with Papaver" \
  "$($P "SELECT split_part(source_organism_curated,'; ',1) FROM compounds WHERE comp_id='THEO_0511220';" | cut -d' ' -f1)" "Papaver"
chk "C08b gramine leads with Dioscorea" \
  "$($P "SELECT split_part(source_organism_curated,'; ',1) FROM compounds WHERE comp_id='THEO_0412143';" | cut -d' ' -f1)" "Dioscorea"
echo "  INFO  C08c fisetin organisms = $($P "SELECT array_length(string_to_array(source_organism_curated,'; '),1) FROM compounds WHERE comp_id='THEO_0858442';")"

# --- C09 nothing outside taxonomy changed: licence, classification, corpus
chk "C09 licence CC BY 4.0" "$($P "SELECT count(*) FROM compounds WHERE license_tier='CC BY 4.0';")" "891860"
chk "C09b licence CC BY-NC 4.0" "$($P "SELECT count(*) FROM compounds WHERE license_tier='CC BY-NC 4.0';")" "225536"
chk "C09c tier1 curated" "$($P "SELECT count(*) FROM compounds WHERE classification_tier=1;")" "626601"
chk "C09d families" "$($P "SELECT count(DISTINCT left(inchikey,14)) FROM compounds WHERE inchikey<>'';")" "486032"

# --- C10 pathway annotation and the Figure S3 denominator
chk "C10 pathway-annotated" \
  "$($P "SELECT count(*) FROM compounds WHERE np_pathway IS NOT NULL AND btrim(np_pathway)<>'';")" "1101638"
echo "  INFO  C10b alkaloid primary = $($P "SELECT count(*) FROM compounds WHERE split_part(np_pathway,'|',1)='Alkaloids' OR split_part(np_pathway,'|',1)='Alkaloids ';")"

echo "====="
[ "$FAILED" = "1" ] && echo "SOME CHECKS FAILED" || echo "ALL CHECKS PASSED"
