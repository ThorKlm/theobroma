#!/usr/bin/env bash
# eval_combinations.sh -- 25 command-line probes of THEOBROMA search combinations:
# multi-filter intersection, result-count monotonicity, duplication, geography
# interaction, and edge cases. Uses the JSON API (/api/search) so counts are exact.
# Run:  bash eval_combinations.sh 2>&1 | tee eval_combinations_out.txt
set +e
B="http://localhost:5000"
# cnt <url> -> prints the "count" field from /api/search JSON
cnt(){ curl -s "$1" | python3 -c "import sys,json;
try:
 d=json.load(sys.stdin); print(d.get('total','ERR'))
except Exception as e: print('ERR:',e)"; }
code(){ curl -s -o /dev/null -w "%{http_code}" "$1"; }
hr(){ echo; echo "----- $1"; }

echo "=================== eval_combinations (25) ==================="
echo "date: $(date -u +%FT%TZ)"

# ---- Baselines (single filters) ----
hr "1. class=Triterpenoids alone (baseline)"
C_TRI=$(cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&limit=1"); echo "  Triterpenoids = $C_TRI"

hr "2. phylum=streptophyta alone (baseline)"
C_STR=$(cnt "$B/api/search?type=phylum&q=streptophyta&limit=1"); echo "  streptophyta = $C_STR"

hr "3. phylum=chlorophyta alone (baseline)"
C_CHL=$(cnt "$B/api/search?type=phylum&q=chlorophyta&limit=1"); echo "  chlorophyta = $C_CHL"

hr "4. region=Europe alone (baseline)"
C_EUR=$(cnt "$B/api/search?type=region&q=Europe&limit=1"); echo "  Europe = $C_EUR"

# ---- 2-way intersections: must be <= each baseline (monotonicity) ----
hr "5. class=Triterpenoids AND phylum=streptophyta (expect <= min(Tri, streptophyta))"
X=$(cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=phylum&extra_q_1=streptophyta&limit=1"); echo "  = $X   (Tri=$C_TRI, str=$C_STR)"

hr "6. class=Triterpenoids AND phylum=chlorophyta"
cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=phylum&extra_q_1=chlorophyta&limit=1" | sed 's/^/  = /'

hr "7. class=Triterpenoids AND region=Europe"
cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=region&extra_q_1=Europe&limit=1" | sed 's/^/  = /'

# ---- Mutually-exclusive taxa: two different phyla -> expect 0 ----
hr "8. phylum=streptophyta AND phylum=chlorophyta (mutually exclusive -> expect 0)"
cnt "$B/api/search?type=phylum&q=streptophyta&extra_type_1=phylum&extra_q_1=chlorophyta&limit=1" | sed 's/^/  = /'

# ---- 3-way and 4-way (the Caulerpa chain) ----
hr "9. Triterpenoids AND phylum=chlorophyta AND genus=Caulerpa (3-way)"
cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=phylum&extra_q_1=chlorophyta&extra_type_2=genus&extra_q_2=Caulerpa&limit=1" | sed 's/^/  = /'

hr "10. + region=Europe (4-way; expect the 2 from the browser test)"
cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=phylum&extra_q_1=chlorophyta&extra_type_2=genus&extra_q_2=Caulerpa&extra_type_3=region&extra_q_3=Europe&limit=1" | sed 's/^/  = /'

# ---- Order independence: same filters, different order -> same count ----
hr "11. order-independence: {genus,phylum} vs {phylum,genus} (should be EQUAL)"
A=$(cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=phylum&extra_q_1=chlorophyta&extra_type_2=genus&extra_q_2=Caulerpa&limit=1")
Bx=$(cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=genus&extra_q_1=Caulerpa&extra_type_2=phylum&extra_q_2=chlorophyta&limit=1")
echo "  phylum-first=$A  genus-first=$Bx  $([ "$A" = "$Bx" ] && echo OK || echo MISMATCH)"

# ---- Redundant/duplicate filter: same filter twice -> same as once ----
hr "12. duplicate filter: genus=Caulerpa given twice (should equal once)"
ONE=$(cnt "$B/api/search?type=genus&q=Caulerpa&limit=1")
DUP=$(cnt "$B/api/search?type=genus&q=Caulerpa&extra_type_1=genus&extra_q_1=Caulerpa&limit=1")
echo "  once=$ONE  twice=$DUP  $([ "$ONE" = "$DUP" ] && echo OK || echo DIFFERENT)"

# ---- Redundant nested taxa: genus + its family -> equals genus alone ----
hr "13. nested: genus=Caulerpa AND family=Caulerpaceae (Caulerpa is in Caulerpaceae -> ~equals genus)"
G=$(cnt "$B/api/search?type=genus&q=Caulerpa&limit=1")
GF=$(cnt "$B/api/search?type=genus&q=Caulerpa&extra_type_1=family&extra_q_1=Caulerpaceae&limit=1")
echo "  genus=$G  genus+family=$GF  $([ "$G" = "$GF" ] && echo CONSISTENT || echo 'differs (some Caulerpa rows lack family=Caulerpaceae)')"

# ---- Kingdom vs phylum consistency ----
hr "14. kingdom=plant AND phylum=chlorophyta (chlorophyta is a plant-kingdom phylum -> should be > 0)"
cnt "$B/api/search?type=kingdom&q=plant&extra_type_1=phylum&extra_q_1=chlorophyta&limit=1" | sed 's/^/  = /'

hr "15. kingdom=animal AND phylum=chlorophyta (contradiction -> expect 0)"
cnt "$B/api/search?type=kingdom&q=animal&extra_type_1=phylum&extra_q_1=chlorophyta&limit=1" | sed 's/^/  = /'

# ---- Property + property (two chemical dimensions) ----
hr "16. class=Triterpenoids AND pathway=Terpenoids (coherent -> should be large)"
cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=pathway&extra_q_1=Terpenoids&limit=1" | sed 's/^/  = /'

hr "17. class=Triterpenoids AND pathway=Alkaloids (incoherent -> expect small/0)"
cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=pathway&extra_q_1=Alkaloids&limit=1" | sed 's/^/  = /'

# ---- License filter interaction ----
hr "18. Triterpenoids commercial-only (CC0+CC BY) vs all (commercial <= all)"
ALL=$(cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&limit=1")
COM=$(cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&license=commercial&limit=1")
echo "  all=$ALL  commercial=$COM  $([ "$COM" -le "$ALL" ] 2>/dev/null && echo OK || echo CHECK)"

# ---- Region interaction: region + region (two regions) ----
hr "19. region=Europe AND region=Africa (a compound can be in both -> may be >0)"
cnt "$B/api/search?type=region&q=Europe&extra_type_1=region&extra_q_1=Africa&limit=1" | sed 's/^/  = /'

# ---- Case-insensitivity of taxon values ----
hr "20. case-insensitivity: phylum=Chlorophyta vs chlorophyta (should be equal)"
U=$(cnt "$B/api/search?type=phylum&q=Chlorophyta&limit=1")
L=$(cnt "$B/api/search?type=phylum&q=chlorophyta&limit=1")
echo "  Chlorophyta=$U  chlorophyta=$L  $([ "$U" = "$L" ] && echo OK || echo CASE-SENSITIVE)"

# ---- Result-set duplication: are comp_ids unique within a result page? ----
hr "21. duplication in results: distinct comp_ids vs returned rows (expect equal)"
curl -s "$B/api/search?type=npclassifier_class&q=Triterpenoids&limit=50" | python3 -c "
import sys,json
d=json.load(sys.stdin); r=d.get('results',[])
ids=[x.get('comp_id') for x in r]
print('  rows=%d distinct=%d %s' % (len(ids), len(set(ids)), 'OK' if len(ids)==len(set(ids)) else 'DUPLICATES'))"

# ---- Pagination consistency: page1 + page2 disjoint ----
hr "22. pagination: page 1 and page 2 have no overlapping comp_ids"
P1=$(curl -s "$B/api/search?type=npclassifier_class&q=Triterpenoids&limit=50&offset=0" | python3 -c "import sys,json;print(' '.join(x['comp_id'] for x in json.load(sys.stdin).get('results',[])))")
P2=$(curl -s "$B/api/search?type=npclassifier_class&q=Triterpenoids&limit=50&offset=50" | python3 -c "import sys,json;print(' '.join(x['comp_id'] for x in json.load(sys.stdin).get('results',[])))")
python3 -c "
p1=set('''$P1'''.split()); p2=set('''$P2'''.split())
ov=p1&p2
print('  page1=%d page2=%d overlap=%d %s' % (len(p1),len(p2),len(ov),'OK' if not ov else 'OVERLAP'))"

# ---- Empty / contradictory / robustness ----
hr "23. contradictory chem: class=Triterpenoids AND class=Alkaloids (single compound can't be both primary+extra same dim? expect small/0)"
cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=npclassifier_class&extra_q_1=Alkaloids&limit=1" | sed 's/^/  = /'

hr "24. nonexistent taxon: genus=Notarealgenus (expect 0, not error)"
echo "  count=$(cnt "$B/api/search?type=genus&q=Notarealgenus&limit=1")  http=$(code "$B/api/search?type=genus&q=Notarealgenus&limit=1")"

hr "25. deep 6-filter chain (the Codium lineage; expect small, no error/0-from-bug)"
echo "  count=$(cnt "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=genus&extra_q_1=Codium&extra_type_2=kingdom&extra_q_2=plant&extra_type_3=phylum&extra_q_3=chlorophyta&extra_type_4=tax_class&extra_q_4=ulvophyceae&extra_type_5=order&extra_q_5=bryopsidales&limit=1")  http=$(code "$B/api/search?type=npclassifier_class&q=Triterpenoids&extra_type_1=genus&extra_q_1=Codium&limit=1")"

echo; echo "=================== DONE ==================="
