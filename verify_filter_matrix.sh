#!/usr/bin/env bash
# Every filter name x every route. Any cell equal to the unfiltered total is a silent drop.
B="https://theobroma.l3s.uni-hannover.de"
tot(){ curl -s "$1" | python3 -c 'import sys,json
try: print(json.load(sys.stdin).get("total") or json.load(sys.stdin).get("total_compounds"))
except: print("ERR")' 2>/dev/null; }
tree(){ curl -s "$1" | python3 -c 'import sys,json
try: print(json.load(sys.stdin).get("total_compounds"))
except: print("ERR")'; }
echo "name                     primary   as_extra   tree_extra"
for pair in "npclassifier_class:Flavonols" "npclassifier_superclass:Flavonoids" \
            "npclassifier_pathway:Terpenoids" "classyfire_superclass:Phenylpropanoids%20and%20polyketides" \
            "pathway:Terpenoids" "classyfire_class:Phenylpropanoids%20and%20polyketides" \
            "class:Flavonols" "genus:Salvia" "family:Lamiaceae" "order:lamiales" \
            "tax_class:liliopsida" "phylum:streptophyta" "region:Middle%20East" "source:LOTUS"; do
  t=${pair%%:*}; q=${pair##*:}
  a=$(tot "$B/api/search?type=$t&q=$q&limit=1")
  b=$(tot "$B/api/search?type=kingdom&q=plant&extra_type_1=$t&extra_q_1=$q&limit=1")
  c=$(tree "$B/api/taxonomy_tree?type=kingdom&q=plant&extra_type_1=$t&extra_q_1=$q")
  printf "%-24s %8s %10s %12s\n" "$t" "$a" "$b" "$c"
done
echo "(plant unfiltered = $(tot "$B/api/search?type=kingdom&q=plant&limit=1"))"
