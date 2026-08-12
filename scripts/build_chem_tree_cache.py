#!/usr/bin/env python3
"""Offline builder for the NPClassifier chemistry-tree cache (static/chem_tree.json).
Run at build/deploy time so the Tree page never builds during a request.
Mirrors update_taxonomy_cache.py: shells out to psql (no psycopg2 import), so it
runs under plain python3. Structure: pathway -> superclass -> class, each node with
an independent DISTINCT-compound count (multi-labels split on ' | '/' $ ')."""
import json, os, subprocess

OUT = "/home/thorben.klamt/theobroma/static/chem_tree.json"

def run_sql(q):
    r = subprocess.run(
        ['sudo','-u','postgres','psql','-d','theobroma','-At','-F','\t','-c', q],
        capture_output=True, text=True, check=True)
    return [line.split('\t') for line in r.stdout.strip().split('\n') if line]

def bulk_counts(col):
    q = ("SELECT trim(v), count(DISTINCT comp_id) FROM compounds, "
         "regexp_split_to_table(%s, ' *[|$] *') AS v "
         "WHERE %s IS NOT NULL AND %s<>'' AND trim(v)<>'' GROUP BY trim(v)" % (col, col, col))
    return {r[0]: int(r[1]) for r in run_sql(q)}

def main():
    super_rows = run_sql("SELECT pathway_name, superclass_name FROM npc_super_parents")
    class_rows = run_sql("SELECT superclass_name, class_name FROM npc_class_parents WHERE parent_rank = 1")
    pw_cnt = bulk_counts("effective_pathway")
    sc_cnt = bulk_counts("effective_superclass")
    cl_cnt = bulk_counts("effective_class")

    classes_by_super = {}
    for sc, cl in class_rows:
        classes_by_super.setdefault(sc, set()).add(cl)
    pathways = {}
    for pw, sc in super_rows:
        pathways.setdefault(pw, {"name": pw, "count": pw_cnt.get(pw, 0), "superclasses": {}})
        if sc not in pathways[pw]["superclasses"]:
            cls = sorted(({"name": c, "count": cl_cnt.get(c, 0)} for c in classes_by_super.get(sc, [])),
                         key=lambda x: -x["count"])
            pathways[pw]["superclasses"][sc] = {"name": sc, "count": sc_cnt.get(sc, 0), "classes": cls}
    out, top_sc, top_n = [], None, -1
    for pw, pd in pathways.items():
        sc_list = sorted(pd["superclasses"].values(), key=lambda x: -x["count"])
        for sd in sc_list:
            if sd["count"] > top_n: top_n, top_sc = sd["count"], sd["name"]
        out.append({"name": pd["name"], "count": pd["count"], "superclasses": sc_list})
    out.sort(key=lambda x: -x["count"])
    result = {"pathways": out, "top_superclass": top_sc}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f: json.dump(result, f)
    print("# wrote %s (%d pathways, top superclass=%s)" % (OUT, len(out), top_sc))

if __name__ == "__main__":
    main()
