"""Export the native/curated training labels for the v1.35 classifier retrain.

Native tier definition (verified against the restored corpus): a compound is a
curated/native training label iff it has a source-inherited np_class AND the
XGBoost model did not also classify it (inferred_class_source blank), excluding
dollar-joined multi-label composites. This reconstructs the manuscript's curated
tier (~397k) source-agnostically and supersedes the stale np_classifications.csv
artifact-membership rule (whose comp_ids no longer align to the current corpus).

Emits one row per native compound: comp_id, np_class, origin=native, and the
NPClassifier-ontology parents / consistency flags from npc_class_parents.
No row_idx here: feature-row assignment happens in the bundle step against the
(regenerated) vectors. Read-only. Run on the production node.
"""
import argparse, csv, json, os
import psycopg2, psycopg2.extras

DSN = os.environ.get("DATABASE_URL", "postgresql://theobroma:theobroma@localhost:5432/theobroma")

def fetch(conn):
    sql = """
        SELECT c.comp_id, c.np_class, c.np_superclass, c.np_pathway,
               array_agg(DISTINCT p.superclass_name) FILTER (WHERE p.superclass_name IS NOT NULL) AS onto_super,
               array_agg(DISTINCT p.pathway_name)    FILTER (WHERE p.pathway_name    IS NOT NULL) AS onto_path
        FROM compounds c
        LEFT JOIN npc_class_parents p ON p.class_name = c.np_class
        WHERE c.np_class IS NOT NULL AND c.np_class <> ''
          AND (c.inferred_class_source IS NULL OR c.inferred_class_source = '')
          AND c.np_class NOT LIKE '%$%'
        GROUP BY c.comp_id, c.np_class, c.np_superclass, c.np_pathway
    """
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        cur.execute(sql)
        return cur.fetchall()

def pathway_consistent(stored, onto):
    if not stored or not onto:
        return None
    members = {m.strip() for m in stored.split("|")}
    return bool(members & set(onto))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="labels_v135b_native.tsv")
    args = ap.parse_args()
    conn = psycopg2.connect(DSN); rows = fetch(conn); conn.close()
    print(f"[db] {len(rows):,} native single-label curated records")
    out, ambiguous = [], 0
    for r in rows:
        onto_path, onto_super = r["onto_path"] or [], r["onto_super"] or []
        if len(onto_path) > 1:
            ambiguous += 1
        out.append({
            "comp_id": r["comp_id"],
            "np_class": r["np_class"],
            "origin": "native",
            "onto_superclass": onto_super[0] if len(onto_super) == 1 else "",
            "onto_pathway": onto_path[0] if len(onto_path) == 1 else "",
            "onto_ambiguous": int(len(onto_path) > 1),
            "in_ontology": int(bool(onto_path)),
            "stored_pathway": r["np_pathway"] or "",
            "pathway_consistent": {True: 1, False: 0, None: -1}[pathway_consistent(r["np_pathway"], onto_path)],
        })
    fields = list(out[0].keys())
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(out)
    summary = {"exported": len(out), "ontology_ambiguous": ambiguous,
               "not_in_ontology": sum(1 for r in out if not r["in_ontology"]),
               "pathway_inconsistent": sum(1 for r in out if r["pathway_consistent"] == 0)}
    print(json.dumps(summary, indent=1))
    with open(args.out.replace(".tsv", "_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

if __name__ == "__main__":
    main()
