"""Export provenance-flagged NPClassifier labels for the v1.35 classifier rebuild.

Emits one row per corpus compound carrying a usable single-label class, annotated
with label origin (native source annotation versus derived from the suspect
np_classifications.csv artifact), ontology consistency against npc_class_parents,
and the row index into the precomputed feature arrays.

Runs on the production node against the live database. Read-only.
"""
import argparse, csv, json, os, sys
import numpy as np
import psycopg2, psycopg2.extras

DSN = os.environ.get("DATABASE_URL", "postgresql://theobroma:theobroma@localhost:5432/theobroma")

def load_suspect_ids(path):
    """comp_ids for which the artifact supplies a non-empty class value."""
    ids = set()
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if (row.get("np_class") or "").strip():
                ids.add(row["comp_id"].strip())
    return ids

def fetch_labels(conn):
    """Single-label curated-tier records with their ontology parents.

    Composite classes (dollar-joined multi-label) are excluded: v1.34 trained
    single-label and excluding them keeps the comparison honest. Classes with
    multiple ontology parents are flagged rather than silently collapsed.
    """
    sql = """
        SELECT c.comp_id, c.np_class, c.np_superclass, c.np_pathway,
               array_agg(DISTINCT p.superclass_name) FILTER (WHERE p.superclass_name IS NOT NULL) AS onto_super,
               array_agg(DISTINCT p.pathway_name)    FILTER (WHERE p.pathway_name    IS NOT NULL) AS onto_path
        FROM compounds c
        LEFT JOIN npc_class_parents p ON p.class_name = c.np_class
        WHERE c.np_class IS NOT NULL AND c.np_class NOT LIKE '%$%'
        GROUP BY c.comp_id, c.np_class, c.np_superclass, c.np_pathway
    """
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        cur.execute(sql)
        return cur.fetchall()

def pathway_consistent(stored, onto):
    """Stored pathway agrees if any ontology parent appears among its pipe-joined members."""
    if not stored or not onto:
        return None
    members = {m.strip() for m in stored.split("|")}
    return bool(members & set(onto))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vectors", default="data/vectors")
    ap.add_argument("--artifact", default="data/vectors/np_classifications.csv")
    ap.add_argument("--out", default="data/vectors/labels_v135.tsv")
    args = ap.parse_args()
    comp_ids = np.load(os.path.join(args.vectors, "comp_ids.npy"), allow_pickle=True).astype(str)
    vi = np.load(os.path.join(args.vectors, "valid_indices.npy"), allow_pickle=True)
    # feature row i corresponds to comp_ids[vi[i]]; the arrays are in ingestion order
    row_of = {c: i for i, c in enumerate(comp_ids[vi.astype(np.int64)])}
    print(f"[feat] {len(comp_ids):,} feature rows")
    suspect = load_suspect_ids(args.artifact)
    print(f"[art ] {len(suspect):,} comp_ids with a class value in the artifact")
    conn = psycopg2.connect(DSN)
    rows = fetch_labels(conn)
    conn.close()
    print(f"[db  ] {len(rows):,} single-label curated records")
    out, skipped_nofeat, ambiguous = [], 0, 0
    for r in rows:
        idx = row_of.get(r["comp_id"])
        if idx is None:
            skipped_nofeat += 1
            continue
        onto_path, onto_super = r["onto_path"] or [], r["onto_super"] or []
        if len(onto_path) > 1:
            ambiguous += 1
        out.append({
            "comp_id": r["comp_id"],
            "row_idx": idx,
            "np_class": r["np_class"],
            "origin": "artifact" if r["comp_id"] in suspect else "native",
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
        w.writeheader()
        w.writerows(out)
    n_native = sum(1 for r in out if r["origin"] == "native")
    summary = {
        "exported": len(out), "native": n_native, "artifact": len(out) - n_native,
        "no_feature_row": skipped_nofeat, "ontology_ambiguous": ambiguous,
        "not_in_ontology": sum(1 for r in out if not r["in_ontology"]),
        "pathway_inconsistent_native": sum(1 for r in out if r["origin"] == "native" and r["pathway_consistent"] == 0),
        "pathway_inconsistent_artifact": sum(1 for r in out if r["origin"] == "artifact" and r["pathway_consistent"] == 0),
    }
    print(json.dumps(summary, indent=1))
    with open(args.out.replace(".tsv", "_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

if __name__ == "__main__":
    main()
