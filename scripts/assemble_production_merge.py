"""Assemble production merge inputs from Phase B bundle.

Combines three sources to produce /tmp/tipdb_novel_rows.csv (472 rows in
the production compounds schema) and /tmp/tipdb_overlap_inchikeys.txt
(7,495 InChIKeys for the enrichment update):

  Source A: data/theobroma_v33_merged.csv (locally produced merge output;
            carries the 472 novel rows with converter-side fields).
  Source B: data/tipdb_phase_b/ (vast.ai bundle: ChemBERTa, Morgan, MACCS,
            scaffolds, np_classifications, admet_predictions, sa_scores).
  Source C: data/theobroma_v33_clean.csv (post-cleanup corpus from L3S;
            carries the 33-column production schema we must match).

Output files:
  /tmp/tipdb_novel_rows.csv         (CSV HEADER, 472 rows, 33 columns)
  /tmp/tipdb_overlap_inchikeys.txt  (7,495 InChIKeys, one per line)

The 472 novel rows go into the live database via tipdb_production_merge.sql
which assigns final THEO_ ids server-side (avoiding the locally-inflated
ids from the 11-phantom-row issue).

Run on the L3S server after the Phase B bundle has been transferred:
  cd ~/theobroma
  venv/bin/python scripts/assemble_production_merge.py
"""
import os, sys
import numpy as np
import pandas as pd

REPO = os.path.expanduser("~/theobroma")
MERGED = os.path.join(REPO, "data", "theobroma_v33_merged.csv")
BUNDLE = os.path.join(REPO, "data", "tipdb_phase_b")
CORPUS = os.path.join(REPO, "data", "theobroma_v33_clean.csv")

OUT_NOVEL = "/tmp/tipdb_novel_rows.csv"
OUT_OVERLAP = "/tmp/tipdb_overlap_inchikeys.txt"

def main():
    if not os.path.exists(MERGED):
        sys.exit(f"merged CSV not found: {MERGED}; run merge_tipdb.py first")
    if not os.path.exists(BUNDLE):
        sys.exit(f"Phase B bundle not found: {BUNDLE}; run phase_b_tipdb.py "
                 f"on vast.ai first and transfer the output back")
    if not os.path.exists(CORPUS):
        sys.exit(f"post-cleanup corpus not found: {CORPUS}")

    print(f"[load] post-cleanup corpus schema from {CORPUS}")
    corpus_cols = pd.read_csv(CORPUS, nrows=0).columns.tolist()
    print(f"  schema: {len(corpus_cols)} columns")

    print(f"[load] merged CSV {MERGED}")
    merged = pd.read_csv(MERGED, low_memory=False)
    novel = merged[merged["source_db"] == "TIPdb"].copy()
    print(f"  novel rows: {len(novel):,}")

    overlap_iks = (merged["all_sources"].astype(str).str.contains("TIPdb",
                   regex=False, na=False) & (merged["source_db"] != "TIPdb"))
    overlap_inchikeys = merged.loc[overlap_iks, "inchikey"].dropna().unique()
    print(f"  overlap InChIKeys: {len(overlap_inchikeys):,}")

    print(f"[load] Phase B bundle from {BUNDLE}")
    pb_compids = np.load(os.path.join(BUNDLE, "comp_ids.npy"),
                         allow_pickle=True).astype(str)
    print(f"  bundle comp_ids: {len(pb_compids):,}")

    npc = pd.read_csv(os.path.join(BUNDLE, "np_classifications.csv"))
    print(f"  np_classifications: {len(npc):,} rows")
    scaff = pd.read_csv(os.path.join(BUNDLE, "scaffolds.csv"))
    print(f"  scaffolds: {len(scaff):,} rows")
    sa = pd.read_csv(os.path.join(BUNDLE, "sa_scores.csv"))
    print(f"  sa_scores: {len(sa):,} rows")
    admet = pd.read_csv(os.path.join(BUNDLE, "admet_predictions.csv"))
    print(f"  admet_predictions: {len(admet):,} rows x {len(admet.columns)} cols")

    # Merge bundle outputs into novel rows by comp_id. The local-side
    # comp_ids in `novel` and bundle-side comp_ids should match since they
    # came from the same merge_tipdb.py run.
    print("[merge] joining bundle into novel rows")
    novel = novel.merge(npc, on="comp_id", how="left", suffixes=("", "_npc"))
    novel = novel.merge(scaff, on="comp_id", how="left", suffixes=("", "_sc"))
    novel = novel.merge(sa, on="comp_id", how="left", suffixes=("", "_sa"))
    # ADMET column names will collide with whatever is in `corpus_cols`; we
    # take the bundle versions, which is why the merge is left-join.
    # admet may carry many columns, only the ones in corpus_cols survive
    # the final select below.

    # Fill the Phase B-derived columns into the production schema.
    if "np_class" in npc.columns and "np_class" in corpus_cols:
        # already merged via the npc merge above
        pass
    if "np_pathway" in npc.columns and "np_pathway" in corpus_cols:
        pass
    if "np_superclass" in npc.columns and "np_superclass" in corpus_cols:
        pass
    # ADMET predictions: the production schema has them in the admet table,
    # not the compounds table, so we don't include them in /tmp/tipdb_novel_rows.csv
    # The admet table is loaded separately via load_data.py post-deploy.

    # Trust score, novelty, sa_score, classyfire_superclass: blank for now,
    # to be filled by L3S-side scripts after the database INSERT.
    for col in ["trust_score", "novelty_morgan", "novelty_maccs",
                "novelty_type", "classyfire_superclass", "inferred_class",
                "inferred_class_source", "inferred_confidence",
                "reference_doi", "trad_medicine", "name_norm"]:
        if col not in novel.columns and col in corpus_cols:
            novel[col] = ""

    # The scaffold field is in the production scaffolds table, not the
    # compounds table either; same handling as ADMET.

    # SA score: production schema has sa_score on compounds; populate from
    # the bundle.
    if "sa_score" in corpus_cols:
        # the merge with `sa` above added a sa_score column; use it
        if "sa_score_sa" in novel.columns:
            novel["sa_score"] = novel["sa_score_sa"]
        novel["sa_score"] = novel["sa_score"].fillna("")

    # Final select: only the columns the production schema has, in order.
    missing = [c for c in corpus_cols if c not in novel.columns]
    if missing:
        print(f"  filling blank for production columns: {missing}")
        for c in missing:
            novel[c] = ""
    novel_out = novel[corpus_cols]
    print(f"[write] {OUT_NOVEL}")
    novel_out.to_csv(OUT_NOVEL, index=False)
    print(f"  {len(novel_out):,} rows x {len(novel_out.columns)} cols, "
          f"{os.path.getsize(OUT_NOVEL)/1024:.1f} KB")

    print(f"[write] {OUT_OVERLAP}")
    with open(OUT_OVERLAP, "w") as f:
        for ik in overlap_inchikeys:
            f.write(f"{ik}\n")
    print(f"  {len(overlap_inchikeys):,} InChIKeys, "
          f"{os.path.getsize(OUT_OVERLAP)/1024:.1f} KB")

    print("\n[summary]")
    print(f"  novel rows for INSERT:    {len(novel_out):,}")
    print(f"  overlap InChIKeys for UPDATE: {len(overlap_inchikeys):,}")
    print("  next step: sudo -u postgres psql -d theobroma -f tipdb_production_merge.sql")

if __name__ == "__main__":
    main()
