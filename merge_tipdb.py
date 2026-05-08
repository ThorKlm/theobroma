"""Merge converted TIPdb CSV into the post-cleanup theobroma corpus.

Append-only merge. Existing THEO_ ids preserved unchanged; TIPdb-novel
compounds get fresh THEO_ ids appended at the tail. Overlap compounds
(TIPdb InChIKey already present in corpus) get TIPdb appended to their
all_sources field but keep all other fields including license tier
unchanged. This implements the most-permissive-license-wins rule by
construction (we never overwrite an existing license).

Targets the post-cleanup CSV (data/theobroma_v33_clean.csv, 1,078,671
rows) as the default input. Override with --input. Writes the augmented
corpus to data/theobroma_v33_merged.csv and the SMILES list of the novel
TIPdb compounds to data/tipdb_novel_smiles.csv (input for Phase B
embedding generation on vast.ai).

Run from the repo root after convert_tipdb.py:
    python merge_tipdb.py
"""
import argparse, os, shutil, warnings
warnings.filterwarnings("ignore")
import pandas as pd

DATADIR = "data"
DEFAULT_INPUT = os.path.join(DATADIR, "theobroma_v33_clean.csv")
DEFAULT_OUTPUT = os.path.join(DATADIR, "theobroma_v33_merged.csv")
TIPDB_CSV = os.path.join(DATADIR, "converted", "tipdb.csv")
NOVEL_SMILES_OUT = os.path.join(DATADIR, "tipdb_novel_smiles.csv")

TIPDB_TAG = "TIPdb"
ID_PREFIX = "THEO_"
ID_WIDTH = 7

def append_source(existing, tag):
    """Append tag to a pipe-separated all_sources string with order
    preserved and no duplication. NaN/None/empty becomes tag alone."""
    if pd.isna(existing) or not str(existing).strip():
        return tag
    parts = [p.strip() for p in str(existing).split("|") if p.strip()]
    if tag not in parts:
        parts.append(tag)
    return "|".join(parts)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=DEFAULT_INPUT,
                    help="canonical corpus CSV; default %(default)s")
    ap.add_argument("--output", default=DEFAULT_OUTPUT,
                    help="merged corpus CSV; default %(default)s")
    args = ap.parse_args()
    if not os.path.exists(TIPDB_CSV):
        raise SystemExit(f"TIPdb CSV not found: {TIPDB_CSV}; run convert_tipdb.py first")
    if not os.path.exists(args.input):
        raise SystemExit(f"corpus CSV not found: {args.input}")

    print(f"[load] TIPdb {TIPDB_CSV}")
    tip = pd.read_csv(TIPDB_CSV, low_memory=False)
    before = len(tip)
    tip = tip[tip["inchikey"].notna() & (tip["inchikey"].astype(str).str.len() > 10)]
    tip = tip.sort_values(["inchikey", "comp_id"]).drop_duplicates(
        subset="inchikey", keep="first").reset_index(drop=True)
    print(f"  {before:,} -> {len(tip):,} after internal InChIKey dedup")

    print(f"[load] corpus {args.input}")
    corpus = pd.read_csv(args.input, low_memory=False)
    print(f"  {len(corpus):,} rows; {len(corpus.columns)} columns")

    corpus_iks = set(corpus["inchikey"].dropna().astype(str))
    overlap_mask = tip["inchikey"].astype(str).isin(corpus_iks)
    overlap = tip[overlap_mask].copy()
    novel = tip[~overlap_mask].copy()
    print(f"[split] overlap (enrichment): {len(overlap):,}")
    print(f"        novel (append):      {len(novel):,}")

    # Enrichment: vectorized in-place update of all_sources for rows whose
    # InChIKey is in the overlap set. Earlier we saw a "rows touched > overlap
    # count" warning when the corpus had un-deduped duplicates; the post-
    # cleanup corpus should report touched == overlap.
    print("[enrich] updating all_sources on overlap rows")
    overlap_iks = set(overlap["inchikey"].astype(str))
    target = corpus["inchikey"].astype(str).isin(overlap_iks)
    new_all_sources = corpus["all_sources"].copy()
    new_all_sources.loc[target] = corpus.loc[target, "all_sources"].apply(
        lambda v: append_source(v, TIPDB_TAG))
    updates = int((new_all_sources != corpus["all_sources"]).sum())
    corpus["all_sources"] = new_all_sources
    print(f"  rows touched: {updates:,} (expected {len(overlap):,})")
    if updates != len(overlap):
        print(f"  WARNING: enrichment update count mismatch; investigate "
              f"residual duplicates in the corpus")

    # Append novel rows. Comp_ids are sequential from len(corpus) onward;
    # this preserves the append-only invariant artifact-extension relies on.
    print("[append] writing novel TIPdb rows with fresh THEO_ ids")
    start_id = len(corpus)
    novel = novel.copy()
    novel["comp_id"] = [f"{ID_PREFIX}{i:0{ID_WIDTH}d}" for i in
                        range(start_id, start_id + len(novel))]
    if "all_sources" not in novel.columns:
        novel["all_sources"] = TIPDB_TAG
    # Fill any column the corpus has but the converter does not (np_class,
    # classyfire_superclass, np_pathway, inferred_class, trust_score,
    # reference_doi, trad_medicine, novelty_*, sa_score, name_norm) with
    # empty strings, to be filled by the per-artifact extension scripts in
    # Phase B.
    for col in corpus.columns:
        if col not in novel.columns:
            novel[col] = ""
    novel = novel[corpus.columns]

    augmented = pd.concat([corpus, novel], ignore_index=True)
    print(f"[write] {args.output}")
    augmented.to_csv(args.output, index=False)
    print(f"  {len(augmented):,} rows ({os.path.getsize(args.output)/1024/1024:.1f} MB)")

    # Phase B input: SMILES list of the novel compounds, with their newly
    # assigned THEO_ ids. This is what gets uploaded to vast.ai for ChemBERTa
    # embedding, ADMET-AI inference, NPClassifier, and XGBoost prediction.
    novel_out = novel[["comp_id", "inchikey", "smiles"]].copy()
    novel_out.to_csv(NOVEL_SMILES_OUT, index=False)
    print(f"[write] {NOVEL_SMILES_OUT}")
    print(f"  {len(novel_out):,} novel SMILES for Phase B "
          f"({os.path.getsize(NOVEL_SMILES_OUT)/1024:.1f} KB)")

    # Final summary block.
    print("\n[summary]")
    print(f"  v33 baseline:           {len(corpus):,}")
    print(f"  TIPdb unique compounds: {len(tip):,}")
    print(f"  enrichment-only:        {len(overlap):,}")
    print(f"  appended:               {len(novel):,}")
    print(f"  v33 augmented corpus:   {len(augmented):,}")
    print(f"  delta:                  +{len(novel):,} "
          f"({100*len(novel)/len(corpus):.4f}%)")
    n_tipdb_anywhere = augmented["all_sources"].astype(str).str.contains(
        TIPDB_TAG, regex=False, na=False).sum()
    print(f"  rows with TIPdb in all_sources: {n_tipdb_anywhere:,} "
          f"(expected {len(tip):,})")

if __name__ == "__main__":
    main()
