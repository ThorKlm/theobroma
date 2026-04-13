"""Pre-pipeline check: verify all source files and flag non-compound files."""
import pandas as pd
import os

CONVERTED = "data/converted"
ENRICHMENT_FILES = {"npass_enrichment.csv", "npass_geo_enrichment.csv"}
EXPECTED_COLS = {"comp_id", "name", "smiles", "source_db"}

print("=== Pre-pipeline file check ===\n")
for f in sorted(os.listdir(CONVERTED)):
    if not f.endswith(".csv"):
        continue
    path = os.path.join(CONVERTED, f)
    try:
        df = pd.read_csv(path, nrows=5, low_memory=False)
        n = len(pd.read_csv(path, usecols=[0], low_memory=False))
        cols = set(df.columns)
        is_source = EXPECTED_COLS.issubset(cols)
        is_enrichment = f in ENRICHMENT_FILES
        has_smiles = "smiles" in cols or "canonical_smiles" in cols
        region = df["region"].iloc[0] if "region" in cols else "N/A"
        tag = "ENRICHMENT" if is_enrichment else ("OK" if is_source else "WRONG COLS")
        print(f"  {f:35s} {n:>8,}  region={str(region):15s}  [{tag}]")
    except Exception as e:
        print(f"  {f:35s}  ERROR: {e}")