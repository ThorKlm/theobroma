"""Apply NPASS organism enrichment to theobroma_final.csv.
Joins 204k organism names by InChIKey into empty source_organism fields."""
import pandas as pd
import os, warnings
warnings.filterwarnings("ignore")

ENRICHMENT = "data/converted/npass_enrichment.csv"
FINAL = "data/theobroma_final.csv"

if not os.path.exists(FINAL):
    print(f"{FINAL} not found -- run dedup_and_merge.py + merge_with_coconut.py first.")
    print("This script runs AFTER the main pipeline produces theobroma_final.csv.")
    exit()

print("Loading enrichment data...")
enrich = pd.read_csv(ENRICHMENT, low_memory=False)
print(f"  Enrichment records: {len(enrich):,}")
print(f"  With organism: {(enrich['source_organism'] != '').sum():,}")
print(f"  With geo location: {(enrich['geo_location'] != '').sum():,}")

# Build InChIKey -> organism mapping
org_map = {}
geo_map = {}
for _, row in enrich.iterrows():
    ik = str(row.get("inchikey", "")).strip()
    if not ik or ik == "nan" or len(ik) < 10:
        continue
    org = str(row.get("source_organism", "")).strip()
    geo = str(row.get("geo_location", "")).strip()
    if org and org != "nan":
        org_map[ik] = org
    if geo and geo != "nan":
        geo_map[ik] = geo

print(f"  Unique InChIKeys with organism: {len(org_map):,}")
print(f"  Unique InChIKeys with location: {len(geo_map):,}")

print("\nLoading theobroma_final.csv...")
df = pd.read_csv(FINAL, low_memory=False)
print(f"  Total compounds: {len(df):,}")

# Count empty organisms before
empty_org_before = (df["source_organism"].fillna("").astype(str).str.strip() == "").sum()
print(f"  Empty source_organism before: {empty_org_before:,}")

# Apply organism enrichment only to empty fields
filled_org = 0
filled_geo = 0
for idx in range(len(df)):
    ik = str(df.at[idx, "inchikey"]).strip() if pd.notna(df.at[idx, "inchikey"]) else ""
    if not ik or len(ik) < 10:
        continue
    # Organism
    current_org = str(df.at[idx, "source_organism"]).strip() if pd.notna(df.at[idx, "source_organism"]) else ""
    if not current_org or current_org == "nan":
        new_org = org_map.get(ik, "")
        if new_org:
            df.at[idx, "source_organism"] = new_org
            filled_org += 1
    # Geographic location (only if currently unresolved)
    current_reg = str(df.at[idx, "region"]).strip() if pd.notna(df.at[idx, "region"]) else ""
    if not current_reg or current_reg in ("", "nan", "global", "unresolved"):
        new_geo = geo_map.get(ik, "")
        if new_geo:
            df.at[idx, "region"] = new_geo
            filled_geo += 1

empty_org_after = (df["source_organism"].fillna("").astype(str).str.strip() == "").sum()
print(f"\n  Filled organism: {filled_org:,}")
print(f"  Filled location: {filled_geo:,}")
print(f"  Empty source_organism after: {empty_org_after:,}")

df.to_csv(FINAL, index=False)
print(f"  Updated {FINAL}")