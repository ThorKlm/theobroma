"""Deduplicate all converted CSVs by InChIKey.
Priority: regional sources over COCONUT. Most permissive license wins."""
import pandas as pd
import os, warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
from rdkit import Chem, RDLogger
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
RDLogger.logger().setLevel(RDLogger.ERROR)

CONVERTED = os.path.join("data", "converted")
OUTFILE = os.path.join("data", "theobroma_merged.csv")

LICENSE_RANK = {"CC0": 0, "CC BY 4.0": 1, "CC BY-NC 4.0": 2, "CC BY-NC-ND 4.0": 3}
SOURCE_PRIORITY = {
    "IMPPAT": 100, "IMPPAT_supp": 100, "Phyto4Health": 95,
    "VIETHERB": 95, "ANPDB": 90, "SANCDB": 90, "AfroDb": 90,
    "CSIRO": 90, "phytochemdb": 85, "MicotoXilico": 85,
    "MeFSAT": 85, "EMNPD": 80, "MIBiG": 80,
    "TM-MC": 75, "HERB": 75, "CMNPD": 70,
    "NPASS": 50, "COCONUT": 10,
}

print("Loading all converted CSVs...")
frames = []
for f in sorted(os.listdir(CONVERTED)):
    if not f.endswith(".csv"): continue
    path = os.path.join(CONVERTED, f)
    df = pd.read_csv(path, low_memory=False)
    if len(df) == 0: continue
    frames.append(df)
    print(f"  {f}: {len(df):,}")

all_df = pd.concat(frames, ignore_index=True)
print(f"\nTotal before dedup: {len(all_df):,}")

# Compute missing InChIKeys from SMILES
missing = all_df["inchikey"].isna() | (all_df["inchikey"].astype(str) == "") | (all_df["inchikey"].astype(str) == "nan") | (all_df["inchikey"].astype(str).str.len() < 10)
print(f"Missing InChIKeys: {missing.sum():,} -- computing from SMILES...")
for idx in tqdm(all_df[missing].index, desc="Computing InChIKeys"):
    smi = str(all_df.at[idx, "smiles"])
    if len(smi) < 3: continue
    mol = Chem.MolFromSmiles(smi)
    if mol:
        try:
            all_df.at[idx, "inchikey"] = InchiToInchiKey(MolToInchi(mol))
        except:
            pass

# Recount
valid_ik = all_df["inchikey"].notna() & (all_df["inchikey"].astype(str).str.len() > 10)
print(f"\nAfter computation:")
print(f"  Valid InChIKeys: {valid_ik.sum():,}")
print(f"  Still missing: {(~valid_ik).sum():,}")
print(f"  Unique InChIKeys: {all_df.loc[valid_ik, 'inchikey'].nunique():,}")

# Assign priority and license rank
all_df["_priority"] = all_df["source_db"].map(SOURCE_PRIORITY).fillna(50)
all_df["_lic_rank"] = all_df["license_tier"].map(LICENSE_RANK).fillna(2)

# Split: compounds with valid InChIKey vs without
has_ik = all_df[valid_ik].copy()
no_ik = all_df[~valid_ik].copy()

print(f"\nWith InChIKey: {len(has_ik):,}")
print(f"Without InChIKey: {len(no_ik):,}")

# Sort: highest priority first, then most permissive license
print("Deduplicating by InChIKey...")
has_ik = has_ik.sort_values(["_priority", "_lic_rank"], ascending=[False, True])

# For license: aggregate most permissive per InChIKey
best_license = has_ik.groupby("inchikey")["_lic_rank"].min().reset_index()
best_license.columns = ["inchikey", "_best_lic_rank"]
lic_rank_inv = {v: k for k, v in LICENSE_RANK.items()}
best_license["best_license"] = best_license["_best_lic_rank"].map(lic_rank_inv)

# Keep first (highest priority) per InChIKey
deduped = has_ik.drop_duplicates(subset="inchikey", keep="first").copy()
deduped = deduped.merge(best_license[["inchikey", "best_license"]], on="inchikey", how="left")
deduped["license_tier"] = deduped["best_license"]

# Collect all source_dbs per InChIKey for provenance
print("Collecting source provenance...")
sources_per_ik = has_ik.groupby("inchikey")["source_db"].apply(
    lambda x: "|".join(sorted(set(x)))).reset_index()
sources_per_ik.columns = ["inchikey", "all_sources"]
deduped = deduped.merge(sources_per_ik, on="inchikey", how="left")

# Collect most specific region (non-"global" preferred)
def best_region(group):
    non_global = [r for r in group if r != "global" and pd.notna(r) and r != "" and str(r) != "nan"]
    return non_global[0] if non_global else "global"

regions_per_ik = has_ik.groupby("inchikey")["region"].apply(best_region).reset_index()
regions_per_ik.columns = ["inchikey", "best_region"]
deduped = deduped.merge(regions_per_ik, on="inchikey", how="left")
deduped["region"] = deduped["best_region"]

# Collect best organism info (non-empty preferred)
def best_organism(group):
    non_empty = [o for o in group if pd.notna(o) and str(o) != "" and str(o) != "nan"]
    return non_empty[0] if non_empty else ""

org_per_ik = has_ik.groupby("inchikey")["source_organism"].apply(best_organism).reset_index()
org_per_ik.columns = ["inchikey", "best_organism"]
deduped = deduped.merge(org_per_ik, on="inchikey", how="left")
deduped["source_organism"] = deduped["best_organism"]

# Clean up temp columns
drop_cols = ["_priority", "_lic_rank", "best_license", "_best_lic_rank", "best_region", "best_organism"]
deduped = deduped.drop(columns=[c for c in drop_cols if c in deduped.columns])

# For no-InChIKey compounds, add all_sources column
no_ik["all_sources"] = no_ik["source_db"]
no_ik = no_ik.drop(columns=[c for c in ["_priority", "_lic_rank"] if c in no_ik.columns])

# Combine
final = pd.concat([deduped, no_ik], ignore_index=True)

# Reassign comp_ids
final["comp_id"] = [f"THEO_{i:07d}" for i in range(len(final))]

# Final column order
out_cols = ["comp_id", "name", "smiles", "inchi", "inchikey", "source_db", "all_sources",
            "kingdom", "region", "source_organism", "mw", "logp", "tpsa", "hba", "hbd",
            "n_rings", "rotatable_bonds", "license_tier"]
for c in out_cols:
    if c not in final.columns:
        final[c] = ""

final[out_cols].to_csv(OUTFILE, index=False)

print(f"\n=== Deduplication complete ===")
print(f"Before: {len(all_df):,}")
print(f"After:  {len(final):,}")
print(f"Removed: {len(all_df) - len(final):,} duplicates")
print(f"\nSaved to {OUTFILE} ({os.path.getsize(OUTFILE) / 1024 / 1024:.1f} MB)")
print(f"\nBy kingdom:")
for k, cnt in final["kingdom"].value_counts().items():
    print(f"  {k}: {cnt:,}")
print(f"\nBy region:")
for r, cnt in final["region"].value_counts().head(15).items():
    print(f"  {r}: {cnt:,}")
print(f"\nBy license:")
for l, cnt in final["license_tier"].value_counts().items():
    print(f"  {l}: {cnt:,}")
print(f"\nBy source (top 20):")
for s, cnt in final["source_db"].value_counts().head(20).items():
    print(f"  {s}: {cnt:,}")
print(f"\nMulti-source compounds: {(final['all_sources'].str.contains('|', regex=False, na=False)).sum():,}")