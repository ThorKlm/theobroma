"""Convert FDP scraper output to THEOBROMA format and split by database."""
import pandas as pd
import os

INPUT = "data/converted/fdp_all.csv"
OUTDIR = "data/converted"

KINGDOM_MAP = {"LMDB": "fungi", "SMDB": "plant", "TMDB": "fungi", "CMDB": "plant"}
REGION_MAP = {"LMDB": "global", "SMDB": "South Asia", "TMDB": "global", "CMDB": "global"}
SOURCE_MAP = {"LMDB": "LMDB_Lichen", "SMDB": "SMDB_Spice", "TMDB": "TMDB_Trichoderma", "CMDB": "CMDB_Cereals"}

COLS = ["comp_id","name","smiles","inchi","inchikey","source_db","kingdom",
        "region","source_organism","mw","logp","tpsa","hba","hbd",
        "n_rings","rotatable_bonds","license_tier"]

print("Loading FDP data...")
df = pd.read_csv(INPUT)
print(f"  Input: {len(df):,} rows")

out = pd.DataFrame({
    "comp_id": [f"{SOURCE_MAP.get(r['database'],r['database'])}_{i:06d}" for i, r in df.iterrows()],
    "name": df["common_name"].fillna(""),
    "smiles": df["canonical_smiles"].fillna(""),
    "inchi": "",
    "inchikey": "",
    "source_db": df["database"].map(SOURCE_MAP).fillna(df["database"]),
    "kingdom": df["database"].map(KINGDOM_MAP).fillna("plant"),
    "region": df["database"].map(REGION_MAP).fillna("global"),
    "source_organism": df["species"].fillna(""),
    "mw": df["molecular_weight"],
    "logp": "",
    "tpsa": "",
    "hba": "",
    "hbd": "",
    "n_rings": "",
    "rotatable_bonds": "",
    "license_tier": "CC BY 4.0",
})

out = out[out["smiles"].astype(str).str.len() > 2]
print(f"  Valid SMILES: {len(out):,}")

for db in sorted(out["source_db"].unique()):
    subset = out[out["source_db"] == db]
    fname = db.lower().replace(" ", "_") + ".csv"
    subset[COLS].to_csv(os.path.join(OUTDIR, fname), index=False)
    print(f"  {fname}: {len(subset):,}")

os.remove(INPUT)
print("Done. Removed fdp_all.csv to avoid double-counting.")