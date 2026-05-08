"""Convert TIPdb-3D archive (per-compound .mol files) to standardized CSV.

Reads tipdb_chemical/3d/Cxxxxxxxx.mol, derives 2D SMILES, InChI, InChIKey via
RDKit, computes physicochemical descriptors, and writes data/converted/tipdb.csv
following the COLS schema consumed by merge_with_coconut.py.

Source-side identifier (Cxxxxxxxx) is stored in comp_id; merge_with_coconut
will reassign THEO_ ids on integration. Name and source_organism are left
empty; the per-compound metadata (plant species, ethnobotanical context) is
not present in tipdb_chemical.zip and would need to come from a separate
TIPdb metadata table that Tung has not yet supplied.

License placeholder is CC BY-NC 4.0 pending Tung's written confirmation.
"""
import os, warnings
warnings.filterwarnings("ignore")
import pandas as pd
from tqdm import tqdm
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, inchi
RDLogger.logger().setLevel(RDLogger.ERROR)

DATADIR = "data"
TIPDB_DIR = os.path.join(DATADIR, "tipdb_extract", "3d")
OUTDIR = os.path.join(DATADIR, "converted")
os.makedirs(OUTDIR, exist_ok=True)

# Schema matches convert_local.py COLS (no all_sources; that is materialized at merge time).
COLS = ["comp_id","name","smiles","inchi","inchikey","source_db","kingdom",
        "region","source_organism","mw","logp","tpsa","hba","hbd",
        "n_rings","rotatable_bonds","license_tier"]

SOURCE_DB = "TIPdb"
KINGDOM = "plant"
REGION = "East Asia"
LICENSE = "CC BY-NC 4.0"  # placeholder pending Tung's written confirmation

def derive(mol):
    """Compute SMILES, InChI, InChIKey, and descriptors from a parsed RDKit Mol.
    Returns a dict; missing values stay None so the merge step can spot them."""
    out = {"smiles": None, "inchi": None, "inchikey": None,
           "mw": None, "logp": None, "tpsa": None, "hba": None, "hbd": None,
           "n_rings": None, "rotatable_bonds": None}
    if mol is None:
        return out
    try:
        out["smiles"] = Chem.MolToSmiles(mol)
    except Exception:
        return out
    try:
        ich = inchi.MolToInchi(mol)
        out["inchi"] = ich
        out["inchikey"] = inchi.InchiToInchiKey(ich) if ich else None
    except Exception:
        pass
    try:
        out["mw"] = round(Descriptors.ExactMolWt(mol), 2)
        out["logp"] = round(Descriptors.MolLogP(mol), 2)
        out["tpsa"] = round(Descriptors.TPSA(mol), 1)
        out["hba"] = rdMolDescriptors.CalcNumHBA(mol)
        out["hbd"] = rdMolDescriptors.CalcNumHBD(mol)
        out["n_rings"] = rdMolDescriptors.CalcNumRings(mol)
        out["rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
    except Exception:
        pass
    return out

def main():
    if not os.path.isdir(TIPDB_DIR):
        raise SystemExit(f"TIPdb directory not found: {TIPDB_DIR}")
    files = sorted(f for f in os.listdir(TIPDB_DIR) if f.endswith(".mol"))
    print(f"[TIPdb] {len(files):,} mol files in {TIPDB_DIR}")
    records, parsed_ok, no_smiles, no_inchikey = [], 0, 0, 0
    for fname in tqdm(files, desc="  TIPdb", leave=False):
        fpath = os.path.join(TIPDB_DIR, fname)
        comp_id_src = os.path.splitext(fname)[0]  # e.g. C00000001
        # removeHs=True normalizes hydrogen handling so SMILES are comparable
        # across the corpus; sanitize=True triggers RDKit valence checks which
        # is what every other source went through.
        mol = Chem.MolFromMolFile(fpath, removeHs=True, sanitize=True)
        d = derive(mol)
        if d["smiles"] is None:
            no_smiles += 1
            continue
        if d["inchikey"] is None:
            no_inchikey += 1
            # We still keep the row; merge_with_coconut routes no-InChIKey rows
            # into the "kept as-is" residual bucket rather than dropping them.
        parsed_ok += 1
        records.append([comp_id_src, "", d["smiles"], d["inchi"] or "",
                        d["inchikey"] or "", SOURCE_DB, KINGDOM, REGION, "",
                        d["mw"], d["logp"], d["tpsa"], d["hba"], d["hbd"],
                        d["n_rings"], d["rotatable_bonds"], LICENSE])
    df = pd.DataFrame(records, columns=COLS)
    df = df[df["smiles"].astype(str).str.len() > 2]
    out_path = os.path.join(OUTDIR, "tipdb.csv")
    df.to_csv(out_path, index=False)
    print(f"  parsed: {parsed_ok:,} / {len(files):,}")
    print(f"  no SMILES (skipped): {no_smiles:,}")
    print(f"  no InChIKey (kept as-is, no dedup): {no_inchikey:,}")
    print(f"  saved {len(df):,} rows to {out_path}")
    print(f"  size: {os.path.getsize(out_path)/1024/1024:.2f} MB")

if __name__ == "__main__":
    main()
