"""
Load compound data from various source CSVs/SDFs into PostgreSQL.
Run once after downloading data from HuggingFace.

Usage: python scripts/load_data.py [--db DATABASE_URL] [--data-dir /path/to/data]
"""
import psycopg2, pandas as pd, argparse, os, sys, glob, json, gzip
from pathlib import Path

SCHEMA = """
DROP TABLE IF EXISTS compounds CASCADE;
CREATE TABLE compounds (
    id SERIAL PRIMARY KEY,
    comp_id TEXT UNIQUE NOT NULL,
    name TEXT,
    smiles TEXT NOT NULL,
    inchi TEXT,
    inchikey TEXT,
    source_db TEXT NOT NULL,
    kingdom TEXT,
    region TEXT,
    source_organism TEXT,
    mw REAL,
    logp REAL,
    tpsa REAL,
    hba INTEGER,
    hbd INTEGER,
    n_rings INTEGER,
    rotatable_bonds INTEGER,
    np_likeness REAL,
    classyfire_superclass TEXT,
    np_class TEXT,
    license_tier TEXT DEFAULT 'CC BY 4.0'
);
CREATE INDEX idx_compounds_inchikey ON compounds(inchikey);
CREATE INDEX idx_compounds_smiles ON compounds(smiles);
CREATE INDEX idx_compounds_name ON compounds USING gin(to_tsvector('english', COALESCE(name, '')));
CREATE INDEX idx_compounds_source ON compounds(source_db);
CREATE INDEX idx_compounds_kingdom ON compounds(kingdom);
CREATE INDEX idx_compounds_region ON compounds(region);
CREATE INDEX idx_compounds_organism ON compounds(source_organism);
"""

# Source DB -> (kingdom, region, license_tier)
SOURCE_META = {
    "COCONUT":       ("multi",    "global",          "CC BY 4.0"),
    "CMAUP":         ("plant",    "East Asia",       "CC BY 4.0"),
    "TCMBank":       ("plant",    "East Asia",       "CC BY 4.0"),
    "HERB":          ("plant",    "East Asia",       "CC BY-NC 4.0"),
    "TM-MC":         ("plant",    "East Asia",       "CC0"),
    "NPASS":         ("multi",    "global",          "CC BY 4.0"),
    "FooDB":         ("food",     "global",          "CC BY-NC 4.0"),
    "CMNPD":         ("marine",   "global",          "CC BY 4.0"),
    "IMPPAT":        ("plant",    "South Asia",      "CC BY-NC 4.0"),
    "Phyto4Health":  ("plant",    "Russia/CIS",      "CC BY-NC 4.0"),
    "StreptomeDB":   ("bacteria", "global",          "CC BY 4.0"),
    "NPAtlas":       ("multi",    "global",          "CC BY 4.0"),
    "MeFSAT":        ("fungi",    "global",          "CC BY 4.0"),
    "MicotoXilico":  ("fungi",    "global",          "CC BY 4.0"),
    "MIBiG":         ("multi",    "global",          "CC BY 4.0"),
    "phytochemdb":   ("plant",    "South Asia",      "CC BY 4.0"),
    "CSIRO":         ("plant",    "Australia",       "CC BY 4.0"),
    "VIETHERB":      ("plant",    "Southeast Asia",  "CC BY-NC 4.0"),
    "ANPDB":         ("multi",    "Africa",          "CC BY 4.0"),
    "SANCDB":        ("plant",    "Africa",          "CC BY 4.0"),
    "AfroDb":        ("multi",    "Africa",          "CC BY 4.0"),
    "ConMedNP":      ("plant",    "Central Africa",  "CC BY 4.0"),
    "CamMedNP":      ("plant",    "West Africa",     "CC BY 4.0"),
    "LANaPDB":       ("plant",    "Latin America",   "CC BY 4.0"),
    "CyanoMetDB":    ("bacteria", "global",          "CC BY 4.0"),
    "SWMD":          ("marine",   "global",          "CC BY 4.0"),
    "YMDB":          ("fungi",    "global",          "CC BY 4.0"),
    "EMNPD":         ("multi",    "global",          "CC BY 4.0"),
}

def insert_batch(cur, records):
    """Batch insert with ON CONFLICT skip."""
    for rec in records:
        try:
            cur.execute("""
                INSERT INTO compounds (comp_id, name, smiles, inchi, inchikey,
                    source_db, kingdom, region, source_organism,
                    mw, logp, tpsa, hba, hbd, n_rings, rotatable_bonds,
                    np_likeness, classyfire_superclass, np_class, license_tier)
                VALUES (%(comp_id)s, %(name)s, %(smiles)s, %(inchi)s, %(inchikey)s,
                    %(source_db)s, %(kingdom)s, %(region)s, %(source_organism)s,
                    %(mw)s, %(logp)s, %(tpsa)s, %(hba)s, %(hbd)s, %(n_rings)s,
                    %(rotatable_bonds)s, %(np_likeness)s, %(classyfire_superclass)s,
                    %(np_class)s, %(license_tier)s)
                ON CONFLICT (comp_id) DO NOTHING
            """, rec)
        except Exception as e:
            print(f"  Skip {rec.get('comp_id','?')}: {e}")
            cur.connection.rollback()

def safe_float(val):
    try: return float(val)
    except: return None

def safe_int(val):
    try: return int(float(val))
    except: return None

def load_csv(filepath, source_db, conn):
    """Generic CSV loader. Expects at minimum a 'smiles' column."""
    meta = SOURCE_META.get(source_db, ("unknown", "unknown", "CC BY 4.0"))
    df = pd.read_csv(filepath, low_memory=False)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
    records = []
    for i, row in df.iterrows():
        smiles = str(row.get("smiles", "")).strip()
        if not smiles or smiles == "nan" or len(smiles) < 2:
            continue
        records.append({
            "comp_id": str(row.get("comp_id", f"{source_db}_{i:06d}")),
            "name": str(row.get("name", row.get("compound_name", ""))) if pd.notna(row.get("name", row.get("compound_name", ""))) else "",
            "smiles": smiles,
            "inchi": str(row.get("inchi", "")) if pd.notna(row.get("inchi", "")) else "",
            "inchikey": str(row.get("inchikey", "")) if pd.notna(row.get("inchikey", "")) else "",
            "source_db": source_db,
            "kingdom": str(row.get("kingdom", meta[0])),
            "region": str(row.get("region", meta[1])),
            "source_organism": str(row.get("source_organism", row.get("plant_sources", ""))) if pd.notna(row.get("source_organism", row.get("plant_sources", ""))) else "",
            "mw": safe_float(row.get("mw", row.get("molecular_weight", None))),
            "logp": safe_float(row.get("logp", row.get("log_p", None))),
            "tpsa": safe_float(row.get("tpsa", None)),
            "hba": safe_int(row.get("hba", row.get("hydrogen_bond_acceptors", None))),
            "hbd": safe_int(row.get("hbd", row.get("hydrogen_bond_donors", None))),
            "n_rings": safe_int(row.get("n_rings", row.get("number_of_rings", None))),
            "rotatable_bonds": safe_int(row.get("rotatable_bonds", None)),
            "np_likeness": safe_float(row.get("np_likeness", None)),
            "classyfire_superclass": str(row.get("classyfire_superclass", "")) if pd.notna(row.get("classyfire_superclass", "")) else "",
            "np_class": str(row.get("np_class", "")) if pd.notna(row.get("np_class", "")) else "",
            "license_tier": meta[2],
        })
    with conn.cursor() as cur:
        insert_batch(cur, records)
    conn.commit()
    print(f"  [{source_db}] Loaded {len(records)} compounds from {os.path.basename(filepath)}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=os.environ.get("DATABASE_URL",
        "postgresql://theobroma:theobroma@localhost:5432/theobroma"))
    ap.add_argument("--data-dir", default="/home/thorben.klamt/theobroma/data")
    ap.add_argument("--skip-schema", action="store_true")
    args = ap.parse_args()
    conn = psycopg2.connect(args.db)
    if not args.skip_schema:
        with conn.cursor() as cur:
            cur.execute(SCHEMA)
        conn.commit()
        print("Schema created.")
    d = args.data_dir
    # Map filenames to (source_db). Add entries as CSVs become available.
    file_map = {
        "imppat2_compounds.csv": "IMPPAT",
        "phyto4health_compounds.csv": "Phyto4Health",
        "mycotoxins_tmap_final.csv": "MicotoXilico",
        # Add more as they are converted to CSV:
        # "coconut_compounds.csv": "COCONUT",
        # "herb_compounds.csv": "HERB",
        # "npass_compounds.csv": "NPASS",
        # etc.
    }
    for fname, source_db in file_map.items():
        path = os.path.join(d, fname)
        if os.path.exists(path):
            load_csv(path, source_db, conn)
        else:
            print(f"  [{source_db}] File not found: {path}")
    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM compounds")
        total = cur.fetchone()[0]
    print(f"\nTotal compounds in database: {total:,}")
    conn.close()

if __name__ == "__main__":
    main()
