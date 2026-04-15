"""Load THEOBROMA SDF files into PostgreSQL using RDKit."""
import sys, os, gzip, psycopg2
sys.path.append("/usr/lib/python3/dist-packages")
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

DB_URI = os.environ.get("DATABASE_URL",
    "postgresql://theobroma:theobroma@localhost:5432/theobroma")

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
CREATE INDEX idx_compounds_source ON compounds(source_db);
CREATE INDEX idx_compounds_kingdom ON compounds(kingdom);
CREATE INDEX idx_compounds_region ON compounds(region);
CREATE INDEX idx_compounds_organism ON compounds(source_organism);
"""

KINGDOM_MAP = {
    "bacteria": "bacteria",
    "fungi": "fungi",
    "marine_algae": "marine",
    "plant_supplementary": "plant",
    "plant": "plant",
    "food": "food",
}

REGION_MAP = {
    "CyanoMetDB_V03_2024": "global",
    "npatlas_bacteria": "global",
    "streptomedb": "global",
    "npatlas_fungi": "global",
    "ymdb": "global",
    "swmd_mol": "global",
    "anpdb": "Africa",
    "cmaup": "East Asia",
    "lanapdb_v2": "Latin America",
    "tcmbank_ingredient_all": "East Asia",
    "COCONUT": "global",
}

LICENSE_MAP = {
    "CyanoMetDB_V03_2024": "CC BY 4.0",
    "npatlas_bacteria": "CC BY 4.0",
    "streptomedb": "CC BY 4.0",
    "npatlas_fungi": "CC BY 4.0",
    "ymdb": "CC BY 4.0",
    "swmd_mol": "CC BY 4.0",
    "anpdb": "CC BY 4.0",
    "cmaup": "CC BY 4.0",
    "lanapdb_v2": "CC BY 4.0",
    "tcmbank_ingredient_all": "CC BY 4.0",
    "COCONUT": "CC BY 4.0",
}

def compute_props(mol):
    try:
        return {
            "mw": round(Descriptors.ExactMolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "tpsa": round(Descriptors.TPSA(mol), 1),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "n_rings": rdMolDescriptors.CalcNumRings(mol),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        }
    except:
        return {"mw":None,"logp":None,"tpsa":None,"hba":None,"hbd":None,"n_rings":None,"rotatable_bonds":None}

def load_sdf(filepath, conn, prefix, is_coconut=False):
    opener = gzip.open if filepath.endswith(".gz") else open
    supplier = Chem.ForwardSDMolSupplier(opener(filepath, "rb"), removeHs=True, sanitize=True)
    batch = []
    loaded = skipped = 0
    cur = conn.cursor()
    for i, mol in enumerate(supplier):
        if mol is None:
            skipped += 1
            continue
        smiles = Chem.MolToSmiles(mol)
        if not smiles:
            skipped += 1
            continue
        src_db = mol.GetProp("source_database") if mol.HasProp("source_database") else prefix
        kingdom_raw = mol.GetProp("source_kingdom") if mol.HasProp("source_kingdom") else prefix
        kingdom = KINGDOM_MAP.get(kingdom_raw, kingdom_raw)
        inchikey = mol.GetProp("InChIKey") if mol.HasProp("InChIKey") else ""
        name = mol.GetProp("original_name") if mol.HasProp("original_name") else ""
        if name == "n.a.":
            name = ""
        comp_id = f"{prefix}_{i:07d}"
        region = REGION_MAP.get(src_db, "global")
        license_tier = LICENSE_MAP.get(src_db, "CC BY 4.0")
        props = compute_props(mol)
        batch.append((comp_id, name, smiles, "", inchikey, src_db, kingdom,
                       region, "", props["mw"], props["logp"], props["tpsa"],
                       props["hba"], props["hbd"], props["n_rings"],
                       props["rotatable_bonds"], None, "", "", license_tier))
        if len(batch) >= 5000:
            cur.executemany("""INSERT INTO compounds
                (comp_id,name,smiles,inchi,inchikey,source_db,kingdom,region,
                 source_organism,mw,logp,tpsa,hba,hbd,n_rings,rotatable_bonds,
                 np_likeness,classyfire_superclass,np_class,license_tier)
                VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
                ON CONFLICT (comp_id) DO NOTHING""", batch)
            conn.commit()
            loaded += len(batch)
            print(f"  {prefix}: {loaded} loaded, {skipped} skipped...", flush=True)
            batch = []
    if batch:
        cur.executemany("""INSERT INTO compounds
            (comp_id,name,smiles,inchi,inchikey,source_db,kingdom,region,
             source_organism,mw,logp,tpsa,hba,hbd,n_rings,rotatable_bonds,
             np_likeness,classyfire_superclass,np_class,license_tier)
            VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
            ON CONFLICT (comp_id) DO NOTHING""", batch)
        conn.commit()
        loaded += len(batch)
    cur.close()
    print(f"  [{prefix}] Done: {loaded} loaded, {skipped} skipped")

if __name__ == "__main__":
    conn = psycopg2.connect(DB_URI)
    with conn.cursor() as cur:
        cur.execute(SCHEMA)
    conn.commit()
    print("Schema created.")
    data = "data"
    files = [
        ("bacteria_unique.sdf", "bact"),
        ("fungi_unique.sdf", "fungi"),
        ("marine_algae_unique.sdf", "marine"),
        ("plant_supplementary_unique.sdf", "plant_sup"),
        ("coconut_prepped.sdf.gz", "coconut"),
    ]
    for fname, prefix in files:
        path = os.path.join(data, fname)
        if os.path.exists(path):
            print(f"Loading {fname}...")
            load_sdf(path, conn, prefix, is_coconut=(prefix=="coconut"))
        else:
            print(f"  Not found: {path}")
    with conn.cursor() as cur:
        cur.execute("SELECT COUNT(*) FROM compounds")
        total = cur.fetchone()[0]
        cur.execute("SELECT kingdom, COUNT(*) FROM compounds GROUP BY kingdom ORDER BY COUNT(*) DESC")
        for row in cur.fetchall():
            print(f"  {row[0]}: {row[1]:,}")
    print(f"\nTotal compounds: {total:,}")
    conn.close()
