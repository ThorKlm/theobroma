"""Convert ALL source files to standardized CSVs. Run from project root."""
import pandas as pd
import os, sys, json, gzip, tarfile, warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, inchi
RDLogger.logger().setLevel(RDLogger.ERROR)

DATADIR = "data"  # os.path.dirname(os.path.abspath(__file__))  # since script is in data/
OUTDIR = os.path.join(DATADIR, "converted")
os.makedirs(OUTDIR, exist_ok=True)

COLS = ["comp_id","name","smiles","inchi","inchikey","source_db","kingdom",
        "region","source_organism","mw","logp","tpsa","hba","hbd",
        "n_rings","rotatable_bonds","license_tier"]

def props(smi):
    mol = Chem.MolFromSmiles(str(smi))
    if mol is None: return {}, None
    try:
        return {"mw": round(Descriptors.ExactMolWt(mol),2),
                "logp": round(Descriptors.MolLogP(mol),2),
                "tpsa": round(Descriptors.TPSA(mol),1),
                "hba": rdMolDescriptors.CalcNumHBA(mol),
                "hbd": rdMolDescriptors.CalcNumHBD(mol),
                "n_rings": rdMolDescriptors.CalcNumRings(mol),
                "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol)}, mol
    except: return {}, mol

def get_ik(mol):
    if mol is None: return ""
    try: return inchi.MolToInchiKey(inchi.MolToInchi(mol))
    except: return ""

def save(records, name):
    df = pd.DataFrame(records, columns=COLS)
    df = df[df["smiles"].astype(str).str.len() > 2]
    path = os.path.join(OUTDIR, f"{name}.csv")
    df.to_csv(path, index=False)
    print(f"  -> {len(df):,} compounds saved to {name}.csv")
    return len(df)

def parse_sdf(filepath, source_db, kingdom, region, license_tier):
    opener = gzip.open if filepath.endswith(".gz") else open
    suppl = Chem.ForwardSDMolSupplier(opener(filepath, "rb"), removeHs=True, sanitize=True)
    recs = []
    for i, mol in enumerate(tqdm(suppl, desc=f"  {source_db}", leave=False)):
        if mol is None: continue
        smi = Chem.MolToSmiles(mol)
        if not smi: continue
        nm = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        p, _ = props(smi)
        recs.append([f"{source_db}_{i:06d}", nm, smi, "", get_ik(mol), source_db,
                      kingdom, region, "",
                      p.get("mw"), p.get("logp"), p.get("tpsa"),
                      p.get("hba"), p.get("hbd"), p.get("n_rings"),
                      p.get("rotatable_bonds"), license_tier])
    return recs

def parse_csv(filepath, source_db, kingdom, region, license_tier,
              name_field=None, organism_field=None, sep=","):
    df = pd.read_csv(filepath, sep=sep, low_memory=False, on_bad_lines="skip")
    df.columns = [c.strip().lower().replace(" ","_") for c in df.columns]
    sc = [c for c in df.columns if "smiles" in c.lower()]
    if not sc:
        print(f"  SKIP {source_db}: no SMILES column. Cols: {list(df.columns)[:10]}")
        return []
    sc = sc[0]
    nf = name_field if name_field and name_field in df.columns else \
         ([c for c in df.columns if "name" in c.lower()] + [""])[0]
    of = organism_field if organism_field and organism_field in df.columns else ""
    recs = []
    for idx, (i, row) in enumerate(tqdm(df.iterrows(), total=len(df), desc=f"  {source_db}", leave=False)):
        smi = str(row.get(sc, ""))
        if len(smi) < 3 or smi == "nan": continue
        p, mol = props(smi)
        ik = get_ik(mol) if mol else (str(row.get("inchikey","")) if "inchikey" in df.columns else "")
        nm = str(row.get(nf, "")) if nf else ""
        org = str(row.get(of, "")) if of else ""
        if nm == "nan": nm = ""
        if org == "nan": org = ""
        recs.append([f"{source_db}_{idx:06d}", nm, smi, "", ik, source_db,
                      kingdom, region, org,
                      p.get("mw"), p.get("logp"), p.get("tpsa"),
                      p.get("hba"), p.get("hbd"), p.get("n_rings"),
                      p.get("rotatable_bonds"), license_tier])
    return recs

def parse_excel(filepath, source_db, kingdom, region, license_tier, name_field=None):
    df = pd.read_excel(filepath, engine="openpyxl")
    df.columns = [c.strip().lower().replace(" ","_") for c in df.columns]
    sc = [c for c in df.columns if "smiles" in c.lower()]
    if not sc:
        print(f"  SKIP {source_db}: no SMILES column. Cols: {list(df.columns)[:10]}")
        return []
    sc = sc[0]
    nf = name_field if name_field and name_field in df.columns else \
         ([c for c in df.columns if "name" in c.lower()] + [""])[0]
    recs = []
    for i, row in tqdm(df.iterrows(), total=len(df), desc=f"  {source_db}", leave=False):
        smi = str(row.get(sc, ""))
        if len(smi) < 3 or smi == "nan": continue
        p, mol = props(smi)
        recs.append([f"{source_db}_{i:06d}", nm if (nm := str(row.get(nf,""))) != "nan" else "",
                      smi, "", get_ik(mol), source_db, kingdom, region, "",
                      p.get("mw"), p.get("logp"), p.get("tpsa"),
                      p.get("hba"), p.get("hbd"), p.get("n_rings"),
                      p.get("rotatable_bonds"), license_tier])
    return recs

total = 0
SOURCES = [
    # (type, filename, source_db, kingdom, region, license, kwargs)
    ("csv", "phyto4health_compounds.csv", "Phyto4Health", "plant", "Russia/CIS", "CC BY-NC 4.0",
     {"organism_field": "plant_sources"}),
    ("csv", "mycotoxins_tmap_final.csv", "MicotoXilico", "fungi", "global", "CC BY 4.0", {}),
    ("sdf", "SDF_Phytochemdb.sdf", "phytochemdb", "plant", "South Asia", "CC BY 4.0", {}),
    ("sdf", "SANCDB_all.sdf", "SANCDB", "plant", "Africa", "CC BY 4.0", {}),
    ("sdf", "pone.0078085.s005.sdf", "AfroDb", "multi", "Africa", "CC BY 4.0", {}),
    ("sdf", "CMNPD_1.0_2d.sdf.gz", "CMNPD", "marine", "global", "CC BY 4.0", {}),
    ("xlsx", "d0ra10322e2.xlsx", "MeFSAT", "fungi", "global", "CC BY 4.0", {}),
    ("xlsx", "ao3c00156_si_001.xlsx", "IMPPAT_supp", "plant", "South Asia", "CC BY-NC 4.0", {}),
    ("xlsx", "chemical_property.xlsx", "TM-MC", "plant", "East Asia", "CC0", {}),
    ("tsv", "HERB_ingredient_info_v2.txt", "HERB", "plant", "East Asia", "CC BY-NC 4.0",
     {"name_field": "ingredient_name"}),
    ("tsv", "NPASS3.0_naturalproducts_structure.txt", "NPASS", "multi", "global", "CC BY 4.0",
     {"organism_field": "organism"}),
    ("csv", "Natural Products General Information.csv", "EMNPD", "multi", "global", "CC BY 4.0", {}),
    ("csv", "Compounds.csv", "CSIRO", "plant", "Australia", "CC BY 4.0",
     {"organism_field": "species"}),
    ("csv", "ANPDB.csv", "ANPDB", "multi", "Africa", "CC BY 4.0", {}),
]

for entry in SOURCES:
    ftype, fname, sdb, king, reg, lic = entry[:6]
    kw = entry[6] if len(entry) > 6 else {}
    fpath = os.path.join(DATADIR, fname)
    if not os.path.exists(fpath):
        print(f"[SKIP] {fname} not found")
        continue
    print(f"[{sdb}] {fname}")
    if ftype == "csv":
        recs = parse_csv(fpath, sdb, king, reg, lic, **kw)
    elif ftype == "tsv":
        recs = parse_csv(fpath, sdb, king, reg, lic, sep="\t", **kw)
    elif ftype == "sdf":
        recs = parse_sdf(fpath, sdb, king, reg, lic)
    elif ftype == "xlsx":
        recs = parse_excel(fpath, sdb, king, reg, lic)
    total += save(recs, sdb.lower())

# ==== VIETHERB (OWL -- special) ====
f = os.path.join(DATADIR, "VHO.owl")
if os.path.exists(f):
    print("[VIETHERB] VHO.owl")
    import rdflib
    g = rdflib.Graph()
    g.parse(f, format="xml")
    recs, i = [], 0
    for s, p_rdf, o in tqdm(g, desc="  VIETHERB", leave=False):
        if "smiles" not in str(p_rdf).lower() and "SMILES" not in str(p_rdf): continue
        smi = str(o).strip()
        if len(smi) < 3: continue
        name = ""
        for _, p2, o2 in g.triples((s, None, None)):
            if "label" in str(p2).lower() or "name" in str(p2).lower():
                name = str(o2); break
        p, mol = props(smi)
        recs.append([f"VHERB_{i:06d}", name, smi, "", get_ik(mol),
                      "VIETHERB", "plant", "Southeast Asia", "",
                      p.get("mw"), p.get("logp"), p.get("tpsa"),
                      p.get("hba"), p.get("hbd"), p.get("n_rings"),
                      p.get("rotatable_bonds"), "CC BY-NC 4.0"])
        i += 1
    total += save(recs, "vietherb")

# ==== MIBiG (tar.gz of JSONs -- special) ====
f = os.path.join(DATADIR, "mibig_json_4.0.tar.gz")
if os.path.exists(f):
    print("[MIBiG] mibig_json_4.0.tar.gz")
    recs, i = [], 0
    with tarfile.open(f, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.name.endswith(".json")]
        for member in tqdm(members, desc="  MIBiG", leave=False):
            fobj = tar.extractfile(member)
            if fobj is None: continue
            try: data = json.loads(fobj.read())
            except: continue
            for comp in data.get("cluster", {}).get("compounds", []):
                smi = comp.get("chem_struct", "")
                if not smi or len(smi) < 3: continue
                p, mol = props(smi)
                recs.append([f"MIBIG_{i:06d}", comp.get("compound",""), smi, "",
                              get_ik(mol), "MIBiG", "multi", "global", "",
                              p.get("mw"), p.get("logp"), p.get("tpsa"),
                              p.get("hba"), p.get("hbd"), p.get("n_rings"),
                              p.get("rotatable_bonds"), "CC BY 4.0"])
                i += 1
    total += save(recs, "mibig")

print(f"\n=== Conversion complete: {total:,} total compounds ===")
for f in sorted(os.listdir(OUTDIR)):
    if not f.endswith(".csv"): continue
    path = os.path.join(OUTDIR, f)
    n = sum(1 for _ in open(path)) - 1
    print(f"  {f}: {n:,} ({os.path.getsize(path)/1024/1024:.1f} MB)")