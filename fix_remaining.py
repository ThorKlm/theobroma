import pandas as pd
import os, json, tarfile, warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, inchi
RDLogger.logger().setLevel(RDLogger.ERROR)

OUTDIR = "data/converted"
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
    try: return inchi.InchiToInchiKey(inchi.MolToInchi(mol))
    except: return ""

def save(records, name):
    df = pd.DataFrame(records, columns=COLS)
    df = df[df["smiles"].astype(str).str.len() > 2]
    path = os.path.join(OUTDIR, f"{name}.csv")
    df.to_csv(path, index=False)
    print(f"  -> {len(df):,} saved to {name}.csv")

# Replace the ANPDB section with:
print("[ANPDB]")
from rdkit import Chem
suppl = Chem.ForwardSDMolSupplier(open("data/ANPDB.sdf", "rb"), removeHs=True, sanitize=True)
recs = []
for idx, mol in enumerate(tqdm(suppl, desc="  ANPDB", leave=False)):
    if mol is None: continue
    smi = Chem.MolToSmiles(mol)
    if not smi: continue
    nm = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
    p, _ = props(smi)
    recs.append([f"ANPDB_{idx:06d}", nm, smi, "", get_ik(mol),
                  "ANPDB", "multi", "Africa", "",
                  p.get("mw"), p.get("logp"), p.get("tpsa"),
                  p.get("hba"), p.get("hbd"), p.get("n_rings"),
                  p.get("rotatable_bonds"), "CC BY 4.0"])
save(recs, "anpdb")

# ==== 2. VIETHERB -- no SMILES available ====
print("\n[VIETHERB] No SMILES in OWL -- only formula/MW/CAS. Skipping.")
print("  Predicates found: city_name, formula, MW, CAS_ID, com_name, iupac_name")
print("  Would need PubChem CID lookup to recover SMILES from CAS/names.")

# Try to recover via CAS/name -> PubChem lookup
import rdflib
g = rdflib.Graph()
g.parse("data/VHO.owl", format="xml")
# Extract what we have
recs = []
idx = 0
for s in set(g.subjects()):
    entry = {}
    for _, p, o in g.triples((s, None, None)):
        pk = str(p).split("/")[-1].split("#")[-1]
        entry[pk] = str(o)
    if "meta:C_ID" not in entry and "meta:com_name" not in entry:
        continue
    # We have compound ID or name but no SMILES
    idx += 1

print(f"  Found {idx} compound entries without SMILES")
print("  TODO: batch lookup CAS IDs via PubChem PUG-REST to get SMILES")

# ==== 3. MIBiG -- scan all JSONs for compounds ====
print("\n[MIBiG]")
recs, idx = [], 0
with tarfile.open("data/mibig_json_4.0.tar.gz", "r:gz") as tar:
    jsons = [m for m in tar.getmembers() if m.name.endswith(".json")]
    has_comps = 0
    for member in tqdm(jsons, desc="  MIBiG", leave=False):
        fobj = tar.extractfile(member)
        if fobj is None: continue
        try:
            data = json.loads(fobj.read())
        except:
            continue
        # Try multiple paths for compounds
        comps = data.get("cluster", {}).get("compounds", [])
        if not comps:
            # try alternate structure
            comps = data.get("compounds", [])
        if comps:
            has_comps += 1
        for comp in comps:
            smi = comp.get("chem_struct", comp.get("smiles", comp.get("structure", "")))
            if not smi or len(smi) < 3: continue
            p, mol = props(smi)
            nm = comp.get("compound", comp.get("name", ""))
            recs.append([f"MIBIG_{idx:06d}", nm, smi, "", get_ik(mol),
                          "MIBiG", "multi", "global", "",
                          p.get("mw"), p.get("logp"), p.get("tpsa"),
                          p.get("hba"), p.get("hbd"), p.get("n_rings"),
                          p.get("rotatable_bonds"), "CC BY 4.0"])
            idx += 1
print(f"  JSONs with compounds field: {has_comps}/{len(jsons)}")
save(recs, "mibig")

# ==== Summary ====
print("\n=== All converted ===")
total = 0
for f in sorted(os.listdir(OUTDIR)):
    if not f.endswith(".csv"): continue
    path = os.path.join(OUTDIR, f)
    df = pd.read_csv(path)
    total += len(df)
    print(f"  {f}: {len(df):,}")
print(f"\n  TOTAL: {total:,}")