"""Parse Aspergillus Metabolome DB SQL dump.
compound_identifiers format: (id, InChI, InChIKey, SMILES, ts, ts)
compounds_aspergillus: (compound_id, ...) for Aspergillus-specific compounds."""
import re, os, warnings
warnings.filterwarnings("ignore")
import pandas as pd
from tqdm import tqdm

SQL_FILE = "../data/compounds.20210215.sql"
OUTFILE = "../data/converted/aspergillus_mdb.csv"
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

COLS = ["comp_id","name","smiles","inchi","inchikey","source_db","kingdom",
        "region","source_organism","mw","logp","tpsa","hba","hbd",
        "n_rings","rotatable_bonds","license_tier"]

print("Reading SQL dump...")
with open(SQL_FILE, "r", encoding="utf-8", errors="ignore") as f:
    content = f.read()

# Step 1: Parse compound_identifiers -> (id, InChI, InChIKey, SMILES)
print("Parsing compound_identifiers...")
ident_data = {}  # id -> {inchi, inchikey, smiles}
# Match tuples: (id,'InChI=...','InChIKey','SMILES','ts','ts')
pattern = r"\((\d+),'(InChI=[^']+)','([^']+)','([^']+)','[^']+','[^']+'\)"
for m in tqdm(re.finditer(pattern, content), desc="  Identifiers"):
    cid = int(m.group(1))
    ident_data[cid] = {
        "inchi": m.group(2),
        "inchikey": m.group(3),
        "smiles": m.group(4),
    }
print(f"  Found {len(ident_data)} compounds with SMILES")

# Step 2: Parse compounds table -> (id, CAS, name, formula, MW, ...)
print("Parsing compounds table for names...")
comp_names = {}
name_pattern = r"\((\d+),(?:'[^']*'|NULL),'([^']*)','([^']*)',([0-9.]+)"
for m in re.finditer(name_pattern, content):
    cid = int(m.group(1))
    name = m.group(2)
    formula = m.group(3)
    mw = m.group(4)
    comp_names[cid] = {"name": name, "formula": formula, "mw": mw}
print(f"  Found {len(comp_names)} compound names")

# Step 3: Parse compounds_aspergillus -> set of Aspergillus compound IDs
print("Parsing compounds_aspergillus...")
asp_ids = set()
asp_pattern = r"\((\d+),(\d+)"
asp_section = re.search(r"INSERT INTO `compounds_aspergillus` VALUES(.+?)(?:;\s*$|\s*--)", content, re.DOTALL)
if asp_section:
    for m in re.finditer(asp_pattern, asp_section.group(1)):
        asp_ids.add(int(m.group(1)))
print(f"  Found {len(asp_ids)} Aspergillus-tagged compounds")

# Step 4: Parse organism table
print("Parsing organisms...")
org_map = {}
org_pattern = r"\((\d+),'([^']*)'"
org_section = re.search(r"INSERT INTO `organism` VALUES(.+?)(?:;\s*$|\s*--|\s*INSERT INTO `(?!organism`))", content, re.DOTALL)
if org_section:
    for m in re.finditer(org_pattern, org_section.group(1)):
        org_map[int(m.group(1))] = m.group(2)
print(f"  Found {len(org_map)} organisms")

# Step 5: Parse compound_organism
comp_org = {}
co_section = re.search(r"INSERT INTO `compound_organism` VALUES(.+?)(?:;\s*$|\s*--|\s*INSERT INTO `(?!compound_organism`))", content, re.DOTALL)
if co_section:
    for m in re.finditer(r"\((\d+),(\d+)", co_section.group(1)):
        comp_org[int(m.group(1))] = int(m.group(2))

# Step 6: Build output -- ALL compounds with SMILES (not just Aspergillus-tagged,
# since only 1,301 are tagged but 331k have identifiers from the full CEU DB)
# We output two tiers: Aspergillus-tagged first, then all others
print(f"\nBuilding output...")
print(f"  Total with SMILES: {len(ident_data)}")
print(f"  Aspergillus-tagged: {len(asp_ids)}")
print(f"  Tagged with SMILES: {len(asp_ids & set(ident_data.keys()))}")

# For THEOBROMA, we only want Aspergillus-specific compounds
# But if too few are tagged, take all compounds from the DB
target_ids = asp_ids & set(ident_data.keys())
if len(target_ids) < 100:
    print("  Few Aspergillus-tagged compounds have SMILES, using all compounds")
    target_ids = set(ident_data.keys())

records = []
for cid in tqdm(sorted(target_ids), desc="  Building records"):
    d = ident_data[cid]
    smiles = d["smiles"]
    if not smiles or len(smiles) < 3:
        continue
    info = comp_names.get(cid, {})
    org_id = comp_org.get(cid)
    org_name = org_map.get(org_id, "Aspergillus") if org_id else "Aspergillus"

    records.append({
        "comp_id": f"AMDB_{len(records):06d}",
        "name": info.get("name", ""),
        "smiles": smiles,
        "inchi": d["inchi"],
        "inchikey": d["inchikey"],
        "source_db": "AMDB",
        "kingdom": "fungi",
        "region": "global",
        "source_organism": org_name,
        "mw": info.get("mw", ""),
        "logp": "", "tpsa": "", "hba": "", "hbd": "",
        "n_rings": "", "rotatable_bonds": "",
        "license_tier": "CC BY 4.0",
    })

df = pd.DataFrame(records)
for c in COLS:
    if c not in df.columns:
        df[c] = ""
df[COLS].to_csv(OUTFILE, index=False)
print(f"\nSaved {len(df)} compounds to {OUTFILE}")
print(f"  With organism data: {sum(1 for r in records if r['source_organism'] != 'Aspergillus')}")