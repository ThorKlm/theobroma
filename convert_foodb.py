"""Unified converter for all THEOBROMA data sources.

Usage:
    python convert.py foodb
    python convert.py fdp
    python convert.py naturar
    python convert.py cmaupv2
    python convert.py aspergillus
    python convert.py enrich_cmaup
    python convert.py enrich_npass
    python convert.py enrich_gbif
    python convert.py normalize_regions
    python convert.py check
    python convert.py all          # run all converters
"""
import argparse, os, sys, re, json, zipfile, warnings
warnings.filterwarnings("ignore")
import pandas as pd
from tqdm import tqdm

OUTDIR = "data/converted"
os.makedirs(OUTDIR, exist_ok=True)

COLS = ["comp_id","name","smiles","inchi","inchikey","source_db","kingdom",
        "region","source_organism","mw","logp","tpsa","hba","hbd",
        "n_rings","rotatable_bonds","license_tier"]

# ============================================================
# FooDB
# ============================================================
def convert_foodb():
    ZIPFILE = "data/databases/foodb_2020_04_07_json.zip"
    OUTFILE = os.path.join(OUTDIR, "foodb.csv")
    if not os.path.exists(ZIPFILE):
        print(f"  SKIP: {ZIPFILE} not found"); return

    KINGDOM_KW = {
        "fungi": ["mushroom","truffle","yeast","fungus","fungi","shiitake",
                  "maitake","reishi","enoki","chanterelle","portobello",
                  "agaricus","pleurotus","lentinula","ganoderma"],
        "bacteria": ["fermented","kefir","yogurt","kombucha","tempeh","miso","natto"],
        "marine": ["seaweed","kelp","algae","nori","wakame","spirulina","chlorella",
                   "fish","shrimp","crab","lobster","mussel","oyster","clam","squid"],
        "animal": ["beef","pork","chicken","turkey","lamb","goat","duck","venison",
                   "egg","milk","cheese","butter","cream","liver","kidney","heart",
                   "honey","beeswax","propolis","lard","tallow","gelatin"],
    }
    def infer_kingdom(food_name):
        fl = food_name.lower()
        for k, kws in KINGDOM_KW.items():
            if any(w in fl for w in kws):
                return k
        return "plant"

    print("=== FooDB ===")
    print("  Loading Food.json...")
    food_map = {}
    with zipfile.ZipFile(ZIPFILE, "r") as zf:
        with zf.open("foodb_2020_04_07_json/Food.json") as f:
            for line in f:
                obj = json.loads(line)
                food_map[obj.get("id")] = {
                    "species": str(obj.get("name_scientific") or obj.get("name","")),
                    "common": str(obj.get("name",""))
                }
    print(f"  Foods: {len(food_map):,}")

    print("  Loading Content.json (compound -> food mapping)...")
    compound_foods = {}
    with zipfile.ZipFile(ZIPFILE, "r") as zf:
        with zf.open("foodb_2020_04_07_json/Content.json") as f:
            for line in tqdm(f, desc="  Content"):
                obj = json.loads(line)
                if obj.get("source_type") == "Compound":
                    cid, fid = obj.get("source_id"), obj.get("food_id")
                    if cid and fid:
                        compound_foods.setdefault(cid, set()).add(fid)
    print(f"  Compounds with food associations: {len(compound_foods):,}")

    print("  Loading Compound.json...")
    from collections import Counter
    records = []
    kingdom_counts = {}
    with zipfile.ZipFile(ZIPFILE, "r") as zf:
        with zf.open("foodb_2020_04_07_json/Compound.json") as f:
            for line in tqdm(f, desc="  Compounds"):
                c = json.loads(line)
                smiles = c.get("moldb_smiles") or ""
                if not smiles or len(smiles) < 3:
                    continue
                cid = c.get("id")
                food_ids = compound_foods.get(cid, set())
                species_list, kingdoms = [], []
                for fid in food_ids:
                    info = food_map.get(fid, {})
                    sp = info.get("species", "")
                    if sp and sp != "None":
                        species_list.append(sp)
                    kingdoms.append(infer_kingdom(info.get("common","") or sp))
                best_k = Counter(kingdoms).most_common(1)[0][0] if kingdoms else "plant"
                kingdom_counts[best_k] = kingdom_counts.get(best_k, 0) + 1
                records.append({
                    "comp_id": f"FOODB_{len(records):06d}",
                    "name": c.get("name",""), "smiles": smiles,
                    "inchi": c.get("moldb_inchi","") or "",
                    "inchikey": c.get("moldb_inchikey","") or "",
                    "source_db": "FooDB", "kingdom": best_k, "region": "global",
                    "source_organism": "; ".join(sorted(set(species_list))[:3]),
                    "mw": c.get("moldb_mono_mass",""), "logp": "", "tpsa": "",
                    "hba": "", "hbd": "", "n_rings": "", "rotatable_bonds": "",
                    "license_tier": "CC BY-NC 4.0",
                })
    print(f"  Extracted: {len(records):,}")
    for k, cnt in sorted(kingdom_counts.items(), key=lambda x: -x[1]):
        print(f"    {k}: {cnt:,}")
    pd.DataFrame(records)[COLS].to_csv(OUTFILE, index=False)
    print(f"  Saved to {OUTFILE}")


# ============================================================
# FDP (LMDB, SMDB, TMDB, CMDB)
# ============================================================
def convert_fdp():
    INPUT = os.path.join(OUTDIR, "fdp_all.csv")
    if not os.path.exists(INPUT):
        INPUT = "fdp_scraper/fdp_all.csv"
    if not os.path.exists(INPUT):
        print("  SKIP: fdp_all.csv not found (already split?)"); return

    print("=== FDP ===")
    KINGDOM_MAP = {"LMDB":"fungi", "SMDB":"plant", "TMDB":"fungi", "CMDB":"plant"}
    SOURCE_MAP = {"LMDB":"LMDB_Lichen", "SMDB":"SMDB_Spice", "TMDB":"TMDB_Trichoderma", "CMDB":"CMDB_Cereals"}

    df = pd.read_csv(INPUT)
    print(f"  Input: {len(df):,}")
    out = pd.DataFrame({
        "comp_id": [f"{SOURCE_MAP.get(r['database'],r['database'])}_{i:06d}" for i,r in df.iterrows()],
        "name": df["common_name"].fillna(""), "smiles": df["canonical_smiles"].fillna(""),
        "inchi": "", "inchikey": "",
        "source_db": df["database"].map(SOURCE_MAP).fillna(df["database"]),
        "kingdom": df["database"].map(KINGDOM_MAP).fillna("plant"),
        "region": "South Asia",
        "source_organism": df["species"].fillna(""),
        "mw": df["molecular_weight"], "logp": "", "tpsa": "", "hba": "", "hbd": "",
        "n_rings": "", "rotatable_bonds": "", "license_tier": "CC BY 4.0",
    })
    out = out[out["smiles"].astype(str).str.len() > 2]
    for db in sorted(out["source_db"].unique()):
        subset = out[out["source_db"] == db]
        fname = db.lower().replace(" ","_") + ".csv"
        subset[COLS].to_csv(os.path.join(OUTDIR, fname), index=False)
        print(f"    {fname}: {len(subset):,}")
    os.remove(INPUT)
    print("  Done.")


# ============================================================
# NaturAr
# ============================================================
def convert_naturar():
    INPUT = "additional_data/NaturAr_query.csv"
    OUTFILE = os.path.join(OUTDIR, "naturar.csv")
    if not os.path.exists(INPUT):
        print("  SKIP: NaturAr not found"); return
    print("=== NaturAr ===")
    df = pd.read_csv(INPUT)
    out = pd.DataFrame({
        "comp_id": [f"NATAR_{i:06d}" for i in range(len(df))],
        "name": df["Synonym"].fillna(""), "smiles": df["SMILES"].fillna(""),
        "inchi": "", "inchikey": "",
        "source_db": "NaturAr", "kingdom": "plant", "region": "Latin America",
        "source_organism": df["Source_species"].fillna(""),
        "mw": "", "logp": "", "tpsa": "", "hba": "", "hbd": "",
        "n_rings": "", "rotatable_bonds": "", "license_tier": "CC BY 4.0",
    })
    out = out[out["smiles"].astype(str).str.len() > 2]
    out[COLS].to_csv(OUTFILE, index=False)
    print(f"  Saved: {len(out):,} to {OUTFILE}")


# ============================================================
# CMAUPv2
# ============================================================
def convert_cmaupv2():
    INPUT = "additional_data/CMAUPv2.0_download_Ingredients_All.txt"
    OUTFILE = os.path.join(OUTDIR, "cmaupv2.csv")
    if not os.path.exists(INPUT):
        print("  SKIP: CMAUPv2 not found"); return
    print("=== CMAUPv2 ===")
    cmaup = pd.read_csv(INPUT, sep="\t", low_memory=False)
    out = pd.DataFrame({
        "comp_id": [f"CMAUP2_{i:06d}" for i in range(len(cmaup))],
        "name": cmaup["pref_name"].fillna(""), "smiles": cmaup["SMILES"].fillna(""),
        "inchi": cmaup["InChI"].fillna(""), "inchikey": cmaup["InChIKey"].fillna(""),
        "source_db": "CMAUPv2", "kingdom": "plant", "region": "East Asia",
        "source_organism": "",
        "mw": cmaup["MW"], "logp": cmaup["LogP"], "tpsa": cmaup["TPSA"],
        "hba": cmaup["nHA"], "hbd": cmaup["nHD"],
        "n_rings": cmaup.get("nRing",""), "rotatable_bonds": cmaup.get("nRot",""),
        "license_tier": "CC BY 4.0",
    })
    out = out[out["smiles"].astype(str).str.len() > 2]
    out[COLS].to_csv(OUTFILE, index=False)
    print(f"  Saved: {len(out):,} to {OUTFILE}")


# ============================================================
# Enrich CMAUPv2 with plant species names
# ============================================================
def enrich_cmaup():
    PLANTS_FILE = "additional_data/CMAUPv2.0_download_Plants.txt"
    ASSOC_FILE = "additional_data/CMAUPv2.0_download_Plant_Ingredient_Associations_allIngredients.txt"
    INGR_FILE = "additional_data/CMAUPv2.0_download_Ingredients_All.txt"
    CMAUP_CSV = os.path.join(OUTDIR, "cmaupv2.csv")
    for f in [PLANTS_FILE, ASSOC_FILE, INGR_FILE, CMAUP_CSV]:
        if not os.path.exists(f):
            print(f"  SKIP: {f} not found"); return

    print("=== Enrich CMAUPv2 with species names ===")
    plants = pd.read_csv(PLANTS_FILE, sep="\t")
    plant_names = {}
    for _, row in plants.iterrows():
        pid = str(row.iloc[0]).strip()
        sp = str(row.get("Species_Name","")).strip()
        if sp == "nan" or not sp:
            sp = str(row.get("Plant_Name","")).strip()
        if sp != "nan" and sp:
            plant_names[pid] = sp
    print(f"  Plant names: {len(plant_names):,}")

    assoc = pd.read_csv(ASSOC_FILE, sep="\t", header=None, names=["plant_id","np_id"])
    np_to_species = {}
    for _, row in assoc.iterrows():
        name = plant_names.get(row["plant_id"])
        if name:
            np_to_species.setdefault(row["np_id"], set()).add(name)
    np_species_str = {k: "; ".join(sorted(v)[:5]) for k,v in np_to_species.items()}

    cmaup_orig = pd.read_csv(INGR_FILE, sep="\t", low_memory=False, usecols=["np_id"])
    cmaup = pd.read_csv(CMAUP_CSV, low_memory=False)
    cmaup["source_organism"] = cmaup_orig["np_id"].iloc[:len(cmaup)].map(np_species_str).fillna("")
    cmaup.to_csv(CMAUP_CSV, index=False)
    with_org = (cmaup["source_organism"] != "").sum()
    print(f"  With species: {with_org:,}/{len(cmaup):,}")


# ============================================================
# Enrich NPASS organisms (resolve NPO codes)
# ============================================================
def enrich_npass():
    SPECIES_FILE = "additional_data/npass_speciesInfo.txt"
    PLANTS_FILE = "additional_data/CMAUPv2.0_download_Plants.txt"
    PAIRS_FILE = "additional_data/NPASS3.0_naturalproducts_species_pair.txt"
    GI_FILE = "additional_data/NPASS3.0_naturalproducts_generalinfo.txt"
    FINAL = "data/theobroma_final.csv"
    if not os.path.exists(FINAL):
        print("  SKIP: theobroma_final.csv not found"); return

    print("=== NPASS Organism Enrichment ===")
    # Build NPO -> species name from both CMAUP plants and NPASS species info
    npo_to_name = {}

    if os.path.exists(PLANTS_FILE):
        plants = pd.read_csv(PLANTS_FILE, sep="\t")
        for _, row in plants.iterrows():
            pid = str(row.iloc[0]).strip()
            sp = str(row.get("Species_Name","")).strip()
            if sp == "nan" or not sp:
                sp = str(row.get("Plant_Name","")).strip()
            if sp != "nan" and sp:
                npo_to_name[pid] = sp
        print(f"  CMAUP plant names: {len(npo_to_name):,}")

    if os.path.exists(SPECIES_FILE):
        sp_df = pd.read_csv(SPECIES_FILE, sep="\t", low_memory=False, encoding="latin-1")
        for _, row in sp_df.iterrows():
            oid = str(row.get("org_id","")).strip()
            name = str(row.get("species_name","")).strip()
            if name == "nan" or not name:
                name = str(row.get("org_name","")).strip()
            if name != "nan" and name and oid not in npo_to_name:
                npo_to_name[oid] = name
        print(f"  Total NPO names (CMAUP+NPASS): {len(npo_to_name):,}")

    # Build np_id -> organism via species pairs
    if os.path.exists(PAIRS_FILE) and os.path.exists(GI_FILE):
        pairs = pd.read_csv(PAIRS_FILE, sep="\t", low_memory=False, usecols=["org_id","np_id"])
        np_to_org = {}
        for _, row in pairs.iterrows():
            name = npo_to_name.get(str(row["org_id"]).strip())
            if name:
                np_to_org.setdefault(str(row["np_id"]).strip(), set()).add(name)
        np_best = {k: "; ".join(sorted(v)[:3]) for k,v in np_to_org.items()}

        gi = pd.read_csv(GI_FILE, sep="\t", low_memory=False, usecols=["np_id","inchikey"])
        ik_to_org = {}
        for _, row in gi.iterrows():
            npid = str(row["np_id"]).strip()
            ik = str(row.get("inchikey","")).strip()
            if ik and len(ik) > 10 and npid in np_best:
                ik_to_org[ik] = np_best[npid]
        print(f"  InChIKeys with organism: {len(ik_to_org):,}")

    # Apply to theobroma_final.csv
    df = pd.read_csv(FINAL, low_memory=False)
    ik_col = df["inchikey"].fillna("").astype(str).str.strip()
    org_col = df["source_organism"].fillna("").astype(str).str.strip()

    # Fix empty organisms via InChIKey lookup
    empty_org = (org_col == "") | (org_col == "nan")
    new_orgs = ik_col.map(ik_to_org).fillna("")
    fill_mask = empty_org & (new_orgs != "")
    df.loc[fill_mask, "source_organism"] = new_orgs[fill_mask]
    print(f"  Filled empty: {fill_mask.sum():,}")

    # Fix NPO codes
    is_npo = org_col.str.contains(r"NPO\d+", na=False)
    npo_new = ik_col.map(ik_to_org).fillna("")
    npo_mask = is_npo & (npo_new != "")
    df.loc[npo_mask, "source_organism"] = npo_new[npo_mask]
    print(f"  Fixed NPO codes: {npo_mask.sum():,}")

    # Direct NPO code resolution for remaining
    still_npo = df["source_organism"].fillna("").astype(str).str.contains(r"NPO\d+", na=False)
    for idx in df[still_npo].index:
        m = re.search(r"(NPO\d+)", str(df.at[idx, "source_organism"]))
        if m:
            name = npo_to_name.get(m.group(1))
            if name:
                df.at[idx, "source_organism"] = name

    org_total = (df["source_organism"].fillna("").astype(str).str.strip() != "").sum()
    npo_left = df["source_organism"].fillna("").astype(str).str.contains(r"NPO\d+", na=False).sum()
    print(f"  Total with organism: {org_total:,} ({100*org_total/len(df):.1f}%)")
    print(f"  Remaining NPO codes: {npo_left:,}")
    df.to_csv(FINAL, index=False)
    print(f"  Saved.")


# ============================================================
# GBIF geographic enrichment
# ============================================================
def enrich_gbif():
    GBIF_FILE = "additional_data/cmaup_plant_countries.csv"
    ASSOC_FILE = "additional_data/CMAUPv2.0_download_Plant_Ingredient_Associations_allIngredients.txt"
    INGR_FILE = "additional_data/CMAUPv2.0_download_Ingredients_All.txt"
    FINAL = "data/theobroma_final.csv"
    for f in [GBIF_FILE, FINAL]:
        if not os.path.exists(f):
            print(f"  SKIP: {f} not found"); return

    print("=== GBIF Geographic Enrichment ===")
    from collections import Counter
    gbif = pd.read_csv(GBIF_FILE)
    resolved = gbif[gbif["region"] != "unresolved"]
    gbif_map = dict(zip(resolved["plant_id"], resolved["region"]))
    print(f"  GBIF resolved plants: {len(gbif_map):,}")

    if os.path.exists(ASSOC_FILE) and os.path.exists(INGR_FILE):
        assoc = pd.read_csv(ASSOC_FILE, sep="\t", header=None, names=["plant_id","np_id"])
        np_regions = {}
        for _, row in assoc.iterrows():
            r = gbif_map.get(row["plant_id"])
            if r:
                np_regions.setdefault(row["np_id"], []).append(r)
        np_best = {k: Counter(v).most_common(1)[0][0] for k,v in np_regions.items()}

        cmaup_ingr = pd.read_csv(INGR_FILE, sep="\t", low_memory=False, usecols=["np_id","InChIKey"])
        ik_to_region = {}
        for _, row in cmaup_ingr.iterrows():
            npid, ik = row["np_id"], str(row.get("InChIKey","")).strip()
            if ik and len(ik) > 10 and npid in np_best:
                ik_to_region[ik] = np_best[npid]
        print(f"  InChIKeys with region: {len(ik_to_region):,}")

        df = pd.read_csv(FINAL, low_memory=False)
        ik_col = df["inchikey"].fillna("").astype(str).str.strip()
        reg_col = df["region"].fillna("").astype(str).str.strip()
        unresolved = (reg_col == "") | (reg_col == "nan") | (reg_col == "global") | (reg_col == "East Asia")
        new_regs = ik_col.map(ik_to_region).fillna("")
        mask = unresolved & (new_regs != "")
        df.loc[mask, "region"] = new_regs[mask]
        print(f"  Filled region: {mask.sum():,}")
        df.to_csv(FINAL, index=False)
        print(f"  Saved.")


# ============================================================
# Normalize regions
# ============================================================
def normalize_regions():
    FINAL = "data/theobroma_final.csv"
    if not os.path.exists(FINAL):
        print("  SKIP: theobroma_final.csv not found"); return

    print("=== Normalize Regions ===")
    df = pd.read_csv(FINAL, low_memory=False)

    VALID = {"East Asia","South Asia","Southeast Asia","Australia","Africa","Europe",
             "Latin America","North America","Russia/CIS","Central Asia","Middle East",
             "Oceania","New Zealand","Caribbean","Antarctic",""}

    LOC_MAP = {
        "Indonesian":"Southeast Asia", "Formosan soft coral":"East Asia",
        "South China Sea":"East Asia", "Okinawan":"East Asia", "Korean":"East Asia",
        "Taiwan":"East Asia", "Tibetan":"East Asia", "Chinese":"East Asia",
        "Vietnamese":"Southeast Asia", "Thai":"Southeast Asia",
        "Philippine":"Southeast Asia", "Philippines":"Southeast Asia",
        "Myanmar":"Southeast Asia", "Palauan":"Oceania",
        "Papua New Guinea":"Oceania", "Brazilian":"Latin America",
        "Mexican":"Latin America", "South African":"Africa",
        "Tanzanian":"Africa", "Madagascar rainforest":"Africa",
        "Mediterranean":"Europe", "Australian":"Australia",
        "Indian Ocean":"Africa", "Suriname rainforest":"Latin America",
    }

    KW_MAP = [
        (["china","yunnan","sichuan","hainan","guangxi","jiangxi","tibet","beijing",
          "shanghai","guangdong","fujian","zhejiang","hubei","hunan","hebei","henan",
          "shandong","shaanxi","guizhou","anhui","jilin","liaoning","heilongjiang",
          "gansu","qinghai","xinjiang","ningxia","inner mongolia","chongqing",
          "weizhou","wenshan","xishuangbanna","dongfang","longlin","pingle",
          "fengqi","zhongdian"], "East Asia"),
        (["japan","okinawa","korea","jeju","formosa","korean","japanese"], "East Asia"),
        (["vietnam","hanoi","thai","myanmar","philippines","indonesia",
          "malaysian","cambodia","laos"], "Southeast Asia"),
        (["india","hyderabad","andhra","bengal","kerala","tamil","karnataka",
          "maharashtra","rajasthan","gujarat","assam","manipur","meghalaya",
          "mizoram","nagaland","sikkim","tripura","uttarakhand"], "South Asia"),
        (["africa","mali","tanzania","madagascar","cameroon","nigeria","kenya",
          "ethiopia","ghana","senegal","congo","uganda","angola","mozambique",
          "zimbabwe","zambia","malawi","rwanda","burkina","benin","togo","guinea",
          "sierra","liberia","namibia","botswana","south africa"], "Africa"),
        (["caribbean","suriname","brazil","mexico","argentina","colombia","peru",
          "chile","venezuela","ecuador","bolivia","costa rica","panama","cuba",
          "guatemala","honduras"], "Latin America"),
        (["antarctic","antarctica"], "Antarctic"),
        (["mediterranean","tenerife","sardinia","spain","france","italy","greece",
          "germany","portugal","poland","hungary","bulgaria","romania","croatia",
          "serbia","turkey"], "Europe"),
        (["australia","queensland","new south wales","victoria","tasmania"], "Australia"),
    ]

    BODY_TERMS = {"blood","serum","milk","body","cytoplasm","fruit body","mycelium",
                  "subtidal","coral","urine","saliva","faeces","cerebrospinal",
                  "breast milk","blood serum","aerial part","root","leaf","stem",
                  "seed","fruit","twig","bark","flower","rhizome","whole plant"}

    reg_col = df["region"].fillna("").astype(str).str.strip()
    df["region"] = reg_col.replace({"global": ""})
    reg_col = df["region"].fillna("").astype(str).str.strip()

    fixed = 0
    for idx in range(len(df)):
        r = reg_col.iloc[idx]
        if r in VALID:
            continue
        mapped = LOC_MAP.get(r)
        if mapped:
            df.at[idx, "region"] = mapped; fixed += 1; continue
        rl = r.lower()
        if rl in BODY_TERMS or any(b in rl for b in BODY_TERMS):
            df.at[idx, "region"] = ""; fixed += 1; continue
        done = False
        for kws, region in KW_MAP:
            if any(k in rl for k in kws):
                df.at[idx, "region"] = region; fixed += 1; done = True; break
        if not done:
            df.at[idx, "region"] = ""
            fixed += 1

    print(f"  Fixed: {fixed:,}")
    regions = df["region"].fillna("unresolved").replace({"":"unresolved"})
    for r, cnt in regions.value_counts().head(15).items():
        print(f"    {r}: {cnt:,}")
    df.to_csv(FINAL, index=False)
    print(f"  Saved.")


# ============================================================
# Check all source files
# ============================================================
def check_sources():
    print("=== Source file check ===")
    ENRICHMENT = {"npass_enrichment.csv","npass_geo_enrichment.csv"}
    for f in sorted(os.listdir(OUTDIR)):
        if not f.endswith(".csv"): continue
        path = os.path.join(OUTDIR, f)
        try:
            df = pd.read_csv(path, nrows=3, low_memory=False)
            n = len(pd.read_csv(path, usecols=[0], low_memory=False))
            has_cols = {"comp_id","name","smiles","source_db"}.issubset(set(df.columns))
            tag = "ENRICHMENT" if f in ENRICHMENT else ("OK" if has_cols else "WRONG COLS")
            region = df["region"].iloc[0] if "region" in df.columns else "N/A"
            print(f"  {f:35s} {n:>8,}  region={str(region):15s}  [{tag}]")
        except Exception as e:
            print(f"  {f:35s}  ERROR: {e}")

    if os.path.exists("data/theobroma_final.csv"):
        df = pd.read_csv("data/theobroma_final.csv", low_memory=False)
        print(f"\n  theobroma_final.csv: {len(df):,} compounds")
        org = (df["source_organism"].fillna("").astype(str).str.strip() != "").sum()
        npo = df["source_organism"].fillna("").astype(str).str.contains("NPO",na=False).sum()
        print(f"  With organism: {org:,} (real: {org-npo:,}, NPO codes: {npo:,})")
        print(f"  Sources: {df['source_db'].nunique()}")
        print(f"  Multi-source: {(df['all_sources'].fillna('').str.contains('|',regex=False,na=False)).sum():,}")


# ============================================================
# Main
# ============================================================
COMMANDS = {
    "foodb": convert_foodb,
    "fdp": convert_fdp,
    "naturar": convert_naturar,
    "cmaupv2": convert_cmaupv2,
    "enrich_cmaup": enrich_cmaup,
    "enrich_npass": enrich_npass,
    "enrich_gbif": enrich_gbif,
    "normalize_regions": normalize_regions,
    "check": check_sources,
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="THEOBROMA unified converter")
    parser.add_argument("command", nargs="?", default="check",
                        choices=list(COMMANDS.keys()) + ["all"],
                        help="Which conversion to run")
    args = parser.parse_args()

    if args.command == "all":
        for name, func in COMMANDS.items():
            if name != "check":
                print(f"\n{'='*60}")
                func()
        print(f"\n{'='*60}")
        check_sources()
    else:
        COMMANDS[args.command]()