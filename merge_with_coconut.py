"""Full pipeline: add IMPPAT, extract COCONUT, cross-dedup, final merge."""
import pandas as pd
import os, warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
RDLogger.logger().setLevel(RDLogger.ERROR)

DATADIR = "data"
CONVERTED = os.path.join(DATADIR, "converted")
MERGED = os.path.join(DATADIR, "theobroma_merged.csv")
FINAL = os.path.join(DATADIR, "theobroma_final.csv")

COLS = ["comp_id","name","smiles","inchi","inchikey","source_db","all_sources",
        "kingdom","region","source_organism","mw","logp","tpsa","hba","hbd",
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
    try: return InchiToInchiKey(MolToInchi(mol))
    except: return ""

# ============================================================
# STEP 1: Convert IMPPAT full scrape if available
# ============================================================
imppat_file = os.path.join(DATADIR, "imppat2_compounds.csv")
if os.path.exists(imppat_file):
    print(f"[IMPPAT full] Loading {imppat_file}...")
    df = pd.read_csv(imppat_file, low_memory=False)
    df.columns = [c.strip().lower().replace(" ","_") for c in df.columns]
    sc = [c for c in df.columns if "smiles" in c.lower()]
    if sc:
        sc = sc[0]
        recs = []
        for idx, (i, row) in enumerate(tqdm(df.iterrows(), total=len(df), desc="  IMPPAT full", leave=False)):
            smi = str(row.get(sc, ""))
            if len(smi) < 3 or smi == "nan": continue
            p, mol = props(smi)
            ik = get_ik(mol)
            nm = str(row.get("name", row.get("comp_id", "")))
            if nm == "nan": nm = ""
            org = str(row.get("source_organism", row.get("plant_sources", "")))
            if org == "nan": org = ""
            recs.append([f"IMPPAT_{idx:06d}", nm, smi, "", ik,
                          "IMPPAT", "IMPPAT", "plant", "South Asia", org,
                          p.get("mw"), p.get("logp"), p.get("tpsa"),
                          p.get("hba"), p.get("hbd"), p.get("n_rings"),
                          p.get("rotatable_bonds"), "CC BY-NC 4.0"])
        out = os.path.join(CONVERTED, "imppat_full.csv")
        outdf = pd.DataFrame(recs, columns=COLS)
        outdf = outdf[outdf["smiles"].astype(str).str.len() > 2]
        outdf.to_csv(out, index=False)
        print(f"  -> {len(outdf):,} saved to imppat_full.csv")
    else:
        print(f"  No SMILES column. Cols: {list(df.columns)[:10]}")
else:
    print("[IMPPAT full] Not found. Skipping.")

# ============================================================
# STEP 2: Load all supplementary CSVs and dedup
# ============================================================
print("\n[Loading all converted CSVs]")
frames = []
for f in sorted(os.listdir(CONVERTED)):
    if not f.endswith(".csv"): continue
    path = os.path.join(CONVERTED, f)
    df = pd.read_csv(path, low_memory=False)
    if len(df) == 0: continue
    frames.append(df)
    print(f"  {f}: {len(df):,}")

supp = pd.concat(frames, ignore_index=True)
print(f"\n  Total supplementary: {len(supp):,}")

# Compute missing InChIKeys
missing = supp["inchikey"].isna() | (supp["inchikey"].astype(str) == "") | (supp["inchikey"].astype(str) == "nan") | (supp["inchikey"].astype(str).str.len() < 10)
if missing.sum() > 0:
    print(f"  Computing {missing.sum():,} missing InChIKeys...")
    for idx in tqdm(supp[missing].index, desc="  InChIKeys", leave=False):
        smi = str(supp.at[idx, "smiles"])
        if len(smi) < 3: continue
        mol = Chem.MolFromSmiles(smi)
        if mol:
            try: supp.at[idx, "inchikey"] = InchiToInchiKey(MolToInchi(mol))
            except: pass

# Dedup supplementary internally
print("  Deduplicating supplementary internally...")
SOURCE_PRIORITY = {
    "IMPPAT": 100, "Phyto4Health": 95, "ANPDB": 90, "SANCDB": 90,
    "AfroDb": 90, "CSIRO": 90, "phytochemdb": 85, "MicotoXilico": 85,
    "MeFSAT": 85, "EMNPD": 80, "MIBiG": 80, "TM-MC": 75, "HERB": 75,
    "CMNPD": 70, "NPASS": 50,
}
supp["_priority"] = supp["source_db"].map(SOURCE_PRIORITY).fillna(50)
valid_supp = supp[supp["inchikey"].notna() & (supp["inchikey"].astype(str).str.len() > 10)].copy()
invalid_supp = supp[~supp.index.isin(valid_supp.index)].copy()
valid_supp = valid_supp.sort_values("_priority", ascending=False)

supp_prov = valid_supp.groupby("inchikey").agg({
    "source_db": lambda x: "|".join(sorted(set(x))),
    "region": lambda x: next((r for r in x if r != "global" and pd.notna(r) and str(r) != "nan"), "global"),
    "source_organism": lambda x: next((o for o in x if pd.notna(o) and str(o) != "" and str(o) != "nan"), ""),
}).reset_index()
supp_prov.columns = ["inchikey", "all_src", "best_region", "best_org"]

supp_deduped = valid_supp.drop_duplicates(subset="inchikey", keep="first").copy()
supp_deduped = supp_deduped.merge(supp_prov, on="inchikey", how="left")
supp_deduped["all_sources"] = supp_deduped["all_src"]
supp_deduped["region"] = supp_deduped["best_region"]
supp_deduped["source_organism"] = supp_deduped["best_org"]
supp_deduped = supp_deduped.drop(columns=["_priority","all_src","best_region","best_org"], errors="ignore")
print(f"  Supplementary after internal dedup: {len(supp_deduped):,}")

# ============================================================
# STEP 3: Extract COCONUT compounds
# ============================================================
print("\n[Extracting COCONUT from coconut_prepped.sdf]")
coconut_sdf = os.path.join(DATADIR, "databases", "coconut_prepped.sdf")
coconut_records = []
if os.path.exists(coconut_sdf):
    suppl = Chem.ForwardSDMolSupplier(open(coconut_sdf, "rb"), removeHs=True, sanitize=True)
    for i, mol in enumerate(tqdm(suppl, desc="  COCONUT SDF")):
        if mol is None: continue
        smi = Chem.MolToSmiles(mol)
        if not smi: continue
        ik = get_ik(mol)
        nm = mol.GetProp("original_name") if mol.HasProp("original_name") else ""
        if nm == "n.a.": nm = ""
        src = mol.GetProp("source_database") if mol.HasProp("source_database") else "COCONUT"
        kingdom_raw = mol.GetProp("source_kingdom") if mol.HasProp("source_kingdom") else "plant"
        km = {"bacteria":"bacteria","fungi":"fungi","marine_algae":"marine",
              "plant_supplementary":"plant","plant":"plant"}.get(kingdom_raw, kingdom_raw)
        p, _ = props(smi)
        coconut_records.append({
            "comp_id": f"CNP_{i:07d}", "name": nm, "smiles": smi, "inchi": "",
            "inchikey": ik, "source_db": "COCONUT",
            "all_sources": f"COCONUT|{src}" if src != "COCONUT" else "COCONUT",
            "kingdom": km, "region": "global", "source_organism": "",
            "mw": p.get("mw"), "logp": p.get("logp"), "tpsa": p.get("tpsa"),
            "hba": p.get("hba"), "hbd": p.get("hbd"),
            "n_rings": p.get("n_rings"), "rotatable_bonds": p.get("rotatable_bonds"),
            "license_tier": "CC BY 4.0",
        })
    coconut_df = pd.DataFrame(coconut_records)
    print(f"  COCONUT extracted: {len(coconut_df):,}")
    before = len(coconut_df)
    coconut_df = coconut_df.drop_duplicates(subset="inchikey", keep="first")
    print(f"  COCONUT after internal dedup: {len(coconut_df):,} (removed {before - len(coconut_df):,})")
else:
    print(f"  NOT FOUND: {coconut_sdf}")
    coconut_df = pd.DataFrame(columns=COLS)

# ============================================================
# STEP 4: Cross-dedup supplementary vs COCONUT
# ============================================================
print("\n[Cross-deduplication]")
coconut_iks = set(coconut_df["inchikey"].dropna().unique())
print(f"  COCONUT unique InChIKeys: {len(coconut_iks):,}")

novel_mask = ~supp_deduped["inchikey"].isin(coconut_iks)
novel_supp = supp_deduped[novel_mask].copy()
overlap_supp = supp_deduped[~novel_mask].copy()

print(f"  Novel (not in COCONUT): {len(novel_supp):,}")
print(f"  Overlap with COCONUT: {len(overlap_supp):,}")
print(f"  No InChIKey (kept as-is): {len(invalid_supp):,}")

# Enrich COCONUT with supplementary provenance
print("  Enriching COCONUT with supplementary metadata...")
if len(overlap_supp) > 0:
    overlap_info = overlap_supp.groupby("inchikey").agg({
        "source_db": lambda x: "|".join(sorted(set(x))),
        "region": lambda x: next((r for r in x if r != "global" and pd.notna(r) and str(r) != "nan"), "global"),
        "source_organism": lambda x: next((o for o in x if pd.notna(o) and str(o) != "" and str(o) != "nan"), ""),
    }).reset_index()
    overlap_info.columns = ["inchikey", "supp_sources", "supp_region", "supp_organism"]

    coconut_enriched = coconut_df.merge(overlap_info, on="inchikey", how="left")
    has_supp = coconut_enriched["supp_sources"].notna()
    coconut_enriched.loc[has_supp, "all_sources"] = (
        coconut_enriched.loc[has_supp, "all_sources"] + "|" + coconut_enriched.loc[has_supp, "supp_sources"])
    has_region = has_supp & (coconut_enriched["supp_region"] != "global")
    coconut_enriched.loc[has_region, "region"] = coconut_enriched.loc[has_region, "supp_region"]
    has_org = has_supp & (coconut_enriched["supp_organism"].astype(str) != "")
    coconut_enriched.loc[has_org, "source_organism"] = coconut_enriched.loc[has_org, "supp_organism"]
    coconut_enriched = coconut_enriched.drop(columns=["supp_sources","supp_region","supp_organism"], errors="ignore")
    enriched_count = has_supp.sum()
else:
    coconut_enriched = coconut_df.copy()
    enriched_count = 0

# ============================================================
# STEP 5: Final merge
# ============================================================
print("\n[Final merge]")
invalid_supp = invalid_supp.drop(columns=["_priority"], errors="ignore")
if "all_sources" not in invalid_supp.columns:
    invalid_supp["all_sources"] = invalid_supp["source_db"]

for c in COLS:
    for df in [coconut_enriched, novel_supp, invalid_supp]:
        if c not in df.columns:
            df[c] = ""

final = pd.concat([
    coconut_enriched[COLS],
    novel_supp[COLS],
    invalid_supp[COLS],
], ignore_index=True)

final["comp_id"] = [f"THEO_{i:07d}" for i in range(len(final))]
final.to_csv(FINAL, index=False)

print(f"\n{'='*60}")
print(f"FINAL THEOBROMA DATABASE")
print(f"{'='*60}")
print(f"Total compounds: {len(final):,}")
print(f"  From COCONUT: {len(coconut_enriched):,}")
print(f"  Novel supplementary: {len(novel_supp):,}")
print(f"  No InChIKey: {len(invalid_supp):,}")
print(f"  COCONUT enriched with supp metadata: {enriched_count:,}")
print(f"\nFile: {FINAL} ({os.path.getsize(FINAL)/1024/1024:.1f} MB)")
print(f"\nBy kingdom:")
for k, cnt in final["kingdom"].value_counts().items():
    print(f"  {k}: {cnt:,}")
print(f"\nBy region (top 15):")
for r, cnt in final["region"].value_counts().head(15).items():
    print(f"  {r}: {cnt:,}")
print(f"\nBy license:")
for l, cnt in final["license_tier"].value_counts().items():
    print(f"  {l}: {cnt:,}")
print(f"\nBy source (top 25):")
for s, cnt in final["source_db"].value_counts().head(25).items():
    print(f"  {s}: {cnt:,}")
print(f"\nMulti-source compounds: {(final['all_sources'].astype(str).str.contains('|', regex=False, na=False)).sum():,}")