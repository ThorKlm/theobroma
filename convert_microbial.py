"""Fix organism fields in npatlas.csv and cyanometdb.csv."""
import pandas as pd

# NP Atlas: combine genus + origin_species
print("=== Fixing NP Atlas organisms ===")
npa = pd.read_csv("final_additional_data/np_atlas_2024_09.tsv", sep="\t", low_memory=False)
out = pd.read_csv("data/converted/npatlas.csv", low_memory=False)
orgs = (npa["genus"].fillna("") + " " + npa["origin_species"].fillna("")).str.strip()
out["source_organism"] = orgs.values[:len(out)]
out.to_csv("data/converted/npatlas.csv", index=False)
with_org = (out["source_organism"].str.strip() != "").sum()
print(f"  With organism: {with_org:,}/{len(out):,}")

# CyanoMetDB: combine Genus + Species
print("=== Fixing CyanoMetDB organisms ===")
cya = pd.read_csv("final_additional_data/CyanoMetDB_V03_2024.csv", low_memory=False)
out2 = pd.read_csv("data/converted/cyanometdb.csv", low_memory=False)
orgs2 = (cya["Genus"].fillna("") + " " + cya["Species"].fillna("")).str.strip()
out2["source_organism"] = orgs2.values[:len(out2)]
out2.to_csv("data/converted/cyanometdb.csv", index=False)
with_org2 = (out2["source_organism"].str.strip() != "").sum()
print(f"  With organism: {with_org2:,}/{len(out2):,}")