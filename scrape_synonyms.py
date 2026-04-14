"""PubChem synonym retrieval -- individual GET, named compounds only."""
import pandas as pd
import requests, time, os, warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm

FINAL = "data/theobroma_final.csv"
OUTFILE = "data/compound_synonyms.csv"
MAX_SYNS = 10

print("=== PubChem Synonym Scraper (GET, named only) ===")
df = pd.read_csv(FINAL, low_memory=False, usecols=["inchikey", "name"])

# Only query compounds that have names (these are the ones users would search)
df = df[df["name"].notna() & (df["name"].astype(str).str.strip() != "")]
df = df[df["inchikey"].notna() & (df["inchikey"].astype(str).str.len() > 10)]
df["inchikey"] = df["inchikey"].astype(str).str.strip()
df = df.drop_duplicates(subset="inchikey")
print(f"  Named compounds with InChIKey: {len(df):,}")

# Resume
records = []
done_iks = set()
if os.path.exists(OUTFILE) and os.path.getsize(OUTFILE) > 10:
    try:
        existing = pd.read_csv(OUTFILE)
        if len(existing) > 0:
            records = existing.to_dict("records")
            done_iks = set(existing["inchikey"].unique())
    except:
        pass
print(f"  Already done: {len(done_iks):,}")

remaining = df[~df["inchikey"].isin(done_iks)]["inchikey"].tolist()
print(f"  Remaining: {len(remaining):,}")

session = requests.Session()
errors = 0
not_found = 0

for i, ik in enumerate(tqdm(remaining, desc="PubChem")):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{ik}/synonyms/JSON"
    try:
        r = session.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            syns = data.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
            for s in syns[:MAX_SYNS]:
                s = s.strip()
                if 2 < len(s) < 200:
                    records.append({"inchikey": ik, "synonym": s})
        elif r.status_code == 404:
            not_found += 1
        elif r.status_code in (429, 503):
            tqdm.write(f"  Rate limited, sleeping 30s...")
            time.sleep(30)
        else:
            errors += 1
    except Exception as e:
        errors += 1
        if errors <= 3:
            tqdm.write(f"  Error: {e}")
        time.sleep(2)

    # 5 requests per second
    if i % 10 == 9:
        time.sleep(1.5)

    # Checkpoint every 500 compounds
    if (i + 1) % 500 == 0:
        pd.DataFrame(records).to_csv(OUTFILE, index=False)
        uq = len(set(r["inchikey"] for r in records))
        tqdm.write(f"  [Checkpoint] {uq:,} compounds, {len(records):,} synonyms, {not_found:,} not in PubChem")

pd.DataFrame(records).to_csv(OUTFILE, index=False)
uq = len(set(r["inchikey"] for r in records))
print(f"\nDone. {uq:,} compounds, {len(records):,} synonyms")
print(f"  Not in PubChem: {not_found:,}")
print(f"  Errors: {errors}")