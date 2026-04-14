"""PubChem synonym scraper -- fixed batch format."""
import pandas as pd
import requests, time, os, warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
import numpy as np

FINAL = "data/theobroma_final.csv"
OUTFILE = "data/compound_synonyms.csv"
MAX_SYNS = 15
BATCH = 10  # smaller batches are more reliable
CHECKPOINT = 1000

SESSION = requests.Session()

print("=== PubChem Synonym Scraper ===")
df = pd.read_csv(FINAL, low_memory=False, usecols=["inchikey"])
df = df[df["inchikey"].notna() & (df["inchikey"].astype(str).str.len() > 10)]
iks = df["inchikey"].astype(str).str.strip().drop_duplicates().tolist()
print(f"  Unique InChIKeys: {len(iks):,}")

# Resume
records = []
done_iks = set()
if os.path.exists(OUTFILE):
    existing = pd.read_csv(OUTFILE)
    if len(existing) > 0:
        records = existing.to_dict("records")
        done_iks = set(existing["inchikey"].unique())
print(f"  Already done: {len(done_iks):,}")

remaining = [ik for ik in iks if ik not in done_iks]
print(f"  Remaining: {len(remaining):,}")

errors = 0
batch_count = 0

for start in tqdm(range(0, len(remaining), BATCH), desc="PubChem"):
    batch = remaining[start:start+BATCH]
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/synonyms/JSON"

    try:
        # POST with newline-separated keys
        payload = "inchikey=" + "\n".join(batch)
        r = SESSION.post(url, data=payload, timeout=30,
                         headers={"Content-Type": "application/x-www-form-urlencoded"})

        if r.status_code == 200:
            data = r.json()
            for info in data.get("InformationList", {}).get("Information", []):
                ik = info.get("InChIKey", "")
                syns = info.get("Synonym", [])
                for s in syns[:MAX_SYNS]:
                    s = s.strip()
                    if 2 < len(s) < 200:
                        records.append({"inchikey": ik, "synonym": s})
        elif r.status_code == 404:
            pass  # none found
        elif r.status_code in (429, 503):
            time.sleep(10)
        else:
            errors += 1

    except Exception as e:
        errors += 1
        time.sleep(2)

    add_delay = np.random.random()*0.1
    time.sleep(0.15+add_delay)  # ~4 req/s
    batch_count += 1

    if batch_count % CHECKPOINT == 0:
        pd.DataFrame(records).to_csv(OUTFILE, index=False)
        uq = len(set(r["inchikey"] for r in records))
        tqdm.write(f"  [Checkpoint] {uq:,} compounds, {len(records):,} synonyms, {errors} errors")

pd.DataFrame(records).to_csv(OUTFILE, index=False)
uq = len(set(r["inchikey"] for r in records))
print(f"\nDone. {uq:,} compounds, {len(records):,} synonyms, {errors} errors")