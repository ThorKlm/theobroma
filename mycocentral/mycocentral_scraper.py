"""Scrape MycoCentral: 904 mycotoxins from individual compound pages."""
import requests, os, time, re, warnings, csv
warnings.filterwarnings("ignore")
from bs4 import BeautifulSoup
from tqdm import tqdm
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE = "https://www.mycocentral.eu"
OUTFILE = os.path.join("..", "data", "converted", "mycocentral.csv")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

session = requests.Session()
retry = Retry(total=3, backoff_factor=2)
session.mount("https://", HTTPAdapter(max_retries=retry))
session.verify = False
HEADERS = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"}

# Step 1: Get all compound IDs from browse pages
print("Collecting compound IDs from browse pages...")
comp_ids = []
page = 1
while True:
    try:
        r = session.get(f"{BASE}/mycotoxins/?page={page}", headers=HEADERS, timeout=30)
        if r.status_code != 200:
            break
        soup = BeautifulSoup(r.text, "html.parser")
        links = soup.find_all("a", href=re.compile(r"/mycotoxins/\d+"))
        if not links:
            break
        for link in links:
            m = re.search(r"/mycotoxins/(\d+)", link.get("href", ""))
            if m:
                comp_ids.append(int(m.group(1)))
        print(f"  Page {page}: found {len(links)} links, total {len(comp_ids)}")
        page += 1
        time.sleep(0.3)
    except Exception as e:
        print(f"  Page {page} error: {e}")
        break

comp_ids = sorted(set(comp_ids))
print(f"Total unique compound IDs: {len(comp_ids)}")

# Step 2: Scrape each compound page
COLS = ["comp_id","name","smiles","inchi","inchikey","source_db","kingdom",
        "region","source_organism","mw","logp","tpsa","hba","hbd",
        "n_rings","rotatable_bonds","license_tier"]

records = []
checkpoint_interval = max(1, len(comp_ids) // 10)

for idx, cid in enumerate(tqdm(comp_ids, desc="Scraping compounds")):
    try:
        r = session.get(f"{BASE}/mycotoxins/{cid}", headers=HEADERS, timeout=30)
        if r.status_code != 200:
            continue
        soup = BeautifulSoup(r.text, "html.parser")
        text = soup.get_text("\n", strip=True)

        # Extract fields
        name = ""
        smiles = ""
        inchi_val = ""
        inchikey = ""
        mw = ""
        cas = ""
        fungi_list = ""

        # Name
        for pattern in [r"Mycotoxin name:\s*(.+)", r"<h2>(.+?)</h2>"]:
            m = re.search(pattern, text)
            if m:
                name = m.group(1).strip()
                break

        # SMILES
        m = re.search(r"Smiles:\s*(\S+)", text)
        if m:
            smiles = m.group(1).strip()

        # InChI
        m = re.search(r"Inchi:\s*(InChI=\S+)", text)
        if m:
            inchi_val = m.group(1).strip()

        # InChIKey
        m = re.search(r"Inchikey:\s*(\S+)", text)
        if m:
            inchikey = m.group(1).strip()

        # MW from formula line
        m = re.search(r"(\d+\.\d+)\s*g/mol", text)
        if m:
            mw = m.group(1)

        # Fungi (source organisms) -- look for fungi section
        fungi_section = soup.find("h3", string=re.compile("Fungi", re.I))
        if fungi_section:
            fungi_items = []
            for sib in fungi_section.find_next_siblings():
                if sib.name == "h3":
                    break
                ft = sib.get_text(strip=True)
                if ft and len(ft) > 2:
                    fungi_items.append(ft)
            fungi_list = "; ".join(fungi_items[:5])

        if smiles and len(smiles) > 3:
            records.append({
                "comp_id": f"MYCOC_{idx:06d}",
                "name": name,
                "smiles": smiles,
                "inchi": inchi_val,
                "inchikey": inchikey,
                "source_db": "MycoCentral",
                "kingdom": "fungi",
                "region": "global",
                "source_organism": fungi_list,
                "mw": mw,
                "logp": "", "tpsa": "", "hba": "", "hbd": "",
                "n_rings": "", "rotatable_bonds": "",
                "license_tier": "CC BY 4.0",
            })

        time.sleep(0.3)

        if (idx + 1) % checkpoint_interval == 0:
            tqdm.write(f"  [Checkpoint] {idx+1}/{len(comp_ids)}, {len(records)} compounds scraped")

    except Exception as e:
        tqdm.write(f"  ID {cid} error: {e}")
        time.sleep(1)

# Save
import pandas as pd
df = pd.DataFrame(records)
for c in COLS:
    if c not in df.columns:
        df[c] = ""
df[COLS].to_csv(OUTFILE, index=False)
print(f"\nSaved {len(df)} compounds to {OUTFILE}")