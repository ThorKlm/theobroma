"""Batch-query GBIF for native country ranges of CMAUP plants. 5% checkpoints."""
import pandas as pd
import requests, time, os, warnings
warnings.filterwarnings("ignore")
from tqdm import tqdm
from collections import Counter

plants = pd.read_csv("additional_data/CMAUPv2.0_download_Plants.txt", sep="\t")
print(f"Total plants: {len(plants):,}")

OUTFILE = "additional_data/cmaup_plant_countries.csv"
SESSION = requests.Session()

COUNTRY_TO_REGION = {
    "CN": "East Asia", "JP": "East Asia", "KR": "East Asia", "TW": "East Asia", "MN": "East Asia",
    "IN": "South Asia", "LK": "South Asia", "NP": "South Asia", "BD": "South Asia", "PK": "South Asia",
    "TH": "Southeast Asia", "VN": "Southeast Asia", "ID": "Southeast Asia", "MY": "Southeast Asia",
    "PH": "Southeast Asia", "MM": "Southeast Asia", "KH": "Southeast Asia", "LA": "Southeast Asia",
    "SG": "Southeast Asia",
    "IR": "Central Asia", "AF": "Central Asia", "UZ": "Central Asia", "KZ": "Central Asia",
    "TR": "Central Asia", "GE": "Central Asia", "AZ": "Central Asia",
    "AU": "Australia", "NZ": "Oceania", "FJ": "Oceania", "PG": "Oceania",
    "BR": "Latin America", "MX": "Latin America", "AR": "Latin America", "CO": "Latin America",
    "PE": "Latin America", "CL": "Latin America", "VE": "Latin America", "EC": "Latin America",
    "BO": "Latin America", "PY": "Latin America", "UY": "Latin America", "CR": "Latin America",
    "PA": "Latin America", "CU": "Latin America", "GT": "Latin America", "HN": "Latin America",
    "NI": "Latin America", "SV": "Latin America",
    "US": "North America", "CA": "North America",
    "NG": "Africa", "ZA": "Africa", "KE": "Africa", "ET": "Africa", "TZ": "Africa",
    "GH": "Africa", "CM": "Africa", "SN": "Africa", "MG": "Africa", "MZ": "Africa",
    "UG": "Africa", "CI": "Africa", "AO": "Africa", "ML": "Africa", "NE": "Africa",
    "BF": "Africa", "ZW": "Africa", "MW": "Africa", "RW": "Africa", "TD": "Africa",
    "SD": "Africa", "CD": "Africa", "CG": "Africa", "GA": "Africa", "BJ": "Africa",
    "TG": "Africa", "GN": "Africa", "SL": "Africa", "LR": "Africa", "NA_": "Africa",
    "BW": "Africa", "ZM": "Africa", "ER": "Africa", "DJ": "Africa", "SO": "Africa",
    "EG": "Africa", "MA": "Africa", "DZ": "Africa", "TN": "Africa", "LY": "Africa",
    "DE": "Europe", "FR": "Europe", "GB": "Europe", "IT": "Europe", "ES": "Europe",
    "PT": "Europe", "GR": "Europe", "NL": "Europe", "BE": "Europe", "AT": "Europe",
    "CH": "Europe", "SE": "Europe", "NO": "Europe", "DK": "Europe", "FI": "Europe",
    "PL": "Europe", "CZ": "Europe", "RO": "Europe", "HU": "Europe", "BG": "Europe",
    "HR": "Europe", "RS": "Europe", "SK": "Europe", "SI": "Europe", "IE": "Europe",
    "RU": "Europe", "UA": "Europe", "BY": "Europe",
    "SA": "Middle East", "AE": "Middle East", "IL": "Middle East", "JO": "Middle East",
    "LB": "Middle East", "IQ": "Middle East", "OM": "Middle East", "YE": "Middle East",
    "KW": "Middle East", "QA": "Middle East", "BH": "Middle East", "SY": "Middle East",
}

def get_gbif_countries(species_name):
    try:
        r = SESSION.get("https://api.gbif.org/v1/species/match",
                        params={"name": species_name, "strict": False}, timeout=8)
        data = r.json()
        if data.get("matchType") == "NONE" or "usageKey" not in data:
            return []
        key = data["usageKey"]
        r2 = SESSION.get("https://api.gbif.org/v1/occurrence/search",
                         params={"taxonKey": key, "limit": 0, "facet": "country", "facetLimit": 10},
                         timeout=8)
        facets = r2.json().get("facets", [])
        for f in facets:
            if f["field"] == "COUNTRY":
                return [c["name"] for c in f["counts"][:10]]
        return []
    except:
        return []

def countries_to_region(countries):
    regions = []
    for c in countries:
        r = COUNTRY_TO_REGION.get(c, "")
        if r:
            regions.append(r)
    if not regions:
        return "unresolved"
    return Counter(regions).most_common(1)[0][0]

# Resume from checkpoint
start_idx = 0
results = []
if os.path.exists(OUTFILE):
    existing = pd.read_csv(OUTFILE)
    start_idx = len(existing)
    results = existing.to_dict("records")
    print(f"Resuming from {start_idx}")

checkpoint_interval = max(1, len(plants) // 20)  # 5% checkpoints

for i in tqdm(range(start_idx, len(plants)), desc="GBIF lookup", initial=start_idx, total=len(plants)):
    row = plants.iloc[i]
    pid = str(row["Plant_ID"])
    species = str(row.get("Species_Name", ""))
    if species == "nan" or not species:
        species = str(row.get("Plant_Name", ""))
    if species == "nan":
        results.append({"plant_id": pid, "species": "", "countries": "", "region": "unresolved"})
        continue
    countries = get_gbif_countries(species)
    region = countries_to_region(countries)
    results.append({
        "plant_id": pid,
        "species": species,
        "countries": "; ".join(countries[:5]),
        "region": region,
    })
    # No sleep -- GBIF rate limit is generous (10 req/s documented)
    if (i + 1) % checkpoint_interval == 0:
        pd.DataFrame(results).to_csv(OUTFILE, index=False)
        resolved = sum(1 for r in results if r["region"] != "unresolved")
        tqdm.write(f"  [Checkpoint {100*(i+1)//len(plants)}%] {resolved:,} resolved")

pd.DataFrame(results).to_csv(OUTFILE, index=False)
resolved = sum(1 for r in results if r["region"] != "unresolved")
print(f"\nDone. {resolved:,}/{len(results):,} plants with geographic region")