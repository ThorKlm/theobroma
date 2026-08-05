"""Aggregate compounds.region into a country-level JSON for the browse-page map.
Each macro-region (East Asia, South Asia, Europe, ...) expands to a list of ISO-3
country codes. Each country gets assigned the count of its region.
"""
import json, os, subprocess

OUT_DIR = "/home/thorben.klamt/theobroma/static"
OUT_DATA = os.path.join(OUT_DIR, "compounds_by_country.json")

# Macro-region -> ISO-3 country codes
REGION_TO_ISO3 = {
    "East Asia":      ["CHN","JPN","KOR","PRK","TWN","MNG","HKG","MAC"],
    "South Asia":     ["IND","PAK","BGD","NPL","LKA","BTN","MDV","AFG"],
    "Southeast Asia": ["IDN","THA","VNM","MYS","PHL","SGP","KHM","LAO","MMR","BRN","TLS"],
    "Europe":         ["DEU","FRA","GBR","ITA","ESP","NLD","BEL","CHE","AUT","SWE","NOR","DNK","FIN","ISL","PRT","GRC","POL","CZE","SVK","HUN","ROU","BGR","HRV","SVN","SRB","BIH","MNE","MKD","ALB","KOS","UKR","BLR","MDA","LTU","LVA","EST","IRL","LUX","MLT","CYP","GRL"],
    "North America":  ["USA","CAN","MEX"],
    "Latin America":  ["BRA","ARG","COL","CHL","PER","VEN","ECU","BOL","PRY","URY","GUY","SUR","GUF","CUB","DOM","HTI","JAM","TTO","PRI","BHS","BRB","BLZ","CRI","GTM","HND","NIC","PAN","SLV"],
    "Africa":         ["NGA","ETH","EGY","COD","TZA","KEN","ZAF","UGA","DZA","SDN","MAR","AGO","MOZ","GHA","MDG","CMR","CIV","NER","BFA","MLI","MWI","ZMB","SEN","SOM","TCD","ZWE","GIN","RWA","BEN","TUN","BDI","SSD","TGO","SLE","LBY","CAF","COG","LBR","MRT","ERI","NAM","BWA","GMB","LSO","GNB","GNQ","MUS","DJI","COM","CPV","STP","SYC"],
    "Middle East":    ["SAU","IRN","IRQ","ISR","JOR","ARE","KWT","LBN","OMN","QAT","SYR","YEM","BHR","PSE","TUR"],
    "Russia/CIS":     ["RUS","KAZ","UZB","TKM","KGZ","TJK","ARM","AZE","GEO"],
    "Central Asia":   ["KAZ","UZB","TKM","KGZ","TJK"],
    "Australia":      ["AUS"],
    "New Zealand":    ["NZL"],
    "Oceania":        ["PNG","FJI","SLB","VUT","NCL","PYF","WSM","TON","KIR","MHL","FSM","PLW","NRU","TUV","COK","ASM","GUM","MNP"],
}

q = "SELECT macro_region, COUNT(DISTINCT comp_id) FROM compound_region_map GROUP BY macro_region ORDER BY 2 DESC"
r = subprocess.run(['psql','-d','theobroma','-U','theobroma','-h','localhost','-At','-F','\t','-c',q],
                   capture_output=True, text=True)

region_counts = {}
for line in r.stdout.strip().split("\n"):
    if not line: continue
    parts = line.split("\t")
    if len(parts) >= 2:
        region_counts[parts[0]] = int(parts[1])

print(f"regions found: {len(region_counts)}")
for region, cnt in sorted(region_counts.items(), key=lambda x: -x[1]):
    print(f"  {region}: {cnt:,}")

# Build per-country aggregate. If multiple regions overlap on a country
# (e.g., Russia/CIS and Central Asia both contain KAZ, UZB, ...), sum them.
country_counts = {}
country_region = {}  # ISO3 -> primary region (the largest contributor)
for region, cnt in region_counts.items():
    iso3_list = REGION_TO_ISO3.get(region, [])
    if not iso3_list:
        print(f"  WARN: no ISO3 mapping for region '{region}'")
        continue
    for iso3 in iso3_list:
        # Dominant-region (max) not sum: a country shared by >1 macro-region
        # (e.g. Central Asian states in both Russia/CIS and Central Asia) takes
        # its strongest region's count, matching the live world-map colouring.
        if cnt > country_counts.get(iso3, 0):
            country_counts[iso3] = cnt
            country_region[iso3] = region

out = {
    "country_counts": [{"iso3": k, "count": v, "region": country_region[k]} for k, v in sorted(country_counts.items(), key=lambda x: -x[1])],
    "region_counts": [{"region": k, "count": v} for k, v in sorted(region_counts.items(), key=lambda x: -x[1])],
    "total_compounds_with_region": int(subprocess.run(
        ['psql','-d','theobroma','-U','theobroma','-h','localhost','-At','-c',
         'SELECT COUNT(DISTINCT comp_id) FROM compound_region_map'],
        capture_output=True, text=True).stdout.strip() or 0),
}

with open(OUT_DATA, "w") as f:
    json.dump(out, f, indent=2)

print(f"wrote {OUT_DATA}: {len(country_counts)} countries, {sum(region_counts.values()):,} compounds total")
