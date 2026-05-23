"""Aggregate access_log IPs into visitors-by-country JSON for /statistics.
Resolves new IPs via ipinfo.io free tier (no key, IPv4 + IPv6, 50k/mo limit).
Caches resolutions in /home/thorben.klamt/.visitor_country_cache.json so re-runs
do not re-query already-known IPs.
"""
import json, os, time, sys, urllib.request, urllib.error, subprocess

CACHE = "/home/thorben.klamt/.visitor_country_cache.json"
OUT_DIR = "/home/thorben.klamt/theobroma/static"
OUT_DATA = os.path.join(OUT_DIR, "visitors_by_country.json")
OUT_META = os.path.join(OUT_DIR, "visitors_meta.json")

# UAs known to be bots, blocked at robots.txt + middleware level.
BOT_UA_SUBSTRINGS = [
    "GPTBot", "CCBot", "ClaudeBot", "anthropic-ai", "PerplexityBot",
    "Bytespider", "AhrefsBot", "SemrushBot", "MJ12bot", "DotBot", "AmazonBot",
]

def load_cache():
    if not os.path.exists(CACHE): return {}
    with open(CACHE) as f: return json.load(f)

def save_cache(d):
    with open(CACHE, "w") as f: json.dump(d, f, indent=2)

def resolve_ip(ip):
    """Returns (country_code, country_name) or (None, None) on failure."""
    try:
        req = urllib.request.Request(f"https://ipinfo.io/{ip}/json", headers={"Accept":"application/json"})
        with urllib.request.urlopen(req, timeout=8) as r:
            data = json.loads(r.read().decode("utf-8"))
            return data.get("country"), data.get("region")
    except Exception:
        return None, None

# Fetch distinct non-bot IPs from access_log
def fetch_ips():
    bot_pattern = "|".join(BOT_UA_SUBSTRINGS)
    q = f"""SELECT ip, COUNT(*) AS hits FROM access_log
            WHERE ip IS NOT NULL AND ip != ''
              AND COALESCE(user_agent, '') !~* '({bot_pattern})'
              AND path NOT LIKE '/static/%'
              AND path NOT LIKE '/api/depict%'
            GROUP BY ip ORDER BY 2 DESC"""
    r = subprocess.run(['sudo','-u','postgres','psql','-d','theobroma','-At','-F','\t','-c',q],
                       capture_output=True, text=True)
    rows = []
    for line in r.stdout.strip().split("\n"):
        if not line: continue
        parts = line.split("\t")
        if len(parts) >= 2:
            rows.append((parts[0], int(parts[1])))
    return rows

print("fetching IPs from access_log...")
ip_rows = fetch_ips()
print(f"  distinct non-bot IPs: {len(ip_rows):,}")

cache = load_cache()
print(f"  cache hits so far: {sum(1 for ip,_ in ip_rows if ip in cache):,}")

needs_lookup = [ip for ip, _ in ip_rows if ip not in cache]
print(f"  needs lookup: {len(needs_lookup):,}")

for i, ip in enumerate(needs_lookup):
    cc, region = resolve_ip(ip)
    cache[ip] = {"cc": cc, "region": region, "looked_up_at": int(time.time())}
    if (i+1) % 25 == 0:
        print(f"    {i+1}/{len(needs_lookup)} resolved")
        save_cache(cache)
    time.sleep(0.2)  # polite to ipinfo

save_cache(cache)
print(f"  cache size: {len(cache):,}")

# Aggregate
country_counts = {}
country_hits = {}
unresolved = 0
for ip, hits in ip_rows:
    entry = cache.get(ip)
    if not entry or not entry.get("cc"):
        unresolved += 1
        continue
    cc = entry["cc"]
    country_counts[cc] = country_counts.get(cc, 0) + 1
    country_hits[cc] = country_hits.get(cc, 0) + hits

# Convert ISO-2 (ipinfo) to ISO-3 (Plotly choropleth standard)
ISO2_TO_3 = {
    "AD":"AND","AE":"ARE","AF":"AFG","AG":"ATG","AI":"AIA","AL":"ALB","AM":"ARM","AO":"AGO",
    "AR":"ARG","AT":"AUT","AU":"AUS","AZ":"AZE","BA":"BIH","BB":"BRB","BD":"BGD","BE":"BEL",
    "BF":"BFA","BG":"BGR","BH":"BHR","BI":"BDI","BJ":"BEN","BN":"BRN","BO":"BOL","BR":"BRA",
    "BS":"BHS","BT":"BTN","BW":"BWA","BY":"BLR","BZ":"BLZ","CA":"CAN","CD":"COD","CF":"CAF",
    "CG":"COG","CH":"CHE","CI":"CIV","CL":"CHL","CM":"CMR","CN":"CHN","CO":"COL","CR":"CRI",
    "CU":"CUB","CV":"CPV","CY":"CYP","CZ":"CZE","DE":"DEU","DJ":"DJI","DK":"DNK","DM":"DMA",
    "DO":"DOM","DZ":"DZA","EC":"ECU","EE":"EST","EG":"EGY","ER":"ERI","ES":"ESP","ET":"ETH",
    "FI":"FIN","FJ":"FJI","FM":"FSM","FR":"FRA","GA":"GAB","GB":"GBR","GD":"GRD","GE":"GEO",
    "GH":"GHA","GM":"GMB","GN":"GIN","GQ":"GNQ","GR":"GRC","GT":"GTM","GW":"GNB","GY":"GUY",
    "HK":"HKG","HN":"HND","HR":"HRV","HT":"HTI","HU":"HUN","ID":"IDN","IE":"IRL","IL":"ISR",
    "IN":"IND","IQ":"IRQ","IR":"IRN","IS":"ISL","IT":"ITA","JM":"JAM","JO":"JOR","JP":"JPN",
    "KE":"KEN","KG":"KGZ","KH":"KHM","KM":"COM","KN":"KNA","KP":"PRK","KR":"KOR","KW":"KWT",
    "KZ":"KAZ","LA":"LAO","LB":"LBN","LC":"LCA","LI":"LIE","LK":"LKA","LR":"LBR","LS":"LSO",
    "LT":"LTU","LU":"LUX","LV":"LVA","LY":"LBY","MA":"MAR","MC":"MCO","MD":"MDA","ME":"MNE",
    "MG":"MDG","MH":"MHL","MK":"MKD","ML":"MLI","MM":"MMR","MN":"MNG","MR":"MRT","MT":"MLT",
    "MU":"MUS","MV":"MDV","MW":"MWI","MX":"MEX","MY":"MYS","MZ":"MOZ","NA":"NAM","NE":"NER",
    "NG":"NGA","NI":"NIC","NL":"NLD","NO":"NOR","NP":"NPL","NR":"NRU","NZ":"NZL","OM":"OMN",
    "PA":"PAN","PE":"PER","PG":"PNG","PH":"PHL","PK":"PAK","PL":"POL","PT":"PRT","PW":"PLW",
    "PY":"PRY","QA":"QAT","RO":"ROU","RS":"SRB","RU":"RUS","RW":"RWA","SA":"SAU","SB":"SLB",
    "SC":"SYC","SD":"SDN","SE":"SWE","SG":"SGP","SI":"SVN","SK":"SVK","SL":"SLE","SM":"SMR",
    "SN":"SEN","SO":"SOM","SR":"SUR","SS":"SSD","ST":"STP","SV":"SLV","SY":"SYR","SZ":"SWZ",
    "TD":"TCD","TG":"TGO","TH":"THA","TJ":"TJK","TL":"TLS","TM":"TKM","TN":"TUN","TO":"TON",
    "TR":"TUR","TT":"TTO","TV":"TUV","TW":"TWN","TZ":"TZA","UA":"UKR","UG":"UGA","US":"USA",
    "UY":"URY","UZ":"UZB","VC":"VCT","VE":"VEN","VN":"VNM","VU":"VUT","YE":"YEM","ZA":"ZAF",
    "ZM":"ZMB","ZW":"ZWE",
}

countries_iso3 = []
for cc2, n in sorted(country_counts.items(), key=lambda x:-x[1]):
    iso3 = ISO2_TO_3.get(cc2, cc2)
    countries_iso3.append({"iso3": iso3, "iso2": cc2, "visitors": n, "hits": country_hits[cc2]})

# Save
os.makedirs(OUT_DIR, exist_ok=True)
with open(OUT_DATA, "w") as f: json.dump(countries_iso3, f, indent=2)
with open(OUT_META, "w") as f: json.dump({
    "last_updated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "distinct_visitor_ips": len(ip_rows),
    "ips_resolved": sum(1 for ip,_ in ip_rows if cache.get(ip,{}).get("cc")),
    "countries": len(country_counts),
    "ips_unresolved": unresolved,
}, f, indent=2)

print(f"\nresolved {len(country_counts):,} countries, {len(ip_rows)-unresolved:,} IPs")
print(f"unresolved: {unresolved}")
print(f"top 10:")
for c in countries_iso3[:10]:
    print(f"  {c['iso2']} ({c['iso3']}): {c['visitors']} visitors, {c['hits']:,} hits")
print(f"\nwrote {OUT_DATA} and {OUT_META}")
