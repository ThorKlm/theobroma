"""Aggregate access_log IPs into visitors-by-country JSON for /statistics.
Resolves new IPs via ipinfo.io free tier (no key, IPv4 + IPv6, 50k/mo limit).
Caches resolutions in /home/thorben.klamt/.visitor_country_cache.json so re-runs
do not re-query already-known IPs.
"""
import json, os, time, sys, urllib.request, urllib.error, psycopg2

CACHE = "/home/thorben.klamt/.visitor_country_cache.json"
OUT_DIR = "/home/thorben.klamt/theobroma/static"
OUT_DATA = os.path.join(OUT_DIR, "visitors_by_country.json")
OUT_META = os.path.join(OUT_DIR, "visitors_meta.json")
DB_URI = "postgresql://theobroma:theobroma@localhost:5432/theobroma"

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

# Aggregate country totals, independent of access_log retention.
with psycopg2.connect(DB_URI) as conn:
    with conn.cursor() as cur:
        cur.execute("""SELECT country, sum(hits)::bigint, sum(tokens)::bigint
                       FROM visitor_country_totals
                       WHERE country IS NOT NULL AND country <> ''
                       GROUP BY 1 ORDER BY 2 DESC""")
        rows = cur.fetchall()
print("countries from visitor_country_totals:", len(rows))

countries_iso3 = []
for cc2, hits, tokens in rows:
    if not cc2:
        continue
    countries_iso3.append({"iso3": ISO2_TO_3.get(cc2, cc2), "iso2": cc2,
                           "visitors": int(tokens), "hits": int(hits)})
countries_iso3.sort(key=lambda c: -c["visitors"])

# Tripwire: refuse to overwrite if the result is empty or implausibly small.
if len(countries_iso3) < 50:
    sys.stderr.write("ERROR: only %d countries; refusing to overwrite\n" % len(countries_iso3))
    sys.exit(1)

os.makedirs(OUT_DIR, exist_ok=True)
with open(OUT_DATA, "w") as f:
    json.dump(countries_iso3, f, indent=2)
with open(OUT_META, "w") as f:
    json.dump({
        "last_updated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source": "visitor_country_totals",
        "countries": len(countries_iso3),
        "visitor_days": sum(c["visitors"] for c in countries_iso3),
        "hits": sum(c["hits"] for c in countries_iso3),
    }, f, indent=2)

print("wrote %s and %s" % (OUT_DATA, OUT_META))
for c in countries_iso3[:10]:
    print("  %s (%s): %s visitor-days, %s hits" % (c["iso2"], c["iso3"], c["visitors"], c["hits"]))
