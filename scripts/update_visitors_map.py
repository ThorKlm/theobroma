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
