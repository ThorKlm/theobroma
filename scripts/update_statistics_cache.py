"""Build the /statistics page cache.
Runs all aggregation queries on the compounds and admet tables, writes
the result to /home/thorben.klamt/theobroma/static/statistics_cache.json.
"""
import json, os, time, subprocess

OUT = "/home/thorben.klamt/theobroma/static/statistics_cache.json"

def run_sql(q):
    r = subprocess.run(
        ['sudo','-u','postgres','psql','-d','theobroma','-At','-F','\t','-c', q],
        capture_output=True, text=True, check=True
    )
    return [line.split('\t') for line in r.stdout.strip().split('\n') if line]

# Hand-written REGION_SQL extract from app.py. Mirrors the canonical region grouping.
REGION_SQL = """CASE
    WHEN region IS NULL OR region='' THEN 'global / unresolved'
    WHEN LOWER(region) IN ('global','world','cosmopolitan','worldwide') THEN 'global / unresolved'
    ELSE region
END"""

print("# building statistics cache...")
t0 = time.time()
cache = {}

total = run_sql("SELECT COUNT(*) FROM compounds")
cache['total'] = int(total[0][0])

cache['kingdoms'] = [
    {'kingdom': r[0], 'cnt': int(r[1])}
    for r in run_sql("SELECT kingdom, COUNT(*) FROM compounds WHERE kingdom IS NOT NULL AND kingdom!='' GROUP BY kingdom ORDER BY 2 DESC")
]

cache['sources'] = [
    {'source_db': r[0], 'cnt': int(r[1])}
    for r in run_sql("SELECT source_db, COUNT(*) FROM compounds GROUP BY source_db ORDER BY 2 DESC")
]

cache['regions'] = [
    {'region': r[0], 'cnt': int(r[1])}
    for r in run_sql(f"SELECT {REGION_SQL} AS region, COUNT(*) FROM compounds GROUP BY 1 ORDER BY 2 DESC LIMIT 30")
]

prop_rows = run_sql("SELECT AVG(mw), AVG(logp), AVG(tpsa), AVG(hba), AVG(hbd) FROM compounds WHERE mw IS NOT NULL")
cache['prop_stats'] = {
    'avg_mw': float(prop_rows[0][0]) if prop_rows[0][0] else None,
    'avg_logp': float(prop_rows[0][1]) if prop_rows[0][1] else None,
    'avg_tpsa': float(prop_rows[0][2]) if prop_rows[0][2] else None,
    'avg_hba': float(prop_rows[0][3]) if prop_rows[0][3] else None,
    'avg_hbd': float(prop_rows[0][4]) if prop_rows[0][4] else None,
}

cache['licenses'] = [
    {'license_tier': r[0], 'cnt': int(r[1])}
    for r in run_sql("SELECT license_tier, COUNT(*) FROM compounds GROUP BY license_tier ORDER BY 2 DESC")
]

multi = run_sql("SELECT COUNT(*) FROM compounds WHERE all_sources LIKE '%|%'")
cache['multi_source'] = int(multi[0][0])

admet_keys = [
    ("hERG_Karim-et-al", "hERG risk"),
    ("AMES_Li-et-al", "AMES mutagenicity"),
    ("BBB_Martins-et-al", "BBB penetration"),
    ("HIA_Hou-et-al", "Intestinal absorption"),
    ("Caco2_Wang-et-al", "Caco-2 perm."),
    ("Solubility_AqSolDB", "Solubility"),
    ("Lipophilicity_AstraZeneca", "Lipophilicity"),
    ("CYP2D6_Inhibition_Veith-et-al", "CYP2D6 inhib."),
    ("CYP3A4_Inhibition_Veith-et-al", "CYP3A4 inhib."),
    ("CYP2C9_Inhibition_Veith-et-al", "CYP2C9 inhib."),
    ("ClinTox_CT_TOX", "Clinical tox."),
    ("DILI", "DILI risk"),
]
cache['admet_stats'] = {}
for col, label in admet_keys:
    try:
        r = run_sql(f'SELECT AVG("{col}"), MIN("{col}"), MAX("{col}"), COUNT(*) FROM admet WHERE "{col}" IS NOT NULL')
        if r and r[0][0]:
            cache['admet_stats'][col] = {
                'label': label,
                'avg': float(r[0][0]),
                'min': float(r[0][1]),
                'max': float(r[0][2]),
                'n': int(r[0][3]),
            }
    except Exception as e:
        pass

cache['_meta'] = {
    'generated_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    'generation_seconds': time.time() - t0,
}

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w') as f:
    json.dump(cache, f, indent=2, default=str)
print(f"# wrote {OUT} ({os.path.getsize(OUT)/1024:.1f} KB)")
print(f"# generation time: {time.time()-t0:.1f}s")
