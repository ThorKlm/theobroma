import requests, time, json, psycopg2

DB = "dbname=theobroma"  # peer auth, run as postgres user
results = []
with psycopg2.connect(DB) as conn, conn.cursor() as cur:
    cur.execute("SELECT comp_id, smiles FROM compounds WHERE source_db = 'TIPdb' AND np_class IS NULL")
    rows = cur.fetchall()

print(f"processing {len(rows)} compounds")
for i, (comp_id, smiles) in enumerate(rows):
    try:
        r = requests.get('https://npclassifier.gnps2.org/classify', params={'smiles': smiles}, timeout=15)
        results.append({"comp_id": comp_id, "smiles": smiles, "response": r.json() if r.status_code == 200 else None, "status": r.status_code})
    except Exception as e:
        results.append({"comp_id": comp_id, "smiles": smiles, "response": None, "error": str(e)})
    if (i + 1) % 50 == 0:
        print(f"{i+1} done")
    time.sleep(0.1)

with open('/tmp/tipdb_npclassifier.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"wrote {len(results)} results to /tmp/tipdb_npclassifier.json")
