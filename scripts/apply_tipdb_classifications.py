import json, psycopg2

with open('/tmp/tipdb_npclassifier.json') as f:
    results = json.load(f)

updates = 0
skipped = 0
with psycopg2.connect("dbname=theobroma") as conn, conn.cursor() as cur:
    for r in results:
        if not r.get('response'):
            skipped += 1
            continue
        resp = r['response']
        comp_id = r['comp_id']
        np_class = ' $ '.join(resp.get('class_results') or []) or None
        np_superclass = ' $ '.join(resp.get('superclass_results') or []) or None
        np_pathway = ' $ '.join(resp.get('pathway_results') or []) or None
        cur.execute("""
            UPDATE compounds
            SET np_class = %s,
                np_superclass = %s,
                np_pathway = %s,
                inferred_class_source = 'npclassifier_direct'
            WHERE comp_id = %s
        """, (np_class, np_superclass, np_pathway, comp_id))
        updates += cur.rowcount
    conn.commit()
print(f"updated {updates} rows, skipped {skipped}")
