"""Build the /taxonomy page cache from resolved_taxonomy with both primary and
secondary kingdoms unioned, so cross-kingdom compounds appear in all their kingdoms.
"""
import json, os, time, subprocess

OUT = "/home/thorben.klamt/theobroma/static/taxonomy_cache.json"

def run_sql(q):
    r = subprocess.run(
        ['sudo','-u','postgres','psql','-d','theobroma','-At','-F','\t','-c', q],
        capture_output=True, text=True, check=True
    )
    return [line.split('\t') for line in r.stdout.strip().split('\n') if line]

LINEAGE_SQL = """
WITH base AS (
    SELECT rt.comp_id,
           COALESCE(pm.lineage_kingdom, rt.kingdom) AS primary_k,
           rt.secondary_kingdoms,
           rt.phylum, rt.taxclass, rt.taxorder, rt.family, rt.genus
    FROM resolved_taxonomy rt
    LEFT JOIN phylum_kingdom_map pm ON pm.phylum = rt.phylum
),
expanded AS (
    SELECT comp_id, primary_k AS k, phylum, taxclass, taxorder, family, genus FROM base
    UNION ALL
    SELECT comp_id, unnest(secondary_kingdoms), NULL, NULL, NULL, NULL, NULL
    FROM base WHERE secondary_kingdoms IS NOT NULL AND secondary_kingdoms <> '{}'
)
SELECT k AS theobroma_kingdom, NULL AS lineage_kingdom, phylum, taxclass, taxorder, family, genus,
       COUNT(DISTINCT comp_id) AS n
FROM expanded
GROUP BY 1,2,3,4,5,6,7
"""

print("# building taxonomy cache...")
t0 = time.time()
rows = run_sql(LINEAGE_SQL)
print(f"# fetched {len(rows):,} distinct lineage paths")

def insert(tree, path, count):
    node = tree
    for level in path:
        if level is None or level == '':
            level = '(unresolved)'
        if 'children' not in node:
            node['children'] = {}
        if level not in node['children']:
            node['children'][level] = {'name': level, 'value': 0}
        node = node['children'][level]
        node['value'] += count

root = {'name': 'THEOBROMA', 'value': 0}
for row in rows:
    theobroma_k, _, phylum, klass, order_, family, genus, n = row
    n = int(n)
    root['value'] += n
    insert(root, [theobroma_k, phylum, klass, order_, family, genus], n)

def to_d3(node):
    out = {'name': node['name'], 'value': node['value']}
    if 'children' in node:
        children = [to_d3(c) for c in node['children'].values()]
        children.sort(key=lambda x: -x['value'])
        out['children'] = children
    return out

tree = to_d3(root)

true_total_rows = run_sql("SELECT COUNT(*) FROM resolved_taxonomy")
total_compounds = int(true_total_rows[0][0])

kingdom_totals = run_sql("""
    SELECT kingdom, COUNT(*) FROM resolved_taxonomy
    WHERE kingdom IS NOT NULL AND kingdom != ''
    GROUP BY kingdom ORDER BY 2 DESC
""")
total_resolved = sum(int(r[7]) for r in rows if r[1] or r[2] or r[3] or r[4] or r[5] or r[6])

cache = {
    'tree': tree,
    'kingdom_totals': [{'kingdom': r[0], 'cnt': int(r[1])} for r in kingdom_totals],
    'total_compounds': total_compounds,
    'total_resolved': total_resolved,
    'distinct_paths': len(rows),
    '_meta': {
        'generated_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'generation_seconds': time.time() - t0,
    }
}

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, 'w') as f:
    json.dump(cache, f, separators=(',', ':'))
print(f"# wrote {OUT} ({os.path.getsize(OUT)/1024:.1f} KB)")
print(f"# total compounds: {total_compounds:,}, resolved: {total_resolved:,}")

KINGDOMS = ["plant", "fungi", "bacteria", "animal", "unresolved"]
for k in KINGDOMS:
    tk = time.time()
    krows = [r for r in rows if r[0] == k]
    if not krows:
        continue
    kroot = {'name': 'THEOBROMA', 'value': 0}
    for row in krows:
        theobroma_k, _, phylum, klass, order_, family, genus, n = row
        n = int(n)
        kroot['value'] += n
        insert(kroot, [theobroma_k, phylum, klass, order_, family, genus], n)
    ktree = to_d3(kroot)
    kresolved = sum(int(r[7]) for r in krows if r[1] or r[2] or r[3] or r[4] or r[5] or r[6])
    ktotal = sum(int(r[7]) for r in krows)
    kcache = {
        'tree': ktree, 'kingdom_filter': k,
        'total_compounds': ktotal, 'total_resolved': kresolved,
        'distinct_paths': len(krows),
        '_meta': {
            'generated_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'generation_seconds': time.time() - tk,
        }
    }
    kpath = OUT.replace('.json', f'_{k}.json')
    with open(kpath, 'w') as f:
        json.dump(kcache, f, separators=(',', ':'))
    print(f"# wrote {kpath} ({os.path.getsize(kpath)/1024:.1f} KB, {ktotal:,} compounds)")

for old_k in ["marine", "multi"]:
    p = OUT.replace('.json', f'_{old_k}.json')
    if os.path.exists(p):
        os.remove(p)
        print(f"# removed stale {p}")

print(f"# generation time: {time.time()-t0:.1f}s")
