"""Unify lower-to-parent rank mappings via majority vote (>=75% threshold).
For each lower rank (family, genus, taxorder, taxclass, phylum), find values
that map to multiple parent ranks. If one parent has >=75% of rows, rewrite
the minority rows to use the majority parent.
"""
import subprocess

LOWER_TO_PARENT = [
    ("genus", "family"),
    ("family", "taxorder"),
    ("taxorder", "taxclass"),
    ("taxclass", "phylum"),
    ("phylum", "kingdom"),
]

def run_sql(sql):
    r = subprocess.run(['sudo','-u','postgres','psql','-d','theobroma','-At','-F','\t','-c',sql],
                       capture_output=True, text=True, check=True)
    return r.stdout.strip()

def report_inconsistencies(lower, parent, threshold=0.75):
    """Find lower-rank values where the dominant parent has < threshold share."""
    sql = f"""
    WITH counts AS (
      SELECT {lower} AS lower_val, {parent} AS parent_val, COUNT(*) AS n
      FROM resolved_taxonomy
      WHERE {lower} IS NOT NULL AND {parent} IS NOT NULL
      GROUP BY 1, 2
    ),
    totals AS (
      SELECT lower_val, SUM(n) AS total
      FROM counts GROUP BY 1
    ),
    inconsistencies AS (
      SELECT c.lower_val, c.parent_val, c.n, t.total, c.n::float / t.total AS share
      FROM counts c JOIN totals t ON t.lower_val = c.lower_val
      WHERE t.total > 1
    )
    SELECT lower_val, parent_val, n, total, share
    FROM inconsistencies
    WHERE lower_val IN (
      SELECT lower_val FROM inconsistencies
      GROUP BY 1 HAVING COUNT(*) > 1
    )
    ORDER BY lower_val, n DESC
    """
    return run_sql(sql)

def fix_via_majority(lower, parent, threshold=0.75):
    """For each lower-rank value with a clear majority parent (>= threshold),
    update minority rows to use the majority parent."""
    # Compute majority parent per lower_val
    sql_find = f"""
    WITH counts AS (
      SELECT {lower} AS lower_val, {parent} AS parent_val, COUNT(*) AS n
      FROM resolved_taxonomy
      WHERE {lower} IS NOT NULL AND {parent} IS NOT NULL
      GROUP BY 1, 2
    ),
    totals AS (
      SELECT lower_val, SUM(n) AS total FROM counts GROUP BY 1
    ),
    majorities AS (
      SELECT c.lower_val, c.parent_val AS majority_parent, c.n AS majority_n, t.total
      FROM counts c JOIN totals t ON t.lower_val = c.lower_val
      WHERE c.n::float / t.total >= {threshold}
    )
    SELECT lower_val, majority_parent, majority_n, total
    FROM majorities
    WHERE total > majority_n
    ORDER BY lower_val
    """
    rows = run_sql(sql_find).split('\n')
    if not rows or rows == ['']:
        print(f"  {lower}->{parent}: no clear majorities to apply")
        return 0

    fix_count = 0
    print(f"  {lower}->{parent}: {len(rows)} values with clear majority")
    for line in rows:
        parts = line.split('\t')
        if len(parts) < 4: continue
        lower_val, majority_parent, majority_n, total = parts
        sql_update = f"""
        UPDATE resolved_taxonomy
        SET {parent} = '{majority_parent.replace("'", "''")}'
        WHERE {lower} = '{lower_val.replace("'", "''")}'
          AND ({parent} IS DISTINCT FROM '{majority_parent.replace("'", "''")}')
        """
        result = subprocess.run(['sudo','-u','postgres','psql','-d','theobroma','-c',sql_update],
                               capture_output=True, text=True)
        # parse UPDATE N
        if 'UPDATE' in result.stdout:
            n_updated = int(result.stdout.strip().split()[-1])
            if n_updated > 0:
                fix_count += n_updated
                if n_updated > 10:
                    print(f"    {lower_val}: {n_updated} rows -> {majority_parent} (was {int(total) - int(majority_n)} minority)")
    return fix_count

print("=== Inconsistency cleanup with >=75% majority threshold ===")
for lower, parent in LOWER_TO_PARENT:
    n = fix_via_majority(lower, parent, threshold=0.75)
    print(f"  total {lower}->{parent} fixed: {n}")
print("Done.")
