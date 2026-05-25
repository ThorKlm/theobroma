"""Rebuild resolved_taxonomy from compound_taxonomy with majority-vote kingdom
selection and post-2022 NCBI kingdom name mapping. Run as theobroma OS user;
uses sudo -u postgres for DB ops.
"""
import subprocess

# Updated kingdom mapping including post-2022 NCBI clade-level names.
KINGDOM_MAP = {
    'Viridiplantae': 'plant', 'Plantae': 'plant',
    'Fungi': 'fungi',
    'Bacteria': 'bacteria', 'Bacillati': 'bacteria', 'Pseudomonadati': 'bacteria',
    'Thermotogati': 'bacteria', 'Methanobacteriati': 'bacteria',
    'Thermoproteati': 'bacteria', 'Archaea': 'bacteria',
    'Metazoa': 'animal', 'Animalia': 'animal',
    'Chromista': 'animal', 'Protozoa': 'animal',
}

def run_sql(sql):
    r = subprocess.run(['sudo','-u','postgres','psql','-d','theobroma','-c',sql],
                       capture_output=True, text=True, check=True)
    return r.stdout

# 1. Temp table holding the mapping, used in CTE join
print('# Setting up kingdom mapping...', flush=True)
mapping_values = ',\n  '.join(f"('{nm}','{km}')" for nm, km in KINGDOM_MAP.items())
sql_setup = f"""
DROP TABLE IF EXISTS _kingdom_mapping_v2;
CREATE TABLE _kingdom_mapping_v2 (ncbi_name text PRIMARY KEY, theobroma_kingdom text);
INSERT INTO _kingdom_mapping_v2 VALUES
  {mapping_values};
"""
run_sql(sql_setup)
print('# kingdom mapping ready')

# 2. Build resolved_taxonomy fresh. Strategy:
#    - For each compound, gather (mapped_kingdom, specificity_score, family, lineage_fields)
#      from compound_taxonomy. Specificity = count of non-null ranks.
#    - Vote: kingdom with the most compound_taxonomy rows wins.
#    - Tie-break: highest specificity, then alphabetic family.
#    - secondary_kingdoms = all other distinct kingdoms with at least one vote.
#    - Lineage taken from the most-specific row of the winning kingdom.
#    - If no kingdom mappable (empty NCBI/WoRMS), kingdom = compounds.kingdom (legacy)
#      or 'unresolved' for retired pseudo-kingdoms (marine, multi).
print('# building new resolved_taxonomy...', flush=True)
sql_build = """
BEGIN;

DROP TABLE IF EXISTS resolved_taxonomy_new;
CREATE TABLE resolved_taxonomy_new (
    comp_id text PRIMARY KEY,
    kingdom text,
    secondary_kingdoms text[],
    phylum text,
    taxclass text,
    taxorder text,
    family text,
    genus text
);

-- Step 1: map every compound_taxonomy row to a theobroma_kingdom and capture
-- its lineage (preferring WCVP's genus/family columns + WoRMS/NCBI lineage JSONB).
WITH mapped_rows AS (
    SELECT
        ct.comp_id,
        COALESCE(km1.theobroma_kingdom, km2.theobroma_kingdom) AS mapped_kingdom,
        COALESCE(ct.ncbi_lineage->>'phylum', ct.worms_lineage->>'phylum') AS phylum,
        COALESCE(ct.ncbi_lineage->>'class', ct.worms_lineage->>'class') AS taxclass,
        COALESCE(ct.ncbi_lineage->>'order', ct.worms_lineage->>'order') AS taxorder,
        COALESCE(ct.family, ct.ncbi_lineage->>'family', ct.worms_lineage->>'family') AS family,
        COALESCE(ct.genus, ct.ncbi_lineage->>'genus', ct.worms_lineage->>'genus') AS genus,
        -- specificity = count of non-null ranks among the 5
        (CASE WHEN COALESCE(ct.ncbi_lineage->>'phylum', ct.worms_lineage->>'phylum') IS NOT NULL THEN 1 ELSE 0 END +
         CASE WHEN COALESCE(ct.ncbi_lineage->>'class', ct.worms_lineage->>'class') IS NOT NULL THEN 1 ELSE 0 END +
         CASE WHEN COALESCE(ct.ncbi_lineage->>'order', ct.worms_lineage->>'order') IS NOT NULL THEN 1 ELSE 0 END +
         CASE WHEN COALESCE(ct.family, ct.ncbi_lineage->>'family', ct.worms_lineage->>'family') IS NOT NULL THEN 1 ELSE 0 END +
         CASE WHEN COALESCE(ct.genus, ct.ncbi_lineage->>'genus', ct.worms_lineage->>'genus') IS NOT NULL THEN 1 ELSE 0 END
        ) AS specificity
    FROM compound_taxonomy ct
    LEFT JOIN _kingdom_mapping_v2 km1 ON km1.ncbi_name = ct.ncbi_lineage->>'kingdom'
    LEFT JOIN _kingdom_mapping_v2 km2 ON km2.ncbi_name = ct.worms_lineage->>'kingdom'
    WHERE COALESCE(km1.theobroma_kingdom, km2.theobroma_kingdom) IS NOT NULL
),
-- Step 2: per compound, vote on kingdom (most rows wins)
kingdom_votes AS (
    SELECT comp_id, mapped_kingdom, COUNT(*) AS votes
    FROM mapped_rows
    GROUP BY 1, 2
),
ranked_kingdoms AS (
    SELECT comp_id, mapped_kingdom, votes,
           ROW_NUMBER() OVER (PARTITION BY comp_id ORDER BY votes DESC, mapped_kingdom ASC) AS rk
    FROM kingdom_votes
),
primary_kingdoms AS (
    SELECT comp_id, mapped_kingdom AS primary_kingdom
    FROM ranked_kingdoms WHERE rk = 1
),
secondary_kingdoms_per_compound AS (
    SELECT comp_id, ARRAY_AGG(mapped_kingdom ORDER BY votes DESC, mapped_kingdom) AS secondaries
    FROM ranked_kingdoms WHERE rk > 1
    GROUP BY comp_id
),
-- Step 3: for each compound's primary kingdom, pick the most-specific lineage row
best_lineage_rows AS (
    SELECT mr.comp_id,
           mr.phylum, mr.taxclass, mr.taxorder, mr.family, mr.genus,
           pk.primary_kingdom,
           ROW_NUMBER() OVER (
             PARTITION BY mr.comp_id
             ORDER BY mr.specificity DESC, mr.family ASC NULLS LAST
           ) AS lr
    FROM mapped_rows mr
    JOIN primary_kingdoms pk ON pk.comp_id = mr.comp_id AND pk.primary_kingdom = mr.mapped_kingdom
),
final_resolved AS (
    SELECT bl.comp_id,
           bl.primary_kingdom AS kingdom,
           COALESCE(sk.secondaries, ARRAY[]::text[]) AS secondary_kingdoms,
           bl.phylum, bl.taxclass, bl.taxorder, bl.family, bl.genus
    FROM best_lineage_rows bl
    LEFT JOIN secondary_kingdoms_per_compound sk ON sk.comp_id = bl.comp_id
    WHERE bl.lr = 1
)
INSERT INTO resolved_taxonomy_new (comp_id, kingdom, secondary_kingdoms, phylum, taxclass, taxorder, family, genus)
SELECT comp_id, kingdom, secondary_kingdoms, phylum, taxclass, taxorder, family, genus
FROM final_resolved;

-- Step 4: fill in compounds with no mappable compound_taxonomy row.
-- Fall back to compounds.kingdom (legacy) when it maps cleanly; otherwise 'unresolved'.
INSERT INTO resolved_taxonomy_new (comp_id, kingdom, secondary_kingdoms, phylum, taxclass, taxorder, family, genus)
SELECT c.comp_id,
       CASE
           WHEN c.kingdom IN ('plant','fungi','bacteria','animal') THEN c.kingdom
           ELSE 'unresolved'
       END AS kingdom,
       ARRAY[]::text[] AS secondary_kingdoms,
       NULL, NULL, NULL, NULL, NULL
FROM compounds c
WHERE c.comp_id NOT IN (SELECT comp_id FROM resolved_taxonomy_new);

-- Step 5: swap. Keep the buggy table for diagnostics; replace authoritative table.
DROP TABLE IF EXISTS resolved_taxonomy_old_bad;
ALTER TABLE resolved_taxonomy RENAME TO resolved_taxonomy_old_bad;
ALTER TABLE resolved_taxonomy_new RENAME TO resolved_taxonomy;

-- Step 6: recreate indexes (the original had per-rank lowercase indexes)
CREATE INDEX IF NOT EXISTS rt_kingdom_lower ON resolved_taxonomy (LOWER(kingdom));
CREATE INDEX IF NOT EXISTS rt_phylum_lower ON resolved_taxonomy (LOWER(phylum));
CREATE INDEX IF NOT EXISTS rt_taxclass_lower ON resolved_taxonomy (LOWER(taxclass));
CREATE INDEX IF NOT EXISTS rt_taxorder_lower ON resolved_taxonomy (LOWER(taxorder));
CREATE INDEX IF NOT EXISTS rt_family_lower ON resolved_taxonomy (LOWER(family));
CREATE INDEX IF NOT EXISTS rt_genus_lower ON resolved_taxonomy (LOWER(genus));
CREATE INDEX IF NOT EXISTS rt_secondary_kingdoms ON resolved_taxonomy USING GIN (secondary_kingdoms);

GRANT SELECT ON resolved_taxonomy TO theobroma;

COMMIT;

-- Verification
SELECT 'Distribution after rebuild:' AS msg;
SELECT kingdom, COUNT(*) FROM resolved_taxonomy GROUP BY kingdom ORDER BY 2 DESC;

SELECT 'Misclassified microbial-family count (should be 0):' AS msg;
SELECT COUNT(*) FROM resolved_taxonomy rt
JOIN compound_taxonomy ct ON ct.comp_id = rt.comp_id
WHERE rt.kingdom = 'plant'
  AND ct.family IN ('Streptomycetaceae','Mycobacteriaceae','Pseudomonadaceae',
                    'Bacillaceae','Enterobacteriaceae','Vibrionaceae');

SELECT 'Plant-with-Streptophyta-phylum bacteria check (should be 0):' AS msg;
SELECT COUNT(*) FROM resolved_taxonomy
WHERE kingdom = 'bacteria' AND phylum IN ('Streptophyta','Chlorophyta','Anthophyta','Marchantiophyta');

SELECT 'Asteraceae+Middle East primary distribution check:' AS msg;
SELECT rt.kingdom, COUNT(DISTINCT c.comp_id) AS n
FROM compounds c
JOIN compound_taxonomy ct ON ct.comp_id = c.comp_id
JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id
WHERE LOWER(ct.family) = 'asteraceae' AND LOWER(c.region) = 'middle east'
GROUP BY 1 ORDER BY 2 DESC;
"""
print('# (this will take 30-60s)', flush=True)
print(run_sql(sql_build))
print('# done')
