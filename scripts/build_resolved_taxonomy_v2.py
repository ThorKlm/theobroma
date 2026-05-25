"""Rebuild resolved_taxonomy with majority-vote kingdom + three-tier lineage
fallback. Preserves WCVP-derived family/genus for compounds whose compound_taxonomy
rows lack mappable NCBI/WoRMS kingdom (typical for plants resolved via WCVP only).
"""
import subprocess

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

print('# building new resolved_taxonomy with three-tier lineage fallback...', flush=True)
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

-- Tier 1: compounds with mappable kingdom from NCBI/WoRMS lineage.
-- Majority-vote on mapped kingdom across compound_taxonomy rows.
-- Lineage taken from the most-specific row of the winning kingdom.
WITH mapped_rows AS (
    SELECT
        ct.comp_id,
        COALESCE(km1.theobroma_kingdom, km2.theobroma_kingdom) AS mapped_kingdom,
        COALESCE(ct.ncbi_lineage->>'phylum', ct.worms_lineage->>'phylum') AS phylum,
        COALESCE(ct.ncbi_lineage->>'class', ct.worms_lineage->>'class') AS taxclass,
        COALESCE(ct.ncbi_lineage->>'order', ct.worms_lineage->>'order') AS taxorder,
        COALESCE(ct.family, ct.ncbi_lineage->>'family', ct.worms_lineage->>'family') AS family,
        COALESCE(ct.genus, ct.ncbi_lineage->>'genus', ct.worms_lineage->>'genus') AS genus,
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
kingdom_votes AS (
    SELECT comp_id, mapped_kingdom, COUNT(*) AS votes
    FROM mapped_rows GROUP BY 1, 2
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
    FROM ranked_kingdoms WHERE rk > 1 GROUP BY comp_id
),
best_lineage_rows AS (
    SELECT mr.comp_id, mr.phylum, mr.taxclass, mr.taxorder, mr.family, mr.genus,
           pk.primary_kingdom,
           ROW_NUMBER() OVER (
             PARTITION BY mr.comp_id
             ORDER BY mr.specificity DESC, mr.family ASC NULLS LAST
           ) AS lr
    FROM mapped_rows mr
    JOIN primary_kingdoms pk ON pk.comp_id = mr.comp_id AND pk.primary_kingdom = mr.mapped_kingdom
)
INSERT INTO resolved_taxonomy_new (comp_id, kingdom, secondary_kingdoms, phylum, taxclass, taxorder, family, genus)
SELECT bl.comp_id, bl.primary_kingdom, COALESCE(sk.secondaries, ARRAY[]::text[]),
       bl.phylum, bl.taxclass, bl.taxorder, bl.family, bl.genus
FROM best_lineage_rows bl
LEFT JOIN secondary_kingdoms_per_compound sk ON sk.comp_id = bl.comp_id
WHERE bl.lr = 1;

-- Tier 2: compounds with compound_taxonomy rows but NO mappable kingdom.
-- Typically WCVP-resolved plants where ncbi_lineage and worms_lineage are NULL
-- but compound_taxonomy.family/genus are populated. Preserve their lineage,
-- use compounds.kingdom as primary.
WITH lineage_only AS (
    SELECT
        ct.comp_id,
        ct.family AS ct_family,
        ct.genus AS ct_genus,
        ct.ncbi_lineage->>'phylum' AS ct_phylum,
        ct.ncbi_lineage->>'class' AS ct_class,
        ct.ncbi_lineage->>'order' AS ct_order,
        (CASE WHEN ct.ncbi_lineage->>'phylum' IS NOT NULL THEN 1 ELSE 0 END +
         CASE WHEN ct.ncbi_lineage->>'class' IS NOT NULL THEN 1 ELSE 0 END +
         CASE WHEN ct.ncbi_lineage->>'order' IS NOT NULL THEN 1 ELSE 0 END +
         CASE WHEN ct.family IS NOT NULL THEN 1 ELSE 0 END +
         CASE WHEN ct.genus IS NOT NULL THEN 1 ELSE 0 END) AS specificity
    FROM compound_taxonomy ct
    WHERE ct.comp_id NOT IN (SELECT comp_id FROM resolved_taxonomy_new)
      AND (ct.family IS NOT NULL OR ct.genus IS NOT NULL
           OR ct.ncbi_lineage->>'family' IS NOT NULL
           OR ct.ncbi_lineage->>'order' IS NOT NULL)
),
best_lineage_only AS (
    SELECT comp_id, ct_phylum, ct_class, ct_order, ct_family, ct_genus,
           ROW_NUMBER() OVER (
             PARTITION BY comp_id
             ORDER BY specificity DESC, ct_family ASC NULLS LAST
           ) AS lr
    FROM lineage_only
)
INSERT INTO resolved_taxonomy_new (comp_id, kingdom, secondary_kingdoms, phylum, taxclass, taxorder, family, genus)
SELECT bl.comp_id,
       CASE WHEN c.kingdom IN ('plant','fungi','bacteria','animal') THEN c.kingdom
            ELSE 'unresolved' END AS kingdom,
       ARRAY[]::text[] AS secondary_kingdoms,
       bl.ct_phylum, bl.ct_class, bl.ct_order, bl.ct_family, bl.ct_genus
FROM best_lineage_only bl
JOIN compounds c ON c.comp_id = bl.comp_id
WHERE bl.lr = 1;

-- Tier 3: compounds with NO compound_taxonomy rows at all (or all-NULL lineage).
-- Use compounds.kingdom, NULL lineage.
INSERT INTO resolved_taxonomy_new (comp_id, kingdom, secondary_kingdoms, phylum, taxclass, taxorder, family, genus)
SELECT c.comp_id,
       CASE WHEN c.kingdom IN ('plant','fungi','bacteria','animal') THEN c.kingdom
            ELSE 'unresolved' END,
       ARRAY[]::text[], NULL, NULL, NULL, NULL, NULL
FROM compounds c
WHERE c.comp_id NOT IN (SELECT comp_id FROM resolved_taxonomy_new);

-- Swap and reindex
DROP INDEX IF EXISTS rt_kingdom_lower;
DROP INDEX IF EXISTS rt_phylum_lower;
DROP INDEX IF EXISTS rt_taxclass_lower;
DROP INDEX IF EXISTS rt_taxorder_lower;
DROP INDEX IF EXISTS rt_family_lower;
DROP INDEX IF EXISTS rt_genus_lower;
DROP INDEX IF EXISTS rt_secondary_kingdoms;
DROP TABLE IF EXISTS resolved_taxonomy_pre_v2;
ALTER TABLE resolved_taxonomy RENAME TO resolved_taxonomy_pre_v2;
ALTER TABLE resolved_taxonomy_new RENAME TO resolved_taxonomy;
CREATE INDEX rt_kingdom_lower ON resolved_taxonomy (LOWER(kingdom));
CREATE INDEX rt_phylum_lower ON resolved_taxonomy (LOWER(phylum));
CREATE INDEX rt_taxclass_lower ON resolved_taxonomy (LOWER(taxclass));
CREATE INDEX rt_taxorder_lower ON resolved_taxonomy (LOWER(taxorder));
CREATE INDEX rt_family_lower ON resolved_taxonomy (LOWER(family));
CREATE INDEX rt_genus_lower ON resolved_taxonomy (LOWER(genus));
CREATE INDEX rt_secondary_kingdoms ON resolved_taxonomy USING GIN (secondary_kingdoms);
GRANT SELECT ON resolved_taxonomy TO theobroma;

COMMIT;

-- Verification
SELECT 'Distribution:' AS msg;
SELECT kingdom, COUNT(*) FROM resolved_taxonomy GROUP BY kingdom ORDER BY 2 DESC;

SELECT 'Lineage coverage:' AS msg,
       COUNT(*) FILTER (WHERE family IS NOT NULL) AS with_family,
       COUNT(*) AS total FROM resolved_taxonomy;

SELECT 'bacteria+Streptophyta (should be 0):' AS msg;
SELECT COUNT(*) FROM resolved_taxonomy
WHERE kingdom = 'bacteria' AND phylum IN ('Streptophyta','Chlorophyta','Anthophyta','Marchantiophyta');

SELECT 'Asteraceae+Middle East primary distribution:' AS msg;
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
