-- =====================================================================
-- THEOBROMA pipeline stage 08: geographic region derivation (multi-valued)
-- =====================================================================
-- Purpose : Build compound_region_map(comp_id, macro_region), multi-valued.
--           Per compound, the SET of:
--             { crosswalk(UNION native_distribution regions across all tokens) }
--             U { fixed-scope regional-source tags from all_sources }
--             U { existing compounds.region if not global }  (superset guarantee)
--           Applies the WCVP 52-subregion -> 13-macro-region crosswalk.
-- Inputs  : compound_taxonomy(native_distribution), compounds(all_sources, region)
-- Output  : compound_region_map (multi-valued, deduped)
-- Run     : after stage 03 (compound_taxonomy) is final.
-- Status  : DONE / canonical. Live: 127,992 single-region, 41,548 two,
--           25,461 three, long tail. Supersedes region_derive_stage2_backfill.sql.
-- NOTE    : original ran interactively with an operator COMMIT after review.
--           Consolidated here as a runnable script; review verification output.
-- Verbatim from region_derive_stage2_persist.sql.
-- =====================================================================

-- Crosswalk as a VALUES table (the 52->13 mapping; unmapped subregions ignored).
CREATE TABLE sub2macro(sub text PRIMARY KEY, macro text) ;
INSERT INTO sub2macro(sub,macro) VALUES
 ('China','East Asia'),('Eastern Asia','East Asia'),('Mongolia','East Asia'),
 ('Indian Subcontinent','South Asia'),
 ('Malesia','Southeast Asia'),('Indo-China','Southeast Asia'),('Papuasia','Southeast Asia'),
 ('Arabian Peninsula','Middle East'),('Western Asia','Middle East'),
 ('Siberia','Russia/CIS'),('Russian Far East','Russia/CIS'),('Caucasus','Russia/CIS'),
 ('Central Asia','Central Asia'),
 ('Middle Europe','Europe'),('Northern Europe','Europe'),('Southwestern Europe','Europe'),
 ('Southeastern Europe','Europe'),('Eastern Europe','Europe'),('Macaronesia','Europe'),
 ('Mexico','Latin America'),('Central America','Latin America'),('Caribbean','Latin America'),
 ('Brazil','Latin America'),('Western South America','Latin America'),
 ('Northern South America','Latin America'),('Southern South America','Latin America'),
 ('Northeastern U.S.A.','North America'),('Southeastern U.S.A.','North America'),
 ('North-Central U.S.A.','North America'),('South-Central U.S.A.','North America'),
 ('Northwestern U.S.A.','North America'),('Southwestern U.S.A.','North America'),
 ('Eastern Canada','North America'),('Western Canada','North America'),('Subarctic America','North America'),
 ('Northern Africa','Africa'),('Northeast Tropical Africa','Africa'),('East Tropical Africa','Africa'),
 ('West Tropical Africa','Africa'),('West-Central Tropical Africa','Africa'),
 ('South Tropical Africa','Africa'),('Southern Africa','Africa'),('Western Indian Ocean','Africa'),
 ('Australia','Australia'),
 ('New Zealand','New Zealand'),
 ('Southwestern Pacific','Oceania'),('South-Central Pacific','Oceania'),
 ('North-Central Pacific','Oceania'),('Northwestern Pacific','Oceania');

-- Fixed-scope regional source databases -> macro (n_regions=1, editorially confirmed).
CREATE TABLE source_region(src text PRIMARY KEY, macro text) ;
INSERT INTO source_region(src,macro) VALUES
 ('ANPDB','Africa'),('AfroDb','Africa'),('SANCDB','Africa'),
 ('CSIRO','Australia'),
 ('IMPPAT','South Asia'),('SMDB_Spice','South Asia'),('CMDB_Cereals','South Asia'),
 ('TMDB_Trichoderma','South Asia'),('phytochemdb','South Asia'),
 ('TIPdb','East Asia'),
 ('NaturAr','Latin America');

DROP TABLE IF EXISTS compound_region_map;
CREATE TABLE compound_region_map(comp_id text NOT NULL, macro_region text NOT NULL);

-- (1) WCVP: union native_distribution 'regions' across ALL tokens per compound, crosswalk to macro.
INSERT INTO compound_region_map(comp_id, macro_region)
SELECT DISTINCT ct.comp_id, s.macro
FROM compound_taxonomy ct
CROSS JOIN LATERAL jsonb_array_elements_text(ct.native_distribution->'regions') AS wr(sub)
JOIN sub2macro s ON s.sub = wr.sub
WHERE ct.native_distribution IS NOT NULL;

-- (2) fixed-scope regional source tags from all_sources (pipe-delimited).
INSERT INTO compound_region_map(comp_id, macro_region)
SELECT DISTINCT c.comp_id, sr.macro
FROM compounds c
CROSS JOIN LATERAL unnest(string_to_array(c.all_sources,'|')) AS src(s)
JOIN source_region sr ON sr.src = trim(src.s)
WHERE c.all_sources IS NOT NULL;

-- (3) existing single region (superset guarantee), unless global/blank.
INSERT INTO compound_region_map(comp_id, macro_region)
SELECT DISTINCT c.comp_id, c.region
FROM compounds c
WHERE c.region IS NOT NULL AND c.region <> '' AND c.region <> 'global';

-- Deduplicate (a compound may get the same macro from several sources).
CREATE TABLE compound_region_map_dedup AS
  SELECT DISTINCT comp_id, macro_region FROM compound_region_map;
DROP TABLE compound_region_map;
ALTER TABLE compound_region_map_dedup RENAME TO compound_region_map;
CREATE INDEX idx_crm_comp ON compound_region_map(comp_id);
CREATE INDEX idx_crm_macro ON compound_region_map(macro_region);

-- ---------- VERIFICATION ----------
\echo '--- total rows + distinct compounds ---'
SELECT count(*) AS rows, count(DISTINCT comp_id) AS distinct_comps FROM compound_region_map;
\echo '--- distinct macro_region values (should be 13) ---'
SELECT macro_region, count(*) AS n FROM compound_region_map GROUP BY 1 ORDER BY 2 DESC;
\echo '--- SUPERSET INVARIANT: old single region not in new map (expect 0) ---'
SELECT count(*) FROM compounds c
WHERE c.region IS NOT NULL AND c.region <> '' AND c.region <> 'global'
  AND NOT EXISTS (SELECT 1 FROM compound_region_map m WHERE m.comp_id=c.comp_id AND m.macro_region=c.region);
\echo '--- multi-region compounds (gained >1) ---'
SELECT count(*) FROM (SELECT comp_id FROM compound_region_map GROUP BY 1 HAVING count(*)>1) t;

-- cleanup helper tables
DROP TABLE IF EXISTS sub2macro;
DROP TABLE IF EXISTS source_region;
