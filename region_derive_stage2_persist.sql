-- Stage 2: build compound_region_map (comp_id, macro_region) multi-valued.
-- Composition per compound (a SET):
--   { crosswalk(UNION of native_distribution 'regions' across ALL organism-tokens) }
--   U { fixed-scope regional-source tags from all_sources }
--   U { existing compounds.region if not global }   -- superset guarantee
-- Transactional: builds into a NEW table, verifies, leaves COMMIT to the operator.
-- Run interactively:  psql -h localhost -U theobroma -d theobroma -f region_derive_stage2_backfill.sql
-- (Review the verification output BEFORE typing COMMIT.)

-- (autocommit: no explicit transaction)

-- Crosswalk as a VALUES table (your drafted 52->13 mapping; DROP set omitted -> unmapped ignored).
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

-- ---------- VERIFICATION (review before COMMIT) ----------
\echo '--- total rows + distinct compounds ---'
SELECT count(*) AS rows, count(DISTINCT comp_id) AS distinct_comps FROM compound_region_map;

\echo '--- distinct macro_region values (should be your 13) ---'
SELECT macro_region, count(*) AS n FROM compound_region_map GROUP BY 1 ORDER BY 2 DESC;

\echo '--- curcumin (expect a broad multi-region set, ~12) ---'
SELECT macro_region FROM compound_region_map WHERE comp_id='THEO_0901197' ORDER BY 1;

\echo '--- SUPERSET INVARIANT: compounds whose old single region is NOT in the new map (expect 0) ---'
SELECT count(*) FROM compounds c
WHERE c.region IS NOT NULL AND c.region <> '' AND c.region <> 'global'
  AND NOT EXISTS (SELECT 1 FROM compound_region_map m WHERE m.comp_id=c.comp_id AND m.macro_region=c.region);

\echo '--- coverage vs old single region (expect new >= old 181,820) ---'
SELECT (SELECT count(DISTINCT comp_id) FROM compound_region_map) AS new_covered,
       (SELECT count(*) FROM compounds WHERE region IS NOT NULL AND region<>'' AND region<>'global') AS old_single;

\echo '--- multi-region compounds (how many gained >1 region) ---'
SELECT count(*) FROM (SELECT comp_id FROM compound_region_map GROUP BY 1 HAVING count(*)>1) t;

\echo '>>> Review the above. If correct, type:  COMMIT;   else:  ROLLBACK;'

-- cleanup helper tables (they were created as regular tables in autocommit mode)
DROP TABLE IF EXISTS sub2macro;
DROP TABLE IF EXISTS source_region;

\echo '=== PERSISTED. Final check: ==='
SELECT count(*) AS rows, count(DISTINCT comp_id) AS comps FROM compound_region_map;
SELECT count(*) FROM compounds c
 WHERE c.region IS NOT NULL AND c.region<>'' AND c.region<>'global'
   AND NOT EXISTS (SELECT 1 FROM compound_region_map m WHERE m.comp_id=c.comp_id AND m.macro_region=c.region)
 ; -- superset violations, expect 0
