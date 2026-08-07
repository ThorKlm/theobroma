-- =====================================================================
-- THEOBROMA pipeline stage 06: NPClassifier load + apply (canonical)
-- =====================================================================
-- Purpose : Load the full-corpus NPClassifier re-run (canonical GNPS2 API)
--           and overwrite np_pathway/np_superclass/np_class on compounds,
--           yielding uniform coherent triples for all compounds (also fixes
--           the v134-exclusive residual classification scramble).
-- Inputs  : npc_output.csv (comp_id, np_pathway, np_superclass, np_class),
--           produced by npc_api_batch.py against the GNPS2 API.
-- Output  : compounds.np_pathway/np_superclass/np_class; npc_reclassified staging.
-- Run     : classification stage, after compounds is loaded. Followed by
--           06b (pathway-straggler ontology fix).
-- NOTE    : original was staged with an operator COMMIT after inspecting
--           curcumin + a cholestane; the trailing COMMIT is left commented
--           exactly as in the source. Multi-value separator normalized ';' -> ' | '.
-- Verbatim from npc_load_apply.sql.
-- =====================================================================
\set ON_ERROR_STOP on
\timing on

-- 0. Load the results into a staging table.
DROP TABLE IF EXISTS npc_reclassified;
CREATE TABLE npc_reclassified (
    comp_id text PRIMARY KEY,
    np_pathway text,
    np_superclass text,
    np_class text
);
\copy npc_reclassified FROM '/home/thorben.klamt/theobroma/npc_output.csv' CSV HEADER

\echo '=== coverage: how many rows loaded, how many classified vs empty vs error ==='
SELECT count(*) AS total,
       count(*) FILTER (WHERE np_pathway <> '' AND np_pathway <> '__ERROR__') AS has_pathway,
       count(*) FILTER (WHERE np_pathway = '') AS empty_pathway,
       count(*) FILTER (WHERE np_pathway = '__ERROR__') AS errored
FROM npc_reclassified;

\echo '=== do all live compounds have a row? ==='
SELECT (SELECT count(*) FROM compounds) AS live_compounds,
       (SELECT count(*) FROM npc_reclassified) AS reclassified_rows,
       (SELECT count(*) FROM compounds c WHERE NOT EXISTS
          (SELECT 1 FROM npc_reclassified r WHERE r.comp_id=c.comp_id)) AS compounds_missing;

-- 1. normalize the multi-value separator to the app convention ' | '.
UPDATE npc_reclassified SET
  np_pathway    = replace(np_pathway, ';', ' | '),
  np_superclass = replace(np_superclass, ';', ' | '),
  np_class      = replace(np_class, ';', ' | ')
WHERE np_pathway LIKE '%;%' OR np_superclass LIKE '%;%' OR np_class LIKE '%;%';

BEGIN;

-- 2. archive current classification for rollback.
DROP TABLE IF EXISTS compounds_class_pre_npcrerun_20260803;
CREATE TABLE compounds_class_pre_npcrerun_20260803 AS
SELECT comp_id, np_pathway, np_superclass, np_class FROM compounds;
\echo '=== archived current classification ==='
SELECT count(*) FROM compounds_class_pre_npcrerun_20260803;

-- 3. BEFORE snapshot of the two check compounds.
\echo '=== BEFORE: curcumin THEO_0854403 + cholestane THEO_1055282 ==='
SELECT comp_id, np_pathway, np_superclass, np_class FROM compounds
WHERE comp_id IN ('THEO_0854403','THEO_1055282') ORDER BY comp_id;

-- 4. apply: overwrite classification from the re-run (only where non-empty).
UPDATE compounds c
SET np_pathway    = r.np_pathway,
    np_superclass = r.np_superclass,
    np_class      = r.np_class
FROM npc_reclassified r
WHERE c.comp_id = r.comp_id
  AND r.np_pathway <> '' AND r.np_pathway <> '__ERROR__';

\echo '=== rows updated ==='
SELECT count(*) FROM npc_reclassified WHERE np_pathway <> '' AND np_pathway <> '__ERROR__';

-- 5. AFTER verification on the two known compounds.
\echo '=== AFTER: curcumin (expect Shikimates and Phenylpropanoids / Diarylheptanoids / Linear diarylheptanoids) ==='
SELECT comp_id, np_pathway, np_superclass, np_class FROM compounds WHERE comp_id='THEO_0854403';
\echo '=== AFTER: THEO_1055282 (expect Terpenoids / Steroids / Cholestane steroids) ==='
SELECT comp_id, np_pathway, np_superclass, np_class FROM compounds WHERE comp_id='THEO_1055282';

-- 6. coherence checks.
\echo '=== Quassinoids pathway/superclass spread AFTER ==='
SELECT np_pathway, np_superclass, count(*) FROM compounds
WHERE np_class LIKE '%Quassinoids%' GROUP BY 1,2 ORDER BY 3 DESC LIMIT 8;
\echo '=== Quassinoids by kingdom AFTER ==='
SELECT kingdom, count(*) FROM compounds WHERE np_class LIKE '%Quassinoids%' GROUP BY 1 ORDER BY 2 DESC;
\echo '=== distinct superclasses per pathway (coherence sanity) ==='
SELECT np_pathway, count(DISTINCT np_superclass) AS distinct_super, count(*) n
FROM compounds WHERE np_pathway<>'' GROUP BY 1 ORDER BY 3 DESC LIMIT 10;

\echo '======================================================================'
\echo 'INSPECT curcumin + THEO_1055282 correct, Quassinoids coherent.'
\echo 'If good: COMMIT;  else ROLLBACK;'
\echo '======================================================================'
-- COMMIT;
