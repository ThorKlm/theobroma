-- =====================================================================
-- FIX Anne's #1 (classification scramble) and #5 (organism truncation)
-- by porting the CORRECT, coherent values from compounds_v135_preswap
-- onto the v1.34 live base, joined on comp_id (543,635 shared).
--
-- #1: v134 np_pathway/superclass/class are mis-associated (Quassinoids
--     paired with 17 different pathway/superclass combos). v135_preswap
--     holds coherent triples (THEO_1055282 -> Terpenoids/Steroids/
--     Cholestane steroids, matching NPClassifier on its SMILES).
-- #5: v134 source_organism truncated at 500 chars (77,660 compounds cut).
--     v135_preswap holds the full untruncated list (Morphine incl. Papaver).
--
-- Scope: only the 543,635 comp_ids present in v135_preswap are corrected.
-- v134-exclusive compounds keep their values (documented limitation).
-- Archives originals; transaction OPEN; inspect; COMMIT / ROLLBACK.
-- =====================================================================
\set ON_ERROR_STOP on
\timing on
BEGIN;

-- 0. Archive current values for rollback/audit.
DROP TABLE IF EXISTS compounds_class_organism_pre_fix_20260803;
CREATE TABLE compounds_class_organism_pre_fix_20260803 AS
SELECT comp_id, np_pathway, np_superclass, np_class, source_organism
FROM compounds;
\echo '=== archived originals ==='
SELECT count(*) FROM compounds_class_organism_pre_fix_20260803;

-- 1. PREVIEW the two target compounds before/after.
\echo '=== BEFORE: THEO_1055282 classification (scrambled) + Morphine organism (truncated) ==='
SELECT comp_id, np_pathway, np_superclass, np_class FROM compounds WHERE comp_id='THEO_1055282';
SELECT comp_id, length(source_organism) AS org_len, source_organism FROM compounds WHERE comp_id='THEO_0511220';

-- 2. Apply classification fix (coherent triple as a unit from v135_preswap).
UPDATE compounds c
SET np_pathway   = p.np_pathway,
    np_superclass= p.np_superclass,
    np_class     = p.np_class
FROM compounds_v135_preswap_20260802 p
WHERE c.comp_id = p.comp_id;
\echo '=== classification rows updated (= shared comp_ids) ==='
SELECT count(*) AS class_updated
FROM compounds c JOIN compounds_v135_preswap_20260802 p ON p.comp_id=c.comp_id;

-- 3. Apply organism fix (full untruncated list from v135_preswap).
UPDATE compounds c
SET source_organism = p.source_organism
FROM compounds_v135_preswap_20260802 p
WHERE c.comp_id = p.comp_id
  AND length(p.source_organism) > length(c.source_organism);  -- only where preswap is fuller
\echo '=== organism rows updated ==='
SELECT count(*) AS organism_updated
FROM compounds c JOIN compounds_v135_preswap_20260802 p ON p.comp_id=c.comp_id
WHERE length(p.source_organism) > 0;

-- 4. VERIFY the two target compounds now correct.
\echo '=== AFTER: THEO_1055282 (expect Terpenoids/Steroids/Cholestane steroids) ==='
SELECT comp_id, np_pathway, np_superclass, np_class FROM compounds WHERE comp_id='THEO_1055282';
\echo '=== AFTER: Morphine (expect full list incl. Papaver) ==='
SELECT comp_id, length(source_organism) AS org_len, source_organism FROM compounds WHERE comp_id='THEO_0511220';

-- 5. VERIFY classification coherence improved: Quassinoids pathway spread.
\echo '=== Quassinoids pathway/superclass spread AFTER (should collapse toward coherent) ==='
SELECT np_pathway, np_superclass, count(*) FROM compounds
WHERE np_class='Quassinoids' GROUP BY 1,2 ORDER BY 3 DESC LIMIT 10;

\echo '=== Quassinoids by kingdom AFTER (should be more plant-dominated) ==='
SELECT kingdom, count(*) FROM compounds WHERE np_class='Quassinoids' GROUP BY 1 ORDER BY 2 DESC;

\echo '======================================================================'
\echo 'INSPECT: THEO_1055282 coherent, Morphine full w/ Papaver. COMMIT / ROLLBACK.'
\echo '======================================================================'
-- COMMIT;
