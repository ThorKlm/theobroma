-- =====================================================================
-- Remove 199 non-natural-product entries (<=3 total atoms: single elements,
-- monatomic ions, noble gases, and simple inorganic molecules like water/CO2/
-- HCN). comp_ids come from small_atoms.csv (produced by find_small_atoms.py
-- via RDKit, which correctly counts implicit H).
--
-- Deletes ONLY from LIVE tables; historical archive/backup tables are left
-- intact for provenance. Full rows archived first for reversibility.
-- Staged: load -> archive -> preview -> delete -> verify -> COMMIT/ROLLBACK.
--
-- Prereq on server:  \copy loads the CSV (run from ~/theobroma where the file is)
-- =====================================================================
\set ON_ERROR_STOP on
\pset pager off

-- 1. Load the removal list.
DROP TABLE IF EXISTS _inorganic_ids;
CREATE TEMP TABLE _inorganic_ids (comp_id text, name text, mw text, smiles text, total_atoms int);
\copy _inorganic_ids FROM 'small_atoms.csv' WITH (FORMAT csv, HEADER true)
\echo '=== loaded removal ids (expect 199) ==='
SELECT count(*) FROM _inorganic_ids;

\echo '=== preview: these will be removed (sample) ==='
SELECT comp_id, name, smiles, total_atoms FROM _inorganic_ids ORDER BY total_atoms, comp_id LIMIT 15;

BEGIN;

-- 2. Archive full live rows for reversibility (only the affected comp_ids).
DROP TABLE IF EXISTS compounds_pre_inorganic_removal_20260804;
CREATE TABLE compounds_pre_inorganic_removal_20260804 AS
SELECT c.* FROM compounds c JOIN _inorganic_ids i ON i.comp_id=c.comp_id;
DROP TABLE IF EXISTS resolved_taxonomy_pre_inorganic_removal_20260804;
CREATE TABLE resolved_taxonomy_pre_inorganic_removal_20260804 AS
SELECT r.* FROM resolved_taxonomy r JOIN _inorganic_ids i ON i.comp_id=r.comp_id;
DROP TABLE IF EXISTS admet_pre_inorganic_removal_20260804;
CREATE TABLE admet_pre_inorganic_removal_20260804 AS
SELECT a.* FROM admet a JOIN _inorganic_ids i ON i.comp_id=a.comp_id;
DROP TABLE IF EXISTS attest_pre_inorganic_removal_20260804;
CREATE TABLE attest_pre_inorganic_removal_20260804 AS
SELECT p.* FROM per_source_license_attestation p JOIN _inorganic_ids i ON i.comp_id=p.comp_id;
DROP TABLE IF EXISTS comptax_pre_inorganic_removal_20260804;
CREATE TABLE comptax_pre_inorganic_removal_20260804 AS
SELECT ct.* FROM compound_taxonomy ct JOIN _inorganic_ids i ON i.comp_id=ct.comp_id;
DROP TABLE IF EXISTS scaffolds_pre_inorganic_removal_20260804;
CREATE TABLE scaffolds_pre_inorganic_removal_20260804 AS
SELECT s.* FROM scaffolds s JOIN _inorganic_ids i ON i.comp_id=s.comp_id;
\echo '=== archived row counts per table ==='
SELECT 'compounds' t, count(*) n FROM compounds_pre_inorganic_removal_20260804
UNION ALL SELECT 'resolved_taxonomy', count(*) FROM resolved_taxonomy_pre_inorganic_removal_20260804
UNION ALL SELECT 'admet', count(*) FROM admet_pre_inorganic_removal_20260804
UNION ALL SELECT 'attestation', count(*) FROM attest_pre_inorganic_removal_20260804
UNION ALL SELECT 'compound_taxonomy', count(*) FROM comptax_pre_inorganic_removal_20260804
UNION ALL SELECT 'scaffolds', count(*) FROM scaffolds_pre_inorganic_removal_20260804;

-- 3. Count before.
\echo '=== compounds count BEFORE ==='
SELECT count(*) AS before_count FROM compounds;

-- 4. Delete from LIVE tables only.
DELETE FROM resolved_taxonomy       WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids);
DELETE FROM admet                   WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids);
DELETE FROM per_source_license_attestation WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids);
DELETE FROM compound_taxonomy       WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids);
DELETE FROM scaffolds               WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids);
DELETE FROM compounds               WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids);

-- 5. Verify.
\echo '=== compounds count AFTER (expect before - 199) ==='
SELECT count(*) AS after_count FROM compounds;
\echo '=== any of the 199 still present anywhere in live tables? (expect 0 each) ==='
SELECT 'compounds' t, count(*) n FROM compounds WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids)
UNION ALL SELECT 'resolved_taxonomy', count(*) FROM resolved_taxonomy WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids)
UNION ALL SELECT 'admet', count(*) FROM admet WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids)
UNION ALL SELECT 'attestation', count(*) FROM per_source_license_attestation WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids)
UNION ALL SELECT 'compound_taxonomy', count(*) FROM compound_taxonomy WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids)
UNION ALL SELECT 'scaffolds', count(*) FROM scaffolds WHERE comp_id IN (SELECT comp_id FROM _inorganic_ids);

\echo '======================================================================'
\echo 'INSPECT: after_count = before - 199, all live tables show 0 remaining.'
\echo 'If good: COMMIT; then REFRESH MATERIALIZED VIEW search_names; else ROLLBACK;'
\echo '======================================================================'
-- COMMIT;
