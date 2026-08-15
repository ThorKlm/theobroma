-- =====================================================================
-- Null out literal placeholder strings in compounds.name so they render
-- as "no name" ( -> em-dash in the UI) instead of misleading literal text
-- like "No Data". Does NOT touch real identifiers (Mol_ IDs, CAS numbers,
-- registry codes) -- only strings that mean "there is no name here".
-- Staged: preview -> archive -> update -> verify -> COMMIT/ROLLBACK.
-- =====================================================================
\set ON_ERROR_STOP on
\pset pager off

-- The placeholder set (case-insensitive), plus trivial junk tokens.
-- Kept intentionally NARROW: only strings that are semantically "no name".
\echo '=== PREVIEW: rows whose name is a placeholder (will be nulled) ==='
SELECT name, count(*) n FROM compounds
WHERE lower(btrim(name)) IN (
  'no data','not named','unknown','n/a','na','none','null','unnamed',
  'no name','not available','not assigned','-','--','---','?','n.a.','nan','.'
)
GROUP BY name ORDER BY n DESC;

\echo '=== total rows affected ==='
SELECT count(*) AS to_null FROM compounds
WHERE lower(btrim(name)) IN (
  'no data','not named','unknown','n/a','na','none','null','unnamed',
  'no name','not available','not assigned','-','--','---','?','n.a.','nan','.'
);

BEGIN;

-- Archive the affected comp_id + original name for provenance/rollback.
DROP TABLE IF EXISTS compounds_name_placeholder_pre_null_20260804;
CREATE TABLE compounds_name_placeholder_pre_null_20260804 AS
SELECT comp_id, name AS name_original FROM compounds
WHERE lower(btrim(name)) IN (
  'no data','not named','unknown','n/a','na','none','null','unnamed',
  'no name','not available','not assigned','-','--','---','?','n.a.','nan','.'
);
\echo '=== archived ==='
SELECT count(*) FROM compounds_name_placeholder_pre_null_20260804;

-- Null them (set to empty string to match the existing "no name" convention;
-- the app already renders empty/NULL name as an em-dash).
UPDATE compounds SET name = NULL
WHERE lower(btrim(name)) IN (
  'no data','not named','unknown','n/a','na','none','null','unnamed',
  'no name','not available','not assigned','-','--','---','?','n.a.','nan','.'
);

\echo '=== verify: none of these placeholders remain as names ==='
SELECT count(*) AS remaining FROM compounds
WHERE lower(btrim(name)) IN (
  'no data','not named','unknown','n/a','na','none','null','unnamed',
  'no name','not available','not assigned','-','--','---','?','n.a.','nan','.'
);
\echo '=== confirm Mol_ IDs and CAS numbers were NOT touched (should still be present) ==='
SELECT count(*) FILTER (WHERE name ~ '^Mol_[0-9]+$') AS mol_ids_intact,
       count(*) FILTER (WHERE name ~ '^[0-9]{2,7}-[0-9]{2}-[0-9]$') AS cas_numbers_intact
FROM compounds WHERE name IS NOT NULL;

\echo '=== INSPECT above. If good: COMMIT; else ROLLBACK; ==='
-- COMMIT;
