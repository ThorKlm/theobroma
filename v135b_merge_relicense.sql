-- =====================================================================
-- THEOBROMA v1.35-B : MERGE v1.34 base + genuinely-new v1.35 compounds,
-- then RE-LICENSE with most-restrictive-wins.  Column-safe construction:
-- build from the v134 archive (its own columns), ADD the 6 missing v135
-- columns, then insert the new-210k rows selecting matching columns.
--
--   Base  = compounds_v134_archive (1,133,004) preserved exactly
--   Add   = 210,357 v135 compounds whose InChIKey not in v134
--           (comp_ids THEO_1135328..THEO_1345684, disjoint - verified)
--   Keep v134 versions of the 543,635 shared InChIKeys.
--   Relicense all via source_license_ref (54 src) most-restrictive-wins.
--   Taxonomy untouched (separate optional follow-on).
--   Expected final: 1,343,361 compounds.
--
-- Staged; transaction OPEN; inspect then COMMIT; or ROLLBACK;.
-- Full dumps verified in two locations already.
-- =====================================================================
\set ON_ERROR_STOP on
\timing on
BEGIN;

-- 0. Safety archive of current live compounds.
DROP TABLE IF EXISTS compounds_pre_merge_20260802;
CREATE TABLE compounds_pre_merge_20260802 AS SELECT * FROM compounds;
\echo '=== archived current live compounds ==='
SELECT count(*) AS pre_merge_rows FROM compounds_pre_merge_20260802;

-- 1. Build merged staging FROM the v134 archive (inherits all its columns).
DROP TABLE IF EXISTS compounds_merged_stage;
CREATE TABLE compounds_merged_stage AS SELECT * FROM compounds_v134_archive;
\echo '=== v134 base copied (expect 1133004) ==='
SELECT count(*) AS base_rows FROM compounds_merged_stage;

-- 1a. Add the 6 v135-only columns (only if absent). Defaults for v134 rows.
ALTER TABLE compounds_merged_stage
  ADD COLUMN IF NOT EXISTS tier_rank int,
  ADD COLUMN IF NOT EXISTS structure_is_fact boolean,
  ADD COLUMN IF NOT EXISTS source_organism_curated text,
  ADD COLUMN IF NOT EXISTS is_matched boolean,
  ADD COLUMN IF NOT EXISTS annotation_source text,
  ADD COLUMN IF NOT EXISTS provenance text;

-- 1b. Populate the new columns for the v134 base rows.
UPDATE compounds_merged_stage
SET source_organism_curated = COALESCE(source_organism_curated, source_organism),
    is_matched        = COALESCE(is_matched, TRUE),
    annotation_source = COALESCE(annotation_source, 'v134'),
    provenance        = COALESCE(provenance, 'v134_archive');

-- 1c. Insert the genuinely-new v135 compounds with EXPLICIT name-matched
--     columns (generated from the live schema) so insertion is by name,
--     never by position -- no misalignment possible.
INSERT INTO compounds_merged_stage (comp_id, name, smiles, inchi, inchikey, source_db, kingdom, region, source_organism, mw, logp, tpsa, hba, hbd, n_rings, rotatable_bonds, license_tier, all_sources, classyfire_superclass, np_class, reference_doi, name_norm, inferred_class, inferred_class_source, inferred_confidence, np_pathway, np_superclass, trad_medicine, trust_score, novelty_morgan, novelty_maccs, novelty_type, sa_score, chebi_id, chebi_name, chebi_definition, chebi_iupac_name, chebi_pmids, chebi_xrefs, annotation_source, provenance, is_matched, tier_rank, structure_is_fact, source_organism_curated)
SELECT comp_id, name, smiles, inchi, inchikey, source_db, kingdom, region, source_organism, mw, logp, tpsa, hba, hbd, n_rings, rotatable_bonds, license_tier, all_sources, classyfire_superclass, np_class, reference_doi, name_norm, inferred_class, inferred_class_source, inferred_confidence, np_pathway, np_superclass, trad_medicine, trust_score, novelty_morgan, novelty_maccs, novelty_type, sa_score, chebi_id, chebi_name, chebi_definition, chebi_iupac_name, chebi_pmids, chebi_xrefs, annotation_source, provenance, is_matched, tier_rank, structure_is_fact, source_organism_curated
FROM compounds c WHERE NOT EXISTS (SELECT 1 FROM compounds_v134_archive a WHERE a.inchikey = c.inchikey);
\echo '=== after adding new-210k (expect 1343361) ==='
SELECT count(*) AS merged_total FROM compounds_merged_stage;

\echo '=== uniqueness check ==='
SELECT count(*) total, count(DISTINCT comp_id) uniq_id, count(DISTINCT inchikey) uniq_ik
FROM compounds_merged_stage;

-- 2. RELICENSE most-restrictive-wins (independent-CC0 exception).
DROP TABLE IF EXISTS merged_license_resolved;
CREATE TEMP TABLE merged_license_resolved AS
WITH att AS (
  SELECT a.comp_id, r.tier_rank
  FROM per_source_license_attestation a
  JOIN source_license_ref r ON r.src = lower(a.source)
  WHERE a.comp_id IN (SELECT comp_id FROM compounds_merged_stage)
),
res AS (
  SELECT comp_id, max(tier_rank) AS worst_rank, bool_or(tier_rank = 0) AS has_cc0
  FROM att GROUP BY comp_id
)
SELECT comp_id,
  CASE WHEN worst_rank >= 2 AND has_cc0 THEN 0 ELSE worst_rank END AS final_rank
FROM res;

UPDATE compounds_merged_stage m
SET tier_rank = lr.final_rank,
    license_tier = CASE lr.final_rank
        WHEN 0 THEN 'CC0' WHEN 1 THEN 'CC BY 4.0' WHEN 2 THEN 'CC BY-NC 4.0'
        WHEN 3 THEN 'CC BY-NC-SA 4.0' WHEN 4 THEN 'CC BY-NC-ND 4.0'
        ELSE 'Unspecified' END
FROM merged_license_resolved lr
WHERE m.comp_id = lr.comp_id;

UPDATE compounds_merged_stage
SET tier_rank = 5, license_tier = 'Unspecified'
WHERE tier_rank IS NULL;

\echo '=== MERGED license distribution ==='
SELECT tier_rank, license_tier, count(*) FROM compounds_merged_stage GROUP BY 1,2 ORDER BY 1;
\echo '=== commercial / non-commercial / unspecified ==='
SELECT CASE WHEN tier_rank<=1 THEN 'commercial' WHEN tier_rank<=4 THEN 'non-commercial'
            ELSE 'unspecified' END AS bucket, count(*)
FROM compounds_merged_stage GROUP BY 1 ORDER BY 2 DESC;

-- 3. Spot-checks.
\echo '=== curcumin family present & licensed ==='
SELECT comp_id, name, license_tier, tier_rank FROM compounds_merged_stage
WHERE inchikey LIKE 'VFLDPWHFBUODDF%' ORDER BY comp_id LIMIT 5;
\echo '=== provenance split + tier coverage ==='
SELECT provenance, count(*), count(*) FILTER (WHERE tier_rank IS NOT NULL) AS have_tier
FROM compounds_merged_stage GROUP BY 1;

-- 4. SWAP (uncomment when previews look right). Reversible via archive + dumps.
-- ALTER TABLE compounds RENAME TO compounds_v135a_replaced_20260802;
-- ALTER TABLE compounds_merged_stage RENAME TO compounds;
-- Recreate indexes on the new compounds table (see companion note).

\echo '======================================================================'
\echo 'INSPECT merged_total=1343361, license dist, spot-checks.'
\echo 'If good: uncomment SWAP (step 4), then COMMIT;  else ROLLBACK;'
\echo 'Transaction OPEN.'
\echo '======================================================================'
-- COMMIT;
