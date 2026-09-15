-- =====================================================================
-- THEOBROMA : RESTORE live database to a clean, consistent v1.34 state.
-- Reset to the co-authors' known 1,133,004-compound baseline, from which a
-- new README is written and this session's corrections are re-applied.
--
-- Satellite layer (search_names matview over compounds_v134_archive,
-- scaffolds, admet, compound_taxonomy) is ALREADY v1.34-scale & comp_id-
-- keyed, so it aligns automatically. Only compounds + resolved_taxonomy
-- (+ attestation) are swapped back. No view depends on compounds (checked).
--
-- Reversible: v1.35 live preserved under *_v135_preswap_20260802 + the
-- verified full v1.35 dump (ac6100d5...) in two locations.
-- Indexes recreated from the ACTUAL captured definitions (v135-only-column
-- indexes omitted: tier_rank / source_organism_curated / apg_clade are not
-- in v1.34).
--
-- Transaction OPEN; inspect; COMMIT / ROLLBACK.
-- =====================================================================
\set ON_ERROR_STOP on
\timing on
BEGIN;

\echo '=== PRE-STATE (v135 live: compounds ~753992) ==='
SELECT 'compounds' t, count(*) FROM compounds
UNION ALL SELECT 'resolved_taxonomy', count(*) FROM resolved_taxonomy
UNION ALL SELECT 'per_source_license_attestation', count(*) FROM per_source_license_attestation;

-- 1. Preserve current v1.35 live tables under preswap names (reversible).
DROP TABLE IF EXISTS compounds_v135_preswap_20260802;
DROP TABLE IF EXISTS resolved_taxonomy_v135_preswap_20260802;
DROP TABLE IF EXISTS per_source_license_attestation_v135_preswap_20260802;
ALTER TABLE compounds                       RENAME TO compounds_v135_preswap_20260802;
ALTER TABLE resolved_taxonomy               RENAME TO resolved_taxonomy_v135_preswap_20260802;
ALTER TABLE per_source_license_attestation  RENAME TO per_source_license_attestation_v135_preswap_20260802;

-- 2. Install v1.34 archives as the live tables (copy; archives stay intact).
CREATE TABLE compounds AS SELECT * FROM compounds_v134_archive;
CREATE TABLE resolved_taxonomy AS SELECT * FROM resolved_taxonomy_v134_archive;
CREATE TABLE per_source_license_attestation AS SELECT * FROM per_source_license_attestation_v134_archive;

\echo '=== POST-SWAP counts (expect compounds=1133004) ==='
SELECT 'compounds' t, count(*) FROM compounds
UNION ALL SELECT 'resolved_taxonomy', count(*) FROM resolved_taxonomy
UNION ALL SELECT 'per_source_license_attestation', count(*) FROM per_source_license_attestation;

-- 3. Recreate indexes with FRESH unique names (_live) so they never collide
--    with index names already held by *_v134_archive / *_preswap tables.
--    v135-only-column indexes omitted (tier_rank/source_organism_curated/apg_clade
--    absent in v1.34).

-- compounds
CREATE UNIQUE INDEX compounds_pkey_live ON compounds USING btree (comp_id);
CREATE INDEX idx_btree_live ON compounds USING btree (name_norm text_pattern_ops);
CREATE INDEX idx_trgm_live ON compounds USING gin (name_norm gin_trgm_ops);
CREATE INDEX idx_compounds_chebi_id_live ON compounds USING btree (chebi_id) WHERE (chebi_id IS NOT NULL);
CREATE INDEX idx_compounds_kingdom_lower_live ON compounds USING btree (lower(kingdom), comp_id);
CREATE INDEX idx_compounds_organism_tokens_live ON compounds USING gin (string_to_array(lower(source_organism), '; '::text)) WHERE (source_organism IS NOT NULL);
CREATE INDEX idx_compounds_region_lower_live ON compounds USING btree (lower(region));
CREATE INDEX idx_compounds_source_db_lower_live ON compounds USING btree (lower(source_db));
CREATE INDEX idx_compounds_lower_organism_live ON compounds USING hash (lower(source_organism));
CREATE INDEX idx_ik_live ON compounds USING btree (inchikey);
CREATE INDEX idx_ik_prefix_live ON compounds USING btree ("substring"(inchikey, 1, 14));
CREATE INDEX idx_king_live ON compounds USING btree (kingdom);
CREATE INDEX idx_name_live ON compounds USING btree (name);
CREATE INDEX idx_novelty_live ON compounds USING btree (novelty_morgan);
CREATE INDEX idx_npc_live ON compounds USING btree (np_class);
CREATE INDEX idx_pathway_live ON compounds USING btree (np_pathway);
CREATE INDEX idx_reg_live ON compounds USING btree (region);
CREATE INDEX idx_src_live ON compounds USING btree (source_db);
CREATE INDEX idx_trust_live ON compounds USING btree (trust_score);

-- resolved_taxonomy
CREATE UNIQUE INDEX resolved_taxonomy_pkey_live ON resolved_taxonomy USING btree (comp_id);
CREATE INDEX rt_kingdom_lower_live  ON resolved_taxonomy USING btree (lower(kingdom));
CREATE INDEX rt_phylum_lower_live   ON resolved_taxonomy USING btree (lower(phylum));
CREATE INDEX rt_taxclass_lower_live ON resolved_taxonomy USING btree (lower(taxclass));
CREATE INDEX rt_taxorder_lower_live ON resolved_taxonomy USING btree (lower(taxorder));
CREATE INDEX rt_family_lower_live   ON resolved_taxonomy USING btree (lower(family));
CREATE INDEX rt_genus_lower_live    ON resolved_taxonomy USING btree (lower(genus));
-- optional columns guarded by existence
DO $$
BEGIN
  IF EXISTS (SELECT 1 FROM information_schema.columns
             WHERE table_name='resolved_taxonomy' AND column_name='secondary_kingdoms') THEN
    EXECUTE 'CREATE INDEX rt_secondary_kingdoms_live ON resolved_taxonomy USING gin (secondary_kingdoms)';
  END IF;
  IF EXISTS (SELECT 1 FROM information_schema.columns
             WHERE table_name='resolved_taxonomy' AND column_name='apg_clade') THEN
    EXECUTE 'CREATE INDEX rt_apg_clade_lower_live ON resolved_taxonomy USING btree (lower(apg_clade))';
  END IF;
END$$;

-- per_source_license_attestation
CREATE INDEX psla_pkey_live ON per_source_license_attestation USING btree (comp_id, source);  -- non-unique: v134 archive may not be dedup'd on this key
CREATE INDEX psla_comp_id_live ON per_source_license_attestation USING btree (comp_id);
CREATE INDEX psla_tier_live ON per_source_license_attestation USING btree (license_tier);

-- 4. Refresh search_names matview (defined over compounds_v134_archive).
REFRESH MATERIALIZED VIEW search_names;
\echo '=== search_names rows after refresh ==='
SELECT count(*) AS sn_rows, count(DISTINCT comp_id) AS sn_distinct_comp FROM search_names;

-- 5. ANALYZE for planner stats on the fresh tables.
ANALYZE compounds;
ANALYZE resolved_taxonomy;
ANALYZE per_source_license_attestation;

-- 6. Verification.
\echo '=== curcumin present (v134)? ==='
SELECT comp_id, name, kingdom, license_tier FROM compounds
WHERE inchikey LIKE 'VFLDPWHFBUODDF%' ORDER BY comp_id LIMIT 5;

\echo '=== integrity: joins across the restored ecosystem ==='
SELECT c.comp_id,
       (SELECT kingdom FROM resolved_taxonomy rt WHERE rt.comp_id=c.comp_id) AS tax_kingdom,
       (SELECT count(*) FROM per_source_license_attestation p WHERE p.comp_id=c.comp_id) AS lic_rows,
       (SELECT count(*) FROM scaffolds s WHERE s.comp_id=c.comp_id) AS scaffold_rows,
       (SELECT count(*) FROM admet a WHERE a.comp_id=c.comp_id) AS admet_rows,
       (SELECT count(*) FROM search_names sn WHERE sn.comp_id=c.comp_id) AS name_rows
FROM compounds c WHERE c.inchikey LIKE 'VFLDPWHFBUODDF%' LIMIT 3;

\echo '=== kingdom distribution (v134 expected) ==='
SELECT kingdom, count(*) FROM compounds GROUP BY 1 ORDER BY 2 DESC LIMIT 10;

\echo '======================================================================'
\echo 'INSPECT: compounds=1133004, curcumin present, joins non-zero.'
\echo 'If good: COMMIT;   else ROLLBACK;'
\echo 'v135 preserved under *_v135_preswap_20260802 (+ full dump ac6100d5).'
\echo '======================================================================'
-- COMMIT;
