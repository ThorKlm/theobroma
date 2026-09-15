-- =====================================================================
-- THEOBROMA validation probes — run on the server:
--   psql -h localhost -U theobroma -d theobroma -f validation_probes.sql
-- Read-only. Each block prints an invariant that SHOULD hold; deviations
-- are the failure points to investigate.
-- =====================================================================
\pset pager off
\timing off

\echo '========== 1. CORE ROW-COUNT INVARIANTS =========='
\echo '-- compounds should be exactly 1,133,004; one row per comp_id --'
SELECT count(*) AS n_compounds, count(DISTINCT comp_id) AS distinct_ids,
       count(*) - count(DISTINCT comp_id) AS dup_rows FROM compounds;
\echo '-- one row per inchikey? (stereo families share 14-char, but full IK should be unique) --'
SELECT count(*) - count(DISTINCT inchikey) AS dup_inchikeys FROM compounds;

\echo '========== 2. LICENSE DISTRIBUTION SUMS TO 100% =========='
SELECT license_tier, count(*),
       round(100.0*count(*)/sum(count(*)) OVER (), 2) AS pct
FROM compounds GROUP BY 1 ORDER BY 2 DESC;
\echo '-- tier_rank present for all? any NULL? --'
SELECT count(*) FILTER (WHERE tier_rank IS NULL) AS null_tier_rank FROM compounds;

\echo '========== 3. CLASSIFICATION SANITY =========='
\echo '-- any __ERROR__ sentinel leaked into the live table? (should be 0) --'
SELECT count(*) FROM compounds WHERE np_pathway = '__ERROR__' OR np_superclass='__ERROR__' OR np_class='__ERROR__';
\echo '-- how many have empty classification (legit unclassifiable + organometallics) --'
SELECT count(*) FILTER (WHERE np_pathway='' OR np_pathway IS NULL) AS empty_pathway,
       count(*) FILTER (WHERE np_class='' OR np_class IS NULL) AS empty_class FROM compounds;
\echo '-- separator sanity: any lingering ; that should be | ? --'
SELECT count(*) FILTER (WHERE np_pathway LIKE '%;%') AS pathway_semicolons,
       count(*) FILTER (WHERE np_class LIKE '%;%') AS class_semicolons FROM compounds;

\echo '========== 4. REFERENTIAL INTEGRITY (satellite tables) =========='
\echo '-- resolved_taxonomy comp_ids that do NOT exist in compounds (orphans) --'
SELECT count(*) AS rt_orphans FROM resolved_taxonomy rt
WHERE NOT EXISTS (SELECT 1 FROM compounds c WHERE c.comp_id=rt.comp_id);
\echo '-- compounds with NO taxonomy row (tree/kingdom will be blank for these) --'
SELECT count(*) AS compounds_without_taxonomy FROM compounds c
WHERE NOT EXISTS (SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id=c.comp_id);
\echo '-- attestation orphans --'
SELECT count(*) AS attest_orphans FROM per_source_license_attestation a
WHERE NOT EXISTS (SELECT 1 FROM compounds c WHERE c.comp_id=a.comp_id);

\echo '========== 5. KNOWN-UNFIXED: taxonomy chimeras / deprecated classes =========='
\echo '-- deprecated taxonomy classes still present (Liliopsida/Magnoliopsida) --'
SELECT taxclass, count(*) FROM resolved_taxonomy
WHERE taxclass IN ('Liliopsida','Magnoliopsida') GROUP BY 1;
\echo '-- curcumin (THEO_0854403) lineage: is it still chimeric (e.g. fern)? --'
SELECT comp_id, kingdom, phylum, taxclass, taxorder, family, genus
FROM resolved_taxonomy WHERE comp_id='THEO_0854403';

\echo '========== 6. KNOWN-UNFIXED: organism truncation on v134-exclusive rows =========='
\echo '-- how many still sit exactly at the 500-char cap (likely truncated) --'
SELECT count(*) AS at_500_cap FROM compounds WHERE length(source_organism)=500;

\echo '========== 7. SEARCH MATVIEW FRESHNESS =========='
\echo '-- search_names row count; does curcumin resolve? --'
SELECT count(*) AS search_names_rows FROM search_names;
SELECT comp_id, original_name, source_kind FROM search_names WHERE name_norm='curcumin' LIMIT 5;
