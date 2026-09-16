-- Measure the licence movement from the September re-audit against the
-- pre-audit snapshot. Run AFTER pipeline/09_license.sql, 09b and
-- scripts/recompute_license_ranks.sql have been committed.
-- Read-only: creates a temp table and reports, changes nothing.
\timing on

CREATE TEMP TABLE pre (comp_id text, license_tier text, tier_rank int, tier_rank_min int);
\copy pre FROM 'backups/reaudit_20260915/compound_tiers_pre.csv' WITH (FORMAT csv, HEADER)

\echo '=== rows compared ==='
SELECT count(*) AS matched FROM compounds c JOIN pre p USING (comp_id);

\echo '=== transition table: resolved tier before and after ==='
SELECT p.license_tier AS was, c.license_tier AS now, count(*) AS n
FROM compounds c JOIN pre p USING (comp_id)
WHERE c.license_tier IS DISTINCT FROM p.license_tier
GROUP BY 1, 2 ORDER BY 3 DESC;

\echo '=== headline movement across the open boundary ==='
SELECT count(*) FILTER (WHERE p.tier_rank <= 1 AND c.tier_rank >  1) AS left_open,
       count(*) FILTER (WHERE p.tier_rank >  1 AND c.tier_rank <= 1) AS became_open,
       count(*) FILTER (WHERE c.tier_rank IS DISTINCT FROM p.tier_rank) AS any_rank_change
FROM compounds c JOIN pre p USING (comp_id);

\echo '=== open fraction before and after ==='
SELECT (SELECT count(*) FROM pre WHERE tier_rank <= 1) AS open_before,
       (SELECT count(*) FROM compounds WHERE tier_rank <= 1) AS open_after,
       (SELECT count(*) FROM compounds) AS corpus;

\echo '=== least-restrictive bound before and after ==='
SELECT (SELECT count(*) FROM pre WHERE tier_rank_min = 0) AS min_cc0_before,
       (SELECT count(*) FROM compounds WHERE tier_rank_min = 0) AS min_cc0_after,
       (SELECT count(*) FROM pre WHERE tier_rank_min = 1) AS min_ccby_before,
       (SELECT count(*) FROM compounds WHERE tier_rank_min = 1) AS min_ccby_after,
       (SELECT count(*) FROM pre WHERE tier_rank_min <> tier_rank) AS nondegen_before,
       (SELECT count(*) FROM compounds WHERE tier_rank_min <> tier_rank) AS nondegen_after;

\echo '=== movement by attesting source, for the response letter ==='
SELECT a.source, count(DISTINCT c.comp_id) AS moved
FROM compounds c JOIN pre p USING (comp_id)
JOIN per_source_license_attestation a ON a.comp_id = c.comp_id
WHERE c.license_tier IS DISTINCT FROM p.license_tier
GROUP BY 1 ORDER BY 2 DESC LIMIT 10;
