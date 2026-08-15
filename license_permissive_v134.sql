-- =====================================================================
-- LICENSE RE-RESOLUTION -- OPTION A: MOST-PERMISSIVE-WINS.
-- Rationale: a molecular structure is a fact independently available from
-- its most permissive source; a user may obtain it there. So each compound
-- takes the LEAST restrictive license among its attested sources.
--
-- Tier ranks: CC0=0, CC BY 4.0=1, CC BY-NC 4.0=2, CC BY-NC-SA 4.0=3,
--             CC BY-NC-ND 4.0=4, Unspecified=5.
-- "Unspecified" (5) means license UNKNOWN, not a real permission, so it is
-- treated as NO SIGNAL: a compound takes the MIN tier_rank among its
-- KNOWN sources (0-4); only if ALL its sources are Unspecified does it
-- resolve to Unspecified.
--
-- Staged; archives current labels; transaction OPEN; inspect; COMMIT/ROLLBACK.
-- =====================================================================
\set ON_ERROR_STOP on
\timing on
BEGIN;

-- 0. Archive current (most-restrictive) labels for rollback/audit.
DROP TABLE IF EXISTS compounds_license_pre_permissive_20260802;
CREATE TABLE compounds_license_pre_permissive_20260802 AS
SELECT comp_id, license_tier, tier_rank FROM compounds;
\echo '=== archived current labels ==='
SELECT count(*) FROM compounds_license_pre_permissive_20260802;

-- 1. Resolve most-permissive-wins.
--    known = min tier_rank among sources with tier_rank <= 4.
--    if no known source (all Unspecified) -> 5.
DROP TABLE IF EXISTS license_permissive_v134;
CREATE TEMP TABLE license_permissive_v134 AS
WITH att AS (
  SELECT a.comp_id, r.tier_rank
  FROM per_source_license_attestation a
  JOIN source_license_ref r ON r.src = lower(a.source)
)
SELECT comp_id,
       COALESCE(min(tier_rank) FILTER (WHERE tier_rank <= 4), 5) AS final_rank
FROM att GROUP BY comp_id;
CREATE INDEX lpv134_cid ON license_permissive_v134 (comp_id);
\echo '=== resolved rows ==='
SELECT count(*) FROM license_permissive_v134;

\echo '=== PROJECTED most-permissive distribution ==='
SELECT final_rank,
  CASE final_rank WHEN 0 THEN 'CC0' WHEN 1 THEN 'CC BY 4.0' WHEN 2 THEN 'CC BY-NC 4.0'
    WHEN 3 THEN 'CC BY-NC-SA 4.0' WHEN 4 THEN 'CC BY-NC-ND 4.0' ELSE 'Unspecified' END AS tier,
  count(*) AS compounds, round(100.0*count(*)/sum(count(*)) over (),2) AS pct
FROM license_permissive_v134 GROUP BY 1 ORDER BY 1;

\echo '=== commercial / non-commercial / unspecified (projected) ==='
SELECT CASE WHEN final_rank<=1 THEN 'commercial' WHEN final_rank<=4 THEN 'non-commercial'
            ELSE 'unspecified' END AS bucket, count(*),
       round(100.0*count(*)/sum(count(*)) over (),2) pct
FROM license_permissive_v134 GROUP BY 1 ORDER BY 2 DESC;

\echo '=== curcumin THEO_0854403 under most-permissive (LOTUS CC0 present -> expect CC0) ==='
SELECT final_rank FROM license_permissive_v134 WHERE comp_id='THEO_0854403';

-- 2. Apply.
UPDATE compounds c
SET tier_rank = lp.final_rank,
    license_tier = CASE lp.final_rank
        WHEN 0 THEN 'CC0' WHEN 1 THEN 'CC BY 4.0' WHEN 2 THEN 'CC BY-NC 4.0'
        WHEN 3 THEN 'CC BY-NC-SA 4.0' WHEN 4 THEN 'CC BY-NC-ND 4.0'
        ELSE 'Unspecified' END
FROM license_permissive_v134 lp
WHERE c.comp_id = lp.comp_id;

\echo '=== NEW live distribution after apply ==='
SELECT tier_rank, license_tier, count(*) FROM compounds GROUP BY 1,2 ORDER BY 1;

\echo '=== curcumin family live (expect CC0) ==='
SELECT comp_id, license_tier, tier_rank FROM compounds WHERE comp_id='THEO_0854403';

\echo '======================================================================'
\echo 'INSPECT projected vs applied + curcumin=CC0. If good: COMMIT; else ROLLBACK;'
\echo '======================================================================'
-- COMMIT;
