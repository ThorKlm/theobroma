-- =====================================================================
-- THEOBROMA pipeline stage 09: license re-resolution (canonical)
-- =====================================================================
-- Purpose : Resolve per-compound license_tier by pure most-restrictive-wins
--           across attested sources (final_rank = max source tier_rank).
--           The erroneous CC0 override is discarded.
-- Inputs  : per_source_license_attestation, source_license_ref
-- Output  : compounds.license_tier + compounds.tier_rank
-- Run     : after compounds + attestations are loaded.
-- Live dist (v1.35): CC0 1,013,320; CC BY-NC 84,956; Unspecified 27,051;
--           CC BY-NC-ND 3,758; CC BY 4.0 3,720. Commercial (CC0+CCBY) 1,017,040.
-- NOTE    : original staged with operator COMMIT; trailing COMMIT left commented
--           as in source. Pair with 09b (reconcile_attestations) which aligns
--           the attestation rows to the ref WITHOUT changing resolved tiers.
-- Verbatim from license_apply_v134.sql.
-- =====================================================================
\set ON_ERROR_STOP on
\timing on
BEGIN;

-- 0. Archive current license labels for rollback / audit.
DROP TABLE IF EXISTS compounds_license_pre_fix_20260802;
CREATE TABLE compounds_license_pre_fix_20260802 AS
SELECT comp_id, license_tier FROM compounds;
\echo '=== archived current license labels ==='
SELECT count(*) FROM compounds_license_pre_fix_20260802;

-- 1. Resolve pure most-restrictive-wins per compound.
DROP TABLE IF EXISTS license_resolved_v134;
CREATE TEMP TABLE license_resolved_v134 AS
WITH att AS (
  SELECT a.comp_id, r.tier_rank
  FROM per_source_license_attestation a
  JOIN source_license_ref r ON r.src = lower(a.source)
)
SELECT comp_id, max(tier_rank) AS final_rank
FROM att GROUP BY comp_id;
CREATE INDEX lrv134_cid ON license_resolved_v134 (comp_id);

\echo '=== resolved rows ==='
SELECT count(*) FROM license_resolved_v134;

-- 2. Apply. Compounds with no attestation -> Unspecified (5). Guard tier_rank col.
DO $$
DECLARE has_tr boolean;
BEGIN
  SELECT EXISTS(SELECT 1 FROM information_schema.columns
                WHERE table_name='compounds' AND column_name='tier_rank') INTO has_tr;
  IF NOT has_tr THEN
    EXECUTE 'ALTER TABLE compounds ADD COLUMN tier_rank int';
  END IF;
END$$;

UPDATE compounds c
SET tier_rank = COALESCE(lr.final_rank, 5),
    license_tier = CASE COALESCE(lr.final_rank, 5)
        WHEN 0 THEN 'CC0' WHEN 1 THEN 'CC BY 4.0' WHEN 2 THEN 'CC BY-NC 4.0'
        WHEN 3 THEN 'CC BY-NC-SA 4.0' WHEN 4 THEN 'CC BY-NC-ND 4.0'
        ELSE 'Unspecified' END
FROM (SELECT c2.comp_id, lr.final_rank
      FROM compounds c2 LEFT JOIN license_resolved_v134 lr ON lr.comp_id=c2.comp_id) lr
WHERE c.comp_id = lr.comp_id;

\echo '=== NEW live license distribution ==='
SELECT tier_rank, license_tier, count(*) FROM compounds GROUP BY 1,2 ORDER BY 1;
\echo '=== commercial / non-commercial / unspecified ==='
SELECT CASE WHEN tier_rank<=1 THEN 'commercial' WHEN tier_rank<=4 THEN 'non-commercial'
            ELSE 'unspecified' END AS bucket, count(*),
       round(100.0*count(*)/sum(count(*)) over (),2) pct
FROM compounds GROUP BY 1 ORDER BY 2 DESC;
\echo '=== curcumin family (should be Unspecified now) ==='
SELECT comp_id, license_tier, tier_rank FROM compounds
WHERE inchikey LIKE 'VFLDPWHFBUODDF%' ORDER BY comp_id LIMIT 6;
\echo '=== changed vs archived ==='
SELECT count(*) AS changed
FROM compounds c JOIN compounds_license_pre_fix_20260802 a ON a.comp_id=c.comp_id
WHERE c.license_tier <> a.license_tier;

\echo '======================================================================'
\echo 'INSPECT distribution + curcumin=Unspecified. If good: COMMIT; else ROLLBACK;'
\echo '======================================================================'
-- COMMIT;
