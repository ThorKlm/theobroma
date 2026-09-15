-- Recompute compounds.tier_rank and tier_rank_min from the attestation table.
-- Run after reconcile_attestations.sql and apply_license_resolver.sql.
-- Interactive: review the invariants, then COMMIT by hand.
\timing on
BEGIN;

CREATE TABLE IF NOT EXISTS compounds_ranks_pre_20260915 AS
  SELECT comp_id, license_tier, tier_rank, tier_rank_min FROM compounds;

UPDATE compounds c
SET tier_rank_min = s.min_rank
FROM (SELECT a.comp_id, min(r.tier_rank) AS min_rank
      FROM per_source_license_attestation a
      JOIN source_license_ref r ON lower(r.src) = lower(a.source)
      GROUP BY 1) s
WHERE s.comp_id = c.comp_id AND c.tier_rank_min IS DISTINCT FROM s.min_rank;

UPDATE compounds c
SET tier_rank = s.max_rank
FROM (SELECT a.comp_id, max(r.tier_rank) AS max_rank
      FROM per_source_license_attestation a
      JOIN source_license_ref r ON lower(r.src) = lower(a.source)
      GROUP BY 1) s
WHERE s.comp_id = c.comp_id AND c.tier_rank IS DISTINCT FROM s.max_rank;

SELECT count(*) FILTER (WHERE tier_rank_min > tier_rank) AS w03_inversions,
       count(*) FILTER (WHERE tier_rank_min IS NULL) AS null_min,
       count(*) FILTER (WHERE tier_rank_min <> tier_rank) AS non_degenerate,
       count(*) FILTER (WHERE tier_rank_min = 0) AS least_public_domain,
       count(*) FILTER (WHERE tier_rank_min = 1) AS least_permissive
FROM compounds;

-- COMMIT;  -- type by hand after reviewing
