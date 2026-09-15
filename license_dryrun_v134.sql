-- =====================================================================
-- LICENSE RE-RESOLUTION DRY-RUN on the v1.34 base (READ-ONLY, no changes).
-- Computes what the corrected most-restrictive-wins distribution WOULD be,
-- so it can be inspected before any UPDATE. Nothing is written.
--
-- tier_rank: CC0=0, CC BY 4.0=1, CC BY-NC 4.0=2, CC BY-NC-SA 4.0=3,
--            CC BY-NC-ND 4.0=4, Unspecified=5. Most restrictive = max rank.
-- Independent-CC0 exception: if a compound has a genuine CC0 source, it may
-- resolve to CC0 even amid NC sources.
-- =====================================================================
\timing on

\echo '=== 0. sanity: attestation + reference shapes ==='
SELECT 'attestation_rows' k, count(*) v FROM per_source_license_attestation
UNION ALL SELECT 'attestation_distinct_comp', count(DISTINCT comp_id) FROM per_source_license_attestation
UNION ALL SELECT 'source_license_ref_rows', count(*) FROM source_license_ref
UNION ALL SELECT 'live_compounds', count(*) FROM compounds;

\echo '=== 1. attestation sources that do NOT match source_license_ref (would default) ==='
SELECT DISTINCT lower(a.source) AS unmatched_source
FROM per_source_license_attestation a
LEFT JOIN source_license_ref r ON r.src = lower(a.source)
WHERE r.src IS NULL
ORDER BY 1 LIMIT 40;

\echo '=== 2. how many attestation rows match the ref vs default ==='
SELECT
  count(*) FILTER (WHERE r.src IS NOT NULL) AS matched_rows,
  count(*) FILTER (WHERE r.src IS NULL)     AS unmatched_rows
FROM per_source_license_attestation a
LEFT JOIN source_license_ref r ON r.src = lower(a.source);

\echo '=== 3. COMPUTE the resolved distribution (most-restrictive-wins + independent-CC0) ==='
WITH att AS (
  SELECT a.comp_id, r.tier_rank
  FROM per_source_license_attestation a
  JOIN source_license_ref r ON r.src = lower(a.source)
),
res AS (
  SELECT comp_id, max(tier_rank) AS worst_rank, bool_or(tier_rank = 0) AS has_cc0
  FROM att GROUP BY comp_id
),
final AS (
  SELECT comp_id,
         CASE WHEN worst_rank >= 2 AND has_cc0 THEN 0 ELSE worst_rank END AS final_rank
  FROM res
)
SELECT final_rank,
       CASE final_rank
         WHEN 0 THEN 'CC0' WHEN 1 THEN 'CC BY 4.0' WHEN 2 THEN 'CC BY-NC 4.0'
         WHEN 3 THEN 'CC BY-NC-SA 4.0' WHEN 4 THEN 'CC BY-NC-ND 4.0'
         ELSE 'Unspecified' END AS license_tier,
       count(*) AS compounds
FROM final GROUP BY 1 ORDER BY 1;

\echo '=== 4. compounds with NO attestation match at all (would become Unspecified) ==='
SELECT count(*) AS compounds_without_attestation
FROM compounds c
WHERE NOT EXISTS (
  SELECT 1 FROM per_source_license_attestation a
  JOIN source_license_ref r ON r.src = lower(a.source)
  WHERE a.comp_id = c.comp_id);

\echo '=== 5. full projected distribution INCLUDING the no-attestation -> Unspecified bucket ==='
WITH att AS (
  SELECT a.comp_id, r.tier_rank
  FROM per_source_license_attestation a
  JOIN source_license_ref r ON r.src = lower(a.source)
),
res AS (
  SELECT comp_id, max(tier_rank) AS worst_rank, bool_or(tier_rank=0) AS has_cc0
  FROM att GROUP BY comp_id
),
resolved AS (
  SELECT c.comp_id,
    COALESCE(
      CASE WHEN r.worst_rank >= 2 AND r.has_cc0 THEN 0 ELSE r.worst_rank END,
      5) AS final_rank
  FROM compounds c LEFT JOIN res r ON r.comp_id = c.comp_id
)
SELECT final_rank,
       CASE final_rank
         WHEN 0 THEN 'CC0' WHEN 1 THEN 'CC BY 4.0' WHEN 2 THEN 'CC BY-NC 4.0'
         WHEN 3 THEN 'CC BY-NC-SA 4.0' WHEN 4 THEN 'CC BY-NC-ND 4.0'
         ELSE 'Unspecified' END AS license_tier,
       count(*) AS compounds,
       round(100.0*count(*)/sum(count(*)) over (), 2) AS pct
FROM resolved GROUP BY 1 ORDER BY 1;

\echo '=== 6. commercial / non-commercial / unspecified split (projected) ==='
WITH att AS (
  SELECT a.comp_id, r.tier_rank FROM per_source_license_attestation a
  JOIN source_license_ref r ON r.src = lower(a.source)),
res AS (SELECT comp_id, max(tier_rank) worst_rank, bool_or(tier_rank=0) has_cc0 FROM att GROUP BY comp_id),
resolved AS (
  SELECT c.comp_id, COALESCE(CASE WHEN r.worst_rank>=2 AND r.has_cc0 THEN 0 ELSE r.worst_rank END,5) AS fr
  FROM compounds c LEFT JOIN res r ON r.comp_id=c.comp_id)
SELECT CASE WHEN fr<=1 THEN 'commercial' WHEN fr<=4 THEN 'non-commercial' ELSE 'unspecified' END AS bucket,
       count(*), round(100.0*count(*)/sum(count(*)) over (),2) AS pct
FROM resolved GROUP BY 1 ORDER BY 2 DESC;

\echo '=== 7. current (OLD) live distribution for comparison ==='
SELECT license_tier, count(*) FROM compounds GROUP BY 1 ORDER BY 2 DESC;

\echo '=== 8. spot-check: curcumin family projected vs current ==='
WITH att AS (
  SELECT a.comp_id, r.tier_rank FROM per_source_license_attestation a
  JOIN source_license_ref r ON r.src=lower(a.source)),
res AS (SELECT comp_id, max(tier_rank) wr, bool_or(tier_rank=0) cc0 FROM att GROUP BY comp_id)
SELECT c.comp_id, c.license_tier AS current_tier,
       CASE COALESCE(CASE WHEN r.wr>=2 AND r.cc0 THEN 0 ELSE r.wr END,5)
         WHEN 0 THEN 'CC0' WHEN 1 THEN 'CC BY 4.0' WHEN 2 THEN 'CC BY-NC 4.0'
         WHEN 3 THEN 'CC BY-NC-SA 4.0' WHEN 4 THEN 'CC BY-NC-ND 4.0' ELSE 'Unspecified' END AS projected_tier
FROM compounds c LEFT JOIN res r ON r.comp_id=c.comp_id
WHERE c.inchikey LIKE 'VFLDPWHFBUODDF%' ORDER BY c.comp_id LIMIT 6;
