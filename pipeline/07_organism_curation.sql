-- =====================================================================
-- THEOBROMA pipeline stage 07: source-organism curation (LOTUS-enriched)
-- =====================================================================
-- Purpose : Build compounds.source_organism_curated from the union of
--           compound_taxonomy tokens and raw source_organism tokens,
--           normalized to clean "Genus species" binomials, deduplicated,
--           ordered resolved-genus-first. Surfaces LOTUS occurrences into
--           the curated display. Compounds with no clean binomial get NULL
--           (the template falls back to raw source_organism).
-- Inputs  : compound_taxonomy(comp_id, token), compounds(comp_id,
--           source_organism), resolved_taxonomy(comp_id, genus)
-- Output  : compounds.source_organism_curated
-- Run     : after stage 04 (resolved_taxonomy) is final.
-- Effect  : curated non-null 163,827 -> 335,530; Gramine 36 -> 49;
--           Morphine leads Papaver; 10,480 nulled-to-fallback, no data loss.
-- Supersedes: build_source_organism_curated_apply.sql (initial),
--             rebuild_source_organism_curated_apply.sql (pre-LOTUS rebuild).
--             Both archived. This is the canonical curation script.
-- Verbatim from rebuild_curated_lotus_merge_apply.sql (validated + applied).
-- =====================================================================
BEGIN;

CREATE TEMP TABLE _binoms AS
WITH raw_tokens AS (
  SELECT comp_id, token FROM compound_taxonomy WHERE token IS NOT NULL AND token<>''
  UNION ALL
  SELECT c.comp_id, trim(t) FROM compounds c
    CROSS JOIN LATERAL unnest(string_to_array(c.source_organism,'; ')) AS t
    WHERE c.source_organism IS NOT NULL AND c.source_organism<>''
),
extracted AS (
  SELECT comp_id, (regexp_match(token, '^([A-Za-z]+)[[:space:]]+([A-Za-z]+)')) AS m
  FROM raw_tokens
)
SELECT DISTINCT comp_id, initcap(m[1]) || ' ' || lower(m[2]) AS binom
FROM extracted
WHERE m IS NOT NULL
  AND lower(m[2]) NOT IN ('sp','spp','var','subsp','cf','aff','f','x','l','ssp','nov','gen')
  AND length(m[2]) > 1;

CREATE TEMP TABLE _curated AS
SELECT b.comp_id,
       string_agg(b.binom, '; ' ORDER BY
           CASE WHEN rt.genus IS NOT NULL AND lower(split_part(b.binom,' ',1))=lower(rt.genus) THEN 0 ELSE 1 END,
           b.binom) AS curated
FROM _binoms b
LEFT JOIN resolved_taxonomy rt ON rt.comp_id=b.comp_id
GROUP BY b.comp_id;

-- apply: set curated from rebuild; compounds absent from _curated get NULL (fall back to raw)
UPDATE compounds c
SET source_organism_curated = cu.curated
FROM _curated cu
WHERE cu.comp_id = c.comp_id
  AND c.source_organism_curated IS DISTINCT FROM cu.curated;

UPDATE compounds c
SET source_organism_curated = NULL
WHERE c.source_organism_curated IS NOT NULL
  AND NOT EXISTS (SELECT 1 FROM _curated cu WHERE cu.comp_id=c.comp_id);

\echo '=== POST-APPLY: Gramine (expect 49) ==='
SELECT (SELECT count(*) FROM unnest(string_to_array(source_organism_curated,'; '))) AS n,
       left(source_organism_curated,200) AS val
FROM compounds WHERE comp_id='THEO_0412143';

\echo '=== POST-APPLY: Morphine leads Papaver? ==='
SELECT split_part(source_organism_curated,'; ',1) AS leads FROM compounds WHERE comp_id='THEO_0511220';

\echo '=== POST-APPLY: curated non-null total ==='
SELECT count(*) FILTER (WHERE source_organism_curated IS NOT NULL) AS has_curated,
       count(*) FILTER (WHERE source_organism_curated IS NULL) AS is_null
FROM compounds;

\echo '=== residual check: any curated value still containing author-cite / sp. / joined junk? (expect 0) ==='
SELECT count(*) AS dirty_curated
FROM compounds
WHERE source_organism_curated ~ '( sp\.| spp\.| L\.| var\.|\$|\|)';

COMMIT;
\echo '=== COMMITTED ==='
