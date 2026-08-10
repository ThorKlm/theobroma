-- LOTUS-enriched source_organism_curated rebuild (DRY RUN, rolls back).
-- Surfaces the fuller organism set from compound_taxonomy (LOTUS occurrences),
-- normalized to clean deduped binomials, resolved-genus-first ordering.
-- Only compounds.source_organism_curated changes. Fallback: compounds with no
-- clean binomial keep NULL (template shows raw source_organism instead).
BEGIN;

-- normalized (comp_id, binomial) from the union of compound_taxonomy tokens + raw source_organism tokens
CREATE TEMP TABLE _binoms AS
WITH raw_tokens AS (
  SELECT comp_id, token FROM compound_taxonomy WHERE token IS NOT NULL AND token<>''
  UNION ALL
  SELECT c.comp_id, trim(t) FROM compounds c
    CROSS JOIN LATERAL unnest(string_to_array(c.source_organism,'; ')) AS t
    WHERE c.source_organism IS NOT NULL AND c.source_organism<>''
),
extracted AS (
  SELECT comp_id,
         (regexp_match(token, '^([A-Za-z]+)[[:space:]]+([A-Za-z-]+)')) AS m
  FROM raw_tokens
)
SELECT DISTINCT comp_id,
       initcap(m[1]) || ' ' || lower(m[2]) AS binom
FROM extracted
WHERE m IS NOT NULL
  AND lower(m[2]) NOT IN ('sp','spp','var','subsp','cf','aff','f','x','l','ssp','nov','gen')
  AND length(m[2]) > 1;

-- ordering: resolved genus first (taxonomically primary), then alphabetical
CREATE TEMP TABLE _curated AS
SELECT b.comp_id,
       string_agg(b.binom, '; ' ORDER BY
           CASE WHEN rt.genus IS NOT NULL AND lower(split_part(b.binom,' ',1))=lower(rt.genus) THEN 0 ELSE 1 END,
           b.binom) AS curated
FROM _binoms b
LEFT JOIN resolved_taxonomy rt ON rt.comp_id=b.comp_id
GROUP BY b.comp_id;

\echo '=== Gramine before/after ==='
SELECT 'BEFORE' w, (SELECT count(*) FROM unnest(string_to_array(source_organism_curated,'; '))) n, left(source_organism_curated,220) val FROM compounds WHERE comp_id='THEO_0412143'
UNION ALL
SELECT 'AFTER', (SELECT count(*) FROM unnest(string_to_array(curated,'; '))), left(curated,220) FROM _curated WHERE comp_id='THEO_0412143';

\echo '=== Morphine: does resolved genus (Papaver) still lead? ==='
SELECT split_part(curated,'; ',1) AS leads FROM _curated WHERE comp_id='THEO_0511220';

\echo '=== magnitude: compounds whose curated organism count grows / shrinks / same ==='
SELECT
  count(*) FILTER (WHERE new_n > old_n) AS grew,
  count(*) FILTER (WHERE new_n < old_n) AS shrank,
  count(*) FILTER (WHERE new_n = old_n) AS same,
  round(avg(new_n - old_n),2) AS avg_delta
FROM (
  SELECT c.comp_id,
    (SELECT count(*) FROM unnest(string_to_array(c.source_organism_curated,'; '))) AS old_n,
    COALESCE((SELECT count(*) FROM unnest(string_to_array(cu.curated,'; '))),0) AS new_n
  FROM compounds c LEFT JOIN _curated cu ON cu.comp_id=c.comp_id
) d;

\echo '=== how many compounds would get NULL curated (no clean binomial -> fall back to raw) ==='
SELECT count(*) AS would_be_null FROM compounds c WHERE NOT EXISTS (SELECT 1 FROM _curated cu WHERE cu.comp_id=c.comp_id);

\echo '=== 5 random before/after samples ==='
SELECT c.comp_id,
       left(c.source_organism_curated,90) AS before,
       left(cu.curated,90) AS after
FROM compounds c JOIN _curated cu ON cu.comp_id=c.comp_id
WHERE c.source_organism_curated IS NOT NULL
ORDER BY random() LIMIT 5;

ROLLBACK;
\echo '=== ROLLED BACK. Review, then APPLY. ==='
