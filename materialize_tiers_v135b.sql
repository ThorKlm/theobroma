-- Three provenance tiers for v1.35 classification (no confidence thresholds).
-- Precedence: curated (source-inherited) > NPClassifier (tool re-run) > inferred (XGBoost).
-- Adds classification_tier (1/2/3, NULL = unclassified) and the effective class columns.
-- Idempotent w.r.t. the tier column; run inside a transaction, verify, then COMMIT.
--
-- Prereqs:
--   * curated np_class already restored (Tier 1 source-inherited labels).
--   * npc_reclassified holds the NPClassifier tool re-run (Tier 2 source).
--   * inferred_v135b loaded into a staging table inferred_new(comp_id, inferred_class,
--     inferred_confidence) from the vast.ai inference output (Tier 3).

BEGIN;

-- 0. stage the new inferred predictions (load via \copy before running this file):
--    CREATE TABLE inferred_new(comp_id text PRIMARY KEY, inferred_class text, inferred_confidence real);
--    \copy inferred_new FROM 'inferred_v135b.tsv' WITH (FORMAT csv, DELIMITER E'\t', HEADER);

-- 1. reset the tier column
ALTER TABLE compounds ADD COLUMN IF NOT EXISTS classification_tier smallint;
UPDATE compounds SET classification_tier = NULL;

-- 2. refresh the inferred columns from the NEW model, but ONLY where we will actually
--    use them (residual). First clear old inferred_* so stale v1.34 predictions do not leak.
UPDATE compounds
SET inferred_class = NULL, inferred_class_source = NULL, inferred_confidence = NULL
WHERE inferred_class IS NOT NULL;

-- 3. Tier 1 (curated): compound already has a source-inherited np_class.
UPDATE compounds
SET classification_tier = 1
WHERE np_class IS NOT NULL AND np_class <> '';

-- 4. Tier 2 (NPClassifier): no curated class, but the NPClassifier re-run assigned one.
UPDATE compounds c
SET np_class = r.np_class,
    np_superclass = r.np_superclass,
    np_pathway = r.np_pathway,
    classification_tier = 2
FROM npc_reclassified r
WHERE c.comp_id = r.comp_id
  AND (c.np_class IS NULL OR c.np_class = '')
  AND r.np_class IS NOT NULL AND r.np_class <> '';

-- 5. Tier 3 (inferred): still no class, but the retrained XGBoost predicted one.
UPDATE compounds c
SET inferred_class = n.inferred_class,
    inferred_confidence = n.inferred_confidence,
    inferred_class_source = 'inferred_xgb_v135b',
    classification_tier = 3
FROM inferred_new n
WHERE c.comp_id = n.comp_id
  AND (c.np_class IS NULL OR c.np_class = '')
  AND n.inferred_class IS NOT NULL AND n.inferred_class <> '';

-- 6. verification (compare against expectations before COMMIT)
SELECT classification_tier,
       count(*) AS n,
       count(*) FILTER (WHERE np_class IS NOT NULL AND np_class<>'') AS has_npclass,
       count(*) FILTER (WHERE inferred_class IS NOT NULL AND inferred_class<>'') AS has_inferred
FROM compounds GROUP BY classification_tier ORDER BY classification_tier NULLS LAST;

-- Expected shape:
--   tier 1 (curated)      ~391k  (all has_npclass)
--   tier 2 (NPClassifier) = npc_reclassified coverage of the non-curated set
--   tier 3 (inferred)     = residual the model could score
--   NULL (unclassified)   = the remainder
-- If correct: COMMIT;  else: ROLLBACK;
