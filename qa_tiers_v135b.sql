-- Post-materialization QA. Run after COMMIT to validate the tier partition.
\echo === tier distribution (every compound in exactly one tier or unclassified) ===
SELECT COALESCE(classification_tier::text,'unclassified') AS tier, count(*)
FROM compounds GROUP BY 1 ORDER BY 1;

\echo === does every classified compound have the right label source per tier? ===
SELECT classification_tier,
       count(*) FILTER (WHERE np_class IS NOT NULL AND np_class<>'') AS has_npclass,
       count(*) FILTER (WHERE inferred_class IS NOT NULL AND inferred_class<>'') AS has_inferred
FROM compounds WHERE classification_tier IS NOT NULL GROUP BY 1 ORDER BY 1;

\echo === placeholder contamination check (MUST be 0 before real deposit) ===
SELECT count(*) AS placeholder_rows FROM compounds WHERE inferred_class_source='PLACEHOLDER_dummy';

\echo === total classified vs corpus (coverage) ===
SELECT (SELECT count(*) FROM compounds) AS corpus,
       count(*) FILTER (WHERE classification_tier IS NOT NULL) AS classified,
       round(100.0*count(*) FILTER (WHERE classification_tier IS NOT NULL)/(SELECT count(*) FROM compounds),1) AS pct
FROM compounds;
