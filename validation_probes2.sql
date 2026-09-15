-- =====================================================================
-- THEOBROMA validation probes, ROUND 2 — plausibility & coherence.
-- Read-only. Run: psql -h localhost -U theobroma -d theobroma -f validation_probes2.sql
-- Round 1 confirmed structural integrity; this round checks whether the
-- data is SENSIBLE (coherent, plausible, completely rendered).
-- =====================================================================
\pset pager off
\timing off

\echo '========== 1. DEPRECATED TAXONOMY CLASSES (case-insensitive, the real count) =========='
SELECT lower(taxclass) AS taxclass_lc, count(*)
FROM resolved_taxonomy
WHERE lower(taxclass) IN ('liliopsida','magnoliopsida')
GROUP BY 1 ORDER BY 2 DESC;
\echo '-- scan for mixed-case / inconsistent kingdom tokens (data-quality smell) --'
SELECT kingdom, count(*) FROM resolved_taxonomy GROUP BY 1 ORDER BY 2 DESC;

\echo '========== 2. CLASSIFICATION NESTING COHERENCE =========='
\echo '-- INCOHERENT: has a class but NO pathway (should be ~0; a class implies a pathway) --'
SELECT count(*) AS class_without_pathway
FROM compounds WHERE np_class <> '' AND np_class IS NOT NULL
  AND (np_pathway = '' OR np_pathway IS NULL);
\echo '-- INCOHERENT: has a superclass but no pathway --'
SELECT count(*) AS superclass_without_pathway
FROM compounds WHERE np_superclass <> '' AND np_superclass IS NOT NULL
  AND (np_pathway = '' OR np_pathway IS NULL);
\echo '-- EXPECTED: has pathway but no class (NPClassifier confident coarsely only) --'
SELECT count(*) AS pathway_but_no_class
FROM compounds WHERE np_pathway <> '' AND np_pathway IS NOT NULL
  AND (np_class = '' OR np_class IS NULL);

\echo '========== 3. CLASS -> PATHWAY CONSISTENCY AT SCALE (residual scramble detector) =========='
\echo '-- classes that map to MANY distinct pathways across compounds (top offenders) --'
\echo '-- multi-value fields inflate this; a single-token class mapping to >3 pathways is suspect --'
SELECT np_class, count(DISTINCT np_pathway) AS distinct_pathways, count(*) AS n
FROM compounds
WHERE np_class <> '' AND np_class NOT LIKE '%|%'   -- single-label classes only
  AND np_pathway <> '' AND np_pathway NOT LIKE '%|%'
GROUP BY np_class
HAVING count(DISTINCT np_pathway) > 2
ORDER BY distinct_pathways DESC, n DESC
LIMIT 15;

\echo '========== 4. KINGDOM-CHEMISTRY PLAUSIBILITY =========='
\echo '-- plant-associated classes: kingdom breakdown (expect plant-dominated) --'
SELECT np_class,
       count(*) FILTER (WHERE rt.kingdom='plant')    AS plant,
       count(*) FILTER (WHERE rt.kingdom='bacteria') AS bacteria,
       count(*) FILTER (WHERE rt.kingdom='fungi')    AS fungi,
       count(*) FILTER (WHERE rt.kingdom='animal')   AS animal,
       count(*) AS total
FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id=c.comp_id
WHERE c.np_class IN ('Quassinoids','Cardenolides','Flavonoids','Iridoids monoterpenoids')
GROUP BY 1 ORDER BY total DESC;

\echo '========== 5. REGION / MAP VOCABULARY =========='
\echo '-- region distribution: how many have a usable region vs global/unresolved/null --'
SELECT COALESCE(NULLIF(region,''),'(null/empty)') AS region, count(*)
FROM compounds GROUP BY 1 ORDER BY 2 DESC LIMIT 20;

\echo '========== 6. STEREOISOMER FAMILY COHERENCE =========='
\echo '-- within a 14-char skeleton family (>1 member), do members share np_class? --'
\echo '-- families where members DISAGREE on class (smell; same molecule modulo stereo) --'
WITH fam AS (
  SELECT substring(inchikey,1,14) AS skel, count(*) AS members,
         count(DISTINCT np_class) AS distinct_classes
  FROM compounds
  WHERE np_class <> '' AND np_class IS NOT NULL
  GROUP BY 1 HAVING count(*) > 1
)
SELECT count(*) FILTER (WHERE distinct_classes > 1) AS families_with_class_disagreement,
       count(*) AS families_multimember
FROM fam;

\echo '========== 7. LICENSE-SOURCE PLAUSIBILITY =========='
\echo '-- compounds labelled CC0 that have NO CC0-tier attesting source (over-assignment?) --'
SELECT count(*) AS cc0_without_cc0_source
FROM compounds c
WHERE c.license_tier='CC0'
  AND NOT EXISTS (
    SELECT 1 FROM per_source_license_attestation a
    JOIN source_license_ref r ON r.src=lower(a.source)
    WHERE a.comp_id=c.comp_id AND r.tier_rank=0);

\echo '========== 8. EMPTY-FIELD DISPLAY POPULATIONS (UX) =========='
\echo '-- how many compound pages would show blank name / class / organism --'
SELECT count(*) FILTER (WHERE name IS NULL OR name='')            AS no_name,
       count(*) FILTER (WHERE np_class IS NULL OR np_class='')    AS no_class,
       count(*) FILTER (WHERE source_organism IS NULL OR source_organism='') AS no_organism
FROM compounds;

\echo '========== 9. FACET / AGGREGATE CONSISTENCY =========='
\echo '-- kingdom facet counts (primary) vs total --'
SELECT rt.kingdom, count(DISTINCT c.comp_id) AS n
FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id=c.comp_id
GROUP BY 1 ORDER BY 2 DESC;
\echo '-- do they sum to <= total compounds? (secondary kingdoms can push >; primary should = total) --'
SELECT count(*) AS total_compounds FROM compounds;

\echo '========== 10. PATHWAY DISTRIBUTION PLAUSIBILITY =========='
\echo '-- overall pathway distribution: sensible NP proportions? (alkaloids/terpenoids/shikimates large) --'
SELECT CASE WHEN np_pathway='' OR np_pathway IS NULL THEN '(none)' ELSE np_pathway END AS pathway,
       count(*)
FROM compounds GROUP BY 1 ORDER BY 2 DESC LIMIT 20;
