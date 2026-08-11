-- Best-available classification per compound, with explicit provenance.
-- Adds four NON-DESTRUCTIVE derived columns (originals np_* / inferred_* untouched):
--   effective_class      best label: curated/NPClassifier np_class (T1/T2) else inferred_class (T3)
--   effective_superclass ontology-derived from effective_class (internally consistent)
--   effective_pathway    ontology-derived from effective_class (internally consistent)
--   class_source         'curated' | 'npclassifier' | 'inferred_xgb' | 'placeholder'
--
-- effective_superclass/pathway are ALWAYS ontology-derived (via npc_class_parents) so the
-- displayed hierarchy is consistent (fixes the terpenoid-labelled-Alkaloid problem) WITHOUT
-- touching np_pathway (which the manuscript's Figure S3 numbers depend on).
--
-- Transactional. Idempotent (re-runnable). Verify, then COMMIT.
-- Single-label classes only for the ontology join; composite ($-joined) classes keep the
-- effective_class verbatim and get a pipe-joined union of parent pathways where resolvable.

BEGIN;

ALTER TABLE compounds ADD COLUMN IF NOT EXISTS effective_class      text;
ALTER TABLE compounds ADD COLUMN IF NOT EXISTS effective_superclass text;
ALTER TABLE compounds ADD COLUMN IF NOT EXISTS effective_pathway    text;
ALTER TABLE compounds ADD COLUMN IF NOT EXISTS class_source         text;

-- reset (idempotent)
UPDATE compounds SET effective_class=NULL, effective_superclass=NULL,
                     effective_pathway=NULL, class_source=NULL;

-- tidy: clear stale inferred_class_source on curated Tier 1 rows (116 leftovers)
UPDATE compounds SET inferred_class_source=NULL
WHERE classification_tier=1 AND inferred_class_source IS NOT NULL;

-- effective_class + class_source by tier
UPDATE compounds SET
  effective_class = CASE WHEN classification_tier IN (1,2) THEN np_class
                         WHEN classification_tier = 3 THEN inferred_class END,
  class_source = CASE WHEN classification_tier = 1 THEN 'curated'
                      WHEN classification_tier = 2 THEN 'npclassifier'
                      WHEN classification_tier = 3 AND inferred_class='PLACEHOLDER' THEN 'placeholder'
                      WHEN classification_tier = 3 THEN 'inferred_xgb' END;

-- ontology-derived superclass + pathway for SINGLE-label effective_class (authoritative,
-- internally consistent). npc_class_parents maps class -> (superclass, pathway).
UPDATE compounds c SET
  effective_superclass = p.superclass_name,
  effective_pathway    = p.pathway_name
FROM npc_class_parents p
WHERE p.class_name = c.effective_class
  AND c.effective_class IS NOT NULL AND c.effective_class NOT LIKE '%$%';

-- composite ($-joined) effective_class: keep the class string; derive a pipe-joined union
-- of parent pathways / superclasses across the atomic classes where resolvable.
WITH atoms AS (
  SELECT c.comp_id,
         trim(unnest(string_to_array(c.effective_class, '$'))) AS atom
  FROM compounds c
  WHERE c.effective_class LIKE '%$%'
), mapped AS (
  SELECT a.comp_id,
         string_agg(DISTINCT p.pathway_name,    ' | ' ORDER BY p.pathway_name)    AS paths,
         string_agg(DISTINCT p.superclass_name, ' | ' ORDER BY p.superclass_name) AS supers
  FROM atoms a LEFT JOIN npc_class_parents p ON p.class_name = a.atom
  GROUP BY a.comp_id
)
UPDATE compounds c SET
  effective_pathway    = COALESCE(m.paths,  c.effective_pathway),
  effective_superclass = COALESCE(m.supers, c.effective_superclass)
FROM mapped m WHERE m.comp_id = c.comp_id;

-- verification
SELECT class_source,
       count(*) AS n,
       count(*) FILTER (WHERE effective_class IS NOT NULL AND effective_class<>'') AS has_class,
       count(*) FILTER (WHERE effective_pathway IS NOT NULL AND effective_pathway<>'') AS has_pathway,
       count(*) FILTER (WHERE effective_superclass IS NOT NULL AND effective_superclass<>'') AS has_superclass
FROM compounds GROUP BY class_source ORDER BY n DESC NULLS LAST;

-- consistency win: how many now have effective_pathway != np_pathway (the corrected ones)
SELECT count(*) AS pathway_corrected
FROM compounds
WHERE effective_pathway IS NOT NULL AND np_pathway IS NOT NULL
  AND effective_pathway <> np_pathway;

-- If correct: COMMIT;  else ROLLBACK;
