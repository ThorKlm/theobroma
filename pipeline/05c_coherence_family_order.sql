-- =====================================================================
-- THEOBROMA pipeline stage 05c: lineage coherence pass 2 (family<->order/class)
-- =====================================================================
-- Purpose : (a) BUILD family_lineage_ref = dominant coherent (order,class,
--           phylum,apg_clade) per family from single-organism compounds.
--           (b) Snap rows whose order/class disagrees with the family's own
--           coherent lineage. Genus & family untouched.
-- Inputs  : resolved_taxonomy, compounds(source_organism)
-- Output  : family_lineage_ref (reference table); corrected resolved_taxonomy
-- Run     : stage 05, after 05b. Canonical BUILDER of family_lineage_ref.
-- Verbatim from v135_coherence_fix2.sql (trailing COMMIT left commented).
-- =====================================================================
\set ON_ERROR_STOP on
BEGIN;

DROP TABLE IF EXISTS family_lineage_ref;
CREATE TABLE family_lineage_ref AS
WITH clean AS (
  SELECT lower(rt.family) AS f,
         rt.taxorder, rt.taxclass, rt.phylum, rt.apg_clade, count(*) AS n
  FROM resolved_taxonomy rt JOIN compounds c ON c.comp_id=rt.comp_id
  WHERE c.source_organism NOT LIKE '%;%'
    AND rt.family IS NOT NULL AND rt.family<>''
    AND rt.taxorder IS NOT NULL AND rt.taxclass IS NOT NULL
  GROUP BY 1,2,3,4,5
)
SELECT DISTINCT ON (f) f AS family_lc, taxorder, taxclass, phylum, apg_clade, n AS support
FROM clean ORDER BY f, n DESC;
ALTER TABLE family_lineage_ref ADD PRIMARY KEY (family_lc);

\echo '=== families with a coherent single-organism lineage ==='
SELECT count(*) AS families_in_ref FROM family_lineage_ref;

DROP TABLE IF EXISTS coherence_fix2_targets;
CREATE TEMP TABLE coherence_fix2_targets AS
SELECT rt.comp_id, rt.family,
       rt.taxorder  AS old_order,  r.taxorder  AS new_order,
       rt.taxclass  AS old_class,  r.taxclass  AS new_class,
       rt.phylum    AS old_phylum, r.phylum    AS new_phylum,
       rt.apg_clade AS old_clade,  r.apg_clade AS new_clade
FROM resolved_taxonomy rt
JOIN family_lineage_ref r ON r.family_lc = lower(rt.family)
WHERE rt.taxorder IS NOT NULL AND rt.taxclass IS NOT NULL
  AND (lower(rt.taxorder) <> lower(r.taxorder)
       OR lower(rt.taxclass) <> lower(r.taxclass));

\echo '=== targets (order/class stitched) ==='
SELECT count(*) AS targets FROM coherence_fix2_targets;

DROP TABLE IF EXISTS resolved_taxonomy_pre_coherence_fix2_20260802;
CREATE TABLE resolved_taxonomy_pre_coherence_fix2_20260802 AS
SELECT rt.* FROM resolved_taxonomy rt
JOIN coherence_fix2_targets t ON t.comp_id = rt.comp_id;

UPDATE resolved_taxonomy rt
SET taxorder  = t.new_order,
    taxclass  = t.new_class,
    phylum    = t.new_phylum,
    apg_clade = t.new_clade
FROM coherence_fix2_targets t
WHERE rt.comp_id = t.comp_id;

\echo '=== POST: family-order & family-class mismatches (target ~0) ==='
WITH clean AS (
  SELECT lower(rt.family) f, lower(rt.taxorder) o, lower(rt.taxclass) cl, count(*) n
  FROM resolved_taxonomy rt JOIN compounds c ON c.comp_id=rt.comp_id
  WHERE c.source_organism NOT LIKE '%;%' AND rt.family IS NOT NULL
    AND rt.taxorder IS NOT NULL AND rt.taxclass IS NOT NULL GROUP BY 1,2,3),
fd AS (SELECT DISTINCT ON (f) f, o AS dom_order, cl AS dom_class FROM clean ORDER BY f, n DESC)
SELECT count(*) FILTER (WHERE lower(rt.taxorder)<>fd.dom_order) AS family_order_mismatch,
       count(*) FILTER (WHERE lower(rt.taxclass)<>fd.dom_class) AS family_class_mismatch
FROM resolved_taxonomy rt JOIN fd ON fd.f=lower(rt.family)
WHERE rt.taxorder IS NOT NULL AND rt.taxclass IS NOT NULL;

\echo '=== review, then COMMIT; or ROLLBACK; ==='
-- COMMIT;
