-- =====================================================================
-- THEOBROMA pipeline stage 05b: lineage coherence pass 1 (genus<->family)
-- =====================================================================
-- Purpose : (a) BUILD genus_lineage_ref = the dominant COHERENT lineage tuple
--           per genus, taken as a whole unit from single-organism compounds
--           (coherent by construction). (b) Snap incoherent genera (family
--           disagrees with the genus's own coherent family) to that lineage.
--           Genus kept; family/order/class/phylum/apg_clade overwritten as a unit.
-- Inputs  : resolved_taxonomy, compounds(source_organism)
-- Output  : genus_lineage_ref (reference table); corrected resolved_taxonomy
-- Run     : stage 05, after backfill (05a) and resolver (04). This is also the
--           canonical BUILDER of genus_lineage_ref used by the app + reports.
-- Verbatim from v135_coherence_fix.sql (trailing COMMIT left commented as in source).
-- =====================================================================
\set ON_ERROR_STOP on
BEGIN;

DROP TABLE IF EXISTS genus_lineage_ref;
CREATE TABLE genus_lineage_ref AS
WITH clean AS (
  SELECT lower(rt.genus) AS g,
         rt.family, rt.taxorder, rt.taxclass, rt.phylum, rt.apg_clade,
         count(*) AS n
  FROM resolved_taxonomy rt
  JOIN compounds c ON c.comp_id = rt.comp_id
  WHERE c.source_organism NOT LIKE '%;%'
    AND rt.genus  IS NOT NULL AND rt.genus  <> ''
    AND rt.family IS NOT NULL AND rt.family <> ''
  GROUP BY 1,2,3,4,5,6
)
SELECT DISTINCT ON (g)
       g AS genus_lc, family, taxorder, taxclass, phylum, apg_clade, n AS support
FROM clean
ORDER BY g, n DESC;
ALTER TABLE genus_lineage_ref ADD PRIMARY KEY (genus_lc);
CREATE INDEX glr_genus ON genus_lineage_ref (genus_lc);

\echo '=== ref size (distinct genera with a coherent single-organism lineage) ==='
SELECT count(*) AS genera_in_ref FROM genus_lineage_ref;

DROP TABLE IF EXISTS coherence_fix_targets;
CREATE TEMP TABLE coherence_fix_targets AS
SELECT rt.comp_id, rt.genus,
       rt.family    AS old_family,    r.family    AS new_family,
       rt.taxorder  AS old_order,     r.taxorder  AS new_order,
       rt.taxclass  AS old_class,     r.taxclass  AS new_class,
       rt.phylum    AS old_phylum,    r.phylum    AS new_phylum,
       rt.apg_clade AS old_clade,     r.apg_clade AS new_clade
FROM resolved_taxonomy rt
JOIN genus_lineage_ref r ON r.genus_lc = lower(rt.genus)
WHERE rt.family IS NOT NULL AND rt.family <> ''
  AND lower(rt.family) <> lower(r.family);

\echo '=== targets to fix ==='
SELECT count(*) AS targets FROM coherence_fix_targets;

DROP TABLE IF EXISTS resolved_taxonomy_pre_coherence_fix_20260802;
CREATE TABLE resolved_taxonomy_pre_coherence_fix_20260802 AS
SELECT rt.* FROM resolved_taxonomy rt
JOIN coherence_fix_targets t ON t.comp_id = rt.comp_id;

UPDATE resolved_taxonomy rt
SET family    = t.new_family,
    taxorder  = t.new_order,
    taxclass  = t.new_class,
    phylum    = t.new_phylum,
    apg_clade = t.new_clade
FROM coherence_fix_targets t
WHERE rt.comp_id = t.comp_id;

\echo '=== POST-FIX incoherent count (target: near 0) ==='
WITH genus_family AS (
  SELECT lower(genus) g, lower(family) f, count(*) n FROM resolved_taxonomy
  WHERE genus IS NOT NULL AND genus<>'' AND family IS NOT NULL AND family<>'' GROUP BY 1,2),
dom AS (SELECT DISTINCT ON (g) g, f FROM genus_family ORDER BY g, n DESC)
SELECT count(*) AS still_incoherent
FROM resolved_taxonomy rt JOIN dom d ON d.g=lower(rt.genus)
WHERE rt.family IS NOT NULL AND rt.family<>'' AND lower(rt.family)<>d.f;

\echo '=== review, then COMMIT; or ROLLBACK; ==='
-- COMMIT;
