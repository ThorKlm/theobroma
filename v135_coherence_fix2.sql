-- =====================================================================
-- THEOBROMA v1.35 : coherence fix PASS 2 -- family <-> order/class stitching.
-- Complements pass 1 (genus<->family). Here genus & family are already
-- coherent, but order/class/phylum/clade were stitched from another organism.
-- Fix: snap order/class/phylum/apg_clade to the family's OWN coherent lineage
-- (dominant tuple among single-organism compounds). Genus & family untouched.
-- Staged in a transaction; previews + verify; COMMIT at end (commented).
-- =====================================================================
\set ON_ERROR_STOP on
BEGIN;

-- ---------------------------------------------------------------------
-- 0. family -> coherent (taxorder, taxclass, phylum, apg_clade) from
--    single-organism compounds (coherent by construction), most-supported.
-- ---------------------------------------------------------------------
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

-- ---------------------------------------------------------------------
-- 1. Targets: order OR class disagrees with the family's coherent lineage.
-- ---------------------------------------------------------------------
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

\echo '=== targets (order/class stitched; expect ~ a few thousand) ==='
SELECT count(*) AS targets FROM coherence_fix2_targets;

\echo '=== PREVIEW: 15 before/after (family kept, order/class corrected) ==='
SELECT comp_id, family,
       old_class||' -> '||coalesce(new_class,'(null)') AS class_change,
       old_order||' -> '||coalesce(new_order,'(null)') AS order_change
FROM coherence_fix2_targets ORDER BY family LIMIT 15;

\echo '=== PREVIEW: the 4 Adiantum stragglers specifically ==='
SELECT comp_id, family, old_class, new_class, old_order, new_order
FROM coherence_fix2_targets WHERE comp_id IN ('THEO_0207201','THEO_0984573','THEO_0984882','THEO_0111153');

-- ---------------------------------------------------------------------
-- 2. Archive pre-fix rows.
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS resolved_taxonomy_pre_coherence_fix2_20260802;
CREATE TABLE resolved_taxonomy_pre_coherence_fix2_20260802 AS
SELECT rt.* FROM resolved_taxonomy rt
JOIN coherence_fix2_targets t ON t.comp_id = rt.comp_id;
\echo '=== archived ==='
SELECT count(*) AS archived FROM resolved_taxonomy_pre_coherence_fix2_20260802;

-- ---------------------------------------------------------------------
-- 3. Apply: correct order/class/phylum/clade to the family's coherent
--    lineage. Genus and family unchanged.
-- ---------------------------------------------------------------------
UPDATE resolved_taxonomy rt
SET taxorder  = t.new_order,
    taxclass  = t.new_class,
    phylum    = t.new_phylum,
    apg_clade = t.new_clade
FROM coherence_fix2_targets t
WHERE rt.comp_id = t.comp_id;
\echo '=== updated ==='
SELECT count(*) AS updated FROM coherence_fix2_targets;

-- ---------------------------------------------------------------------
-- 4. VERIFY both dimensions now near zero.
-- ---------------------------------------------------------------------
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

\echo '=== POST: genus-family mismatch (pass-1 dimension, should still be ~163) ==='
WITH gf AS (SELECT lower(genus) g, lower(family) f, count(*) n FROM resolved_taxonomy
  WHERE genus IS NOT NULL AND genus<>'' AND family IS NOT NULL AND family<>'' GROUP BY 1,2),
dom AS (SELECT DISTINCT ON (g) g, f FROM gf ORDER BY g, n DESC)
SELECT count(*) AS genus_family_mismatch FROM resolved_taxonomy rt JOIN dom d ON d.g=lower(rt.genus)
WHERE rt.family IS NOT NULL AND lower(rt.family)<>d.f;

\echo '=== VERIFY Adiantum now fully coherent ==='
SELECT family, taxorder, taxclass, apg_clade, count(*) n
FROM resolved_taxonomy WHERE lower(genus)='adiantum' GROUP BY 1,2,3,4 ORDER BY 5 DESC;

\echo '=== kingdom counts unchanged ==='
SELECT kingdom, count(*) FROM resolved_taxonomy GROUP BY 1 ORDER BY 2 DESC;

\echo '=== apg_clade distribution (final) ==='
SELECT coalesce(apg_clade,'(null)') apg_clade, count(*) FROM resolved_taxonomy GROUP BY 1 ORDER BY 2 DESC;

\echo '================ review, then COMMIT; or ROLLBACK; ================'
-- COMMIT;
