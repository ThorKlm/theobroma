-- =====================================================================
-- THEOBROMA v1.35 : multi-organism lineage-stitching coherence fix
-- Option 1 (coherence-only): keep resolved genus, overwrite family/order/
-- class/phylum/apg_clade with that genus's COHERENT lineage taken as a unit
-- from the clean single-organism reference (coherent by construction).
-- Raw attestation lists (source_organism, all_attested_*) are NOT touched.
-- Staged: build ref -> preview -> archive -> update -> verify. Wrapped in a
-- transaction; inspect the previews, then COMMIT (or ROLLBACK) at the end.
-- =====================================================================
\set ON_ERROR_STOP on
BEGIN;

-- ---------------------------------------------------------------------
-- 0. Reference: coherent genus -> (family, taxorder, taxclass, phylum,
--    apg_clade) taken as a WHOLE TUPLE from single-organism compounds
--    (these cannot be stitched: one organism => one coherent lineage).
--    We pick, per genus, the most frequent complete lineage tuple.
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS genus_lineage_ref;
CREATE TABLE genus_lineage_ref AS
WITH clean AS (
  SELECT lower(rt.genus) AS g,
         rt.family, rt.taxorder, rt.taxclass, rt.phylum, rt.apg_clade,
         count(*) AS n
  FROM resolved_taxonomy rt
  JOIN compounds c ON c.comp_id = rt.comp_id
  WHERE c.source_organism NOT LIKE '%;%'          -- single-organism only
    AND rt.genus  IS NOT NULL AND rt.genus  <> ''
    AND rt.family IS NOT NULL AND rt.family <> ''
  GROUP BY 1,2,3,4,5,6
)
SELECT DISTINCT ON (g)
       g AS genus_lc, family, taxorder, taxclass, phylum, apg_clade, n AS support
FROM clean
ORDER BY g, n DESC;                               -- most-supported lineage per genus
ALTER TABLE genus_lineage_ref ADD PRIMARY KEY (genus_lc);
CREATE INDEX glr_genus ON genus_lineage_ref (genus_lc);

\echo '=== ref size (distinct genera with a coherent single-organism lineage) ==='
SELECT count(*) AS genera_in_ref FROM genus_lineage_ref;

-- ---------------------------------------------------------------------
-- 1. Identify the incoherent set: genus assigned a family other than the
--    genus's own coherent (single-organism) family. This is the set to fix.
-- ---------------------------------------------------------------------
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
  AND lower(rt.family) <> lower(r.family);        -- family disagrees => stitched

\echo '=== targets to fix (should be ~17,456) ==='
SELECT count(*) AS targets FROM coherence_fix_targets;

\echo '=== PREVIEW: 15 before/after lineages ==='
SELECT comp_id, genus,
       old_family||' -> '||new_family      AS family_change,
       old_order ||' -> '||coalesce(new_order,'(null)')  AS order_change,
       old_class ||' -> '||coalesce(new_class,'(null)')  AS class_change
FROM coherence_fix_targets ORDER BY genus LIMIT 15;

\echo '=== PREVIEW: the two curcumin entries + THEO_0901197 specifically ==='
SELECT comp_id, genus, old_family, new_family, old_order, new_order, old_class, new_class
FROM coherence_fix_targets WHERE comp_id IN ('THEO_0901197','THEO_0854403');

-- ---------------------------------------------------------------------
-- 2. Archive the pre-fix rows for rollback.
-- ---------------------------------------------------------------------
DROP TABLE IF EXISTS resolved_taxonomy_pre_coherence_fix_20260802;
CREATE TABLE resolved_taxonomy_pre_coherence_fix_20260802 AS
SELECT rt.*
FROM resolved_taxonomy rt
JOIN coherence_fix_targets t ON t.comp_id = rt.comp_id;

\echo '=== archived pre-fix rows (should equal targets) ==='
SELECT count(*) AS archived FROM resolved_taxonomy_pre_coherence_fix_20260802;

-- ---------------------------------------------------------------------
-- 3. Apply the coherent lineage. Keep genus; overwrite the rest as a UNIT.
-- ---------------------------------------------------------------------
UPDATE resolved_taxonomy rt
SET family    = t.new_family,
    taxorder  = t.new_order,
    taxclass  = t.new_class,
    phylum    = t.new_phylum,
    apg_clade = t.new_clade
FROM coherence_fix_targets t
WHERE rt.comp_id = t.comp_id;

\echo '=== rows updated ==='
SELECT count(*) AS updated FROM coherence_fix_targets;

-- ---------------------------------------------------------------------
-- 4. VERIFY: recompute incoherence AFTER the fix. Should drop to ~0
--    (residual = genera absent from the single-organism reference).
-- ---------------------------------------------------------------------
\echo '=== POST-FIX incoherent count (target: near 0) ==='
WITH genus_family AS (
  SELECT lower(genus) g, lower(family) f, count(*) n FROM resolved_taxonomy
  WHERE genus IS NOT NULL AND genus<>'' AND family IS NOT NULL AND family<>'' GROUP BY 1,2),
dom AS (SELECT DISTINCT ON (g) g, f FROM genus_family ORDER BY g, n DESC)
SELECT count(*) AS still_incoherent
FROM resolved_taxonomy rt JOIN dom d ON d.g=lower(rt.genus)
WHERE rt.family IS NOT NULL AND rt.family<>'' AND lower(rt.family)<>d.f;

\echo '=== VERIFY the two curcumin entries now coherent ==='
SELECT comp_id, genus, family, taxorder, taxclass, apg_clade
FROM resolved_taxonomy WHERE comp_id IN ('THEO_0901197','THEO_0854403');

\echo '=== VERIFY: kingdom counts UNCHANGED (kingdom not touched) ==='
SELECT kingdom, count(*) FROM resolved_taxonomy GROUP BY 1 ORDER BY 2 DESC;

\echo '=== apg_clade distribution AFTER (family/clade may shift for corrected plants) ==='
SELECT coalesce(apg_clade,'(null)') apg_clade, count(*) FROM resolved_taxonomy GROUP BY 1 ORDER BY 2 DESC;

\echo '======================================================================'
\echo 'Review previews + POST-FIX incoherent count. If good: COMMIT;  else: ROLLBACK;'
\echo '(Transaction left OPEN on purpose.)'
\echo '======================================================================'
-- COMMIT;   -- <-- uncomment / type COMMIT; to apply, or ROLLBACK; to abort
