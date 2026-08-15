-- v1.35 rebuild, Stage 5: taxonomy class correction (monocot recovery, surgical).
-- Runs on theobroma production DB. Dry-run by default (ROLLBACK; UPDATE commented).
--
-- DEFECT: an untracked ncbi_lineage re-derivation misclassed EVERY monocot-order
-- compound as 'magnoliopsida' instead of 'liliopsida' (verified: 0/13,654 monocot
-- rows were correctly classed pre-fix). The curated family is intact; only the
-- derived class is wrong.
--
-- FIX (surgical): for rows under a monocot ORDER whose FAMILY is a monocot family
-- (per the backup family->class map), set taxclass='liliopsida'. Scoped to monocot
-- orders so it cannot leak across kingdoms or into dicot orders — an earlier broad
-- family-remap wrongly touched 83 classes (animals, fungi, bacteria) and put
-- liliopsida under dicot orders; this scoping eliminates that.
--
-- Genuine dicot families under monocot orders (branch-B, ~853: Asteraceae, Fabaceae,
-- etc. — multi-organism resolution artifacts) correctly stay magnoliopsida and are
-- counted for Anne's review, not auto-resolved.

\set ON_ERROR_STOP on
BEGIN;

-- family -> class map from the backup (single-valued; used only to identify
-- which families are monocot).
CREATE TEMP TABLE fcm AS
SELECT family, lower(MAX(taxclass)) AS taxclass
FROM resolved_taxonomy_pre_streptophyta_fix_20260610
WHERE family IS NOT NULL AND taxclass IS NOT NULL
GROUP BY family;

-- supplementary monocot family (APG IV) absent from backup
INSERT INTO fcm(family, taxclass)
SELECT 'Taccaceae','liliopsida'
WHERE NOT EXISTS (SELECT 1 FROM fcm WHERE family='Taccaceae');

CREATE TEMP TABLE rt_fixed AS SELECT * FROM resolved_taxonomy;

-- surgical: monocot ORDER + monocot FAMILY -> liliopsida
UPDATE rt_fixed rt
SET taxclass='liliopsida'
FROM fcm
WHERE rt.family=fcm.family AND fcm.taxclass='liliopsida'
  AND lower(rt.taxorder) IN ('liliales','asparagales','poales','zingiberales',
    'commelinales','alismatales','dioscoreales','pandanales','acorales','arecales','petrosaviales')
  AND rt.taxclass IS DISTINCT FROM 'liliopsida';

-- assertions
DO $$
DECLARE changed int; lil int; leak int; branchb int; cannabis_cls text;
BEGIN
  SELECT count(*) INTO changed FROM rt_fixed rt JOIN resolved_taxonomy live USING(comp_id)
  WHERE rt.taxclass IS DISTINCT FROM live.taxclass;
  IF changed < 12000 OR changed > 15000 THEN
    RAISE EXCEPTION 'change count out of expected band: %', changed; END IF;

  -- every changed row must be under a monocot order (no leakage)
  SELECT count(*) INTO leak FROM rt_fixed rt JOIN resolved_taxonomy live USING(comp_id)
  WHERE rt.taxclass IS DISTINCT FROM live.taxclass
    AND lower(rt.taxorder) NOT IN ('liliales','asparagales','poales','zingiberales',
      'commelinales','alismatales','dioscoreales','pandanales','acorales','arecales','petrosaviales');
  IF leak <> 0 THEN RAISE EXCEPTION 'leakage: % changed rows outside monocot orders', leak; END IF;

  -- Cannabis must not be monocot
  SELECT taxclass INTO cannabis_cls FROM rt_fixed WHERE family='Cannabaceae' LIMIT 1;
  IF cannabis_cls='liliopsida' THEN RAISE EXCEPTION 'Cannabaceae wrongly monocot'; END IF;

  SELECT count(*) INTO lil FROM rt_fixed WHERE taxclass='liliopsida';
  SELECT count(*) INTO branchb FROM rt_fixed rt
  WHERE lower(rt.taxorder) IN ('liliales','asparagales','poales','zingiberales','commelinales',
    'alismatales','dioscoreales','pandanales','acorales','arecales','petrosaviales')
    AND rt.taxclass='magnoliopsida';

  RAISE NOTICE 'taxclass changes staged: %', changed;
  RAISE NOTICE 'Liliopsida after fix: % (backup snapshot: 12,532; higher = corpus growth)', lil;
  RAISE NOTICE 'no-leakage check passed (0 changes outside monocot orders)';
  RAISE NOTICE 'branch-B residue (dicot family under monocot order, for Anne): %', branchb;
END $$;

SELECT 'AFTER' AS state, taxclass, count(*)
FROM rt_fixed
WHERE lower(taxorder) IN ('liliales','asparagales','poales','zingiberales','commelinales')
GROUP BY taxclass ORDER BY count(*) DESC;

-- APPLY: uncomment next 2 lines + change ROLLBACK->COMMIT.
UPDATE resolved_taxonomy live SET taxclass=rt.taxclass
FROM rt_fixed rt WHERE live.comp_id=rt.comp_id AND live.taxclass IS DISTINCT FROM rt.taxclass;

COMMIT;
