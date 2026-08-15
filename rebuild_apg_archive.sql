-- Rebuild the pre-APG archive from the current state + the deterministic mapping.
-- Every remapped row was originally 'liliopsida' (if its order is a monocot order)
-- or 'magnoliopsida' (otherwise). This faithfully reconstructs the original
-- taxclass for rollback/provenance, since the remap was deterministic by order.
\set ON_ERROR_STOP on

DROP TABLE IF EXISTS resolved_taxonomy_taxclass_pre_apg_20260804;
CREATE TABLE resolved_taxonomy_taxclass_pre_apg_20260804 AS
SELECT comp_id,
       CASE
         WHEN lower(taxorder) IN ('acorales','alismatales','arecales','asparagales',
              'commelinales','dioscoreales','liliales','pandanales','poales',
              'zingiberales','petrosaviales')
           THEN 'liliopsida'
         ELSE 'magnoliopsida'
       END AS taxclass_original,
       taxorder,
       taxclass AS taxclass_new
FROM resolved_taxonomy
WHERE taxclass IN ('eudicots','monocots','magnoliids','basal angiosperms',
                   'gymnosperms','ferns','green algae','red algae','bryophytes',
                   'lycophytes','hornworts');

\echo '=== rebuilt archive row count (expect 174,411) ==='
SELECT count(*) FROM resolved_taxonomy_taxclass_pre_apg_20260804;
\echo '=== original-class split (magnoliopsida vs liliopsida) ==='
SELECT taxclass_original, count(*) FROM resolved_taxonomy_taxclass_pre_apg_20260804 GROUP BY 1 ORDER BY 2 DESC;
\echo '=== spot: curcumin was liliopsida (zingiberales), now monocots ==='
SELECT * FROM resolved_taxonomy_taxclass_pre_apg_20260804 WHERE comp_id='THEO_0854403';
