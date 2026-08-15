-- Scope B (APPLY): rebuild genus/family + full lineage JOINTLY from compound_taxonomy
-- so each compound's ranks come from ONE consistent organism lineage. Commits.
-- kingdom + secondary_kingdoms are left untouched. Pre-fix backup already taken.
BEGIN;

CREATE TEMP TABLE _joint AS
WITH lineages AS (
  SELECT comp_id,
         ct.phylum_any AS phylum, ct.class_any AS taxclass, ct.order_any AS taxorder,
         ct.family AS family, ct.genus AS genus,
         count(*) AS support,
         ( (ct.phylum_any IS NOT NULL)::int + (ct.class_any IS NOT NULL)::int
         + (ct.order_any IS NOT NULL)::int + (ct.family IS NOT NULL)::int
         + (ct.genus IS NOT NULL)::int ) AS specificity
  FROM compound_taxonomy ct
  WHERE ct.genus IS NOT NULL AND ct.genus<>'' AND ct.family IS NOT NULL AND ct.family<>''
  GROUP BY comp_id, ct.phylum_any, ct.class_any, ct.order_any, ct.family, ct.genus
),
pick AS (
  SELECT DISTINCT ON (comp_id)
         comp_id, phylum, taxclass, taxorder, family, genus
  FROM lineages
  ORDER BY comp_id, support DESC, specificity DESC, genus ASC
)
SELECT * FROM pick;

-- update the full consistent lineage for every compound that has joint evidence
UPDATE resolved_taxonomy rt
SET genus = j.genus, family = j.family,
    taxorder = j.taxorder, taxclass = j.taxclass, phylum = j.phylum
FROM _joint j
WHERE j.comp_id = rt.comp_id
  AND ( COALESCE(rt.genus,'')<>COALESCE(j.genus,'')
     OR COALESCE(rt.family,'')<>COALESCE(j.family,'')
     OR COALESCE(rt.taxorder,'')<>COALESCE(j.taxorder,'')
     OR COALESCE(rt.taxclass,'')<>COALESCE(j.taxclass,'')
     OR COALESCE(rt.phylum,'')<>COALESCE(j.phylum,'') );

\echo '=== post-apply: residual impossible pairings (must be 0) ==='
SELECT count(*) AS still_inconsistent
FROM resolved_taxonomy rt
WHERE rt.genus IS NOT NULL AND rt.genus<>'' AND rt.family IS NOT NULL AND rt.family<>''
  AND EXISTS (SELECT 1 FROM compound_taxonomy ct WHERE ct.comp_id=rt.comp_id AND lower(ct.genus)=lower(rt.genus))
  AND NOT EXISTS (SELECT 1 FROM compound_taxonomy ct WHERE ct.comp_id=rt.comp_id AND lower(ct.genus)=lower(rt.genus) AND lower(ct.family)=lower(rt.family));

\echo '=== Gramine after ==='
SELECT genus, family, taxorder, taxclass, phylum FROM resolved_taxonomy WHERE comp_id='THEO_0412143';

COMMIT;
