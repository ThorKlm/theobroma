-- #5 fix (APPLY): add source_organism_curated column (if missing) and populate it
-- with the taxonomy-prioritized ordering. Commits. Run only after the dry-run looks right.
BEGIN;
ALTER TABLE compounds ADD COLUMN IF NOT EXISTS source_organism_curated text;

WITH exploded AS (
  SELECT c.comp_id, trim(o.org) AS org, split_part(trim(o.org),' ',1) AS org_genus
  FROM compounds c
  CROSS JOIN LATERAL unnest(string_to_array(c.source_organism,'; ')) AS o(org)
  WHERE c.source_organism LIKE '%; %'
),
anchor AS (
  SELECT comp_id, genus AS r_genus, family AS r_family
  FROM resolved_taxonomy WHERE genus IS NOT NULL AND genus <> ''
),
famgen AS (
  SELECT DISTINCT comp_id, genus AS ct_genus, family AS ct_family
  FROM compound_taxonomy WHERE genus IS NOT NULL AND genus <> ''
),
tiered AS (
  SELECT e.comp_id, e.org,
         CASE
           WHEN a.r_genus IS NOT NULL AND lower(e.org_genus)=lower(a.r_genus) THEN 0
           WHEN a.r_family IS NOT NULL AND EXISTS (
                  SELECT 1 FROM famgen f WHERE f.comp_id=e.comp_id
                    AND lower(f.ct_genus)=lower(e.org_genus) AND lower(f.ct_family)=lower(a.r_family)
                ) THEN 1
           ELSE 2 END AS tier
  FROM exploded e LEFT JOIN anchor a ON a.comp_id=e.comp_id
),
agg AS (
  SELECT comp_id, string_agg(org,'; ' ORDER BY tier, org) AS curated
  FROM tiered GROUP BY comp_id
)
UPDATE compounds c SET source_organism_curated = a.curated
FROM agg a WHERE a.comp_id = c.comp_id;

SELECT count(*) AS populated FROM compounds WHERE source_organism_curated IS NOT NULL;
SELECT split_part(source_organism_curated,'; ',1) AS morphine_first_org
FROM compounds WHERE comp_id='THEO_0511220';
COMMIT;
