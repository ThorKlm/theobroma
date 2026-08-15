-- ISSUE 11 (APPLY): correct the scrambled pathway on the 27 Triterpenoid-superclass compounds.
BEGIN;
UPDATE compounds SET np_pathway='Terpenoids'
WHERE lower(np_pathway)='alkaloids' AND np_superclass ILIKE '%Triterpenoid%';
SELECT count(*) AS residual_impossible FROM compounds WHERE lower(np_pathway)='alkaloids' AND np_superclass ILIKE '%Triterpenoid%';
COMMIT;
