-- =====================================================================
-- THEOBROMA pipeline stage 06b: NPClassifier pathway-straggler correction
-- =====================================================================
-- Purpose : Correct the scrambled pathway on the 27 compounds carrying a
--           Triterpenoid superclass but an Alkaloids pathway (21 Fusidane,
--           6 Quassinoids). Per the NPClassifier ontology a Triterpenoid
--           superclass implies the Terpenoids pathway.
-- Inputs  : compounds(np_pathway, np_superclass)
-- Output  : corrected compounds.np_pathway; residual impossible = 0
-- Run     : after stage 06 (npc_load_apply) as an ontology-consistency pass.
-- Verbatim from fix_pathway_stragglers_apply.sql (validated + applied).
-- =====================================================================
BEGIN;
UPDATE compounds SET np_pathway='Terpenoids'
WHERE lower(np_pathway)='alkaloids' AND np_superclass ILIKE '%Triterpenoid%';
SELECT count(*) AS residual_impossible FROM compounds WHERE lower(np_pathway)='alkaloids' AND np_superclass ILIKE '%Triterpenoid%';
COMMIT;
