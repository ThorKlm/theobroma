-- Residual = compounds with NO curated np_class AND NOT covered by the NPClassifier
-- re-run. These are the only compounds the inferred (XGBoost) tier scores.
-- Export to residual_comp_ids.tsv for the inference step.
\copy (SELECT c.comp_id FROM compounds c LEFT JOIN npc_reclassified r ON r.comp_id=c.comp_id WHERE (c.np_class IS NULL OR c.np_class='') AND (r.np_class IS NULL OR r.np_class='')) TO 'residual_comp_ids.tsv' WITH (FORMAT csv, HEADER, DELIMITER E'\t')
