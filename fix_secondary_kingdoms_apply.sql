-- FAIL 35 fix (APPLY): commits the change. Run only after reviewing the dry-run.
BEGIN;
UPDATE resolved_taxonomy
SET secondary_kingdoms = array_remove(secondary_kingdoms, kingdom)
WHERE kingdom = ANY(secondary_kingdoms);
SELECT count(*) AS residual_violations FROM resolved_taxonomy WHERE kingdom = ANY(secondary_kingdoms);
COMMIT;
