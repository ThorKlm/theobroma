-- =====================================================================
-- THEOBROMA pipeline stage 09b: reconcile attestation tiers to the ref
-- =====================================================================
-- Purpose : Align per_source_license_attestation.license_tier to the
--           authoritative source_license_ref. Provenance hygiene only:
--           compounds.license_tier is computed from the ref, so resolved
--           per-compound tiers are INVARIANT under this (verified in-script).
-- Inputs/Output : per_source_license_attestation (in place)
-- Run     : after stage 09; interactive (operator COMMIT).
-- Verbatim from reconcile_attestations.sql.
-- =====================================================================
\timing on
BEGIN;

CREATE TABLE IF NOT EXISTS per_source_license_attestation_pre_reconcile_20260805 AS
  SELECT * FROM per_source_license_attestation;

SELECT count(*) AS before_disagree
FROM per_source_license_attestation a
JOIN source_license_ref r ON LOWER(r.src)=LOWER(a.source)
WHERE a.license_tier <> r.license_tier;

CREATE TEMP TABLE _compound_tier_before AS
  SELECT comp_id, license_tier FROM compounds;

UPDATE per_source_license_attestation a
SET license_tier = r.license_tier
FROM source_license_ref r
WHERE LOWER(r.src)=LOWER(a.source) AND a.license_tier <> r.license_tier;

SELECT count(*) AS after_disagree
FROM per_source_license_attestation a
JOIN source_license_ref r ON LOWER(r.src)=LOWER(a.source)
WHERE a.license_tier <> r.license_tier;

-- INVARIANT: no compound's resolved license_tier changed (must be 0).
SELECT count(*) AS compound_tiers_changed
FROM compounds c
JOIN _compound_tier_before b ON b.comp_id=c.comp_id
WHERE c.license_tier <> b.license_tier;

SELECT a.comp_id, c.license_tier AS compound_tier, a.source,
       a.license_tier AS attested_now, r.license_tier AS ref
FROM per_source_license_attestation a
JOIN compounds c ON c.comp_id=a.comp_id
JOIN source_license_ref r ON LOWER(r.src)=LOWER(a.source)
WHERE a.comp_id IN ('THEO_0102053','THEO_0895633','THEO_0854403')
ORDER BY a.comp_id, a.source;

-- before_disagree=1264216, after_disagree=0, compound_tiers_changed=0 -> COMMIT;
