-- Reconcile per_source_license_attestation.license_tier to source_license_ref.
-- Run INTERACTIVELY: psql ... then \i reconcile_attestations.sql
-- The transaction is left OPEN at the end; review the output, then type COMMIT;
-- (or ROLLBACK;) by hand. Do NOT run with -f (it would exit and roll back).

\timing on
BEGIN;

-- In-DB convenience archive (the authoritative archive is the file dump pulled
-- to Downloads/theobroma_pre_license_state; this is a fast rollback aid).
CREATE TABLE IF NOT EXISTS per_source_license_attestation_pre_reconcile_20260805 AS
  SELECT * FROM per_source_license_attestation;

-- BEFORE: rows whose attested tier disagrees with the authoritative ref (expect 1,264,216)
SELECT count(*) AS before_disagree
FROM per_source_license_attestation a
JOIN source_license_ref r ON LOWER(r.src)=LOWER(a.source)
WHERE a.license_tier <> r.license_tier;

-- Snapshot every compound's resolved tier BEFORE, to prove it does NOT change.
CREATE TEMP TABLE _compound_tier_before AS
  SELECT comp_id, license_tier FROM compounds;

-- THE UPDATE: set each attestation's tier to its source's authoritative ref tier.
UPDATE per_source_license_attestation a
SET license_tier = r.license_tier
FROM source_license_ref r
WHERE LOWER(r.src)=LOWER(a.source) AND a.license_tier <> r.license_tier;

-- AFTER: disagreeing rows (must be 0)
SELECT count(*) AS after_disagree
FROM per_source_license_attestation a
JOIN source_license_ref r ON LOWER(r.src)=LOWER(a.source)
WHERE a.license_tier <> r.license_tier;

-- INVARIANT: no compound's resolved license_tier changed (must be 0 rows).
-- (compounds.license_tier is computed from the ref, not from attestations, so
--  reconciling attestations must leave every compound tier identical.)
SELECT count(*) AS compound_tiers_changed
FROM compounds c
JOIN _compound_tier_before b ON b.comp_id=c.comp_id
WHERE c.license_tier <> b.license_tier;

-- Spot-check the report examples: compound_tier unchanged, source rows now match ref.
SELECT a.comp_id, c.license_tier AS compound_tier, a.source,
       a.license_tier AS attested_now, r.license_tier AS ref
FROM per_source_license_attestation a
JOIN compounds c ON c.comp_id=a.comp_id
JOIN source_license_ref r ON LOWER(r.src)=LOWER(a.source)
WHERE a.comp_id IN ('THEO_0102053','THEO_0895633','THEO_0854403')
ORDER BY a.comp_id, a.source;

-- If before_disagree=1264216, after_disagree=0, compound_tiers_changed=0, and the
-- spot-checks read coherently: type   COMMIT;
-- Otherwise:                          ROLLBACK;
