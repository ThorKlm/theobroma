-- THEOBROMA per-compound license resolver
-- Applies most-restrictive-wins resolution across all attesting sources
-- per compound, writing the resolved tier to compounds.license_tier.
--
-- Precedence ordering (most restrictive first):
--
--   1. copyleft         -- share-alike obligations propagate to derivatives
--   2. Unspecified      -- absence of documented license; treated conservatively
--   3. CC BY-NC 4.0     -- attribution + non-commercial-only
--   4. CC BY 4.0        -- attribution required
--   5. CC0              -- public domain dedication; no obligations
--
-- Reasoning for CC BY 4.0 above CC0: a downstream user redistributing a
-- compound co-attested by a CC0 source AND a CC BY 4.0 source must comply
-- with all attestations simultaneously. CC0 imposes no obligations; CC BY 4.0
-- imposes attribution. The combined obligation is "attribute", so the
-- resolved tier must be CC BY 4.0 (the more restrictive). If a compound is
-- attested only by CC0 sources, the resolved tier is CC0 because no
-- attribution-requiring source claims it.
--
-- Prerequisite: per_source_license_attestation must be populated for every
-- compound (via the population script that joins all_sources with the
-- license map derived from sources.yaml).

BEGIN;

-- Precedence ranking table.
CREATE TEMP TABLE tier_rank (
    tier text PRIMARY KEY,
    rank integer NOT NULL  -- lower rank = more restrictive
);
INSERT INTO tier_rank VALUES
    ('copyleft',     1),
    ('Unspecified',  2),
    ('CC BY-NC 4.0', 3),
    ('CC BY 4.0',    4),
    ('CC0',          5);

-- Per-compound resolved tier: the lowest-rank (most restrictive) license
-- among all attesting sources for that compound.
CREATE TEMP TABLE resolved_tier_new AS
SELECT
    psla.comp_id,
    (SELECT psla2.license_tier
     FROM per_source_license_attestation psla2
     JOIN tier_rank tr ON tr.tier = psla2.license_tier
     WHERE psla2.comp_id = psla.comp_id
     ORDER BY tr.rank ASC
     LIMIT 1) AS license_tier_resolved
FROM per_source_license_attestation psla
GROUP BY psla.comp_id;

-- Diagnostic: tier transitions before applying the UPDATE.
SELECT
    c.license_tier AS current_tier,
    r.license_tier_resolved AS new_tier,
    COUNT(*) AS n
FROM compounds c
JOIN resolved_tier_new r ON r.comp_id = c.comp_id
WHERE c.license_tier IS DISTINCT FROM r.license_tier_resolved
GROUP BY c.license_tier, r.license_tier_resolved
ORDER BY n DESC;

-- Apply the resolved tier.
UPDATE compounds c
SET license_tier = r.license_tier_resolved
FROM resolved_tier_new r
WHERE r.comp_id = c.comp_id
  AND c.license_tier IS DISTINCT FROM r.license_tier_resolved;

-- Final tier distribution corpus-wide.
SELECT license_tier, COUNT(*) AS n
FROM compounds
GROUP BY license_tier
ORDER BY n DESC;

COMMIT;
