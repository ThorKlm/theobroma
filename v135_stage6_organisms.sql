-- v135_stage6_organisms.sql
-- Produces v135_coconut_organisms: resolved organism rows for the COCONUT-only
-- population (v1.35 compounds the v1.34 bridge does not annotate but COCONUT does).
-- Bridge-covered compounds are NOT touched here; their organisms already live,
-- fully resolved and untruncated, in compound_taxonomy and are reused by comp_id
-- at assembly. Supplementary compounds are organism-blank at source and produce
-- no rows. This script writes one new table and nothing live.
--
-- Inputs on theobroma: v135(inchikey), compounds(inchikey,comp_id),
--   compound_taxonomy(comp_id,token,accepted_taxon,genus,family,kingdom_any,...),
--   coconut_organisms(inchikey,organisms)  -- organisms pipe-delimited, raw occurrence text
-- Run: psql -h localhost -U theobroma -d theobroma -f v135_stage6_organisms.sql

BEGIN;

-- 1. COCONUT-only compounds: in v135, COCONUT carries an organism, bridge silent.
CREATE TEMP TABLE coc_only ON COMMIT DROP AS
SELECT co.inchikey, co.organisms
FROM coconut_organisms co
JOIN v135 v ON v.inchikey = co.inchikey
WHERE NOT EXISTS (
    SELECT 1 FROM compounds c
    JOIN compound_taxonomy ct ON ct.comp_id = c.comp_id
    WHERE c.inchikey = co.inchikey);

-- 2. Explode pipe-delimited occurrence text to one row per (inchikey, organism name).
CREATE TEMP TABLE coc_pairs ON COMMIT DROP AS
SELECT inchikey, trim(name) AS organism_name
FROM coc_only, unnest(string_to_array(organisms,'|')) AS name
WHERE trim(name) <> '';

-- 3. Resolution reference: reuse existing compound_taxonomy resolutions by name.
--    Prefer an accepted_taxon (species-level) match, fall back to genus match.
--    Majority kingdom per name breaks any within-name disagreement deterministically.
CREATE TEMP TABLE name_resolution ON COMMIT DROP AS
WITH distinct_names AS (
    SELECT DISTINCT organism_name, lower(organism_name) AS key FROM coc_pairs),
species_hit AS (
    SELECT d.organism_name,
           mode() WITHIN GROUP (ORDER BY ct.kingdom_any) AS kingdom_any,
           mode() WITHIN GROUP (ORDER BY ct.family)      AS family,
           mode() WITHIN GROUP (ORDER BY ct.genus)       AS genus
    FROM distinct_names d
    JOIN compound_taxonomy ct ON lower(ct.accepted_taxon) = d.key
    GROUP BY d.organism_name),
genus_hit AS (
    SELECT d.organism_name,
           mode() WITHIN GROUP (ORDER BY ct.kingdom_any) AS kingdom_any,
           mode() WITHIN GROUP (ORDER BY ct.family)      AS family,
           mode() WITHIN GROUP (ORDER BY ct.genus)       AS genus
    FROM distinct_names d
    JOIN compound_taxonomy ct ON lower(ct.genus) = d.key
    GROUP BY d.organism_name),
-- Fallback for names an exact match missed: strip parenthetical subgenus, drop
-- 'sp.'/'strain ...' qualifiers, take the leading genus token, retry on genus.
-- Recovers subgenus and indeterminate-species forms (Petrosia (Petrosia) ficiformis,
-- Agelas sp., Aspergillus sp. strain 05545) without guessing at true species.
genus_fallback_key AS (
    SELECT d.organism_name,
           lower(split_part(
               regexp_replace(d.organism_name, '\s*\([^)]*\)', '', 'g'),
               ' ', 1)) AS gkey
    FROM distinct_names d
    LEFT JOIN species_hit s USING (organism_name)
    LEFT JOIN genus_hit   g USING (organism_name)
    WHERE s.organism_name IS NULL AND g.organism_name IS NULL),
genus_fallback AS (
    SELECT f.organism_name,
           mode() WITHIN GROUP (ORDER BY ct.kingdom_any) AS kingdom_any,
           mode() WITHIN GROUP (ORDER BY ct.family)      AS family,
           mode() WITHIN GROUP (ORDER BY ct.genus)       AS genus
    FROM genus_fallback_key f
    JOIN compound_taxonomy ct ON lower(ct.genus) = f.gkey
    WHERE f.gkey <> ''
    GROUP BY f.organism_name),
-- Explicit non-identifications keep a kingdom hint, genus stays null.
unknown_hint AS (
    SELECT d.organism_name,
           CASE WHEN d.organism_name ILIKE 'unknown-fungus%'    THEN 'Fungi'
                WHEN d.organism_name ILIKE 'unknown-bacterium%' THEN 'Bacteria'
                WHEN d.organism_name ILIKE 'unknown-plant%'     THEN 'Plantae'
                WHEN d.organism_name ILIKE 'unknown-animal%'    THEN 'Animalia'
           END AS kingdom_any
    FROM distinct_names d
    WHERE d.organism_name ILIKE 'unknown-%')
SELECT d.organism_name,
       coalesce(s.kingdom_any, g.kingdom_any, gf.kingdom_any, u.kingdom_any) AS kingdom_any,
       coalesce(s.family,  g.family,  gf.family)          AS family,
       coalesce(s.genus,   g.genus,   gf.genus, d.organism_name) AS genus,
       (s.organism_name IS NOT NULL OR g.organism_name IS NOT NULL
        OR gf.organism_name IS NOT NULL) AS resolved,
       CASE WHEN s.organism_name IS NOT NULL THEN 'exact_species'
            WHEN g.organism_name IS NOT NULL THEN 'exact_genus'
            WHEN gf.organism_name IS NOT NULL THEN 'genus_fallback'
            WHEN u.organism_name IS NOT NULL THEN 'kingdom_hint'
            ELSE 'unresolved' END AS resolution_method
FROM distinct_names d
LEFT JOIN species_hit    s  USING (organism_name)
LEFT JOIN genus_hit      g  USING (organism_name)
LEFT JOIN genus_fallback gf USING (organism_name)
LEFT JOIN unknown_hint   u  USING (organism_name);

-- 4. Output table: per-organism rows, inchikey-keyed, tagged attested_coconut.
DROP TABLE IF EXISTS v135_coconut_organisms;
CREATE TABLE v135_coconut_organisms (
    inchikey            text NOT NULL,
    organism_name       text NOT NULL,
    genus               text,
    family              text,
    kingdom_any         text,
    resolved            boolean NOT NULL,
    resolution_method   text NOT NULL,
    organism_provenance text NOT NULL DEFAULT 'attested_coconut',
    PRIMARY KEY (inchikey, organism_name));

INSERT INTO v135_coconut_organisms
    (inchikey, organism_name, genus, family, kingdom_any, resolved, resolution_method)
SELECT p.inchikey, p.organism_name, r.genus, r.family, r.kingdom_any,
       r.resolved, r.resolution_method
FROM coc_pairs p
JOIN name_resolution r USING (organism_name)
ON CONFLICT (inchikey, organism_name) DO NOTHING;

CREATE INDEX ON v135_coconut_organisms (inchikey);
CREATE INDEX ON v135_coconut_organisms (lower(genus)) WHERE genus IS NOT NULL;

-- 5. Assertions and report.
DO $$
DECLARE
    n_compounds   int; n_rows int; n_names int; n_resolved int; n_kingdom int;
    n_double      int;
BEGIN
    SELECT count(DISTINCT inchikey) INTO n_compounds FROM v135_coconut_organisms;
    SELECT count(*)                 INTO n_rows      FROM v135_coconut_organisms;
    SELECT count(DISTINCT organism_name) INTO n_names FROM v135_coconut_organisms;
    SELECT count(*) FILTER (WHERE resolved)    INTO n_resolved FROM name_resolution;
    SELECT count(*) FILTER (WHERE kingdom_any IS NOT NULL) INTO n_kingdom
        FROM v135_coconut_organisms;

    -- No compound may be annotated by both sources: coc_only was built bridge-silent,
    -- so any overlap with a bridge compound is a construction error.
    SELECT count(*) INTO n_double
    FROM v135_coconut_organisms x
    JOIN compounds c ON c.inchikey = x.inchikey
    JOIN compound_taxonomy ct ON ct.comp_id = c.comp_id;
    IF n_double > 0 THEN
        RAISE EXCEPTION 'source overlap: % rows annotated by both bridge and COCONUT', n_double;
    END IF;

    RAISE NOTICE 'coconut-only compounds annotated: %', n_compounds;
    RAISE NOTICE 'organism rows written:            %', n_rows;
    RAISE NOTICE 'distinct organism names:          %', n_names;
    RAISE NOTICE 'names resolved (any method):      % of %', n_resolved, n_names;
    RAISE NOTICE 'rows carrying a kingdom:          %', n_kingdom;
    RAISE NOTICE 'source-overlap check passed (0 double-annotated)';
END $$;

-- Resolution profile, printed separately so it is easy to read.
DO $$
DECLARE r record;
BEGIN
    RAISE NOTICE '--- resolution method ---';
    FOR r IN
        SELECT resolution_method, count(DISTINCT organism_name) AS names,
               count(*) AS rows
        FROM v135_coconut_organisms GROUP BY resolution_method ORDER BY 2 DESC
    LOOP
        RAISE NOTICE '  % : names=% rows=%', r.resolution_method, r.names, r.rows;
    END LOOP;
END $$;

-- Inspect before trusting. Flip to COMMIT once the NOTICEs look right.
COMMIT;
-- COMMIT;
