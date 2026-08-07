-- 05d: APG taxclass relabelling (verbatim reference).
-- Consolidated from apg_taxclass_fix.sql; the embedded order-to-clade map is
-- kept byte-for-byte and must not be paraphrased. Runs after 05c.

-- =====================================================================
-- APG IV taxclass correction: replace the deprecated 'magnoliopsida' /
-- 'liliopsida' class names with modern APG IV clades, mapping by taxorder.
-- Corrects ~174k rows AND fixes genuine misfilings (conifers/ferns/algae
-- that were wrongly dumped under magnoliopsida).
--
-- Field choice 2a: overwrite taxclass directly (tree/facets read taxclass,
-- so the corrected clades show immediately, no app changes).
-- Only rows whose current taxclass is one of the two deprecated names are
-- touched, and only where their taxorder is in the mapping. All other rows
-- (aves, agaricomycetes, ...) are untouched.
--
-- Staged: build map -> archive -> preview -> update -> verify -> COMMIT/ROLLBACK
-- =====================================================================
\set ON_ERROR_STOP on
\timing on

-- 1. Build the order -> clade lookup as a temp table.
DROP TABLE IF EXISTS apg_order_clade;
CREATE TEMP TABLE apg_order_clade (taxorder text PRIMARY KEY, clade text);
INSERT INTO apg_order_clade (taxorder, clade) VALUES
-- monocots
('acorales','monocots'),('alismatales','monocots'),('arecales','monocots'),
('asparagales','monocots'),('commelinales','monocots'),('dioscoreales','monocots'),
('liliales','monocots'),('pandanales','monocots'),('poales','monocots'),
('zingiberales','monocots'),('petrosaviales','monocots'),
-- eudicots
('apiales','eudicots'),('aquifoliales','eudicots'),('asterales','eudicots'),
('boraginales','eudicots'),('brassicales','eudicots'),('bruniales','eudicots'),
('buxales','eudicots'),('cardiopteridales','eudicots'),('caryophyllales','eudicots'),
('celastrales','eudicots'),('cornales','eudicots'),('crossosomatales','eudicots'),
('cucurbitales','eudicots'),('dilleniales','eudicots'),('dipsacales','eudicots'),
('ericales','eudicots'),('escalloniales','eudicots'),('fabales','eudicots'),
('fagales','eudicots'),('garryales','eudicots'),('gentianales','eudicots'),
('geraniales','eudicots'),('gunnerales','eudicots'),('icacinales','eudicots'),
('lamiales','eudicots'),('malpighiales','eudicots'),('malvales','eudicots'),
('metteniusales','eudicots'),('myrtales','eudicots'),('oxalidales','eudicots'),
('paracryphiales','eudicots'),('picramniales','eudicots'),('proteales','eudicots'),
('ranunculales','eudicots'),('rosales','eudicots'),('santalales','eudicots'),
('sapindales','eudicots'),('saxifragales','eudicots'),('solanales','eudicots'),
('trochodendrales','eudicots'),('vahliales','eudicots'),('vitales','eudicots'),
('zygophyllales','eudicots'),
-- magnoliids
('canellales','magnoliids'),('laurales','magnoliids'),('magnoliales','magnoliids'),
('piperales','magnoliids'),
-- basal angiosperms (ANA grade + Chloranthales + Ceratophyllales)
('austrobaileyales','basal angiosperms'),('nymphaeales','basal angiosperms'),
('chloranthales','basal angiosperms'),('ceratophyllales','basal angiosperms'),
-- gymnosperms
('pinales','gymnosperms'),('cupressales','gymnosperms'),('araucariales','gymnosperms'),
('ephedrales','gymnosperms'),('gnetales','gymnosperms'),('cycadales','gymnosperms'),
('ginkgoales','gymnosperms'),
-- ferns & fern allies
('polypodiales','ferns'),('cyatheales','ferns'),('gleicheniales','ferns'),
('marattiales','ferns'),('ophioglossales','ferns'),('osmundales','ferns'),
('schizaeales','ferns'),('equisetales','ferns'),
-- lycophytes
('lycopodiales','lycophytes'),('selaginellales','lycophytes'),
-- bryophytes / hornworts
('bryales','bryophytes'),('jungermanniales','bryophytes'),('marchantiales','bryophytes'),
('pelliales','bryophytes'),('blasiales','bryophytes'),('anthocerotales','hornworts'),
-- green algae
('charales','green algae'),('chlamydomonadales','green algae'),('chlorellales','green algae'),
('cladophorales','green algae'),('dasycladales','green algae'),('bryopsidales','green algae'),
('mamiellales','green algae'),('ulvales','green algae'),
-- red algae
('bangiales','red algae'),('bonnemaisoniales','red algae'),('ceramiales','red algae'),
('corallinales','red algae'),('gigartinales','red algae'),('gracilariales','red algae');

\echo '=== map size (expect ~100) ==='
SELECT count(*) FROM apg_order_clade;

-- 2. Coverage check: any deprecated-class row whose order is NOT in the map?
\echo '=== deprecated rows with an UNMAPPED order (should be 0) ==='
SELECT rt.taxorder, count(*) n
FROM resolved_taxonomy rt
WHERE lower(rt.taxclass) IN ('magnoliopsida','liliopsida')
  AND NOT EXISTS (SELECT 1 FROM apg_order_clade m WHERE m.taxorder = lower(rt.taxorder))
GROUP BY 1 ORDER BY 2 DESC;

BEGIN;

-- 3. Archive current taxclass for the affected rows.
DROP TABLE IF EXISTS resolved_taxonomy_taxclass_pre_apg_20260804;
CREATE TABLE resolved_taxonomy_taxclass_pre_apg_20260804 AS
SELECT comp_id, taxclass, taxorder FROM resolved_taxonomy
WHERE lower(taxclass) IN ('magnoliopsida','liliopsida');
\echo '=== archived affected rows ==='
SELECT count(*) FROM resolved_taxonomy_taxclass_pre_apg_20260804;

-- 4. Preview the remap distribution.
\echo '=== preview: new clade distribution for the affected rows ==='
SELECT m.clade, count(*) n
FROM resolved_taxonomy rt JOIN apg_order_clade m ON m.taxorder = lower(rt.taxorder)
WHERE lower(rt.taxclass) IN ('magnoliopsida','liliopsida')
GROUP BY 1 ORDER BY 2 DESC;

-- 5. Apply: overwrite taxclass with the mapped clade.
UPDATE resolved_taxonomy rt
SET taxclass = m.clade
FROM apg_order_clade m
WHERE m.taxorder = lower(rt.taxorder)
  AND lower(rt.taxclass) IN ('magnoliopsida','liliopsida');
\echo '=== rows updated ==='

-- 6. Verify: no deprecated names remain; spot-check known compounds.
\echo '=== any magnoliopsida/liliopsida left? (expect 0) ==='
SELECT count(*) FROM resolved_taxonomy WHERE lower(taxclass) IN ('magnoliopsida','liliopsida');
\echo '=== curcumin THEO_0854403 (zingiberales -> monocots) ==='
SELECT comp_id, kingdom, phylum, taxclass, taxorder, genus FROM resolved_taxonomy WHERE comp_id='THEO_0854403';
\echo '=== a eudicot spot-check (asterales) ==='
SELECT taxclass, count(*) FROM resolved_taxonomy WHERE lower(taxorder)='asterales' GROUP BY 1;
\echo '=== a stray spot-check (pinales -> gymnosperms) ==='
SELECT taxclass, count(*) FROM resolved_taxonomy WHERE lower(taxorder)='pinales' GROUP BY 1;
\echo '=== full taxclass vocab now: confirm modern clades present, deprecated gone ==='
SELECT taxclass, count(*) FROM resolved_taxonomy
WHERE taxclass IN ('eudicots','monocots','magnoliids','basal angiosperms','gymnosperms','ferns','green algae','red algae','bryophytes','lycophytes','hornworts','magnoliopsida','liliopsida')
GROUP BY 1 ORDER BY 2 DESC;

\echo '======================================================================'
\echo 'INSPECT: 0 deprecated left, curcumin=monocots, asterales=eudicots,'
\echo 'pinales=gymnosperms. If good: COMMIT; else ROLLBACK;'
\echo '======================================================================'
-- COMMIT;
