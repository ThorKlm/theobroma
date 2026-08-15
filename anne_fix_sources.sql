\pset pager off
\echo '=== Do we have raw per-source organism data to re-derive source_organism correctly? ==='
SELECT table_name FROM information_schema.tables WHERE table_schema='public'
  AND (table_name ILIKE '%organism%' OR table_name ILIKE '%coconut_org%' OR table_name ILIKE '%occurrence%')
ORDER BY table_name;

\echo '=== Do we have original NPClassifier output (per-compound coherent triples) to re-join #1? ==='
SELECT table_name FROM information_schema.tables WHERE table_schema='public'
  AND (table_name ILIKE '%classif%' OR table_name ILIKE '%npclass%' OR table_name ILIKE '%class_ik%' OR table_name ILIKE '%classified%')
ORDER BY table_name;

\echo '=== classified_ik: does it hold coherent per-inchikey classification? (candidate re-join source) ==='
SELECT column_name FROM information_schema.columns WHERE table_name='classified_ik' ORDER BY ordinal_position;
SELECT count(*) FROM classified_ik;

\echo '=== sample classified_ik to see if pathway/superclass/class are coherent there ==='
SELECT * FROM classified_ik ORDER BY random() LIMIT 5;

\echo '=== does classified_ik have THEO_1055282 (the scrambled compound) with a COHERENT classification? ==='
SELECT ci.* FROM classified_ik ci
JOIN compounds c ON c.inchikey = ci.inchikey
WHERE c.comp_id='THEO_1055282';

\echo '=== v135_coconut_organisms / coconut_org: do these hold correct per-inchikey organisms? ==='
SELECT column_name FROM information_schema.columns WHERE table_name='v135_coconut_organisms' ORDER BY ordinal_position;
SELECT count(*) FROM v135_coconut_organisms;

\echo '=== check: is there a source_organism_curated anywhere (the v135 fix) we can port? ==='
SELECT table_name, column_name FROM information_schema.columns
WHERE column_name ILIKE '%organism_curated%' OR column_name ILIKE '%curated%organism%';

\echo '=== morphine: what does the v135 preswap (curated) version show vs current v134? ==='
SELECT 'v134_live' src, comp_id, source_organism FROM compounds WHERE name ILIKE 'morphine' LIMIT 3;
SELECT 'v135_preswap' src, comp_id, source_organism,
       (SELECT string_agg(column_name,',') FROM information_schema.columns
        WHERE table_name='compounds_v135_preswap_20260802' AND column_name ILIKE '%organism%') AS organism_cols
FROM compounds_v135_preswap_20260802 WHERE name ILIKE 'morphine' LIMIT 3;
