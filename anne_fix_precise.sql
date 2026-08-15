\pset pager off
\echo '=== ORGANISM: is v134 source_organism truly truncated in storage? char length check ==='
SELECT comp_id, length(source_organism) AS v134_len,
       (SELECT length(source_organism) FROM compounds_v135_preswap_20260802 p WHERE p.comp_id=c.comp_id) AS v135_raw_len,
       (SELECT length(source_organism_curated) FROM compounds_v135_preswap_20260802 p WHERE p.comp_id=c.comp_id) AS v135_curated_len
FROM compounds c WHERE comp_id='THEO_0511220';

\echo '=== is there a systematic truncation? compounds where v134 organism is SHORTER than v135 ==='
SELECT count(*) AS v134_shorter_than_v135
FROM compounds c JOIN compounds_v135_preswap_20260802 p ON p.comp_id=c.comp_id
WHERE length(c.source_organism) < length(p.source_organism);

\echo '=== distribution of v134 organism lengths (is there a hard cap like 60/64/255 chars?) ==='
SELECT max(length(source_organism)) AS max_len,
       count(*) FILTER (WHERE length(source_organism) BETWEEN 60 AND 64) AS near_64,
       count(*) FILTER (WHERE length(source_organism) BETWEEN 250 AND 256) AS near_255
FROM compounds WHERE source_organism IS NOT NULL;

\echo '=== morphine curated organism (host-contamination removed?) ==='
SELECT comp_id, source_organism_curated FROM compounds_v135_preswap_20260802 WHERE comp_id='THEO_0511220';

\echo '=== CLASSIFICATION: where do np_pathway/superclass/class values live besides compounds? ==='
SELECT table_name, string_agg(column_name, ', ') AS class_cols
FROM information_schema.columns
WHERE column_name IN ('np_pathway','np_superclass','np_class')
GROUP BY table_name ORDER BY table_name;

\echo '=== does v135_preswap have DIFFERENT (correct?) classification for THEO_1055282? ==='
SELECT 'v134_live' src, np_pathway, np_superclass, np_class FROM compounds WHERE comp_id='THEO_1055282'
UNION ALL
SELECT 'v135_preswap', np_pathway, np_superclass, np_class FROM compounds_v135_preswap_20260802 WHERE comp_id='THEO_1055282';

\echo '=== NPClassifier coherence: does the SAME inchikey get consistent class across compounds? ==='
\echo '    (if np_class Quassinoids always pairs with the same pathway, it is coherent) ==='
SELECT np_class, np_pathway, np_superclass, count(*)
FROM compounds WHERE np_class='Quassinoids'
GROUP BY 1,2,3 ORDER BY 4 DESC;

\echo '=== is np_class assignment kingdom-plausible? Quassinoids by kingdom (should be plant) ==='
SELECT kingdom, count(*) FROM compounds WHERE np_class='Quassinoids' GROUP BY 1 ORDER BY 2 DESC;
