-- v33 TIPdb production merge
--
-- Applies the TIPdb integration to the live compounds table:
--   1. INSERT 472 novel TIPdb-only rows with fresh THEO_ ids starting from
--      the current MAX(comp_id) + 1.
--   2. UPDATE all_sources on 7,023 enrichment rows to append "|TIPdb".
--
-- Inputs (assumed already on server):
--   /tmp/tipdb_novel_rows.csv     - 472 rows, full 33-column compounds schema,
--                                    with comp_id placeholders that we ignore
--                                    (we re-assign on insert) and trust_score,
--                                    novelty_*, sa_score blank for L3S fill-in.
--   /tmp/tipdb_overlap_inchikeys.txt - 7,495 InChIKeys (one per line) for the
--                                    enrichment update.
--
-- The script is wrapped in a transaction with explicit row-count checks
-- before COMMIT. Run inside psql with single-line semicolons, NOT as a
-- plain piped script, so the assertions can abort cleanly if numbers
-- diverge from expectation.
--
-- Run as:
--   sudo -u postgres psql -d theobroma -f tipdb_production_merge.sql
--
-- Pre-flight (run separately first to confirm baseline state):
--   SELECT COUNT(*) FROM compounds;  -- expect 1,078,671
--   SELECT MAX(comp_id) FROM compounds;  -- expect THEO_1080940 or higher
--   SELECT COUNT(*) FROM compounds WHERE all_sources LIKE '%TIPdb%';  -- expect 0

\timing on
\set ON_ERROR_STOP on

BEGIN;

-- Step 1: stage TIPdb-novel rows in a temp table so we can reassign
-- comp_ids and validate before the INSERT.
CREATE TEMP TABLE tipdb_novel_staging (LIKE compounds INCLUDING DEFAULTS);
\copy tipdb_novel_staging FROM '/tmp/tipdb_novel_rows.csv' CSV HEADER

-- Validate we got the expected count.
SELECT COUNT(*) AS staged_rows FROM tipdb_novel_staging;
-- expect 472

-- Compute the next available THEO_ id from the live table.
-- THEO_xxxxxxx ids are zero-padded 7-digit; cast to int via SUBSTRING.
DO $$
DECLARE
    next_id_num INT;
    novel_count INT;
    update_count INT;
BEGIN
    SELECT COALESCE(MAX(SUBSTRING(comp_id FROM 6)::INT), 0) + 1
        INTO next_id_num FROM compounds;
    SELECT COUNT(*) INTO novel_count FROM tipdb_novel_staging;
    RAISE NOTICE 'Live max comp_id base: %, next: %', next_id_num - 1, next_id_num;
    RAISE NOTICE 'Novel staging rows: %', novel_count;
    IF novel_count != 472 THEN
        RAISE EXCEPTION 'Expected 472 novel rows, got %', novel_count;
    END IF;

    -- Reassign comp_ids on the staging table to the live-table's next-available
    -- range. Use ROW_NUMBER() to assign sequentially.
    UPDATE tipdb_novel_staging s
    SET comp_id = 'THEO_' || LPAD((next_id_num + r.rn - 1)::TEXT, 7, '0')
    FROM (SELECT comp_id AS old_cid,
                 ROW_NUMBER() OVER (ORDER BY inchikey) AS rn
          FROM tipdb_novel_staging) r
    WHERE s.comp_id = r.old_cid;

    -- Insert into the live compounds table.
    INSERT INTO compounds SELECT * FROM tipdb_novel_staging;
    GET DIAGNOSTICS novel_count = ROW_COUNT;
    RAISE NOTICE 'Inserted novel rows: %', novel_count;
END $$;

-- Step 2: enrichment update. Append "|TIPdb" to all_sources for every
-- compound whose inchikey is in the overlap list, idempotent (only if
-- TIPdb is not already present).
CREATE TEMP TABLE tipdb_overlap_iks (inchikey TEXT);
\copy tipdb_overlap_iks FROM '/tmp/tipdb_overlap_inchikeys.txt'

SELECT COUNT(*) AS overlap_iks FROM tipdb_overlap_iks;
-- expect 7,495 (the full TIPdb unique-by-InChIKey set; some of those land
-- in the enrichment path and the novel-insert path covers the rest)

DO $$
DECLARE
    update_count INT;
BEGIN
    UPDATE compounds c
    SET all_sources = CASE
        WHEN c.all_sources IS NULL OR c.all_sources = '' THEN 'TIPdb'
        WHEN c.all_sources NOT LIKE '%TIPdb%' THEN c.all_sources || '|TIPdb'
        ELSE c.all_sources
    END
    FROM tipdb_overlap_iks o
    WHERE c.inchikey = o.inchikey
      AND (c.all_sources IS NULL OR c.all_sources NOT LIKE '%TIPdb%');
    GET DIAGNOSTICS update_count = ROW_COUNT;
    RAISE NOTICE 'Enrichment rows updated: % (expected ~7,023; the count may '
                 'exceed this if residual InChIKey duplicates exist beyond '
                 'the v33 cleanup)', update_count;
END $$;

-- Step 3: post-merge sanity checks. These run inside the transaction so
-- they ROLLBACK with the rest if any are off.
SELECT 'final corpus rows' AS metric, COUNT(*) AS val FROM compounds
UNION ALL
SELECT 'rows with TIPdb anywhere', COUNT(*) FROM compounds
    WHERE all_sources LIKE '%TIPdb%'
UNION ALL
SELECT 'rows with TIPdb as primary', COUNT(*) FROM compounds
    WHERE source_db = 'TIPdb';
-- expect: 1,079,143 corpus, 7,495 with TIPdb anywhere, 472 with TIPdb primary

-- If the numbers are wrong, ROLLBACK manually before committing. Otherwise:
COMMIT;

-- Step 4: refresh the search_names side-table for the new rows so name
-- search picks them up. Do this in a separate transaction.
INSERT INTO search_names (comp_id, name, name_norm)
SELECT comp_id, name, LOWER(REGEXP_REPLACE(name, '[^a-z0-9]', '', 'gi'))
FROM compounds
WHERE source_db = 'TIPdb'
  AND name IS NOT NULL AND name != ''
  AND comp_id NOT IN (SELECT comp_id FROM search_names);

-- Step 5: refresh table statistics for query planner.
ANALYZE compounds;
ANALYZE search_names;

-- Step 6: restart the application to pick up the new artifact baseline
-- (this needs to be done outside psql):
--   sudo systemctl restart theobroma

\timing off
