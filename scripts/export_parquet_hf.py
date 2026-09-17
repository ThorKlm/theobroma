"""Export the principal tables to Parquet for the HuggingFace mirror.

row_group_size and write_page_index are load-bearing. The first export wrote two
row groups per file, with admet at 382 MB and compounds at 342 MB, which exceeded
the dataset viewer's scan limit and left is-valid false on every subset. At 50,000
rows each file has 23 to 31 groups with the largest at 26 MB, and the viewer loads.
"""
import os
import pandas as pd, pyarrow as pa, pyarrow.parquet as pq, psycopg2

# source_licenses carries every source assessed during the audit. in_corpus marks
# the 29 that contribute compounds; the rest are ingestion-file aliases, sources
# evaluated and excluded, and sources considered but never ingested.
TABLES = [
    ("compounds", "compounds"),
    ("compound_taxonomy", "compound_taxonomy"),
    ("resolved_taxonomy", "resolved_taxonomy"),
    ("license_attestations", "per_source_license_attestation"),
    ("source_licenses",
     "(SELECT r.*, (r.src IN (SELECT DISTINCT lower(source_db) FROM compounds)) AS in_corpus FROM source_license_ref r) s"),
    ("synonyms", "compound_synonyms"),
    ("regions", "compound_region_map"),
    ("admet", "admet"),
    ("scaffolds", "scaffolds"),
]

def main(out="data", dsn="host=localhost dbname=theobroma user=theobroma"):
    os.makedirs(out, exist_ok=True)
    con = psycopg2.connect(dsn)
    for name, tbl in TABLES:
        df = pd.read_sql("SELECT * FROM %s" % tbl, con)
        path = os.path.join(out, name + ".parquet")
        pq.write_table(pa.Table.from_pandas(df, preserve_index=False), path,
                       compression="zstd", row_group_size=50_000, write_page_index=True)
        print("%-24s %9d rows  %6.1f MB" % (name, len(df), os.path.getsize(path) / 1e6))
    con.close()

if __name__ == "__main__":
    main()
