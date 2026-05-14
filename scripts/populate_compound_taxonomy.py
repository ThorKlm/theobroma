"""Populate compound_taxonomy from source_organism via WCVP resolution.

Downloads the latest WCVP names dump from Kew (sftp.kew.org), loads it
into a wcvp_names table, then tokenizes each compound's source_organism
on '; ' and resolves each token to genus and family via WCVP.

Resolution pattern:
    1. Look up the token in wcvp_names by lowercased taxon_name.
    2. If the row is Accepted, return its genus and family directly.
    3. If the row is a Synonym, follow accepted_plant_name_id to the
       canonical row and return its genus and family.

Coverage on the v33 corpus (May 2026 WCVP release): 67% of clean plant
binomials, dropping to 28% overall when non-plant kingdoms and NPO/NPC
identifier pollution are included. Non-plant kingdoms are not resolved
because WCVP covers only vascular plants.

Outputs:
    wcvp_names         (~1.4M rows; SELECT GRANT for theobroma role)
    compound_taxonomy  (~500K rows; SELECT GRANT for theobroma role)

NPO/NPC pollution: tokens beginning with NPO or NPC (NPAtlas internal
organism/compound IDs) are filtered out at populate time. They should
ideally not appear in source_organism at all; see open Issue 16 for
the upstream ingest fix.

Run as the database superuser (postgres) with INSERT/CREATE rights:
    sudo -u postgres ~/theobroma/venv/bin/python scripts/populate_compound_taxonomy.py

Requires:
    - urllib (stdlib) for download
    - psycopg2 for database access
    - The wcvp_names_stage import uses the |-delimited Kew dump format
      with a backspace-quote trick to disable CSV quoting handling.
"""
import os, urllib.request, zipfile, subprocess, sys

REPO = os.path.expanduser("~/theobroma")
EXT_DIR = os.path.join(REPO, "data", "external", "wcvp")
WCVP_URL = "https://sftp.kew.org/pub/data-repositories/WCVP/wcvp.zip"
ZIP_PATH = os.path.join(EXT_DIR, "wcvp.zip")
CSV_PATH = os.path.join(EXT_DIR, "wcvp_names.csv")

def fetch_wcvp():
    os.makedirs(EXT_DIR, exist_ok=True)
    if not os.path.exists(CSV_PATH):
        print(f"downloading {WCVP_URL}")
        urllib.request.urlretrieve(WCVP_URL, ZIP_PATH)
        with zipfile.ZipFile(ZIP_PATH) as z:
            z.extract("wcvp_names.csv", EXT_DIR)
        print(f"extracted to {CSV_PATH}")
    else:
        print(f"reusing existing {CSV_PATH}")

LOAD_SQL = r"""
DROP TABLE IF EXISTS wcvp_names;
CREATE TABLE wcvp_names (
    plant_name_id BIGINT PRIMARY KEY,
    ipni_id TEXT, taxon_rank TEXT, taxon_status TEXT,
    family TEXT, genus_hybrid TEXT, genus TEXT,
    species_hybrid TEXT, species TEXT,
    infraspecific_rank TEXT, infraspecies TEXT,
    parenthetical_author TEXT, primary_author TEXT, publication_author TEXT,
    place_of_publication TEXT, volume_and_page TEXT, first_published TEXT,
    nomenclatural_remarks TEXT, geographic_area TEXT,
    lifeform_description TEXT, climate_description TEXT,
    taxon_name TEXT, taxon_authors TEXT,
    accepted_plant_name_id BIGINT, basionym_plant_name_id BIGINT,
    replaced_synonym_author TEXT, homotypic_synonym TEXT,
    parent_plant_name_id BIGINT, powo_id TEXT,
    hybrid_formula TEXT, reviewed TEXT
);
\copy wcvp_names FROM '{csv}' WITH (FORMAT csv, DELIMITER '|', HEADER true, QUOTE E'\b');
CREATE INDEX wcvp_taxon_name_lower ON wcvp_names (LOWER(taxon_name));
CREATE INDEX wcvp_accepted_fk ON wcvp_names (accepted_plant_name_id);

DROP TABLE IF EXISTS compound_taxonomy;
CREATE TABLE compound_taxonomy (
    comp_id TEXT NOT NULL,
    token TEXT NOT NULL,
    accepted_taxon TEXT,
    genus TEXT,
    family TEXT,
    plant_name_id BIGINT,
    PRIMARY KEY (comp_id, token)
);

INSERT INTO compound_taxonomy (comp_id, token, accepted_taxon, genus, family, plant_name_id)
SELECT
    c.comp_id,
    TRIM(t) AS token,
    a.taxon_name, a.genus, a.family, a.plant_name_id
FROM compounds c,
     regexp_split_to_table(c.source_organism, '; ') AS t
LEFT JOIN wcvp_names s ON LOWER(s.taxon_name) = LOWER(TRIM(t))
LEFT JOIN wcvp_names a ON a.plant_name_id = COALESCE(s.accepted_plant_name_id, s.plant_name_id)
WHERE c.source_organism IS NOT NULL AND c.source_organism != ''
  AND TRIM(t) NOT LIKE 'NPO%' AND TRIM(t) NOT LIKE 'NPC%'
ON CONFLICT (comp_id, token) DO NOTHING;

CREATE INDEX ct_comp_id ON compound_taxonomy (comp_id);
CREATE INDEX ct_genus ON compound_taxonomy (genus) WHERE genus IS NOT NULL;
CREATE INDEX ct_family ON compound_taxonomy (family) WHERE family IS NOT NULL;

GRANT SELECT ON wcvp_names TO theobroma;
GRANT SELECT ON compound_taxonomy TO theobroma;

SELECT
    COUNT(*) AS total_token_rows,
    COUNT(DISTINCT comp_id) AS distinct_compounds_resolved,
    COUNT(*) FILTER (WHERE genus IS NOT NULL) AS rows_with_genus,
    COUNT(DISTINCT genus) FILTER (WHERE genus IS NOT NULL) AS distinct_genera,
    COUNT(DISTINCT family) FILTER (WHERE family IS NOT NULL) AS distinct_families
FROM compound_taxonomy;
"""

def main():
    fetch_wcvp()
    sql_path = "/tmp/_compound_taxonomy_load.sql"
    with open(sql_path, "w") as f:
        f.write(LOAD_SQL.format(csv=CSV_PATH))
    print("running psql load (this takes ~30s for the WCVP COPY)")
    rc = subprocess.call(["sudo", "-u", "postgres", "psql", "-d", "theobroma", "-f", sql_path])
    if rc != 0:
        print(f"psql exited with code {rc}")
        sys.exit(rc)
    print("done.")

if __name__ == "__main__":
    main()
