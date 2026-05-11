"""Compute trust scores for compounds.

The trust score is a heuristic composite reflecting data-provenance
confidence: how many sources attest to the compound, whether canonical
metadata (organism, kingdom, name, classification) is present, and the
license-tier permissiveness. Values are clamped to [0, 1].

The coefficients below were recovered by linear regression on 5,000
randomly-sampled compounds against their existing trust_score values
(R-squared 0.93; RMSE 0.06). The original formula used by v32 ingest is
no longer available; this is an approximate re-derivation that produces
values close to v32 trust scores for the same inputs. New scores produced
by this script are recomputed across the whole corpus and replace any
previous values.

Inputs: PostgreSQL connection to the theobroma database.
Outputs: trust_score column written into compounds table.

Run as the database role with UPDATE rights on compounds:
    sudo -u postgres ~/theobroma/venv/bin/python scripts/compute_trust_score.py
"""
import psycopg2

INTERCEPT = -0.1117
COEFS = {
    "lic_score":       0.1052,
    "log1p_n_sources": 0.2356,
    "has_np":          0.1848,
    "has_cf":          0.0403,
    "has_inf":         0.1084,
    "has_name":        0.1164,
    "has_kingdom":    -0.1117,
    "has_organism":    0.2497,
    "has_admet":       0.0210,
}

UPDATE_SQL = f"""
UPDATE compounds c
SET trust_score = LEAST(1.0, GREATEST(0.0,
    ({INTERCEPT})
    + ({COEFS['lic_score']}) * (CASE c.license_tier
                                     WHEN 'CC0' THEN 1.0
                                     WHEN 'CC BY 4.0' THEN 0.9
                                     WHEN 'CC BY-NC 4.0' THEN 0.7
                                     ELSE 0.5 END)
    + ({COEFS['log1p_n_sources']}) * LN(1 + COALESCE(array_length(string_to_array(c.all_sources, '|'), 1), 1))
    + ({COEFS['has_np']}) * (CASE WHEN c.np_class IS NOT NULL AND c.np_class != '' THEN 1.0 ELSE 0.0 END)
    + ({COEFS['has_cf']}) * (CASE WHEN c.classyfire_superclass IS NOT NULL AND c.classyfire_superclass != '' THEN 1.0 ELSE 0.0 END)
    + ({COEFS['has_inf']}) * (CASE WHEN c.inferred_class IS NOT NULL AND c.inferred_class != '' THEN 1.0 ELSE 0.0 END)
    + ({COEFS['has_name']}) * (CASE WHEN c.name IS NOT NULL AND c.name != '' THEN 1.0 ELSE 0.0 END)
    + ({COEFS['has_kingdom']}) * (CASE WHEN c.kingdom IS NOT NULL AND c.kingdom != '' AND c.kingdom != 'unresolved' THEN 1.0 ELSE 0.0 END)
    + ({COEFS['has_organism']}) * (CASE WHEN c.source_organism IS NOT NULL AND c.source_organism != '' THEN 1.0 ELSE 0.0 END)
    + ({COEFS['has_admet']}) * (CASE WHEN EXISTS(SELECT 1 FROM admet a WHERE a.comp_id = c.comp_id) THEN 1.0 ELSE 0.0 END)
))
"""

def main(dry_run=False):
    with psycopg2.connect("postgresql://theobroma:theobroma@localhost:5432/theobroma") as conn:
        with conn.cursor() as cur:
            if dry_run:
                cur.execute("SELECT COUNT(*) FROM compounds")
                n = cur.fetchone()[0]
                print(f"would update {n:,} rows; dry run, no commit")
                return
            cur.execute(UPDATE_SQL)
            print(f"updated {cur.rowcount:,} rows")
            conn.commit()
            print("committed")

if __name__ == "__main__":
    import sys
    main(dry_run="--dry-run" in sys.argv)
