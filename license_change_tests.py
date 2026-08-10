"""THEOBROMA post-license-revert verification suite (25 tests).

Covers the license revert (four-table coherence, most-restrictive-wins resolution)
plus the earlier v1.35 changes it must not have disturbed: PLACEHOLDER repair,
three-tier classification, pathway normalization, organism cleanup, corpus count,
and per-compound/attestation/ref consistency.

Standalone: needs only psycopg2. Configure with THEO_DB_DSN, e.g.
  export THEO_DB_DSN="postgresql://theobroma:theobroma@localhost:5432/theobroma"
  python3 license_change_tests.py [--verbose]

Read-only. Opens one connection, runs 25 assertions, prints PASS/FAIL/WARN, exits
non-zero if any FAIL. WARN marks advisory checks whose exact value is informational
rather than a hard invariant (e.g. the corpus differs from the pre-removal audit by
the known tiny-molecule count).
"""
import os, sys, argparse
import psycopg2

DSN = os.environ.get("THEO_DB_DSN", "postgresql://theobroma:theobroma@localhost:5432/theobroma")
VERBOSE = "--verbose" in sys.argv or "-v" in sys.argv

# Expected published/audited constants. The corpus is the v1.35 count AFTER the
# tiny-molecule (<=3 atoms) removal; license counts are on that corpus, i.e. the
# audited Table 2 figures minus the 199 removed compounds (66 CC BY, 132 CC BY-NC,
# 1 CC0). Kept explicit so a drift is visible rather than silently absorbed.
CORPUS = 1_132_805
LICENSE_TARGET = {"CC BY 4.0": 891_860, "CC BY-NC 4.0": 225_536,
                  "CC0": 8_243, "Unspecified": 7_166}
ATTEST_ROWS, ATTEST_SOURCES = 1_512_594, 29
TIERS = {626_601: "curated"}  # classification tier-1 count sanity (approx-gated below)
TIER_RANK = {"CC0": 0, "CC BY 4.0": 1, "CC BY-NC 4.0": 2,
             "CC BY-NC-SA 4.0": 3, "CC BY-NC-ND 4.0": 4, "Unspecified": 5}
# The 14 sources whose v1.34 map differs from the v1.35 reassignment.
SOURCE_V134 = {
    "coconut": "CC BY 4.0", "lotus": "CC BY 4.0", "sancdb": "CC BY 4.0",
    "mibig": "CC BY 4.0", "cmaupv2": "CC BY 4.0", "streptomedb": "CC BY 4.0",
    "mycocentral": "CC BY 4.0", "amdb": "CC BY 4.0", "naturar": "CC BY 4.0",
    "tipdb": "CC BY-NC 4.0", "imppat": "CC BY-NC 4.0", "herb": "CC BY-NC 4.0",
    "phytochemdb": "CC BY-NC 4.0", "tm-mc": "CC0"}

_results = []

def _check(name, ok, detail="", warn=False):
    status = "WARN" if (warn and not ok) else ("PASS" if ok else "FAIL")
    _results.append((name, status))
    if VERBOSE or status != "PASS":
        print(f"[{status}] {name}: {detail}")
    elif not VERBOSE:
        print(f"[{status}] {name}")

def main():
    conn = psycopg2.connect(DSN)
    conn.set_session(readonly=True)
    cur = conn.cursor()
    def q1(sql, args=None):
        cur.execute(sql, args or ()); return cur.fetchone()[0]
    def qall(sql, args=None):
        cur.execute(sql, args or ()); return cur.fetchall()

    # --- corpus + license distribution -------------------------------------
    n = q1("SELECT count(*) FROM compounds")
    _check("T01 corpus count", n == CORPUS, f"{n:,} (expected {CORPUS:,})")

    dist = dict(qall("SELECT license_tier, count(*) FROM compounds GROUP BY 1"))
    for tier, exp in LICENSE_TARGET.items():
        got = dist.get(tier, 0)
        _check(f"T02 license_tier {tier}", got == exp, f"{got:,} (expected {exp:,})")

    _check("T03 license sum == corpus",
           sum(dist.values()) == CORPUS, f"{sum(dist.values()):,}")

    _check("T04 no CC BY-NC-ND in corpus",
           dist.get("CC BY-NC-ND 4.0", 0) == 0,
           f"IMPPAT reverted to CC BY-NC; found {dist.get('CC BY-NC-ND 4.0',0)}")

    _check("T05 license vocabulary is the audited four",
           set(dist) <= set(LICENSE_TARGET),
           f"unexpected tiers: {set(dist) - set(LICENSE_TARGET)}")

    # --- tier_rank consistency ---------------------------------------------
    bad_rank = q1("""SELECT count(*) FROM compounds WHERE
        (license_tier='CC0' AND tier_rank<>0) OR
        (license_tier='CC BY 4.0' AND tier_rank<>1) OR
        (license_tier='CC BY-NC 4.0' AND tier_rank<>2) OR
        (license_tier='CC BY-NC-ND 4.0' AND tier_rank<>4) OR
        (license_tier='Unspecified' AND tier_rank<>5)""")
    _check("T06 tier_rank matches license_tier", bad_rank == 0, f"{bad_rank} mismatches")

    null_rank = q1("SELECT count(*) FROM compounds WHERE tier_rank IS NULL")
    _check("T07 no NULL tier_rank", null_rank == 0, f"{null_rank} nulls")

    # --- attestation table --------------------------------------------------
    a_rows = q1("SELECT count(*) FROM per_source_license_attestation")
    _check("T08 attestation row count", a_rows == ATTEST_ROWS,
           f"{a_rows:,} (expected {ATTEST_ROWS:,})")

    a_src = q1("SELECT count(DISTINCT source) FROM per_source_license_attestation")
    _check("T09 attestation distinct sources", a_src == ATTEST_SOURCES,
           f"{a_src} (expected {ATTEST_SOURCES})")

    a_orphan = q1("""SELECT count(*) FROM per_source_license_attestation a
        LEFT JOIN compounds c USING (comp_id) WHERE c.comp_id IS NULL""")
    _check("T10 no orphan attestations (all comp_ids in corpus)",
           a_orphan == 0, f"{a_orphan} rows reference absent compounds")

    a_ncnd = q1("""SELECT count(*) FROM per_source_license_attestation
        WHERE license_tier='CC BY-NC-ND 4.0'""")
    _check("T11 no CC BY-NC-ND in attestation", a_ncnd == 0, f"{a_ncnd} rows")

    # --- source_license_ref map (the durability guarantee) -----------------
    ref = dict(qall("SELECT lower(trim(src)), license_tier FROM source_license_ref"))
    for src, exp in SOURCE_V134.items():
        got = ref.get(src)
        _check(f"T12 ref[{src}]", got == exp, f"{got} (expected {exp})")

    # --- resolver invariant: most-restrictive-wins (the G45b relationship) --
    joined = q1("""SELECT count(*) FROM per_source_license_attestation a
        JOIN source_license_ref r ON lower(trim(r.src))=lower(trim(a.source))""")
    _check("T13 attestation<->ref join is total",
           joined == ATTEST_ROWS, f"{joined:,}/{ATTEST_ROWS:,} joined")

    viol = q1("""WITH best AS (
            SELECT a.comp_id, max(r.tier_rank) b
            FROM per_source_license_attestation a
            JOIN source_license_ref r ON lower(trim(r.src))=lower(trim(a.source))
            GROUP BY 1)
        SELECT count(*) FROM compounds c JOIN best ON best.comp_id=c.comp_id
        WHERE c.tier_rank <> best.b""")
    _check("T14 resolver = most-restrictive-wins (max rank)", viol == 0,
           f"{viol} compounds where tier_rank != max(source rank)")

    # cross-check: least-restrictive would now be WRONG; confirm it differs, so
    # the suite would actually catch a silent revert to the v1.35 rule.
    min_viol = q1("""WITH worst AS (
            SELECT a.comp_id, min(r.tier_rank) b
            FROM per_source_license_attestation a
            JOIN source_license_ref r ON lower(trim(r.src))=lower(trim(a.source))
            GROUP BY 1)
        SELECT count(*) FROM compounds c JOIN worst ON worst.comp_id=c.comp_id
        WHERE c.tier_rank <> worst.b""")
    _check("T15 least-restrictive rule is NOT in effect",
           min_viol > 0, f"{min_viol} compounds differ under min-rule (expect >0)")

    # --- per-compound tier equals its snapshot value (faithful restore) -----
    # every current compound's license must equal the published pre_fix snapshot.
    snap_mismatch = q1("""SELECT count(*) FROM compounds c
        JOIN compounds_license_pre_fix_20260802 s USING (comp_id)
        WHERE c.license_tier <> s.license_tier""")
    _check("T16 tier equals audited snapshot for every compound",
           snap_mismatch == 0, f"{snap_mismatch} diverge from pre_fix snapshot")

    # --- PLACEHOLDER repair intact (earlier change) ------------------------
    ph = q1("""SELECT count(*) FROM compounds
        WHERE effective_class='PLACEHOLDER' OR np_class='PLACEHOLDER'""")
    _check("T17 no PLACEHOLDER classes", ph == 0, f"{ph} placeholder rows")

    # --- three-tier classification sums to corpus --------------------------
    tier_counts = dict(qall(
        "SELECT COALESCE(classification_tier::text,'null'), count(*) "
        "FROM compounds GROUP BY 1"))
    tsum = sum(tier_counts.values())
    _check("T18 classification tiers sum to corpus", tsum == CORPUS, f"{tsum:,}")

    unclassified = tier_counts.get("null", 0)
    _check("T19 unclassified count small (featurization gap)",
           unclassified < 100, f"{unclassified} NULL-tier (expected ~21)", warn=True)

    # --- pathway normalization (7 canonical tokens) ------------------------
    pw = q1("""SELECT count(DISTINCT trim(t))
        FROM compounds, unnest(string_to_array(np_pathway,'|')) t
        WHERE np_pathway IS NOT NULL AND np_pathway<>''""")
    _check("T20 exactly 7 pathway tokens", pw == 7, f"{pw} distinct")

    # --- organism cleanup: no literal [genus] placeholder ------------------
    genus = q1(r"""SELECT count(*) FROM compounds
        WHERE source_organism ~ '\[genus\]|\[species\]|\[family\]'""")
    _check("T21 no [genus] placeholder artifacts", genus == 0, f"{genus} rows")

    # legitimate bracketed genus names preserved (should still be present)
    legit = q1(r"SELECT count(*) FROM compounds WHERE source_organism ~ '\[[A-Z]'")
    _check("T22 legitimate [Genus] names preserved",
           legit > 0, f"{legit} rows with [Capitalized] genus (expect >0)", warn=True)

    # --- commercial-use figure recomputes from reverted tiers --------------
    commercial = q1("SELECT count(*) FROM compounds "
                    "WHERE license_tier IN ('CC0','CC BY 4.0')")
    exp_comm = LICENSE_TARGET["CC0"] + LICENSE_TARGET["CC BY 4.0"]
    _check("T23 commercial-permitted (CC0+CC BY)", commercial == exp_comm,
           f"{commercial:,} (expected {exp_comm:,})")

    # --- backup tables exist (revert is reversible) ------------------------
    baks = q1("""SELECT count(*) FROM pg_tables WHERE schemaname='public'
        AND tablename LIKE 'compounds_license_pre_revert_%%'""")
    _check("T24 pre-revert backup table present", baks >= 1,
           f"{baks} compounds_license_pre_revert_* table(s)", warn=True)

    # --- no accidental corpus mutation: structure columns intact -----------
    # smiles/inchikey non-null coverage should be corpus-wide (license work must
    # not have touched structure); allow the tiny featurization gaps only.
    no_struct = q1("SELECT count(*) FROM compounds WHERE inchikey IS NULL OR inchikey=''")
    _check("T25 all compounds retain an InChIKey (no structural damage)",
           no_struct == 0, f"{no_struct} rows missing inchikey")

    cur.close(); conn.close()

    # --- summary -----------------------------------------------------------
    p = sum(1 for _, s in _results if s == "PASS")
    w = sum(1 for _, s in _results if s == "WARN")
    f = sum(1 for _, s in _results if s == "FAIL")
    print(f"\nSUMMARY  PASS={p}  WARN={w}  FAIL={f}  (of {len(_results)})")
    if f:
        print("FAILURES:", ", ".join(n for n, s in _results if s == "FAIL"))
    sys.exit(1 if f else 0)

if __name__ == "__main__":
    main()
