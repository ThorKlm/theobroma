#!/usr/bin/env python3
"""THEOBROMA manuscript reconciliation, 25 checks.

Reports every printed figure in the manuscript against the value the live
database returns, with the query that produced it. This suite makes no
assertions about which value is correct: several printed figures are measured
against intermediate pipeline states that the current database no longer holds,
and replacing them would silently change their meaning. The output is evidence
for a decision, not a decision.

  MATCH      live value equals the printed one
  DIFFERS    live value differs; the row shows both and the delta
  CONTEXT    reported for interpretation, nothing printed to compare against
  ERROR      the check could not run

Usage:  python3 manuscript_reconcile.py [--verbose]
Exit 0 always. This suite reports, it does not gate.
"""
import argparse, sys, traceback

import psycopg2

DSN = "host=localhost dbname=theobroma user=theobroma"

CHECKS = []
COUNT = {"MATCH": 0, "DIFFERS": 0, "CONTEXT": 0, "ERROR": 0}


def check(tag, printed, where):
    """printed: the value in the manuscript, or None for CONTEXT rows.
    where:   where in the manuscript it appears."""
    def deco(fn):
        CHECKS.append((tag, printed, where, fn))
        return fn
    return deco


def fmt(v):
    if isinstance(v, int):
        return format(v, ",")
    if isinstance(v, float):
        return "%.4g" % v
    return str(v)


def emit(tag, code, where, printed, live, note=None):
    COUNT[code] += 1
    if code == "CONTEXT":
        print("[%-7s] %-4s %-34s live %s" % (code, tag, where, fmt(live)))
    elif code == "ERROR":
        print("[%-7s] %-4s %-34s %s" % (code, tag, where, live))
    else:
        delta = ""
        if isinstance(printed, int) and isinstance(live, int):
            d = live - printed
            delta = "  delta %+d" % d if d else ""
        print("[%-7s] %-4s %-34s printed %-14s live %-14s%s"
              % (code, tag, where, fmt(printed), fmt(live), delta))
    if note:
        print("            %s" % note)


class Ctx:
    def __init__(self, con):
        self.con = con
        self.cur = con.cursor()

    def q1(self, sql, args=None):
        self.cur.execute(sql, args or ())
        r = self.cur.fetchone()
        return r[0] if r else None

    def q(self, sql, args=None):
        self.cur.execute(sql, args or ())
        return self.cur.fetchall()


# ------------------------------------------------------------------ A. corpus

@check("A01", 1132805, "abstract, throughout")
def a01(c):
    return c.q1("SELECT count(*) FROM compounds")


@check("A02", 486032, "abstract, Section 2, Table 1")
def a02(c):
    return c.q1("SELECT count(DISTINCT left(inchikey,14)) FROM compounds")


@check("A04", 29, "Section 2, Table 2")
def a04(c):
    return c.q1("SELECT count(DISTINCT source_db) FROM compounds")


@check("A05", 998699, "Section 2 P1, aggregator group")
def a05(c):
    return c.q1("""SELECT count(*) FROM compounds
                   WHERE source_db IN ('COCONUT','LOTUS','NPASS','NPAtlas')""")


@check("A06", 58775, "Section 2 P1, food and metabolomics group")
def a06(c):
    return c.q1("""SELECT count(*) FROM compounds
                   WHERE source_db IN ('FooDB','LMDB_Lichen','CMDB_Cereals')""")


# --------------------------------------------------------------- B. licensing

@check("B01", 891860, "Section 3, abstract, Figure S1")
def b01(c):
    return c.q1("SELECT count(*) FROM compounds WHERE license_tier='CC BY 4.0'")


@check("B02", 8243, "Section 3, Figure S1")
def b02(c):
    return c.q1("SELECT count(*) FROM compounds WHERE license_tier='CC0'")


@check("B03", 225536, "Section 3, Figure S1")
def b03(c):
    return c.q1("SELECT count(*) FROM compounds WHERE license_tier='CC BY-NC 4.0'")


@check("B04", 7166, "Section 3, Figure S1")
def b04(c):
    return c.q1("SELECT count(*) FROM compounds WHERE license_tier='Unspecified'")


@check("B05", 900103, "abstract, Figure 1, persona scenario")
def b05(c):
    return c.q1("SELECT count(*) FROM compounds WHERE tier_rank <= 1")


@check("B06", 1512594, "S3, attestation table")
def b06(c):
    return c.q1("SELECT count(*) FROM per_source_license_attestation")


@check("B07", 137149, "Section 3, licence interval")
def b07(c):
    return c.q1("SELECT count(*) FROM compounds WHERE tier_rank_min <> tier_rank")


@check("B08", 1014300, "Section 3, least restrictive permissive")
def b08(c):
    return c.q1("SELECT count(*) FROM compounds WHERE tier_rank_min = 1")


@check("B09", 18925, "Section 3, least restrictive public domain")
def b09(c):
    return c.q1("SELECT count(*) FROM compounds WHERE tier_rank_min = 0")


@check("B10", 108769, "S3, transition ledger total")
def b10(c):
    live = c.q1("""SELECT count(*) FROM compounds c
                   JOIN source_license_ref r ON lower(r.src) = lower(c.source_db)
                   WHERE c.tier_rank > r.tier_rank""")
    return live, ("the printed figure is measured against the pre-correction "
                  "first-stage map that S3 itself flags as unpublished; the live "
                  "value uses the current map, so the two are different objects")


@check("B11", None, "S3 ledger, live transition breakdown")
def b11(c):
    rows = c.q("""SELECT rp.tier_rank, c.tier_rank, count(*)
                  FROM compounds c
                  JOIN source_license_ref rp ON lower(rp.src) = lower(c.source_db)
                  WHERE c.tier_rank > rp.tier_rank
                  GROUP BY 1,2 ORDER BY 3 DESC""")
    return "; ".join("%d->%d %s" % (a, b, format(n, ",")) for a, b, n in rows)


@check("B12", None, "Table 2, first-stage tier totals")
def b12(c):
    rows = c.q("""SELECT r.license_tier, count(*)
                  FROM compounds c JOIN source_license_ref r
                    ON lower(r.src) = lower(c.source_db)
                  GROUP BY 1 ORDER BY 2 DESC""")
    return "; ".join("%s %s" % (t, format(n, ",")) for t, n in rows)


# ----------------------------------------------------------- C. classification

@check("C01", 626601, "Section 4.1, curated tier")
def c01(c):
    return c.q1("SELECT count(*) FROM compounds WHERE class_source='curated'")


@check("C02", 206042, "Section 4.1, tool tier")
def c02(c):
    return c.q1("SELECT count(*) FROM compounds WHERE class_source='npclassifier'")


@check("C03", 300141, "Section 4.1, model tier")
def c03(c):
    return c.q1("SELECT count(*) FROM compounds WHERE class_source='inferred_xgb'")


@check("C04", 1101638, "Section 4.2, Figure S3 denominator")
def c04(c):
    return c.q1("SELECT count(*) FROM compounds WHERE coalesce(np_pathway,'') <> ''")


@check("C05", 426042, "Section 4.1, Figure S3 alkaloid wedge")
def c05(c):
    return c.q1("""SELECT count(*) FROM compounds
                   WHERE split_part(np_pathway,'|',1) = 'Alkaloids'""")


# ------------------------------------------------------------- D. annotations

@check("D01", 1131184, "Section 4.2, Figure S4")
def d01(c):
    return c.q1("SELECT count(*) FROM compounds WHERE sa_score IS NOT NULL")



# ---------------------------------------------------------------- E. taxonomy

@check("E01", None, "Table S1, kingdom distribution")
def e01(c):
    rows = c.q("""SELECT kingdom, count(*) FROM resolved_taxonomy
                  WHERE coalesce(kingdom,'') <> '' GROUP BY 1 ORDER BY 2 DESC""")
    return "; ".join("%s %s" % (k, format(n, ",")) for k, n in rows)


@check("E02", 264039, "Discussion, geographic coverage")
def e02(c):
    return c.q1("SELECT count(DISTINCT comp_id) FROM compound_region_map")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    con = psycopg2.connect(DSN)
    con.set_session(readonly=True, autocommit=True)
    ctx = Ctx(con)
    print("THEOBROMA manuscript reconciliation  %d checks" % len(CHECKS))
    print("-" * 96)
    for tag, printed, where, fn in CHECKS:
        try:
            r = fn(ctx)
            note = None
            if isinstance(r, tuple):
                r, note = r
            if printed is None:
                emit(tag, "CONTEXT", where, None, r, note)
            elif r == printed:
                emit(tag, "MATCH", where, printed, r, note)
            else:
                emit(tag, "DIFFERS", where, printed, r, note)
        except Exception as e:
            emit(tag, "ERROR", where, None, "%s: %s" % (type(e).__name__, str(e)[:60]))
            if a.verbose:
                traceback.print_exc()
            try:
                ctx.cur.close()
                ctx.cur = con.cursor()
            except Exception:
                pass
    print("-" * 96)
    print("SUMMARY  MATCH=%d DIFFERS=%d CONTEXT=%d ERROR=%d"
          % (COUNT["MATCH"], COUNT["DIFFERS"], COUNT["CONTEXT"], COUNT["ERROR"]))
    con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
