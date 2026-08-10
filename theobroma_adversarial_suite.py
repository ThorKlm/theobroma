#!/usr/bin/env python3
"""THEOBROMA adversarial suite, 25 checks.

Targets the places where the implementation is known or suspected to be weak,
rather than the happy path. Every documented limitation is asserted explicitly,
so it cannot drift without the suite noticing.

  PASS   invariant holds
  FAIL   regression, investigate before the deposit
  KNOWN  documented limitation, still exactly as documented
  ERROR  the check itself could not run (missing column, route down)

Usage:  python3 theobroma_adversarial_suite.py [--base http://localhost:5000]
Exit 1 if any FAIL or ERROR.
"""
import argparse, json, sys, time, traceback
import urllib.error, urllib.parse, urllib.request

import psycopg2

DSN = "host=localhost dbname=theobroma user=theobroma"
CORPUS = 1132805
MS = {"cc_by": 891860, "cc_by_nc": 225536, "tier1": 626601,
      "families": 486032, "pathway": 1101638, "alkaloid": 426042}

EXEMPLAR = ("type=tax_class&q=liliopsida&extra_type_1=classification"
            "&extra_prop_type_1=npclassifier_class&extra_q_1="
            + urllib.parse.quote("Germacrane sesquiterpenoids")
            + "&extra_type_2=region&extra_q_2=" + urllib.parse.quote("Middle East"))

TESTS = []
COUNT = {"PASS": 0, "FAIL": 0, "KNOWN": 0, "ERROR": 0}


def test(tag, note):
    def deco(fn):
        TESTS.append((tag, note, fn))
        return fn
    return deco


def emit(tag, code, headline, detail=None):
    COUNT[code] += 1
    print("[%-5s] %-4s %s" % (code, tag, headline))
    for line in (detail or [])[:2]:
        if line:
            print("             %s" % line)


class Ctx:
    def __init__(self, con, base):
        self.con, self.base = con, base
        self.cur = con.cursor()

    def q(self, sql, args=None):
        self.cur.execute(sql, args or ())
        return self.cur.fetchall()

    def q1(self, sql, args=None):
        r = self.q(sql, args)
        return r[0][0] if r else None

    def get(self, path, timeout=30):
        try:
            t0 = time.time()
            with urllib.request.urlopen(self.base.rstrip("/") + path, timeout=timeout) as r:
                return r.status, r.read(), time.time() - t0
        except urllib.error.HTTPError as e:
            return e.code, b"", 0.0
        except Exception:
            return 0, b"", 0.0

    def total(self, query):
        st, body, _ = self.get("/api/search?%s&limit=1" % query)
        if st != 200:
            return None
        try:
            return json.loads(body).get("total")
        except Exception:
            return None


# ---------------------------------------------------------------- A integrity

@test("A01", "corpus count")
def a01(c):
    n = c.q1("SELECT count(*) FROM compounds")
    return n == CORPUS, "corpus %s" % format(n, ","), \
        "manuscript prints %s" % format(CORPUS, ",")


@test("A02", "organism binomial regex")
def a02(c):
    n = c.q1("""SELECT count(*) FROM (
        SELECT 1 FROM compounds c
        CROSS JOIN LATERAL unnest(string_to_array(c.source_organism,'; ')) AS r(tok)
        CROSS JOIN LATERAL unnest(string_to_array(c.source_organism_curated,'; ')) AS k(tok)
        WHERE c.source_organism ~ '[a-z]+-[a-z]+'
          AND r.tok LIKE k.tok||'-%%' AND r.tok <> k.tok) x""")
    return n == 0, "hyphen-truncated binomials %s" % n, \
        "nux-vomica and foenum-graecum were truncated across 3,860 compounds"


@test("A03", "classification column whitespace")
def a03(c):
    r = c.q("""SELECT
        count(*) FILTER (WHERE np_pathway ~ '(^ | $|  )'),
        count(*) FILTER (WHERE np_class ~ '(^ | $|  )'),
        count(*) FILTER (WHERE np_superclass ~ '(^ | $|  )'),
        count(*) FILTER (WHERE effective_class ~ '(^ | $|  )') FROM compounds""")[0]
    return sum(r) == 0, "artifacts per column %s" % (tuple(r),), \
        "this once split every class and pathway count into two buckets"


@test("A04", "post-dedup synonym orphans")
def a04(c):
    n = c.q1("""SELECT count(*) FROM compound_synonyms s
        LEFT JOIN compounds x ON x.inchikey = s.inchikey WHERE x.inchikey IS NULL""")
    return ("KNOWN" if n is not None and n < 3000 else "FAIL"), \
        "synonym rows with no matching compound %s" % n, \
        "harmless for queries, prune before the deposit"


# ---------------------------------------------------------------- B taxonomy

@test("B01", "phylum to display-kingdom coverage")
def b01(c):
    r = c.q("""SELECT DISTINCT rt.phylum FROM resolved_taxonomy rt
        LEFT JOIN phylum_kingdom_map m ON m.phylum = rt.phylum
        WHERE rt.phylum IS NOT NULL AND rt.phylum <> '' AND m.phylum IS NULL""")
    return not r, "unmapped phyla %s" % len(r), ", ".join(x[0] for x in r[:5]) or "none"


@test("B02", "taxonomy case consistency")
def b02(c):
    r = c.q("""SELECT count(*) FILTER (WHERE taxclass <> lower(taxclass)),
                      count(*) FILTER (WHERE taxorder <> lower(taxorder)),
                      count(*) FILTER (WHERE phylum   <> lower(phylum))
               FROM resolved_taxonomy""")[0]
    return sum(r) == 0, "mixed-case class/order/phylum %s" % (tuple(r),), \
        "the APG IV table is capitalised, the corpus is not, the SET must lower()"


@test("B03", "monocot class assignment")
def b03(c):
    n = c.q1("""SELECT count(*) FROM resolved_taxonomy
        WHERE taxorder IN ('acorales','alismatales','arecales','asparagales','commelinales',
                           'dioscoreales','liliales','pandanales','poales','zingiberales',
                           'petrosaviales')
          AND taxclass IS DISTINCT FROM 'liliopsida'""")
    return n == 0, "monocot rows not under liliopsida %s" % n, \
        "a gap-only backfill guard let the resolver overwrite the curated table"


@test("B04", "single phylum name for land plants")
def b04(c):
    n = c.q1("SELECT count(*) FROM resolved_taxonomy WHERE phylum='tracheophyta'")
    return n == 0, "tracheophyta rows %s" % n, \
        "two names for one clade drew parallel branches for the same family"


@test("B05", "orders spanning more than one class")
def b05(c):
    n = c.q1("""SELECT count(*) FROM (SELECT taxorder FROM resolved_taxonomy
        WHERE coalesce(taxorder,'')<>'' AND coalesce(taxclass,'')<>''
        GROUP BY 1 HAVING count(DISTINCT taxclass) > 1) x""")
    return ("KNOWN" if n is not None and n <= 26 else "FAIL"), \
        "orders under several classes %s" % n, \
        "eudicots and magnoliids against magnoliopsida, deliberately not decided"


@test("B06", "lineage kingdom against voted kingdom")
def b06(c):
    n = c.q1("""SELECT count(*) FROM resolved_taxonomy rt
        JOIN phylum_kingdom_map m ON m.phylum = rt.phylum
        WHERE m.lineage_kingdom IS DISTINCT FROM rt.kingdom""")
    return ("KNOWN" if n is not None and 4000 < n < 9000 else "FAIL"), \
        "compounds where the two disagree %s" % n, \
        "the tree groups by lineage, Table S1 reports the vote, S2 must say so"


@test("B07", "plant rows still missing a class")
def b07(c):
    n = c.q1("""SELECT count(*) FROM resolved_taxonomy
        WHERE taxclass IS NULL AND coalesce(family,'') <> ''
          AND (phylum IS NULL OR phylum IN ('streptophyta','tracheophyta'))""")
    return ("KNOWN" if n is not None and n < 1000 else "FAIL"), \
        "gap-fill candidates remaining %s" % n, \
        "a large number means the narrowed backfill guard has not been run"


# ---------------------------------------------------------- C classification

@test("C01", "materialized class column integrity")
def c01(c):
    n = c.q1("SELECT count(*) FROM compounds WHERE effective_class = 'PLACEHOLDER'")
    return n == 0, "PLACEHOLDER rows %s" % n, \
        "the validation gate checked np_class while the write went elsewhere"


@test("C02", "ontology primary parent is designated")
def c02(c):
    a = c.q1("SELECT count(*) FROM npc_class_parents WHERE parent_rank IS NULL")
    b = c.q1("""SELECT count(*) FROM (SELECT class_name FROM npc_class_parents
                WHERE parent_rank=1 GROUP BY 1 HAVING count(*)>1) x""")
    return a == 0 and b == 0, "unranked %s, duplicate rank 1 %s" % (a, b), \
        "without parent_rank a VACUUM FULL could silently reclassify compounds"


@test("C03", "model tier stays inside its training vocabulary")
def c03(c):
    n = c.q1("""SELECT count(DISTINCT c.effective_class) FROM compounds c
        WHERE c.class_source='inferred_xgb' AND coalesce(c.effective_class,'') <> ''
          AND NOT EXISTS (SELECT 1 FROM compounds t WHERE t.class_source='curated'
                          AND t.effective_class = c.effective_class)""")
    return n == 0, "model classes absent from the curated set %s" % n


@test("C04", "two pathway denominators, one in the paper")
def c04(c):
    a = c.q1("SELECT count(*) FROM compounds WHERE coalesce(np_pathway,'') <> ''")
    b = c.q1("SELECT count(*) FROM compounds WHERE coalesce(effective_pathway,'') <> ''")
    return ("KNOWN" if a == MS["pathway"] else "FAIL"), \
        "np_pathway %s, effective_pathway %s" % (format(a, ","), format(b, ",")), \
        "Figure S3 quotes the former, regression guard H53 the latter"


# ---------------------------------------------------------------- D licensing

@test("D01", "most-restrictive-wins resolver")
def d01(c):
    n = c.q1("""SELECT count(*) FROM compounds c
        JOIN (SELECT a.comp_id, max(r.tier_rank) AS worst
              FROM per_source_license_attestation a
              JOIN source_license_ref r ON lower(r.src) = lower(a.source)
              GROUP BY 1) s ON s.comp_id = c.comp_id
        WHERE c.tier_rank IS DISTINCT FROM s.worst""")
    return n == 0, "compounds not at the most restrictive attested tier %s" % n, \
        "a least-restrictive resolver would invert the paper's central result"


@test("D02", "licence tier totals")
def d02(c):
    t = dict(c.q("SELECT license_tier, count(*) FROM compounds GROUP BY 1"))
    ok = (t.get("CC BY 4.0") == MS["cc_by"] and t.get("CC BY-NC 4.0") == MS["cc_by_nc"]
          and sum(t.values()) == CORPUS)
    return ok, "tiers sum %s" % format(sum(t.values()), ","), \
        "CC BY %s, CC BY-NC %s" % (t.get("CC BY 4.0"), t.get("CC BY-NC 4.0"))


@test("D03", "attestation join must not match vacuously")
def d03(c):
    n = c.q1("""SELECT count(*) FROM per_source_license_attestation a
        LEFT JOIN source_license_ref r ON lower(r.src) = lower(a.source)
        WHERE r.src IS NULL""")
    return n == 0, "attestations with no licence-map match %s" % n, \
        "a case-sensitive join once made the licence test pass on zero rows"


# ------------------------------------------------------- E routes and filters

@test("E01", "routes a reviewer might open")
def e01(c):
    paths = ["/", "/search?type=name&q=fisetin", "/browse", "/tree", "/statistics",
             "/annotate", "/api/stats", "/api/docs", "/compound/THEO_0858442",
             "/api/compound/THEO_0858442",
             "/api/compound/THEO_0858442/license-provenance",
             "/api/stereoisomers/THEO_0858442",
             "/api/taxonomy_tree?type=tax_class&q=liliopsida"]
    bad = [(p, s) for p, s in ((p, c.get(p)[0]) for p in paths) if s != 200]
    return not bad, "non-200 routes %s of %s" % (len(bad), len(paths)), \
        "; ".join("%s %s" % (p, s) for p, s in bad[:2]) or "all 200"


@test("E02", "classification filter as an extra must filter")
def e02(c):
    p = c.total("type=npclassifier_class&q=Flavonols")
    x = c.total("type=kingdom&q=plant&extra_type_1=classification"
                "&extra_prop_type_1=npclassifier_class&extra_q_1=Flavonols")
    k = c.total("type=kingdom&q=plant")
    ok = bool(p and x and k) and x <= p and x < k
    return ok, "primary %s, extra under plant %s, plant alone %s" % (p, x, k), \
        "the extra was silently dropped in three of five filter implementations"


@test("E03", "two filter names unusable as a primary type")
def e03(c):
    a = c.total("type=npclassifier_pathway&q=Alkaloids")
    b = c.total("type=classyfire_superclass&q="
                + urllib.parse.quote("Alkaloids and derivatives"))
    return ("KNOWN" if not a and not b else "FAIL"), \
        "npclassifier_pathway %s, classyfire_superclass %s" % (a, b), \
        "primary maps key these as pathway and classyfire_class; aliases not added"


@test("E04", "source semantics differ between API and tree")
def e04(c):
    api = c.total("type=source&q=LOTUS")
    any_ = c.q1("SELECT count(*) FROM compounds WHERE all_sources ILIKE '%%LOTUS%%'")
    return ("KNOWN" if api is not None and any_ and api < any_ else "FAIL"), \
        "api primary-attribution %s vs any-attestation %s" % (api, any_), \
        "the tree's reading matches the paper's provenance framing; the API should change"


@test("E05", "figure exemplar reproduces from its own caption")
def e05(c):
    n = c.total(EXEMPLAR)
    return n == 30, "exemplar returns %s compounds" % n, \
        "the published exemplar returned 5,030 against a caption claiming 25"


@test("E06", "exemplar tree payload carries compound leaves")
def e06(c):
    st, body, dt = c.get("/api/taxonomy_tree?" + EXEMPLAR)
    n = body.count(b'"is_compound": true') + body.count(b'"is_compound":true')
    return st == 200 and n >= 30, "%s compound leaves in %.2fs" % (n, dt), \
        "if this passes while the page is blank, the render fault is client-side"


# ---------------------------------------------------- F manuscript regression

@test("F01", "figures printed in the manuscript")
def f01(c):
    fams = c.q1("SELECT count(DISTINCT left(inchikey,14)) FROM compounds")
    t1 = c.q1("SELECT count(*) FROM compounds WHERE class_source='curated'")
    alk = c.q1("SELECT count(*) FROM compounds WHERE split_part(np_pathway,'|',1)='Alkaloids'")
    ok = fams == MS["families"] and t1 == MS["tier1"] and alk == MS["alkaloid"]
    return ok, "families %s, curated %s, alkaloid primary %s" % (
        format(fams or 0, ","), format(t1 or 0, ","), format(alk or 0, ",")), \
        "expected %s / %s / %s" % (format(MS["families"], ","),
                                   format(MS["tier1"], ","), format(MS["alkaloid"], ","))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://localhost:5000")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    con = psycopg2.connect(DSN)
    con.set_session(readonly=True, autocommit=True)
    ctx = Ctx(con, a.base)
    print("THEOBROMA adversarial suite  %d checks  base=%s" % (len(TESTS), a.base))
    print("-" * 74)
    for tag, note, fn in TESTS:
        try:
            r = fn(ctx)
            code = r[0] if isinstance(r[0], str) else ("PASS" if r[0] else "FAIL")
            emit(tag, code, r[1], list(r[2:]))
        except Exception as e:
            emit(tag, "ERROR", "%s: %s" % (type(e).__name__, str(e).strip()[:70]), [note])
            if a.verbose:
                traceback.print_exc()
            try:
                ctx.cur.close()
                ctx.cur = con.cursor()
            except Exception:
                pass
    print("-" * 74)
    print("SUMMARY  PASS=%d FAIL=%d KNOWN=%d ERROR=%d" % (
        COUNT["PASS"], COUNT["FAIL"], COUNT["KNOWN"], COUNT["ERROR"]))
    con.close()
    return 1 if COUNT["FAIL"] or COUNT["ERROR"] else 0


if __name__ == "__main__":
    sys.exit(main())
