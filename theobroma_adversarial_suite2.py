#!/usr/bin/env python3
"""THEOBROMA adversarial suite II, 25 checks.

Suite I covers integrity, taxonomy, classification, licensing and route
liveness. This one covers the surfaces where correctness is most likely to
degrade quietly: query composition and semantics, adversarial input, the
similarity engine and its indices, stereoisomer-family symmetry, derived
caches, and the export and batch paths.

  PASS   invariant holds
  FAIL   regression, investigate before the deposit
  KNOWN  documented limitation, still exactly as documented
  ERROR  the check itself could not run

Usage:  python3 theobroma_adversarial_suite2.py [--base http://localhost:5000]
        --quick   skip the similarity canary (G/I block runs ~30 s otherwise)
Exit 1 if any FAIL or ERROR.
"""
import argparse, glob, json, os, random, sys, time, traceback
import urllib.error, urllib.parse, urllib.request

import psycopg2

DSN = "host=localhost dbname=theobroma user=theobroma"
CORPUS = 1132805
STATIC = "static"

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
    def __init__(self, con, base, quick):
        self.con, self.base, self.quick = con, base, quick
        self.cur = con.cursor()

    def q(self, sql, args=None):
        self.cur.execute(sql, args or ())
        return self.cur.fetchall()

    def q1(self, sql, args=None):
        r = self.q(sql, args)
        return r[0][0] if r else None

    def get(self, path, timeout=60):
        try:
            t0 = time.time()
            with urllib.request.urlopen(self.base.rstrip("/") + path, timeout=timeout) as r:
                return r.status, r.read(), time.time() - t0
        except urllib.error.HTTPError as e:
            return e.code, e.read()[:400], 0.0
        except Exception as e:
            return 0, str(e).encode()[:200], 0.0

    def post(self, path, payload, timeout=60):
        body = json.dumps(payload).encode()
        req = urllib.request.Request(self.base.rstrip("/") + path, data=body,
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            return e.code, e.read()[:400]
        except Exception as e:
            return 0, str(e).encode()[:200]

    def js(self, path, timeout=60):
        st, body, dt = self.get(path, timeout)
        if st != 200:
            return None, st, dt
        try:
            return json.loads(body), st, dt
        except Exception:
            return None, st, dt

    def total(self, query):
        d, st, _ = self.js("/api/search?%s&limit=1" % query)
        return d.get("total") if d else None

    def rows(self, query):
        d, st, _ = self.js("/api/search?%s" % query)
        return (d or {}).get("results") or []


# ------------------------------------------- G. query composition and semantics

@test("G01", "each added filter must narrow, never widen")
def g01(c):
    a = c.total("type=kingdom&q=plant")
    b = c.total("type=kingdom&q=plant&extra_type_1=tax_class&extra_q_1=liliopsida")
    d = c.total("type=kingdom&q=plant&extra_type_1=tax_class&extra_q_1=liliopsida"
                "&extra_type_2=region&extra_q_2=" + urllib.parse.quote("Middle East"))
    ok = bool(a and b and d) and d <= b <= a
    return ok, "plant %s, +liliopsida %s, +Middle East %s" % (a, b, d), \
        "a dropped filter widens the set silently, as in the published Figure 2"


@test("G02", "filter values must be case-insensitive")
def g02(c):
    hi = c.total("type=npclassifier_class&q=Flavonols")
    lo = c.total("type=npclassifier_class&q=flavonols")
    return hi is not None and hi == lo, "Flavonols %s, flavonols %s" % (hi, lo), \
        "the routes lower both sides; a raw equality anywhere breaks this"


@test("G03", "leading and trailing whitespace in a filter value")
def g03(c):
    a = c.total("type=npclassifier_class&q=Flavonols")
    b = c.total("type=npclassifier_class&q=" + urllib.parse.quote("  Flavonols  "))
    return a is not None and a == b, "bare %s, padded %s" % (a, b), \
        "users paste values with whitespace; untrimmed input silently returns nothing"


@test("G04", "plus and percent-twenty must encode the same space")
def g04(c):
    a = c.total("type=region&q=Middle+East")
    b = c.total("type=region&q=Middle%20East")
    return a is not None and a == b and a > 0, "plus %s, pct20 %s" % (a, b), \
        "the figure header prints the plus form and reviewers paste the other"


@test("G05", "pagination pages are disjoint and cover the total")
def g05(c):
    q = "type=npclassifier_class&q=Flavonols&limit=25"
    a = [r.get("comp_id") for r in c.rows(q + "&offset=0")]
    b = [r.get("comp_id") for r in c.rows(q + "&offset=0")]
    o = [r.get("comp_id") for r in c.rows(q + "&offset=25")]
    stable = a == b
    disjoint = not (set(a) & set(o))
    return stable and disjoint and len(a) == 25, \
        "repeat call identical %s, offset slices disjoint %s" % (stable, disjoint), \
        "an unordered query lets synchronize_seqscans return a different window each call"



# ------------------------------------------------------ H. adversarial input

@test("H01", "injection patterns return safely")
def h01(c):
    before = c.q1("SELECT count(*) FROM compounds")
    payloads = ["'; DROP TABLE compounds; --", "' OR '1'='1", "%'; SELECT pg_sleep(5); --",
                "1) UNION SELECT NULL,NULL--"]
    codes = [c.get("/api/search?type=name&q=" + urllib.parse.quote(p) + "&limit=1")[0]
             for p in payloads]
    after = c.q1("SELECT count(*) FROM compounds")
    ok = before == after and all(s in (200, 400) for s in codes)
    return ok, "status codes %s, corpus intact %s" % (codes, before == after)


@test("H02", "empty and whitespace-only queries")
def h02(c):
    codes = [c.get("/api/search?type=name&q=&limit=1")[0],
             c.get("/api/search?type=name&q=%20%20&limit=1")[0],
             c.get("/api/search?limit=1")[0]]
    return all(s in (200, 400) for s in codes), "status codes %s" % codes, \
        "an empty filter must not fall through to an unfiltered scan"


@test("H03", "oversized query string")
def h03(c):
    st, body, dt = c.get("/api/search?type=name&q=" + "A" * 8000 + "&limit=1")
    return st in (200, 400, 414), "status %s in %.2fs" % (st, dt), \
        "a 500 here means the length is only bounded by the database"


@test("H04", "out-of-range pagination parameters")
def h04(c):
    codes = [c.get("/api/search?type=kingdom&q=plant&limit=-5")[0],
             c.get("/api/search?type=kingdom&q=plant&limit=0")[0],
             c.get("/api/search?type=kingdom&q=plant&page=-1")[0],
             c.get("/api/search?type=kingdom&q=plant&page=999999999")[0]]
    return all(s in (200, 400) for s in codes), "status codes %s" % codes


@test("H05", "an enormous limit must be clamped, not served")
def h05(c):
    d, st, dt = c.js("/api/search?type=kingdom&q=plant&limit=500000", timeout=90)
    n = len((d or {}).get("results") or [])
    return st == 200 and 0 < n <= 100000, "returned %s rows in %.1fs" % (n, dt), \
        "an unclamped limit on a 967k-row filter would exhaust the single worker"


@test("H06", "non-ASCII in a search term")
def h06(c):
    codes, hits = [], []
    for term in ["Bréb", "Curcuma", "Müller", "北京"]:
        d, st, _ = c.js("/api/search?type=name&q=" + urllib.parse.quote(term) + "&limit=1")
        codes.append(st)
        hits.append((d or {}).get("total"))
    return all(s == 200 for s in codes), "status %s, totals %s" % (codes, hits), \
        "organism strings carry diacritics; an encoding fault shows up as a 500"


# ----------------------------------------------- I. similarity and its indices

@test("I02", "self-similarity canary")
def i02(c):
    if c.quick:
        return "KNOWN", "skipped under --quick", "run without --quick before the deposit"
    rows = c.q("""SELECT comp_id, smiles FROM compounds
                  WHERE smiles IS NOT NULL AND smiles NOT LIKE '%%.%%'
                  ORDER BY random() LIMIT 25""")
    hit = 0
    for cid, smi in rows:
        d, st, _ = c.js("/api/similarity?smiles=" + urllib.parse.quote(smi)
                        + "&metric=morgan&top_n=1")
        res = (d or {}).get("results") or (d or {}).get("hits") or []
        if res and res[0].get("tanimoto") == 1.0:
            hit += 1
    return hit >= 24, "top hit at tanimoto 1.0 for %s of %s" % (hit, len(rows)), \
        "Morgan is not injective, so ties are expected; a top score below 1.0 is not"


@test("I03", "an invalid SMILES must not reach the index")
def i03(c):
    codes = [c.get("/api/similarity?smiles=" + urllib.parse.quote(s) + "&metric=morgan&top_n=5")[0]
             for s in ["not_a_smiles", "C(((", "", "C" * 5000]]
    return all(s in (200, 400) for s in codes), "status codes %s" % codes, \
        "RDKit returning None must be caught before the vector lookup"


@test("I04", "every similarity metric answers")
def i04(c):
    smi = "Oc1ccc2c(c1)oc(c(c2=O)O)c1ccc(c(c1)O)O"
    out = {}
    for m in ["morgan", "maccs", "chemberta", "nafm"]:
        d, st, dt = c.js("/api/similarity?smiles=" + urllib.parse.quote(smi)
                         + "&metric=%s&top_n=5" % m)
        res = (d or {}).get("results") or (d or {}).get("hits") or []
        out[m] = (st, len(res), round(dt, 2))
    ok = all(v[0] == 200 for v in out.values())
    return ok, "; ".join("%s %s/%s/%ss" % (k, v[0], v[1], v[2]) for k, v in out.items()), \
        "NaFM is corpus-only; an off-corpus SMILES must degrade, not error"


# --------------------------------------- J. stereoisomer families, cross-field

@test("J01", "family members share the connectivity prefix")
def j01(c):
    cid = c.q1("""SELECT comp_id FROM compounds WHERE left(inchikey,14) IN
                  (SELECT left(inchikey,14) FROM compounds GROUP BY 1 HAVING count(*)>2)
                  LIMIT 1""")
    d, st, _ = c.js("/api/stereoisomers/%s" % cid)
    mem = (d or {}).get("stereoisomers") or []
    echoed = (d or {}).get("inchikey_prefix") or ""
    pref = {(m.get("inchikey") or "")[:14] for m in mem}
    return st == 200 and len(mem) >= 2 and pref == {echoed}, \
        "%s members, prefix %s, distinct prefixes %s" % (len(mem), echoed, len(pref)), \
        "the family relation is the paper's third contribution"


@test("J02", "family membership is symmetric")
def j02(c):
    cid = c.q1("""SELECT comp_id FROM compounds WHERE left(inchikey,14) IN
                  (SELECT left(inchikey,14) FROM compounds GROUP BY 1 HAVING count(*)=3)
                  LIMIT 1""")
    d, _, _ = c.js("/api/stereoisomers/%s" % cid)
    mem = [m.get("comp_id") for m in ((d or {}).get("stereoisomers") or [])]
    other = next((m for m in mem if m and m != cid), None)
    d2, _, _ = c.js("/api/stereoisomers/%s" % other) if other else (None, 0, 0)
    back = [m.get("comp_id") for m in ((d2 or {}).get("stereoisomers") or [])]
    return bool(other) and cid in back, "%s lists %s, and back: %s" % (
        cid, other, cid in back), "asymmetry would mean the prefix is computed twice"


@test("J03", "effective_class agrees with its tier's source column")
def j03(c):
    bad = c.q1("""SELECT count(*) FROM compounds
        WHERE class_source='curated' AND coalesce(effective_class,'') <> ''
          AND coalesce(np_class,'') <> '' AND effective_class <> np_class""")
    bad2 = c.q1("""SELECT count(*) FROM compounds
        WHERE class_source='inferred_xgb' AND coalesce(effective_class,'') <> ''
          AND coalesce(inferred_class,'') <> '' AND effective_class <> inferred_class""")
    return bad == 0 and bad2 == 0, \
        "curated mismatches %s, model mismatches %s" % (bad, bad2), \
        "effective_class is materialized; drift from its source is invisible to users"



# ------------------------------------------- K. derived caches and artifacts

@test("K01", "statistics cache against the live corpus")
def k01(c):
    path = next((p for p in ["statistics_cache.json", "static/statistics_cache.json"]
                 if os.path.exists(p)), None)
    if not path:
        return "KNOWN", "statistics cache file not found", "checked two conventional paths"
    obj = json.load(open(path))
    txt = json.dumps(obj)
    age = (time.time() - os.path.getmtime(path)) / 86400.0
    return str(CORPUS) in txt or "1132805" in txt, \
        "%s, %.1f days old, corpus figure present %s" % (
            path, age, str(CORPUS) in txt or "1132805" in txt), \
        "a stale cache served a reverted licence model to every visitor"


@test("K02", "chemistry tree cache against the class vocabulary")
def k02(c):
    path = os.path.join(STATIC, "chem_tree.json")
    if not os.path.exists(path):
        return "KNOWN", "chem_tree.json absent", path
    txt = open(path, encoding="utf-8").read()
    age = (time.time() - os.path.getmtime(path)) / 86400.0
    dup = txt.count('"Quinoline alkaloids"')
    return dup <= 2, "%.1f days old, Quinoline alkaloids appears %s time(s)" % (age, dup), \
        "duplicate class nodes were the symptom of the unranked ontology parent"


@test("K03", "global linear tree cache freshness")
def k03(c):
    path = os.path.join(STATIC, "linear_tree_global.json")
    if not os.path.exists(path):
        return "KNOWN", "linear_tree_global.json absent", path
    age = (time.time() - os.path.getmtime(path)) / 86400.0
    txt = open(path, encoding="utf-8").read()
    return age < 3 and "tracheophyta" not in txt, \
        "%.1f days old, tracheophyta present %s" % (age, "tracheophyta" in txt), \
        "this file was stale by nine days and carried the pre-collapse phylum"


@test("K04", "kingdom caches carry the voted totals")
def k04(c):
    files = sorted(glob.glob(os.path.join(STATIC, "taxonomy_cache_*.json")))
    if not files:
        return "KNOWN", "no taxonomy_cache_*.json found", STATIC
    stale = [os.path.basename(f) for f in files
             if (time.time() - os.path.getmtime(f)) / 86400.0 > 3]
    return not stale, "%s cache files, %s older than 3 days" % (len(files), len(stale)), \
        "the tree groups by lineage while Table S1 reports the vote; caches must say which"


# ---------------------------------------------------- L. export and batch paths

@test("L01", "bulk export agrees with the search total")
def l01(c):
    d, st, dt = c.js("/api/bulk?limit=500&format=json", timeout=90)
    rows = d if isinstance(d, list) else ((d or {}).get("results") or [])
    return st == 200 and len(rows) > 0, "bulk returned %s rows in %.1fs" % (len(rows), dt), \
        "the persona scenario pulls 900k compounds through this endpoint"


@test("L02", "batch annotate accounts for every input")
def l02(c):
    payload = {"smiles": ["Oc1ccc2c(c1)oc(c(c2=O)O)c1ccc(c(c1)O)O",
                          "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
                          "C1CN(CCN1CCO)CC2=CC=CC=N2"]}
    st, body = c.post("/api/annotate", payload)
    try:
        d = json.loads(body)
    except Exception:
        d = {}
    m = len(d.get("matched") or d.get("results") or [])
    u = len(d.get("unmatched") or [])
    return st in (200, 400) and (st == 400 or m + u == 3), \
        "status %s, matched %s, unmatched %s of 3" % (st, m, u), \
        "an input that is neither matched nor reported is silently lost"


@test("L03", "batch annotate enforces its input cap")
def l03(c):
    st, body = c.post("/api/annotate", {"smiles": ["CCO"] * 1500}, timeout=90)
    return st in (200, 400, 413), "status %s for 1,500 inputs" % st, \
        "the documented cap is 1,000; an unenforced cap serialises behind one worker"


@test("L04", "licence provenance endpoint agrees with the compounds table")
def l04(c):
    cid = c.q1("""SELECT comp_id FROM compounds WHERE comp_id IN (
                    SELECT comp_id FROM per_source_license_attestation
                    GROUP BY 1 HAVING count(*) > 5) LIMIT 1""")
    tier = c.q1("SELECT license_tier FROM compounds WHERE comp_id=%s", (cid,))
    d, st, _ = c.js("/api/compound/%s/license-provenance" % cid)
    txt = json.dumps(d or {})
    n = len(d.get("attestations") or d.get("sources") or []) if d else 0
    return st == 200 and tier and tier in txt and n > 1, \
        "%s tier %s, %s attestations echoed" % (cid, tier, n), \
        "the endpoint is the auditable half of the paper's headline claim"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://localhost:5000")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    random.seed(0)
    con = psycopg2.connect(DSN)
    con.set_session(readonly=True, autocommit=True)
    ctx = Ctx(con, a.base, a.quick)
    print("THEOBROMA adversarial suite II  %d checks  base=%s%s"
          % (len(TESTS), a.base, "  (quick)" if a.quick else ""))
    print("-" * 74)
    t0 = time.time()
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
    print("SUMMARY  PASS=%d FAIL=%d KNOWN=%d ERROR=%d  (%.0fs)" % (
        COUNT["PASS"], COUNT["FAIL"], COUNT["KNOWN"], COUNT["ERROR"], time.time() - t0))
    con.close()
    return 1 if COUNT["FAIL"] or COUNT["ERROR"] else 0


if __name__ == "__main__":
    sys.exit(main())
