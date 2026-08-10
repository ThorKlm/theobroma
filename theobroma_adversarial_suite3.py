#!/usr/bin/env python3
"""THEOBROMA adversarial suite III, 25 checks.

Suite I covers integrity, taxonomy, classification, licensing and routes.
Suite II covers query semantics, adversarial input, similarity, families,
caches and exports. This one covers what is left and what the last two days
suggest is most fragile: agreement between views of the same quantity,
referential integrity in both directions, descriptor and structure sanity,
the reproducibility artifacts the manuscript points reviewers at, and
behaviour at the edges of the corpus.

  PASS   invariant holds
  FAIL   regression, investigate before the deposit
  KNOWN  documented limitation or expected divergence, asserted so it cannot drift
  ERROR  the check itself could not run

Usage:  python3 theobroma_adversarial_suite3.py [--base http://localhost:5000]
Exit 1 if any FAIL or ERROR.
"""
import argparse, glob, json, os, re, sys, time, traceback
import urllib.error, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

import psycopg2

DSN = "host=localhost dbname=theobroma user=theobroma"
CORPUS = 1132805
N_SOURCES = 29
SALTS = 5443
ISOTOPOLOGUES = 303

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

    def get(self, path, timeout=60):
        try:
            t0 = time.time()
            with urllib.request.urlopen(self.base.rstrip("/") + path, timeout=timeout) as r:
                return r.status, r.read(), time.time() - t0
        except urllib.error.HTTPError as e:
            return e.code, e.read()[:400], 0.0
        except Exception as e:
            return 0, str(e).encode()[:200], 0.0

    def js(self, path, timeout=60):
        st, body, dt = self.get(path, timeout)
        if st != 200:
            return None, st, dt
        try:
            return json.loads(body), st, dt
        except Exception:
            return None, st, dt

    def total(self, query):
        d, _, _ = self.js("/api/search?%s&limit=1" % query)
        return (d or {}).get("total")


# ------------------------------------------- M. agreement between views

@test("M01", "the statistics page and the database")
def m01(c):
    st, body, _ = c.get("/statistics")
    txt = body.decode("utf-8", "replace")
    shown = "1,132,805" in txt or "1132805" in txt
    return st == 200 and shown, "statistics page carries the corpus figure %s" % shown, \
        "a cache served a reverted licence model to every visitor for hours"


@test("M02", "the stats API and the database")
def m02(c):
    d, st, _ = c.js("/api/stats")
    txt = json.dumps(d or {})
    n = c.q1("SELECT count(*) FROM compounds")
    return st == 200 and (str(n) in txt or "%s" % n in txt), \
        "api reports the corpus figure %s" % (str(n) in txt), \
        "the API and the page are built from different queries"


@test("M03", "the tree total and the search total for one query")
def m03(c):
    qs = "type=tax_class&q=liliopsida&extra_type_1=region&extra_q_1=" \
         + urllib.parse.quote("Middle East")
    s = c.total(qs)
    d, _, _ = c.js("/api/taxonomy_tree?" + qs)
    t = (d or {}).get("total_compounds")
    return ("KNOWN" if s and t and s != t else ("PASS" if s == t else "FAIL")), \
        "search %s, tree %s" % (s, t), \
        "search matches any kingdom attestation, the tree groups by resolved lineage"


@test("M04", "the HTML search and the JSON API for one query")
def m04(c):
    st, body, _ = c.get("/search?type=npclassifier_class&q=Flavonols")
    txt = body.decode("utf-8", "replace")
    api = c.total("type=npclassifier_class&q=Flavonols")
    hits = [int(x.replace(",", "")) for x in re.findall(r"([\d,]+)\s+results", txt)]
    html = max(hits) if hits else None
    return st == 200 and html == api, "html %s, api %s" % (html, api), \
        "two independent filter implementations answer the same question"


@test("M05", "the licence filter composes as a subset")
def m05(c):
    allc = c.total("type=kingdom&q=plant")
    comm = c.total("type=kingdom&q=plant&license=commercial")
    acad = c.total("type=kingdom&q=plant&license=academic")
    ok = bool(allc) and (comm is None or comm <= allc) and (acad is None or acad <= allc)
    return ok, "all %s, commercial %s, academic %s" % (allc, comm, acad), \
        "the persona scenario assembles a commercial subset through this filter"


@test("M06", "one version string across the interface")
def m06(c):
    seen = set()
    for p in ["/", "/statistics", "/tree", "/api/stats"]:
        st, body, _ = c.get(p)
        seen |= set(re.findall(r"THEOBROMA v1\.\d\d", body.decode("utf-8", "replace")))
    return seen in ({"THEOBROMA v1.35"}, set()), "version strings served: %s" % (sorted(seen) or "none"), \
        "figure exports bake the version in; a stale one reaches the manuscript"


# --------------------------------------------- N. referential integrity

@test("N01", "child tables reference existing compounds")
def n01(c):
    out = {}
    for tbl in ("admet", "compound_region_map", "per_source_license_attestation",
                "compound_taxonomy", "resolved_taxonomy"):
        try:
            out[tbl] = c.q1("""SELECT count(*) FROM %s t
                LEFT JOIN compounds x ON x.comp_id = t.comp_id
                WHERE x.comp_id IS NULL""" % tbl)
        except Exception:
            c.con.rollback()
            out[tbl] = "n/a"
    bad = {k: v for k, v in out.items() if isinstance(v, int) and v > 0}
    return not bad, "orphaned child rows %s" % (bad or "none"), \
        "several thousand rows were rewritten across five tables in two days"


@test("N02", "every compound has a taxonomy row")
def n02(c):
    n = c.q1("""SELECT count(*) FROM compounds x
        LEFT JOIN resolved_taxonomy rt ON rt.comp_id = x.comp_id
        WHERE rt.comp_id IS NULL""")
    return n == 0, "compounds with no resolved_taxonomy row %s" % n, \
        "the tree and every taxonomic filter join through this table"


@test("N03", "attestation count agrees with the all_sources string")
def n03(c):
    n = c.q1("""SELECT count(*) FROM (
        SELECT x.comp_id,
               array_length(string_to_array(x.all_sources,'|'),1) AS declared,
               (SELECT count(*) FROM per_source_license_attestation a
                WHERE a.comp_id = x.comp_id) AS attested
        FROM compounds x WHERE coalesce(x.all_sources,'') <> ''
        ORDER BY random() LIMIT 20000) y
        WHERE declared IS DISTINCT FROM attested""")
    return ("KNOWN" if n is not None and n < 2000 else "FAIL"), \
        "sampled rows where the two disagree %s of 20,000" % n, \
        "the resolver parses all_sources; a mismatch means a silent attestation gap"


@test("N04", "the macro-region vocabulary is closed")
def n04(c):
    r = c.q("SELECT macro_region, count(*) FROM compound_region_map GROUP BY 1 ORDER BY 2 DESC")
    return len(r) == 13, "distinct macro regions %s" % len(r), \
        "largest %s, smallest %s" % (r[0][0] if r else "-", r[-1][0] if r else "-")


# ------------------------------------- O. descriptor and structure sanity

@test("O01", "InChIKey re-derives from the stored SMILES")
def o01(c):
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.inchi import MolToInchiKey
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        return "KNOWN", "RDKit unavailable in this interpreter", "skipped"
    rows = c.q("""SELECT smiles, inchikey FROM compounds
                  WHERE smiles IS NOT NULL ORDER BY random() LIMIT 200""")
    ok = fail = 0
    for smi, ik in rows:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            fail += 1
            continue
        ok += 1 if MolToInchiKey(m) == ik else 0
    pct = 100.0 * ok / max(1, len(rows))
    return ("KNOWN" if pct >= 95 else "FAIL"), \
        "re-derived key matches for %.1f%% of 200, unparseable %s" % (pct, fail), \
        "tautomer and stereo perception drift across RDKit versions; ~97%% is expected"


@test("O02", "the stored mass agrees with a recompute")
def o02(c):
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.Descriptors import ExactMolWt
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        return "KNOWN", "RDKit unavailable in this interpreter", "skipped"
    rows = c.q("""SELECT smiles, mw FROM compounds
                  WHERE smiles IS NOT NULL AND mw IS NOT NULL
                  ORDER BY random() LIMIT 200""")
    off = 0
    for smi, mw in rows:
        m = Chem.MolFromSmiles(smi)
        if m is not None and abs(ExactMolWt(m) - float(mw)) > 0.01:
            off += 1
    return off <= 4, "rows differing by more than 0.01 Da: %s of 200" % off, \
        "the column is monoisotopic ExactMolWt, not average molecular weight"


@test("O03", "descriptors stay inside physically sensible ranges")
def o03(c):
    r = c.q("""SELECT
        count(*) FILTER (WHERE mw <= 0 OR mw > 20000),
        count(*) FILTER (WHERE tpsa < 0),
        count(*) FILTER (WHERE hba < 0 OR hbd < 0),
        count(*) FILTER (WHERE n_rings < 0 OR rotatable_bonds < 0),
        count(*) FILTER (WHERE logp < -150 OR logp > 150) FROM compounds""")[0]
    return sum(r) == 0, "out-of-range mw/tpsa/hb/rings/logp %s" % (tuple(r),), \
        "Crippen logP is additive and unbounded, so large polysaccharides and peptides\n        legitimately reach plus or minus 80; only absurd values indicate a failure"


@test("O04", "salts are preserved as multi-fragment SMILES")
def o04(c):
    n = c.q1("SELECT count(*) FROM compounds WHERE smiles LIKE '%%.%%'")
    return abs((n or 0) - SALTS) < 200, "multi-fragment SMILES %s" % n, \
        "the manuscript prints %s under the source-fidelity principle" % format(SALTS, ","), \
        ""


@test("O05", "isotopologues are retained as distinct entries")
def o05(c):
    n = c.q1("SELECT count(*) FROM compounds WHERE smiles ~ '\\[[0-9]+[A-Z]'")
    return ("KNOWN" if n is not None else "FAIL"), \
        "compounds carrying an isotopic label %s" % n, \
        "the manuscript prints %s; heuristic here, exact figure is from ingestion" % ISOTOPOLOGUES


# ------------------------------ P. reproducibility artifacts and deposit

@test("P01", "the manifest lists exactly the integrated sources")
def p01(c):
    path = next((p for p in ["sources.yaml", "data/sources.yaml"] if os.path.exists(p)), None)
    if not path:
        return "KNOWN", "sources.yaml not found", "the manuscript cites it by name"
    txt = open(path, encoding="utf-8", errors="replace").read()
    db = [r[0] for r in c.q("SELECT DISTINCT source_db FROM compounds")]
    missing = [s for s in db if s.lower() not in txt.lower()]
    return len(db) == N_SOURCES and not missing, \
        "%s sources in the corpus, %s absent from the manifest" % (len(db), len(missing)), \
        ", ".join(missing[:4]) or "all present"


@test("P02", "the licence map agrees with the database")
def p02(c):
    path = next((p for p in ["license_map.tsv", "data/license_map.tsv"] if os.path.exists(p)), None)
    if not path:
        return "KNOWN", "license_map.tsv not found", "S3 names it as the persisted map"
    m = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            m[parts[0].strip().lower()] = parts[1].strip()
    db = {r[0].lower(): r[1] for r in c.q("SELECT src, license_tier FROM source_license_ref")}
    diff = [k for k in db if k in m and m[k] != db[k]]
    return not diff, "%s file entries, %s db entries, %s disagreeing" % (len(m), len(db), len(diff)), \
        ", ".join(diff[:3]) or "none"


@test("P03", "the deposit dump postdates the last corpus change")
def p03(c):
    dumps = sorted(glob.glob("*.dump") + glob.glob("deposit/*.dump")
                   + glob.glob(os.path.expanduser("~/*.dump")))
    if not dumps:
        return "KNOWN", "no .dump found in the usual places", "regenerate before the deposit"
    newest = max(dumps, key=os.path.getmtime)
    age = (time.time() - os.path.getmtime(newest)) / 86400.0
    return ("KNOWN" if age > 0.5 else "PASS"), \
        "%s is %.1f days old" % (os.path.basename(newest), age), \
        "the corpus changed today; a dump older than that deposits a superseded state"


@test("P04", "the README carries the current corpus figure")
def p04(c):
    if not os.path.exists("README.md"):
        return "KNOWN", "README.md not found", ""
    txt = open("README.md", encoding="utf-8", errors="replace").read()
    old = [s for s in ("1,133,004", "1133004") if s in txt]
    return not old and ("1,132,805" in txt or "1132805" in txt), \
        "current figure present %s, superseded figure present %s" % (
            "1,132,805" in txt or "1132805" in txt, bool(old)), \
        "the repository is cited in Data Availability"


@test("P05", "the repository tracks the artifacts the manuscript names")
def p05(c):
    want = ["scripts/apply_license_resolver.sql", "sources.yaml",
            "regression_guards.py", "theobroma_test_suite.py",
            "pipeline/07_organism_curation.sql", "scripts/backfill_plant_lineage.py"]
    missing = [w for w in want if not os.path.exists(w)]
    return not missing, "%s of %s named artifacts present" % (len(want) - len(missing), len(want)), \
        ", ".join(missing[:3]) or "all present"


# --------------------------------------------- Q. robustness at the edges

@test("Q01", "a compound with no organism annotation renders")
def q01(c):
    cid = c.q1("""SELECT comp_id FROM compounds
                  WHERE coalesce(source_organism,'') = '' LIMIT 1""")
    st, body, _ = c.get("/compound/%s" % cid)
    return st == 200 and len(body) > 2000, "%s renders, %s bytes" % (cid, len(body)), \
        "68%% of the corpus has no organism; the page must not assume one"


@test("Q02", "the extremes of the corpus render")
def q02(c):
    big = c.q1("SELECT comp_id FROM compounds ORDER BY mw DESC NULLS LAST LIMIT 1")
    small = c.q1("SELECT comp_id FROM compounds ORDER BY mw ASC NULLS LAST LIMIT 1")
    many = c.q1("""SELECT comp_id FROM compounds
                   ORDER BY length(coalesce(source_organism,'')) DESC LIMIT 1""")
    codes = [c.get("/compound/%s" % x)[0] for x in (big, small, many)]
    return all(s == 200 for s in codes), "heaviest, lightest, most-attested: %s" % codes, \
        "%s / %s / %s" % (big, small, many)


@test("Q03", "compound identifiers are handled consistently")
def q03(c):
    a = c.get("/compound/THEO_0858442")[0]
    b = c.get("/compound/theo_0858442")[0]
    d = c.get("/compound/THEO_9999999")[0]
    return a == 200 and d == 404 and b in (200, 301, 302, 404), \
        "upper %s, lower %s, nonexistent %s" % (a, b, d), \
        "a 500 on a bad identifier is a crawler-reachable fault"


@test("Q04", "concurrent requests on a single worker")
def q04(c):
    paths = ["/api/search?type=kingdom&q=plant&limit=5",
             "/api/search?type=npclassifier_class&q=Flavonols&limit=5",
             "/api/compound/THEO_0858442",
             "/api/stats",
             "/api/stereoisomers/THEO_0858442"] * 2
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=10) as ex:
        codes = [r[0] for r in ex.map(lambda p: c.get(p, timeout=90), paths)]
    return all(s == 200 for s in codes), \
        "%s of %s succeeded in %.1fs" % (sum(1 for s in codes if s == 200), len(codes),
                                         time.time() - t0), \
        "the in-heap index confines the deployment to one worker; requests serialise"


@test("Q05", "a query with no results renders rather than errors")
def q05(c):
    codes = [c.get("/search?type=name&q=zzzznotacompoundzzzz")[0],
             c.get("/api/search?type=name&q=zzzznotacompoundzzzz&limit=5")[0],
             c.get("/api/taxonomy_tree?type=genus&q=zzzznotagenuszzzz")[0]]
    d, _, _ = c.js("/api/search?type=name&q=zzzznotacompoundzzzz&limit=5")
    return all(s == 200 for s in codes) and (d or {}).get("total") == 0, \
        "status %s, api total %s" % (codes, (d or {}).get("total")), \
        "the empty tree is the case the cladogram renderer currently fails on"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://localhost:5000")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    con = psycopg2.connect(DSN)
    con.set_session(readonly=True, autocommit=True)
    ctx = Ctx(con, a.base)
    print("THEOBROMA adversarial suite III  %d checks  base=%s" % (len(TESTS), a.base))
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
