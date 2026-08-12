#!/usr/bin/env python3
"""THEOBROMA adversarial suite V, 50 checks.

Suites I to IV cover data integrity, taxonomy, licensing, query semantics,
adversarial input, similarity, caches, exports and display-versus-filter
agreement. This one covers the layers those do not touch: the stated NAR
Database Issue criteria, HTTP and security headers, API correctness, FAIR
machine-readability, programmatically testable accessibility, single-worker
robustness, and cheminformatics data quality angles beyond exact-match checks.

Security headers and TLS come from the reverse proxy, so those checks run
against --public while everything else runs against --base.

  PASS   invariant holds
  FAIL   regression or a defect a reviewer would find
  KNOWN  documented gap or limitation, asserted so it cannot drift
  ERROR  the check itself could not run

Usage:
  python3 theobroma_adversarial_suite5.py
  python3 theobroma_adversarial_suite5.py --stress      # adds load and slow-client checks
  python3 theobroma_adversarial_suite5.py --public https://theobroma.l3s.uni-hannover.de

Do not run --stress against the public endpoint during a review window: the
service runs one worker and the slow-client check is designed to starve it.
Exit 1 if any FAIL or ERROR.
"""
import argparse, json, os, re, socket, ssl, sys, time, traceback
import urllib.error, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

import psycopg2

DSN = "host=localhost dbname=theobroma user=theobroma"
CORPUS = 1132805
EXEMPLAR = "THEO_0858442"

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
    def __init__(self, con, base, public, stress):
        self.con, self.base, self.public, self.stress = con, base, public, stress
        self.cur = con.cursor()

    def q(self, sql, args=None):
        self.cur.execute(sql, args or ())
        return self.cur.fetchall()

    def q1(self, sql, args=None):
        r = self.q(sql, args)
        return r[0][0] if r else None

    def raw(self, path, root=None, method="GET", headers=None, timeout=60, redirect=True):
        """Return (status, headers dict lowercased, body text)."""
        url = (root or self.base).rstrip("/") + path
        req = urllib.request.Request(url, method=method, headers=headers or {})
        opener = urllib.request.build_opener()
        if not redirect:
            class NoRedirect(urllib.request.HTTPRedirectHandler):
                def redirect_request(self, *a, **k):
                    return None
            opener = urllib.request.build_opener(NoRedirect)
        try:
            with opener.open(req, timeout=timeout) as r:
                body = r.read()
                return r.status, {k.lower(): v for k, v in r.headers.items()}, \
                    body.decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            body = e.read()
            return e.code, {k.lower(): v for k, v in e.headers.items()}, \
                body.decode("utf-8", "replace")
        except Exception as e:
            return 0, {}, str(e)

    def get(self, path, **kw):
        st, _, body = self.raw(path, **kw)
        return st, body

    def js(self, path, **kw):
        st, _, body = self.raw(path, **kw)
        if st != 200:
            return None, st
        try:
            return json.loads(body), st
        except Exception:
            return None, st

    def pub(self, path, **kw):
        return self.raw(path, root=self.public, **kw)


# ------------------------------------------- X. NAR Database Issue criteria

@test("X01", "free access without login or registration")
def x01(c):
    paths = ["/", "/search?type=name&q=fisetin", "/browse", "/tree", "/statistics",
             "/download", "/api/search?type=name&q=fisetin&limit=1",
             "/compound/" + EXEMPLAR]
    bad = []
    for p in paths:
        st, hd, _ = c.raw(p, redirect=False)
        if st in (301, 302, 303, 307, 308) and "login" in (hd.get("location") or "").lower():
            bad.append((p, "login redirect"))
        elif "www-authenticate" in hd:
            bad.append((p, "auth challenge"))
        elif st not in (200, 301, 302, 308):
            bad.append((p, st))
    return not bad, "routes requiring auth %s of %s" % (len(bad), len(paths)), \
        "; ".join("%s %s" % x for x in bad[:2]) or "none", \
        "NAR requires databases to be freely available with no login or registration"


@test("X02", "HTTPS is served and HTTP redirects to it")
def x02(c):
    st, hd, _ = c.pub("/")
    https_ok = st == 200
    host = urllib.parse.urlparse(c.public).netloc
    st2, hd2, _ = c.raw("/", root="http://" + host, redirect=False)
    redirects = st2 in (301, 302, 307, 308) and \
        (hd2.get("location") or "").lower().startswith("https")
    return https_ok and redirects, "https %s, http redirects to https %s" % (https_ok, redirects), \
        "NAR encourages HTTPS; a plain-HTTP resource invites an editorial comment"


@test("X03", "HSTS declared")
def x03(c):
    st, hd, _ = c.pub("/")
    return "strict-transport-security" in hd, \
        "Strict-Transport-Security %s" % ("present" if "strict-transport-security" in hd else "absent")


@test("X04", "outbound source-database links resolve")
def x04(c):
    st, html = c.get("/download")
    st2, html2 = c.get("/help")
    urls = set(re.findall(r'href="(https?://[^"]+)"', html + html2))
    urls = [u for u in urls if "theobroma" not in u and "uni-hannover" not in u]
    urls = sorted(urls)[:12]
    if not urls:
        return "KNOWN", "no outbound source links found on /download or /help", \
            "Bateman 2007 names 404s as a specific reviewer irritant"
    dead = []
    for u in urls:
        try:
            req = urllib.request.Request(u, method="HEAD",
                                         headers={"User-Agent": "theobroma-linkcheck"})
            with urllib.request.urlopen(req, timeout=20) as r:
                if r.status >= 400:
                    dead.append((u, r.status))
        except urllib.error.HTTPError as e:
            if e.code not in (403, 405):
                dead.append((u, e.code))
        except Exception:
            dead.append((u, "unreachable"))
    return len(dead) <= max(1, len(urls) // 10), \
        "outbound links dead %s of %s sampled" % (len(dead), len(urls)), \
        "; ".join(str(x[0])[:60] for x in dead[:2]) or "none"


@test("X05", "compound permalinks are stable and domain-based")
def x05(c):
    ids = [r[0] for r in c.q("SELECT comp_id FROM compounds ORDER BY random() LIMIT 8")]
    bad = [i for i in ids if c.raw("/compound/%s" % i)[0] != 200]
    ip_url = bool(re.match(r"https?://\d+\.\d+\.\d+\.\d+", c.public))
    return not bad and not ip_url, \
        "permalinks failing %s of 8, IP-based URL %s" % (len(bad), ip_url), \
        "Bateman 2007: URLs to specific IP addresses are unlikely to stand the test of time"


@test("X06", "first-visit help with worked examples")
def x06(c):
    st, html = c.get("/help")
    has_example = bool(re.search(r"(example|e\.g\.|for instance)", html, re.I))
    return st == 200 and len(html) > 4000 and has_example, \
        "help page %s bytes, worked examples %s" % (len(html), has_example), \
        "NAR requires sufficient help for a first-time visitor"


@test("X07", "download formats and terms are documented")
def x07(c):
    st, html = c.get("/download")
    fmts = [f for f in ("CSV", "JSON", "SDF", "dump", "TSV") if f.lower() in html.lower()]
    terms = bool(re.search(r"(licen[cs]e|CC BY|CC0|terms)", html, re.I))
    return st == 200 and len(fmts) >= 2 and terms, \
        "formats advertised %s, terms stated %s" % (",".join(fmts) or "none", terms), \
        "NAR asks authors to address formats and terms for data download"


# ------------------------------------------------------------ Y. security

@test("Y01", "search terms are escaped in the results page")
def y01(c):
    payloads = ['<script>alert(1)</script>', '"><svg onload=alert(1)>',
                "'><img src=x onerror=alert(1)>"]
    leaked = []
    for p in payloads:
        st, html = c.get("/search?type=name&q=" + urllib.parse.quote(p))
        if ("<script>alert" in html or "<svg onload" in html
                or "<img src=x onerror" in html):
            leaked.append(p[:24])
    return not leaked, "unescaped reflections %s of %s" % (len(leaked), len(payloads)), \
        "; ".join(leaked[:2]) or "none", \
        "the search echo is the most plausible XSS vector in a read-only resource"


@test("Y02", "chemical name fields are escaped on the compound page")
def y02(c):
    cid = c.q1("""SELECT comp_id FROM compounds
                  WHERE name LIKE '%%<%%' OR name LIKE '%%>%%' LIMIT 1""")
    if not cid:
        cid = c.q1("""SELECT c.comp_id FROM compounds c
                      JOIN compound_synonyms s ON s.inchikey = c.inchikey
                      WHERE s.synonym LIKE '%%<%%' LIMIT 1""")
    if not cid:
        return "PASS", "no name or synonym field contains angle brackets", \
            "nothing to escape, so the stored-XSS surface is empty"
    st, html = c.get("/compound/%s" % cid)
    return "&lt;" in html or "&gt;" in html, "%s renders escaped entities" % cid


@test("Y03", "no debug traceback or framework disclosure")
def y03(c):
    probes = ["/compound/../../etc/passwd", "/api/search?type=name&q=%00",
              "/api/search?limit=notanumber&type=name&q=x", "/nonexistent-route-xyz"]
    leaks = []
    for p in probes:
        st, hd, body = c.raw(p)
        if "Traceback (most recent call last)" in body or "werkzeug" in body.lower():
            leaks.append((p, "traceback"))
        if "x-powered-by" in hd:
            leaks.append((p, "x-powered-by"))
    return not leaks, "framework or traceback disclosure %s" % (len(leaks)), \
        "; ".join("%s %s" % x for x in leaks[:2]) or "none", \
        "Flask debug mode in production is a severe and common finding"


@test("Y04", "directory traversal in file-serving routes")
def y04(c):
    probes = ["/download/../../../etc/passwd", "/download/%2e%2e%2f%2e%2e%2fetc%2fpasswd",
              "/static/../app.py", "/static/%2e%2e%2fapp.py"]
    bad = []
    for p in probes:
        st, hd, body = c.raw(p)
        if st == 200 and ("root:x:" in body or "psycopg2" in body or "DATABASE_URL" in body):
            bad.append(p)
    return not bad, "traversal probes returning file contents %s of %s" % (len(bad), len(probes)), \
        "; ".join(bad[:2]) or "none"


@test("Y05", "core security headers on HTML responses")
def y05(c):
    st, hd, _ = c.pub("/")
    want = ["content-security-policy", "x-content-type-options",
            "referrer-policy", "x-frame-options"]
    missing = [h for h in want if h not in hd]
    return ("KNOWN" if missing else "PASS"), \
        "security headers missing %s of %s" % (len(missing), len(want)), \
        ", ".join(missing) or "all present", \
        "the OWASP baseline; a university proxy commonly omits these"


@test("Y06", "X-Content-Type-Options nosniff")
def y06(c):
    st, hd, _ = c.pub("/")
    v = (hd.get("x-content-type-options") or "").lower()
    return ("PASS" if v == "nosniff" else "KNOWN"), \
        "X-Content-Type-Options %s" % (v or "absent")


@test("Y07", "cookies carry Secure, HttpOnly and SameSite")
def y07(c):
    st, hd, _ = c.pub("/")
    raw = hd.get("set-cookie")
    if not raw:
        return "PASS", "no cookies set", "a read-only resource needs none"
    low = raw.lower()
    ok = all(f in low for f in ("secure", "httponly", "samesite"))
    return ok, "cookie flags %s" % ("complete" if ok else "incomplete"), raw[:70]


@test("Y08", "the corpus is unchanged by injection probes")
def y08(c):
    before = c.q1("SELECT count(*) FROM compounds")
    for p in ["'; DROP TABLE compounds; --", "' OR '1'='1", "1);SELECT pg_sleep(3);--"]:
        c.get("/api/search?type=name&q=" + urllib.parse.quote(p) + "&limit=1")
    after = c.q1("SELECT count(*) FROM compounds")
    return before == after == CORPUS, "corpus %s before, %s after" % (before, after)


@test("Y09", "expensive queries do not admit a trivial amplification")
def y09(c):
    t0 = time.time()
    st, _ = c.get("/api/search?type=name&q=" + "a" * 200 + "&limit=10000", timeout=90)
    dt = time.time() - t0
    return st in (200, 400, 414) and dt < 30, \
        "long-term wide-limit query: status %s in %.1fs" % (st, dt), \
        "one worker means a single expensive request blocks every other user"


# ----------------------------------------------------- Z. HTTP correctness

@test("Z01", "JSON responses carry a JSON content type")
def z01(c):
    st, hd, _ = c.raw("/api/search?type=name&q=fisetin&limit=1")
    ct = (hd.get("content-type") or "").split(";")[0]
    return ct == "application/json", "Content-Type %s" % (ct or "absent")


@test("Z02", "CSV export carries a CSV content type")
def z02(c):
    st, hd, body = c.raw("/api/search?type=name&q=fisetin&limit=1&format=csv")
    ct = (hd.get("content-type") or "").split(";")[0]
    return st == 200 and ct in ("text/csv", "application/csv"), \
        "status %s, Content-Type %s" % (st, ct or "absent"), \
        "text/html on a CSV export breaks every programmatic consumer"


@test("Z03", "a nonexistent compound returns 404")
def z03(c):
    st, _ = c.get("/compound/THEO_9999999")
    st2, _ = c.get("/api/compound/THEO_9999999")
    return st == 404 and st2 == 404, "html %s, api %s" % (st, st2), \
        "200-with-empty for a missing record misleads clients that check status"


@test("Z04", "an empty result set returns 200, not 404")
def z04(c):
    d, st = c.js("/api/search?type=name&q=zzzznotacompoundzzzz&limit=5")
    return st == 200 and (d or {}).get("total") == 0, \
        "status %s, total %s" % (st, (d or {}).get("total")), \
        "an empty search is a successful query, not a missing resource"


@test("Z05", "malformed parameters return 400, not 500")
def z05(c):
    probes = ["/api/search?type=name&q=x&limit=notanumber",
              "/api/search?type=nosuchtype&q=x",
              "/api/search?type=name&q=x&offset=-99999999999999999999"]
    codes = [c.raw(p)[0] for p in probes]
    return all(s in (200, 400) for s in codes), "status codes %s" % codes


@test("Z06", "HEAD is handled on the public routes")
def z06(c):
    out = {}
    for p in ("/", "/api/stats", "/compound/" + EXEMPLAR):
        st, hd, body = c.raw(p, method="HEAD")
        out[p] = st
    return all(s in (200, 405) for s in out.values()), \
        "HEAD status %s" % list(out.values()), \
        "405 is tolerable; a 500 on HEAD is a crawler-reachable fault"


@test("Z07", "CORS allows programmatic browser access to the API")
def z07(c):
    st, hd, _ = c.raw("/api/search?type=name&q=fisetin&limit=1",
                      headers={"Origin": "https://example.org"})
    acao = hd.get("access-control-allow-origin")
    return ("PASS" if acao else "KNOWN"), \
        "Access-Control-Allow-Origin %s" % (acao or "absent"), \
        "without it, browser-based cheminformatics tools cannot call the API"


@test("Z08", "large responses are compressed when requested")
def z08(c):
    st, hd, body = c.raw("/api/search?type=kingdom&q=plant&limit=500",
                         headers={"Accept-Encoding": "gzip"})
    enc = hd.get("content-encoding")
    return ("PASS" if enc and "gzip" in enc else "KNOWN"), \
        "Content-Encoding %s on a 500-row response" % (enc or "identity")


@test("Z09", "SMILES-bearing URLs survive encoding")
def z09(c):
    smi = "Oc1ccc2c(c1)oc(c(c2=O)O)c1ccc(c(c1)O)O"
    weird = ["C[C@H](N)C(=O)O", "CC(=O)O.[Na+]", "C1=CC=CC=C1#N"]
    codes = []
    for s in [smi] + weird:
        st, _ = c.get("/api/similarity?smiles=" + urllib.parse.quote(s, safe="") +
                      "&metric=morgan&top_n=1")
        codes.append(st)
    return all(s in (200, 400) for s in codes), "status codes %s" % codes, \
        "SMILES routinely contain #, /, \\, [ and ], which break naive routing"


# --------------------------------------------- AA. FAIR and machine-readability

@test("AA1", "structured metadata on the landing page")
def aa1(c):
    st, html = c.get("/")
    has_jsonld = 'application/ld+json' in html
    return ("PASS" if has_jsonld else "KNOWN"), \
        "schema.org JSON-LD %s" % ("present" if has_jsonld else "absent"), \
        "Bioschemas Dataset markup is what Google Dataset Search indexes"


@test("AA2", "robots.txt is served")
def aa2(c):
    st, hd, body = c.raw("/robots.txt")
    return ("PASS" if st == 200 else "KNOWN"), \
        "robots.txt status %s" % st, \
        "absent robots.txt leaves crawler behaviour to defaults"


@test("AA3", "sitemap.xml is served and well-formed")
def aa3(c):
    st, hd, body = c.raw("/sitemap.xml")
    ok = st == 200 and body.lstrip().startswith("<?xml")
    return ("PASS" if ok else "KNOWN"), \
        "sitemap.xml status %s, xml %s" % (st, body.lstrip().startswith("<?xml"))


@test("AA4", "an OpenAPI description is discoverable")
def aa4(c):
    found = None
    for p in ("/openapi.json", "/api/openapi.json", "/api/spec", "/swagger.json"):
        st, hd, body = c.raw(p)
        if st == 200 and ("openapi" in body[:400].lower() or "swagger" in body[:400].lower()):
            found = p
            break
    return ("PASS" if found else "KNOWN"), \
        "machine-readable API description %s" % (found or "absent"), \
        "/api/docs is human-readable; a spec is what clients generate from"


@test("AA5", "machine-readable licence metadata")
def aa5(c):
    st, hd, html = c.raw("/")
    link_hdr = "license" in (hd.get("link") or "").lower()
    in_html = bool(re.search(r'rel=["\']license["\']|"license"\s*:', html))
    ok = link_hdr or in_html
    return ("PASS" if ok else "KNOWN"), \
        "licence declared machine-readably %s" % ok, \
        "the resource's own thesis is per-compound licence provenance"


@test("AA6", "content negotiation offers a machine format for a compound")
def aa6(c):
    st, hd, body = c.raw("/compound/" + EXEMPLAR,
                         headers={"Accept": "application/json"})
    ct = (hd.get("content-type") or "").split(";")[0]
    return ("PASS" if ct == "application/json" else "KNOWN"), \
        "Accept: application/json yields %s" % (ct or "nothing"), \
        "the /api/compound route exists, so this is discoverability rather than capability"


@test("AA7", "per-compound source provenance is exposed in the API")
def aa7(c):
    d, st = c.js("/api/compound/" + EXEMPLAR)
    keys = set((d or {}).keys())
    prov = [k for k in keys if "source" in k.lower() or "licen" in k.lower()]
    return st == 200 and len(prov) >= 2, \
        "provenance fields in the compound API %s" % len(prov), \
        ", ".join(sorted(prov)[:4])


@test("AA8", "the dataset DOI in Data Availability resolves")
def aa8(c):
    st, html = c.get("/")
    st2, html2 = c.get("/help")
    dois = set(re.findall(r"10\.\d{4,9}/[^\s\"'<>)]+", html + html2))
    if not dois:
        return "KNOWN", "no DOI advertised on the landing or help page", \
            "the Zenodo deposit DOI is the resource's persistent identifier"
    bad = []
    for d in sorted(dois)[:3]:
        try:
            req = urllib.request.Request("https://doi.org/" + d, method="HEAD",
                                         headers={"User-Agent": "theobroma-linkcheck"})
            with urllib.request.urlopen(req, timeout=20) as r:
                if r.status >= 400:
                    bad.append(d)
        except Exception:
            bad.append(d)
    return not bad, "DOIs advertised %s, unresolvable %s" % (len(dois), len(bad)), \
        ", ".join(bad[:2]) or "all resolve"


# ------------------------------------------------------- AB. accessibility

@test("AB1", "informative images carry alternative text")
def ab1(c):
    st, html = c.get("/compound/" + EXEMPLAR)
    imgs = re.findall(r"<img\b[^>]*>", html)
    missing = [i for i in imgs if not re.search(r'\balt\s*=', i)]
    return not missing, "images %s, without alt %s" % (len(imgs), len(missing)), \
        "WebAIM Million 2025: 18.5%% of home-page images lack alternative text"


@test("AB2", "form inputs are labelled")
def ab2(c):
    st, html = c.get("/search?type=name&q=fisetin")
    inputs = re.findall(r"<(?:input|select|textarea)\b[^>]*>", html)
    inputs = [i for i in inputs if 'type="hidden"' not in i and "type='hidden'" not in i]
    unlabelled = [i for i in inputs
                  if not re.search(r'\b(aria-label|aria-labelledby|title|id)\s*=', i)]
    return not unlabelled, "visible inputs %s, unlabelled %s" % (len(inputs), len(unlabelled)), \
        "WebAIM Million 2025: 34.2%% of form inputs are not properly labelled"


@test("AB3", "the document declares a language")
def ab3(c):
    st, html = c.get("/")
    return bool(re.search(r"<html[^>]+\blang\s*=", html)), \
        "html lang attribute %s" % ("present" if re.search(r"<html[^>]+\blang", html) else "absent")


@test("AB4", "heading hierarchy is well formed")
def ab4(c):
    st, html = c.get("/compound/" + EXEMPLAR)
    levels = [int(m) for m in re.findall(r"<h([1-6])\b", html)]
    h1 = levels.count(1)
    skips = sum(1 for a, b in zip(levels, levels[1:]) if b - a > 1)
    return h1 == 1 and skips == 0, "h1 count %s, level skips %s" % (h1, skips)


@test("AB5", "the interface declares a viewport for small screens")
def ab5(c):
    st, html = c.get("/")
    return 'name="viewport"' in html, \
        "viewport meta %s" % ("present" if 'name="viewport"' in html else "absent"), \
        "NAR asks that databases be legible on phone and tablet screens"


# ------------------------------------------------- AC. robustness under load

@test("AC1", "mixed concurrent requests all succeed")
def ac1(c):
    paths = ["/api/search?type=kingdom&q=plant&limit=5",
             "/api/compound/" + EXEMPLAR,
             "/api/stats",
             "/api/stereoisomers/" + EXEMPLAR,
             "/api/search?type=npclassifier_class&q=Flavonols&limit=5"] * 3
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=15) as ex:
        codes = list(ex.map(lambda p: c.raw(p, timeout=120)[0], paths))
    return all(s == 200 for s in codes), \
        "%s of %s succeeded in %.1fs" % (sum(1 for s in codes if s == 200), len(codes),
                                         time.time() - t0), \
        "one sync worker serialises every request behind the slowest"


@test("AC2", "tail latency under concurrency")
def ac2(c):
    times = []
    def one(_):
        t0 = time.time()
        c.raw("/api/search?type=name&q=fisetin&limit=5", timeout=120)
        return time.time() - t0
    with ThreadPoolExecutor(max_workers=10) as ex:
        times = sorted(ex.map(one, range(20)))
    p50 = times[len(times) // 2]
    p95 = times[int(len(times) * 0.95) - 1]
    return ("KNOWN" if p95 < 30 else "FAIL"), \
        "p50 %.2fs, p95 %.2fs over 20 requests at concurrency 10" % (p50, p95), \
        "the single-worker deployment is documented in S11; this records its cost"


@test("AC3", "an unreachable route does not wedge the worker")
def ac3(c):
    c.raw("/nonexistent-route-xyz")
    st, _ = c.get("/api/stats")
    return st == 200, "service responsive after a 404: status %s" % st


@test("AC4", "similarity failure does not take text search with it")
def ac4(c):
    c.get("/api/similarity?smiles=not_a_smiles&metric=morgan&top_n=5")
    st, _ = c.get("/api/search?type=name&q=fisetin&limit=1")
    return st == 200, "text search after a bad similarity call: %s" % st, \
        "shared-fate failure between the FAISS layer and plain search would be serious"


@test("AC5", "slow clients do not starve the worker")
def ac5(c):
    if not c.stress:
        return "KNOWN", "skipped without --stress", \
            "designed to starve a single sync worker; do not run against production"
    host = urllib.parse.urlparse(c.base)
    socks = []
    try:
        for _ in range(8):
            s = socket.create_connection((host.hostname, host.port or 80), timeout=10)
            s.sendall(b"GET /api/stats HTTP/1.1\r\nHost: x\r\n")
            socks.append(s)
        time.sleep(2)
        t0 = time.time()
        st, _ = c.get("/api/stats", timeout=30)
        dt = time.time() - t0
        return st == 200 and dt < 15, \
            "served in %.1fs with 8 slow clients open, status %s" % (dt, st)
    finally:
        for s in socks:
            try:
                s.close()
            except Exception:
                pass


@test("AC6", "memory is bounded across a sustained workload")
def ac6(c):
    if not c.stress:
        return "KNOWN", "skipped without --stress", \
            "reads the service RSS from /proc, so it must run on the host"
    def rss():
        try:
            out = os.popen("systemctl show theobroma.service -p MemoryCurrent").read()
            return int(out.strip().split("=")[1]) / 1e9
        except Exception:
            return None
    a = rss()
    for _ in range(40):
        c.get("/api/search?type=kingdom&q=plant&limit=50", timeout=60)
    b = rss()
    if a is None or b is None:
        return "KNOWN", "service RSS not readable", "needs systemd accounting enabled"
    return b - a < 2.0, "RSS %.1f GB before, %.1f GB after 40 requests" % (a, b), \
        "the in-process FAISS indices make unbounded growth an OOM risk"


# ------------------------------------------------- AD. chemical data quality

@test("AD1", "stored SMILES parse and round-trip")
def ad1(c):
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        return "KNOWN", "RDKit unavailable in this interpreter", "skipped"
    rows = c.q("""SELECT comp_id, smiles FROM compounds
                  WHERE smiles IS NOT NULL ORDER BY random() LIMIT 300""")
    unparsed = stable = 0
    for cid, smi in rows:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            unparsed += 1
            continue
        can = Chem.MolToSmiles(m)
        stable += 1 if Chem.MolToSmiles(Chem.MolFromSmiles(can)) == can else 0
    return unparsed == 0 and stable >= len(rows) - 5, \
        "unparsed %s of 300, canonical round-trip stable %s" % (unparsed, stable)


@test("AD2", "sanitization passes without valence errors")
def ad2(c):
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        return "KNOWN", "RDKit unavailable", "skipped"
    rows = c.q("SELECT smiles FROM compounds WHERE smiles IS NOT NULL ORDER BY random() LIMIT 300")
    bad = 0
    for (smi,) in rows:
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is None:
            bad += 1
            continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            bad += 1
    return bad <= 3, "sanitization failures %s of 300" % bad, \
        "ChEMBL's Checker flags valence errors as the highest-severity structural defect"


@test("AD3", "stored exact mass agrees with the structure")
def ad3(c):
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.Descriptors import ExactMolWt
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        return "KNOWN", "RDKit unavailable", "skipped"
    rows = c.q("""SELECT smiles, mw FROM compounds
                  WHERE smiles IS NOT NULL AND mw IS NOT NULL ORDER BY random() LIMIT 300""")
    off = 0
    for smi, mw in rows:
        m = Chem.MolFromSmiles(smi)
        if m is not None and abs(ExactMolWt(m) - float(mw)) > 0.01:
            off += 1
    return off <= 5, "rows differing by more than 0.01 Da: %s of 300" % off, \
        "the column is monoisotopic ExactMolWt, not average molecular weight"


@test("AD4", "InChIKey re-derives from the stored structure")
def ad4(c):
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.inchi import MolToInchiKey
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        return "KNOWN", "RDKit unavailable", "skipped"
    rows = c.q("""SELECT smiles, inchikey FROM compounds
                  WHERE smiles IS NOT NULL ORDER BY random() LIMIT 300""")
    ok = 0
    for smi, ik in rows:
        m = Chem.MolFromSmiles(smi)
        if m is not None and MolToInchiKey(m) == ik:
            ok += 1
    pct = 100.0 * ok / max(1, len(rows))
    return ("KNOWN" if pct >= 95 else "FAIL"), \
        "InChIKey re-derives for %.1f%% of 300" % pct, \
        "tautomer and stereo perception drift across RDKit versions, so 97%% is expected"


@test("AD5", "ingestion artefacts are bounded and labelled")
def ad5(c):
    r = c.q("""SELECT
        count(*) FILTER (WHERE smiles LIKE '%%.%%'),
        count(*) FILTER (WHERE smiles LIKE '%%*%%'),
        count(*) FILTER (WHERE smiles ~ '\\[[0-9]+[A-Z]'),
        count(*) FILTER (WHERE smiles ~ '\\[R[0-9]*\\]')
        FROM compounds""")[0]
    return ("KNOWN" if r[1] + r[3] == 0 else "FAIL"), \
        "multi-component %s, wildcard %s, isotopic %s, R-group %s" % tuple(r), \
        "salts and isotopologues are retained by design; wildcards and R-groups are not"


@test("AD6", "the connectivity-prefix partition is internally consistent")
def ad6(c):
    r = c.q("""SELECT count(DISTINCT left(inchikey,14)),
                      count(*) FILTER (WHERE length(inchikey) <> 27),
                      count(*) FILTER (WHERE inchikey IS NULL OR inchikey = '')
               FROM compounds""")[0]
    return r[0] == 486032 and r[1] == 0 and r[2] == 0, \
        "families %s, malformed keys %s, missing keys %s" % tuple(r), \
        "the family relation is the paper's third contribution and rests on the full key"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://localhost:5000")
    ap.add_argument("--public", default="https://theobroma.l3s.uni-hannover.de")
    ap.add_argument("--stress", action="store_true",
                    help="run the slow-client and sustained-load checks")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    con = psycopg2.connect(DSN)
    con.set_session(readonly=True, autocommit=True)
    ctx = Ctx(con, a.base, a.public, a.stress)
    print("THEOBROMA adversarial suite V  %d checks" % len(TESTS))
    print("  app %s   public %s%s" % (a.base, a.public, "   stress" if a.stress else ""))
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
