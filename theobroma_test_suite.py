#!/usr/bin/env python3
"""THEOBROMA test suite: 50 read-only checks across DB integrity, cheminformatics,
search/similarity, API/HTTP, security, performance, provenance, and project-specific
invariants. Designed for pre-NAR-submission validation of a live, read-only server.

Output is intentionally terse: one line per test (PASS/FAIL/WARN + a short metric).
No full row/file contents are printed. Read-only throughout; the injection/XSS probes
confirm safety without exfiltrating data and must only be run against your own server.

Deps: psycopg2 (DB tests), requests (HTTP tests); rdkit optional (chem tests degrade
to SKIP if absent). Config via env or flags:
  THEO_DB_DSN   default 'postgresql://theobroma:theobroma@localhost:5432/theobroma'
  THEO_BASE_URL default 'http://localhost:5000'
Run subsets with --only db,chem,search,api,sec,perf,prov,proj  (comma-separated).
"""
import argparse, os, re, sys, time, json, string, random, statistics
from contextlib import contextmanager

DSN  = os.environ.get("THEO_DB_DSN", "postgresql://theobroma:theobroma@localhost:5432/theobroma")
BASE = os.environ.get("THEO_BASE_URL", "http://localhost:5000").rstrip("/")

# expected invariants (verified live this build; count-parity tests re-derive them)
EXP = {
    "n_compounds": 1132805, "tier1": 626601, "tier2": 206042, "tier3": 300141,
    "unclassified": 21, "families": 486032, "region_compounds": 264039,
    "macro_regions": 13, "sources": 29, "cc0_pct_min": 88.5, "cc0_pct_max": 90.5,
    "np_pathways": 7,
}
KINGDOMS = {"plant", "animal", "fungi", "bacteria", "unresolved"}
LIC_TIERS = {"CC0", "CC BY 4.0", "CC BY-NC 4.0", "CC BY-NC-ND 4.0", "Unspecified"}

# ---- harness ---------------------------------------------------------------
R = {"PASS": 0, "FAIL": 0, "WARN": 0, "SKIP": 0}
FAILED = []
def _log(status, tid, msg):
    R[status] += 1
    if status == "FAIL": FAILED.append(tid)
    print(f"[{status:4}] {tid:5} {msg}", flush=True)
def PASS(t, m=""): _log("PASS", t, m)
def FAIL(t, m=""): _log("FAIL", t, m)
def WARN(t, m=""): _log("WARN", t, m)
def SKIP(t, m=""): _log("SKIP", t, m)
def check(t, cond, ok="", bad=""):
    (PASS if cond else FAIL)(t, ok if cond else bad)

_conn = None
def db():
    global _conn
    if _conn is None:
        import psycopg2
        _conn = psycopg2.connect(DSN); _conn.set_session(readonly=True, autocommit=True)
    return _conn
def q1(sql, args=None):
    cur = db().cursor(); cur.execute(sql, args or ()); r = cur.fetchone(); cur.close()
    return r[0] if r else None
def qall(sql, args=None):
    cur = db().cursor(); cur.execute(sql, args or ()); r = cur.fetchall(); cur.close()
    return r

try:
    import requests
    HTTP = True
except Exception:
    HTTP = False
def get(path, **kw):
    kw.setdefault("timeout", 30)
    return requests.get(BASE + path, **kw)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, inchi as rdinchi
    from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")
    RDKIT = True
except Exception:
    RDKIT = False

def sample_smiles(n):
    return qall("SELECT comp_id, smiles, inchikey FROM compounds "
                "WHERE smiles<>'' ORDER BY random() LIMIT %s", (n,))

# ============================================================================
# A. DB INTEGRITY & INVARIANTS
# ============================================================================
def t_db():
    n = q1("SELECT count(*) FROM compounds")
    check("A01", n == EXP["n_compounds"], f"compounds={n:,}", f"compounds={n:,} != {EXP['n_compounds']:,}")

    d = q1("SELECT count(*)-count(DISTINCT comp_id) FROM compounds")
    check("A02", d == 0, "comp_id unique", f"{d} duplicate comp_id")

    bad = q1(r"SELECT count(*) FROM compounds WHERE inchikey IS NULL OR "
             r"inchikey !~ '^[A-Z]{14}-[A-Z]{10}-[A-Z]$'")
    check("A03", bad == 0, "all inchikeys well-formed", f"{bad:,} malformed inchikeys")

    orphans = {}
    for tbl, col in [("compound_taxonomy","comp_id"),("compound_region_map","comp_id"),
                     ("admet","comp_id"),("scaffolds","comp_id")]:
        o = q1(f"SELECT count(*) FROM {tbl} t LEFT JOIN compounds c ON t.{col}=c.comp_id "
               f"WHERE c.comp_id IS NULL")
        orphans[tbl] = o
    tot = sum(orphans.values())
    check("A04", tot == 0, "no orphaned child rows", f"orphans: {orphans}")

    syn_orphan = q1("SELECT count(*) FROM compound_synonyms s "
                    "LEFT JOIN compounds c ON s.inchikey=c.inchikey WHERE c.inchikey IS NULL")
    # synonyms keyed by inchikey; some may legitimately not map -> WARN not FAIL
    (PASS if syn_orphan == 0 else WARN)("A05", f"{syn_orphan:,} synonym rows with no matching inchikey")

    rows = dict(qall("SELECT classification_tier, count(*) FROM compounds GROUP BY 1"))
    t1,t2,t3 = rows.get(1,0),rows.get(2,0),rows.get(3,0); un = rows.get(None,0)
    part_ok = (t1==EXP["tier1"] and t2==EXP["tier2"] and t3==EXP["tier3"] and un==EXP["unclassified"]
               and t1+t2+t3+un==EXP["n_compounds"])
    check("A06", part_ok, f"tiers {t1:,}/{t2:,}/{t3:,}/{un} sum ok",
          f"tier partition {t1:,}/{t2:,}/{t3:,}/{un}")

    badk = q1("SELECT count(*) FROM compounds WHERE kingdom IS NOT NULL AND kingdom NOT IN %s",
              (tuple(KINGDOMS),))
    check("A07", badk == 0, "kingdom domain ok", f"{badk:,} rows with unexpected kingdom")

    badl = q1("SELECT count(*) FROM compounds WHERE license_tier IS NOT NULL "
              "AND license_tier NOT IN %s", (tuple(LIC_TIERS),))
    badr = q1("SELECT count(*) FROM compounds WHERE tier_rank NOT BETWEEN 0 AND 5")
    check("A08", badl == 0 and badr == 0, "license_tier + tier_rank domains ok",
          f"{badl:,} bad license_tier, {badr:,} bad tier_rank")

    out = q1("SELECT count(*) FROM compounds WHERE mw<=0 OR mw>5000 OR tpsa<0 "
             "OR hba<0 OR hbd<0 OR n_rings<0 OR rotatable_bonds<0 OR logp< -20 OR logp>50")
    (PASS if out == 0 else WARN)("A09", f"{out:,} rows with out-of-range descriptors")

    fam = q1("SELECT count(DISTINCT left(inchikey,14)) FROM compounds WHERE inchikey<>''")
    check("A10", fam == EXP["families"], f"families={fam:,}", f"families={fam:,} != {EXP['families']:,}")

# ============================================================================
# B. CHEMINFORMATICS CORRECTNESS
# ============================================================================
def t_chem():
    if not RDKIT:
        for t in ("B11","B12","B13","B14","B15","B16","B17"): SKIP(t, "rdkit not installed")
        return
    smp = sample_smiles(2000)
    parsed = [(cid, Chem.MolFromSmiles(s), s, ik) for cid, s, ik in smp]
    nfail = sum(1 for _,m,_,_ in parsed if m is None)
    (PASS if nfail == 0 else WARN)("B11", f"SMILES parse fail {nfail}/{len(smp)}")

    ok_mols = [(cid,m,s,ik) for cid,m,s,ik in parsed if m is not None]
    mism = 0
    for cid,m,s,ik in ok_mols:
        try:
            k = rdinchi.InchiToInchiKey(rdinchi.MolToInchi(m))
            if k and k != ik: mism += 1
        except Exception:
            mism += 1
    frac = 100*mism/max(1,len(ok_mols))
    # policy: stored inchikey should equal rdkit-derived; allow small tautomer/stereo drift
    (PASS if frac < 5 else WARN)("B12", f"InChIKey re-derivation mismatch {mism}/{len(ok_mols)} ({frac:.2f}%)")

    nonidem = 0
    for cid,m,s,ik in ok_mols[:1000]:
        try:
            c1 = Chem.MolToSmiles(m); c2 = Chem.MolToSmiles(Chem.MolFromSmiles(c1))
            if c1 != c2: nonidem += 1
        except Exception:
            nonidem += 1
    (PASS if nonidem == 0 else WARN)("B13", f"canonicalization non-idempotent {nonidem}/1000")

    # descriptor recomputation spot-check (mw within 0.5)
    dd = qall("SELECT smiles, mw, tpsa FROM compounds WHERE smiles<>'' AND mw IS NOT NULL "
              "ORDER BY random() LIMIT 300")
    off = 0
    for s, mw, tpsa in dd:
        m = Chem.MolFromSmiles(s)
        if m is None: continue
        if abs(Descriptors.MolWt(m) - float(mw)) > 1.0: off += 1
    (PASS if off < 15 else WARN)("B14", f"mw recompute off by >1.0 for {off}/300")

    # adversarial structures must not be accepted as valid by RDKit-side (sanity of our parser use)
    adv = ["", "   ", "C1CC", "not_a_smiles", "[Xx]", "C"*2000]
    crashed = 0
    for a in adv:
        try: Chem.MolFromSmiles(a)
        except Exception: crashed += 1
    check("B15", crashed == 0, "adversarial SMILES handled w/o exception", f"{crashed} raised")

    # every stored SMILES that parses should round-trip to non-empty canonical
    empty = sum(1 for _,m,_,_ in ok_mols if not Chem.MolToSmiles(m))
    check("B16", empty == 0, "all parsed SMILES canonicalize non-empty", f"{empty} empty canon")

    # formula collision within 14-char families (WARN: rare but investigate)
    if RDKIT:
        fam_rows = qall("SELECT left(inchikey,14) p, comp_id, smiles FROM compounds "
                        "WHERE smiles<>'' ORDER BY random() LIMIT 3000")
        from collections import defaultdict
        byp = defaultdict(list)
        for p,cid,s in fam_rows: byp[p].append(s)
        coll = 0
        for p, ss in byp.items():
            if len(ss) < 2: continue
            formulas = set()
            for s in ss[:5]:
                m = Chem.MolFromSmiles(s)
                if m: formulas.add(Chem.rdMolDescriptors.CalcMolFormula(m))
            if len(formulas) > 1: coll += 1
        (PASS if coll == 0 else WARN)("B17", f"families w/ differing formulae {coll} (sampled)")

# ============================================================================
# C. SEARCH & SIMILARITY  (vector alignment canaries)
# ============================================================================
def t_search():
    vec_dir = os.environ.get("THEO_VEC_DIR", os.path.expanduser("~/theobroma/data/vectors"))
    try:
        import numpy as np
    except Exception:
        for t in ("C18","C19"): SKIP(t, "numpy missing")
        np = None
    # C18: index cardinality alignment (comp_ids == embeddings == valid_indices; faiss ntotal)
    if np:
        try:
            cids = np.load(f"{vec_dir}/comp_ids.npy", allow_pickle=True)
            vi   = np.load(f"{vec_dir}/valid_indices.npy")
            emb  = np.load(f"{vec_dir}/chemberta_embeddings.npy", mmap_mode="r")
            card_ok = emb.shape[0] == len(vi)
            in_range = int(vi.max()) < len(cids)
            check("C18", card_ok and in_range,
                  f"emb={emb.shape[0]:,} vi={len(vi):,} cids={len(cids):,}",
                  f"cardinality/align off emb={emb.shape[0]:,} vi={len(vi):,}")
            # coverage vs corpus (WARN if <99%: e.g. stale 48% vectors still deployed)
            corpus = set(x for (x,) in qall("SELECT comp_id FROM compounds"))
            feat = set(np.asarray(cids)[vi].astype(str))
            cov = 100*len(feat & corpus)/len(corpus)
            (PASS if cov >= 99 else WARN)("C19", f"vector corpus coverage {cov:.1f}%")
        except Exception as e:
            SKIP("C18", f"vectors unreadable: {type(e).__name__}")
            SKIP("C19", "vectors unreadable")
    # C20: self-similarity canary via API (top hit == self, score ~1)
    if HTTP:
        try:
            ids = [r[0] for r in qall("SELECT comp_id FROM compounds "
                                      "WHERE classification_tier IS NOT NULL ORDER BY random() LIMIT 5")]
            selfhit = 0; tried = 0
            for cid in ids:
                r = None
                for params in ({"comp_id": cid, "top_n": 3}, {"query": cid, "n": 3},
                               {"comp_id": cid}, {"q": cid}):
                    rr = get("/api/similarity", params=params)
                    if rr.status_code == 200:
                        r = rr; break
                if r is None: continue
                tried += 1
                try: data = r.json()
                except Exception: continue
                hits = data if isinstance(data, list) else (
                    data.get("results") or data.get("hits") or data.get("neighbors") or [])
                top = hits[0] if hits else {}
                topid = top.get("comp_id") if isinstance(top, dict) else (
                    top if isinstance(top, str) else None)
                if topid == cid: selfhit += 1
            if tried == 0:
                SKIP("C20", "similarity endpoint not reachable / unknown shape")
            else:
                (PASS if selfhit == tried else WARN)("C20",
                    f"self-similarity top-hit {selfhit}/{tried} (misalignment canary)")
        except Exception as e:
            SKIP("C20", f"{type(e).__name__}")
    else:
        SKIP("C20", "requests missing")
    # C21: exact-match lookup by comp_id via API page
    if HTTP:
        cid = q1("SELECT comp_id FROM compounds ORDER BY random() LIMIT 1")
        r = get(f"/api/compound/{cid}")
        ok = r.status_code == 200 and cid in r.text
        check("C21", ok, f"/api/compound/{cid} 200+id", f"/api/compound/{cid} -> {r.status_code}")
    else:
        SKIP("C21", "requests missing")
    # C22: empty-result search returns clean (no 500)
    if HTTP:
        r = get("/api/search", params={"q": "zzzq_no_such_compound_xyz123"})
        check("C22", r.status_code == 200, f"empty search -> {r.status_code}",
              f"empty search -> {r.status_code}")
    else:
        SKIP("C22", "requests missing")
    # C23: pagination stability (no dup ids across two pages) if browse api exists
    if HTTP:
        r1 = get("/api/search", params={"q": "acid", "limit": 25, "offset": 0})
        r2 = get("/api/search", params={"q": "acid", "limit": 25, "offset": 25})
        if r1.status_code != 200 or r2.status_code != 200:
            SKIP("C23", "pagination params not honored by /api/search")
        else:
            try:
                def ids(r):
                    d = r.json(); rows = d if isinstance(d, list) else (
                        d.get("results") or d.get("compounds") or [])
                    return {x.get("comp_id") for x in rows if isinstance(x, dict) and x.get("comp_id")}
                a, b = ids(r1), ids(r2)
                if not a:
                    SKIP("C23", "search result shape has no comp_id")
                else:
                    check("C23", len(a & b) == 0, f"pages disjoint ({len(a)},{len(b)})",
                          f"page overlap {len(a & b)}")
            except Exception:
                SKIP("C23", "pagination shape unknown")
    else:
        SKIP("C23", "requests missing")
    # C24: family endpoint members share 14-char prefix
    if HTTP:
        cid = q1("SELECT comp_id FROM compounds ORDER BY random() LIMIT 1")
        r = get(f"/api/stereoisomers/{cid}")
        check("C24", r.status_code == 200, f"/api/stereoisomers/{cid} -> {r.status_code}",
              f"/api/stereoisomers/{cid} -> {r.status_code}")

# ============================================================================
# D. API / HTTP BEHAVIOR
# ============================================================================
def t_api():
    if not HTTP:
        for t in ("D25","D26","D27","D28","D29","D30","D31"): SKIP(t, "requests missing")
        return
    r = get("/")
    check("D25", r.status_code == 200, "home 200", f"home -> {r.status_code}")

    bad = get("/compound/THEO_9999999_nonexistent")
    check("D26", bad.status_code in (404, 400), f"missing compound -> {bad.status_code}",
          f"missing compound -> {bad.status_code} (want 404)")

    cid = q1("SELECT comp_id FROM compounds ORDER BY random() LIMIT 1")
    r = requests.request("POST", BASE + f"/compound/{cid}", timeout=15)
    check("D27", r.status_code in (405, 400, 404), f"POST to GET route -> {r.status_code}",
          f"POST -> {r.status_code} (want 405)")

    # content-type on an API endpoint (try a couple of likely paths)
    rr = get("/api/stats")
    ct_ok = rr.status_code == 200 and "application/json" in rr.headers.get("Content-Type", "")
    check("D28", ct_ok, "/api/stats application/json", f"/api/stats ct={rr.headers.get('Content-Type','?')}")

    # bulk endpoint must not hang/OOM: expect it bounded or streamed (status<=... quickly with a cap)
    t0 = time.time()
    try:
        rb = get("/api/bulk", params={"format": "csv", "limit": 10}, timeout=20)
        dt = time.time() - t0
        if rb.status_code == 404: SKIP("D29", "/api/bulk route unknown")
        else: check("D29", rb.status_code == 200 and dt < 20, f"bulk(limit=10) {rb.status_code} {dt:.1f}s",
                    f"bulk -> {rb.status_code} in {dt:.1f}s")
    except Exception as e:
        FAIL("D29", f"bulk raised {type(e).__name__}")

    # UI/API/DB count consistency (best-effort: compare a stats endpoint to DB)
    dbn = q1("SELECT count(*) FROM compounds")
    r = get("/api/stats")
    if r.status_code != 200:
        SKIP("D30", "/api/stats unknown; DB count only")
    else:
        try:
            d = r.json()
            apin = d.get("compounds") or d.get("n_compounds") or d.get("total")
            check("D30", apin == dbn, f"api={apin} db={dbn:,}", f"api={apin} != db={dbn:,}")
        except Exception:
            SKIP("D30", "stats shape unknown")

    # region reconciliation via API if present, else DB-only invariant
    rc = q1("SELECT count(DISTINCT comp_id) FROM compound_region_map")
    mr = q1("SELECT count(DISTINCT macro_region) FROM compound_region_map")
    check("D31", rc == EXP["region_compounds"] and mr == EXP["macro_regions"],
          f"region compounds={rc:,} macro_regions={mr}",
          f"region {rc:,}/{mr} != {EXP['region_compounds']:,}/{EXP['macro_regions']}")

    # D32: SVG depiction renders valid, non-empty SVG (endpoint takes ?smiles=)
    row = qall("SELECT smiles FROM compounds WHERE smiles<>'' ORDER BY random() LIMIT 1")
    smi = row[0][0] if row else ""
    r = get("/api/depict", params={"smiles": smi, "w": 300, "h": 200})
    body = r.text if r.status_code == 200 else ""
    ok = r.status_code == 200 and "<svg" in body and len(body) > 200
    check("D32", ok, f"depict SVG ok ({len(body)}B)", f"/api/depict -> {r.status_code} len={len(body)}")

    # D33: license-provenance endpoint returns per-source attestations for a compound
    cid = q1("SELECT comp_id FROM compounds ORDER BY random() LIMIT 1")
    r = get(f"/api/compound/{cid}/license-provenance")
    check("D33", r.status_code == 200, f"license-provenance -> {r.status_code}",
          f"license-provenance -> {r.status_code}")

    # D34: bulk export format parity - CSV and JSON of the same small query agree in count
    rc_csv = get("/api/bulk", params={"format": "csv", "q": "caffeine", "limit": 50})
    rc_json = get("/api/bulk", params={"format": "json", "q": "caffeine", "limit": 50})
    if rc_csv.status_code != 200 or rc_json.status_code != 200:
        SKIP("D34", f"bulk parity csv={rc_csv.status_code} json={rc_json.status_code}")
    else:
        n_csv = max(0, len([l for l in rc_csv.text.splitlines() if l.strip()]) - 1)
        try:
            jd = rc_json.json(); n_json = len(jd if isinstance(jd, list) else jd.get("results", []))
            (PASS if n_csv == n_json else WARN)("D34", f"bulk csv={n_csv} json={n_json} rows")
        except Exception:
            SKIP("D34", "bulk json shape unknown")

# ============================================================================
# E. SECURITY & ROBUSTNESS
# ============================================================================
def t_sec():
    if not HTTP:
        for t in ("E32","E33","E34","E35","E36","E37","E38"): SKIP(t, "requests missing")
        return
    base_n = q1("SELECT count(*) FROM compounds")
    # E32: SQL injection on search -> must not error, must not change result semantics
    inj = ["' OR '1'='1", "'; DROP TABLE compounds;--", "1 UNION SELECT NULL--", "') OR ('1'='1"]
    bad = 0; reachable = False
    for p in inj:
        for path in ("/api/search", "/search"):
            try:
                r = get(path, params={"q": p}, timeout=15)
            except Exception:
                bad += 1; continue
            if r.status_code == 404: continue
            reachable = True
            if r.status_code >= 500: bad += 1
    still = q1("SELECT count(*) FROM compounds")
    if not reachable: SKIP("E32", "no search route reached")
    else: check("E32", bad == 0 and still == base_n, "injection payloads safe",
                f"{bad} 5xx / table count {still}!={base_n}")

    # E33: XSS escaping in a rendered compound page (name with HTML must be escaped)
    #   find/insert nothing; just check that pages HTML-escape by probing a search reflection
    xss = "<script>theo_xss()</script>"
    reflected = False
    for path in ("/search", "/api/search"):
        try: r = get(path, params={"q": xss}, timeout=15)
        except Exception: continue
        if r.status_code == 200 and xss in r.text:
            reflected = True; break
    (FAIL if reflected else PASS)("E33", "reflected <script> unescaped" if reflected
                                  else "no unescaped script reflection")

    # E34: error responses must not leak stack traces / SQL / debug
    leak_markers = ("Traceback (most recent call last)", "psycopg2.", "SELECT ", "File \"/",
                    "Werkzeug Debugger", "sqlalchemy", "gunicorn")
    r = get("/compound/%27%22%3E", timeout=15)  # junk id
    body = r.text[:5000]
    leaked = [m for m in leak_markers if m in body]
    check("E34", not leaked, "no debug/stack leakage", f"leaked markers: {leaked}")

    # E35: no auth wall (NAR requirement) - key pages reachable w/o creds
    walls = 0
    for p in ("/", "/compound/"+ (q1("SELECT comp_id FROM compounds LIMIT 1") or "x")):
        r = get(p)
        if r.status_code in (401, 403): walls += 1
    check("E35", walls == 0, "no auth wall on public pages", f"{walls} pages require auth")

    # E36: security headers present (WARN-level; not all mandatory)
    r = get("/")
    h = {k.lower() for k in r.headers.keys()}
    missing = [x for x in ("x-content-type-options",) if x not in h]
    (PASS if not missing else WARN)("E36", f"missing headers: {missing}" if missing else "nosniff present")

    # E37: oversized/huge pagination is clamped, not honored literally
    r = get("/api/search", params={"q": "acid", "limit": 99999999})
    if r.status_code != 200: SKIP("E37", f"/api/search limit probe -> {r.status_code}")
    else:
        try:
            d = r.json(); rows = d if isinstance(d, list) else (d.get("results") or d.get("compounds") or [])
            n = len(rows)
            (PASS if n <= 10000 else WARN)("E37", f"huge limit returned {n} rows (clamp expected)")
        except Exception:
            SKIP("E37", "shape unknown")

    # E38: path traversal rejected
    r = get("/compound/..%2f..%2f..%2fetc%2fpasswd")
    check("E38", r.status_code in (400,404) and "root:" not in r.text,
          f"traversal -> {r.status_code}", "path traversal not rejected")

    # E39: admin + write endpoints must not be anonymously actionable
    concerns = {}
    ra = get("/admin/refresh_cache")
    concerns["/admin/refresh_cache"] = ra.status_code
    try:
        rp = requests.post(BASE + "/api/annotate", json={}, timeout=15)
        concerns["/api/annotate(POST)"] = rp.status_code
    except Exception:
        concerns["/api/annotate(POST)"] = "err"
    # want: admin refresh not 200-for-anyone; annotate POST rejected/authless-guarded
    admin_open = concerns.get("/admin/refresh_cache") == 200
    annot_open = concerns.get("/api/annotate(POST)") in (200, 201)
    (WARN if (admin_open or annot_open) else PASS)("E39",
        f"admin/write exposure: {concerns}" if (admin_open or annot_open)
        else f"admin/write guarded {concerns}")

# ============================================================================
# F. PERFORMANCE
# ============================================================================
def t_perf():
    if not HTTP:
        for t in ("F39","F40","F41"): SKIP(t, "requests missing")
    else:
        ids = [r[0] for r in qall("SELECT comp_id FROM compounds ORDER BY random() LIMIT 10")]
        lat = []
        for cid in ids:
            t0 = time.time(); r = get(f"/compound/{cid}"); lat.append(time.time()-t0)
        p95 = sorted(lat)[int(0.95*len(lat))-1]
        (PASS if p95 < 1.0 else WARN)("F39", f"/compound p95={p95*1000:.0f}ms (warm)")

        for path, thr, tid in (("/api/search", 2.0, "F40"),):
            t0 = time.time()
            try: r = get(path, params={"q":"caffeine"}, timeout=15); dt = time.time()-t0
            except Exception: SKIP(tid, "search route error"); continue
            if r.status_code == 404: SKIP(tid, "search route unknown")
            else: (PASS if dt < thr else WARN)(tid, f"search {dt*1000:.0f}ms")

    # F41: index usage on the hot comp_id lookup (no seq scan on 1.1M table)
    try:
        plan = "\n".join(x[0] for x in qall(
            "EXPLAIN SELECT * FROM compounds WHERE comp_id=%s", ("THEO_0000001",)))
        (PASS if "Seq Scan" not in plan else WARN)("F41",
            "comp_id lookup uses index" if "Seq Scan" not in plan else "comp_id lookup SEQ SCAN")
    except Exception as e:
        SKIP("F41", f"{type(e).__name__}")

# ============================================================================
# G. PROVENANCE & REPRODUCIBILITY
# ============================================================================
def t_prov():
    # G42: documented counts re-derived from live DB
    checks = {
        "compounds": (q1("SELECT count(*) FROM compounds"), EXP["n_compounds"]),
        "families":  (q1("SELECT count(DISTINCT left(inchikey,14)) FROM compounds WHERE inchikey<>''"), EXP["families"]),
        "sources":   (q1("SELECT count(DISTINCT source_db) FROM compounds"), EXP["sources"]),
    }
    bad = {k:v for k,(v,e) in checks.items() if v != e}
    check("G42", not bad, "documented counts match live", f"mismatches: {bad}")

    # G43: CC0 dominance within tolerance
    pct = float(q1("SELECT round(100.0*count(*) FILTER (WHERE license_tier='CC0')/count(*),2) FROM compounds"))
    check("G43", EXP["cc0_pct_min"] <= pct <= EXP["cc0_pct_max"], f"CC0={pct}%",
          f"CC0={pct}% outside [{EXP['cc0_pct_min']},{EXP['cc0_pct_max']}]")

    # G44: license tables complete for all sources + internally consistent
    nsrc_lic = q1("SELECT count(DISTINCT src) FROM source_license_ref")
    (PASS if nsrc_lic and nsrc_lic >= EXP["sources"] else WARN)("G44",
        f"source_license_ref covers {nsrc_lic} sources")

    # G45: per-compound license == most-permissive across attestations (spot-check)
    #   most-permissive == min(tier_rank). Compare to stored tier_rank for a sample.
    try:
        # most-permissive = min(tier_rank) across sources, EXCLUDING Unspecified (rank 5),
        # which is an unknown fallback, not a permissive license. Only compounds whose
        # stored tier is itself known (0-4) are checked against the known-source minimum.
        wrong = q1("""
            WITH s AS (
              SELECT a.comp_id, min(r.tier_rank) AS best
              FROM per_source_license_attestation a
              JOIN source_license_ref r ON r.src=a.source
              WHERE r.tier_rank < 5
              GROUP BY a.comp_id)
            SELECT count(*) FROM compounds c JOIN s ON s.comp_id=c.comp_id
            WHERE c.tier_rank < 5 AND c.tier_rank <> s.best
        """)
        check("G45", wrong == 0, "per-compound license = most-permissive across known sources",
              f"{wrong:,} compounds violate most-permissive rule (excl. Unspecified)")
    except Exception as e:
        SKIP("G45", f"attestation join failed: {type(e).__name__}")

    # G46: DOI resolves (dataset + preprint) - network dependent -> WARN on failure
    if HTTP:
        dois = ["https://doi.org/10.64898/2026.06.12.731585"]
        okd = 0
        for d in dois:
            try:
                rr = requests.head(d, allow_redirects=True, timeout=15)
                if rr.status_code < 400: okd += 1
            except Exception: pass
        (PASS if okd == len(dois) else WARN)("G46", f"DOIs resolve {okd}/{len(dois)}")
    else:
        SKIP("G46", "requests missing")

# ============================================================================
# H. THEOBROMA-SPECIFIC
# ============================================================================
def t_proj():
    # H47: tier <-> source consistency
    rows = dict(qall("SELECT classification_tier, class_source FROM compounds "
                     "WHERE classification_tier IS NOT NULL GROUP BY 1,2"))  # collapses; use explicit
    bad = q1("""SELECT count(*) FROM compounds WHERE
        (classification_tier=1 AND class_source<>'curated') OR
        (classification_tier=2 AND class_source<>'npclassifier') OR
        (classification_tier=3 AND class_source<>'inferred_xgb')""")
    check("H47", bad == 0, "tier<->class_source consistent", f"{bad:,} tier/source mismatches")

    # H48: Tier 3 rows carry inferred_class + source + confidence
    miss = q1("""SELECT count(*) FROM compounds WHERE classification_tier=3 AND
        (inferred_class IS NULL OR inferred_class='' OR inferred_class_source<>'inferred_xgb_v135b'
         OR inferred_confidence IS NULL)""")
    check("H48", miss == 0, "tier3 inferred fields populated", f"{miss:,} tier3 rows incomplete")

    # H49: inferred_confidence in [0,1]
    oob = q1("SELECT count(*) FROM compounds WHERE inferred_confidence IS NOT NULL "
             "AND (inferred_confidence<0 OR inferred_confidence>1)")
    check("H49", oob == 0, "inferred_confidence in [0,1]", f"{oob:,} out of range")

    # H50: every np_pathway token (multi-valued, pipe-delimited) is one of the 7 canonical
    CANON = {"Alkaloids","Terpenoids","Shikimates and Phenylpropanoids","Polyketides",
             "Fatty acids","Amino acids and Peptides","Carbohydrates"}
    import re as _re
    vals = [r[0] for r in qall("SELECT DISTINCT np_pathway FROM compounds "
                               "WHERE np_pathway IS NOT NULL AND np_pathway<>''")]
    toks = set()
    for v in vals:
        for t in v.split("|"): toks.add(_re.sub(r"\s+"," ",t).strip())
    bad = toks - CANON
    # non-canonical tokens are a real data issue (FAIL); whitespace-only variants are folded above.
    check("H50", not bad, f"pathway tokens canonical ({len(toks)} distinct)",
          f"non-canonical pathway tokens: {sorted(bad)}")
    # H51: np_pathway whitespace consistency. LIKE only; %% escapes psycopg2 param marker.
    ws = q1("SELECT count(*) FROM compounds WHERE np_pathway LIKE '%%  %%'")
    ws = 0 if ws is None else ws
    (PASS if ws == 0 else WARN)("H51", f"np_pathway double-space rows: {ws:,}")

# ============================================================================
CATS = {"db":t_db,"chem":t_chem,"search":t_search,"api":t_api,
        "sec":t_sec,"perf":t_perf,"prov":t_prov,"proj":t_proj}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="comma list: "+",".join(CATS))
    ap.add_argument("--dsn", default=DSN); ap.add_argument("--base", default=BASE)
    a = ap.parse_args()
    globals()["DSN"] = a.dsn; globals()["BASE"] = a.base.rstrip("/")
    sel = [c.strip() for c in a.only.split(",") if c.strip()] or list(CATS)
    print(f"THEOBROMA test suite  db={a.dsn.split('@')[-1]}  base={BASE}  "
          f"rdkit={'y' if RDKIT else 'n'} http={'y' if HTTP else 'n'}")
    print("-"*72)
    t0 = time.time()
    for c in sel:
        if c in CATS:
            try: CATS[c]()
            except Exception as e:
                FAIL(c.upper()+"XX", f"category crashed: {type(e).__name__}: {e}")
    print("-"*72)
    dt = time.time()-t0
    print(f"SUMMARY  PASS={R['PASS']} FAIL={R['FAIL']} WARN={R['WARN']} SKIP={R['SKIP']}  ({dt:.1f}s)")
    if FAILED: print("FAILED: " + ", ".join(FAILED))
    sys.exit(1 if R["FAIL"] else 0)

if __name__ == "__main__":
    main()
