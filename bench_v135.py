#!/usr/bin/env python3
"""THEOBROMA v1.35 latency benchmark and self-similarity canary.

Replaces the v1.34 latency figures in the manuscript (Section 5, Supplementary S11),
which were measured against a FAISS index serving ~48 percent coverage and are
therefore not valid for the redeployed full-coverage index.

Two phases, both read-only:
  bench   sequential latency for textual search, compound page, and each similarity
          metric; reports n, p50, p95, p99, mean, and the node spec.
  canary  enlarged self-similarity check (default n=500) replacing the 5-compound
          C20 test, which cannot distinguish sampling variance from index misalignment.

Usage on the production node:
  cd ~/theobroma && source venv/bin/activate
  python3 bench_v135.py --phase all --out bench_v135.json
  python3 bench_v135.py --phase bench --n-text 10000 --n-sim 1000
  python3 bench_v135.py --phase canary --n-canary 500

Env: THEO_DB_DSN, THEO_BASE_URL. Runtime at defaults is roughly 20 min textual,
20 min per similarity metric at ~1.3 s, 12 min canary. Use --n-sim 200 for a fast pass.
"""
import argparse, json, os, platform, random, statistics, subprocess, sys, time
import psycopg2, requests

DSN = os.environ.get("THEO_DB_DSN", "postgresql://theobroma:theobroma@localhost:5432/theobroma")
BASE = os.environ.get("THEO_BASE_URL", "http://localhost:5000").rstrip("/")
METRICS = ["morgan", "maccs", "chemberta", "nafm"]

_conn = None
def db():
    global _conn
    if _conn is None:
        _conn = psycopg2.connect(DSN)
        _conn.set_session(readonly=True, autocommit=True)
    return _conn

def qall(sql, args=None):
    cur = db().cursor(); cur.execute(sql, args or ()); r = cur.fetchall(); cur.close()
    return r

def pct(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(round(p * (len(xs) - 1))))]

def summarize(label, lat, extra=None):
    """Latencies in seconds; report in ms to match the manuscript convention."""
    d = {"label": label, "n": len(lat),
         "p50_ms": round(1000 * pct(lat, 0.50), 1), "p95_ms": round(1000 * pct(lat, 0.95), 1),
         "p99_ms": round(1000 * pct(lat, 0.99), 1), "mean_ms": round(1000 * statistics.fmean(lat), 1),
         "min_ms": round(1000 * min(lat), 1), "max_ms": round(1000 * max(lat), 1)}
    if extra: d.update(extra)
    print(f"  {label:28s} n={d['n']:6d}  p50={d['p50_ms']:8.1f}ms  p95={d['p95_ms']:8.1f}ms  "
          f"p99={d['p99_ms']:8.1f}ms  mean={d['mean_ms']:8.1f}ms", flush=True)
    return d

def node_spec():
    """Captured so the manuscript can state the measurement environment exactly."""
    def sh(cmd):
        try: return subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.DEVNULL).strip()
        except Exception: return "unavailable"
    return {"cpu_model": sh("grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs"),
            "cores": os.cpu_count(),
            "mem_total": sh("grep MemTotal /proc/meminfo | awk '{print $2, $3}'"),
            "kernel": platform.release(),
            "postgres": qall("SHOW server_version")[0][0],
            "python": platform.python_version(),
            "shared_buffers": qall("SHOW shared_buffers")[0][0],
            "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

def warm(n=200):
    """Warm the PostgreSQL cache and the FAISS index before timing, matching the
    manuscript's 'warm PostgreSQL cache' condition."""
    print("warming caches", flush=True)
    ids = [r[0] for r in qall("SELECT comp_id FROM compounds ORDER BY random() LIMIT %s", (n,))]
    for cid in ids:
        try: requests.get(f"{BASE}/api/compound/{cid}", timeout=30)
        except Exception: pass
    for q in ("caffeine", "curcumin", "morphine", "fisetin", "quinine"):
        try: requests.get(f"{BASE}/api/search", params={"type": "name", "q": q, "limit": 25}, timeout=30)
        except Exception: pass

def probe_similarity_shape(smiles, metric):
    """The similarity endpoint's parameter names are not fixed across revisions; resolve
    once so the timed loop does not pay for failed shapes."""
    shapes = [{"smiles": smiles, "metric": metric, "top_n": 25},
              {"smiles": smiles, "metric": metric, "n": 25},
              {"q": smiles, "metric": metric, "top_n": 25},
              {"smiles": smiles, "top_n": 25}]
    for s in shapes:
        try:
            r = requests.get(f"{BASE}/api/similarity", params=s, timeout=120)
            if r.status_code == 200: return s
        except Exception: pass
    return None

def extract_hits(payload):
    if isinstance(payload, list): return payload
    for k in ("results", "hits", "neighbors", "matches"):
        v = payload.get(k)
        if isinstance(v, list): return v
    return []

def phase_bench(n_text, n_sim, n_page, seed):
    random.seed(seed)
    out = {"node": node_spec(), "runs": []}
    warm()
    print("\ntextual search", flush=True)
    names = [r[0] for r in qall(
        "SELECT name FROM compounds WHERE name IS NOT NULL AND name<>'' "
        "ORDER BY random() LIMIT %s", (n_text,))]
    lat, errs = [], 0
    for nm in names:
        t0 = time.perf_counter()
        try:
            r = requests.get(f"{BASE}/api/search", params={"type": "name", "q": nm, "limit": 25}, timeout=60)
            if r.status_code != 200: errs += 1
        except Exception:
            errs += 1; continue
        lat.append(time.perf_counter() - t0)
    out["runs"].append(summarize("api_search_textual", lat, {"errors": errs}))
    print("\ncompound page", flush=True)
    ids = [r[0] for r in qall("SELECT comp_id FROM compounds ORDER BY random() LIMIT %s", (n_page,))]
    lat, errs = [], 0
    for cid in ids:
        t0 = time.perf_counter()
        try:
            r = requests.get(f"{BASE}/compound/{cid}", timeout=60)
            if r.status_code != 200: errs += 1
        except Exception:
            errs += 1; continue
        lat.append(time.perf_counter() - t0)
    out["runs"].append(summarize("compound_page", lat, {"errors": errs}))
    print("\nstructural similarity", flush=True)
    smis = [r[0] for r in qall(
        "SELECT smiles FROM compounds WHERE smiles<>'' ORDER BY random() LIMIT %s", (n_sim,))]
    for metric in METRICS:
        shape = probe_similarity_shape(smis[0], metric)
        if shape is None:
            print(f"  {metric:28s} endpoint shape not resolved, skipped", flush=True)
            out["runs"].append({"label": f"similarity_{metric}", "n": 0, "skipped": "shape unresolved"})
            continue
        key = "smiles" if "smiles" in shape else "q"
        lat, errs, empty = [], 0, 0
        for s in smis:
            p = dict(shape); p[key] = s
            t0 = time.perf_counter()
            try:
                r = requests.get(f"{BASE}/api/similarity", params=p, timeout=300)
                if r.status_code != 200: errs += 1; continue
            except Exception:
                errs += 1; continue
            lat.append(time.perf_counter() - t0)
            try:
                if not extract_hits(r.json()): empty += 1
            except Exception: pass
        if lat:
            out["runs"].append(summarize(f"similarity_{metric}", lat,
                                         {"errors": errs, "empty_results": empty, "params": shape}))
    return out

def phase_canary(n_canary, seed):
    """A compound must retrieve itself as top hit. Any miss rate above chance indicates
    comp_id-to-vector misalignment, which is what the 48 percent index exhibited."""
    random.seed(seed)
    rows = qall("SELECT comp_id, smiles FROM compounds WHERE smiles<>'' "
                "ORDER BY random() LIMIT %s", (n_canary,))
    res = {}
    for metric in METRICS:
        shape = probe_similarity_shape(rows[0][1], metric)
        if shape is None:
            res[metric] = {"skipped": "shape unresolved"}; continue
        key = "smiles" if "smiles" in shape else "q"
        hit, top5, tried, misses = 0, 0, 0, []
        for cid, smi in rows:
            p = dict(shape); p[key] = smi; p[[k for k in ("top_n", "n") if k in p][0]] = 5
            try:
                r = requests.get(f"{BASE}/api/similarity", params=p, timeout=300)
                if r.status_code != 200: continue
                hits = extract_hits(r.json())
            except Exception:
                continue
            tried += 1
            ids = [h.get("comp_id") if isinstance(h, dict) else h for h in hits]
            if ids and ids[0] == cid: hit += 1
            elif cid in ids: top5 += 1; misses.append((cid, ids[0]))
            else: misses.append((cid, ids[0] if ids else None))
        rate = 100.0 * hit / tried if tried else 0.0
        res[metric] = {"n": tried, "self_top1": hit, "self_in_top5_only": top5,
                       "top1_rate_pct": round(rate, 2), "example_misses": misses[:10]}
        print(f"  canary {metric:12s} n={tried:4d}  top1={hit:4d} ({rate:6.2f}%)  "
              f"in_top5_only={top5}", flush=True)
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="all", choices=["all", "bench", "canary"])
    ap.add_argument("--n-text", type=int, default=10000)
    ap.add_argument("--n-sim", type=int, default=1000)
    ap.add_argument("--n-page", type=int, default=1000)
    ap.add_argument("--n-canary", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="bench_v135.json")
    a = ap.parse_args()
    print(f"THEOBROMA benchmark  base={BASE}  db={DSN.split('@')[-1]}", flush=True)
    t0 = time.time()
    result = {"config": vars(a)}
    if a.phase in ("all", "bench"):
        result["bench"] = phase_bench(a.n_text, a.n_sim, a.n_page, a.seed)
    if a.phase in ("all", "canary"):
        print("\nself-similarity canary", flush=True)
        result["canary"] = phase_canary(a.n_canary, a.seed)
    result["elapsed_s"] = round(time.time() - t0, 1)
    with open(a.out, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nwrote {a.out} in {result['elapsed_s']}s", flush=True)

if __name__ == "__main__":
    main()
