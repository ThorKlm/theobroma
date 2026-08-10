"""FAISS/retrieval latency re-benchmark for THEOBROMA S11.

Times the deployed retrieval path through the running gunicorn service, so each
index is measured in its actual resident-memory state rather than a freshly
loaded second copy. This matters on the 15 GiB host (~286 MiB free): loading a
standalone FAISS copy would evict Postgres page cache or swap, and the numbers
would reflect paging, not retrieval. Hitting localhost:5000 avoids that.

Memory is guarded and reported: MemAvailable before/after and, decisively, swap
pages in/out across the run (from /proc/vmstat). If swap activity occurred the
run is flagged untrustworthy. This is the "watch memory headroom before trusting
the result" requirement made mechanical.

The published S11 figures (1,308 / 1,578 ms) were the two embedding metrics on a
~48%-coverage index; this re-benchmarks chemberta and nafm at full coverage over
the same 10,000-request workload, plus the textual comparator.

Usage:
  cd ~/theobroma && source venv/bin/activate
  python3 bench_latency_v135.py                 # 10000 reqs, all metrics + textual
  python3 bench_latency_v135.py --n 2000        # shorter smoke run
Reads nothing sensitive; localhost only; writes bench_latency_v135.json.
"""
import time, statistics, json, urllib.request, urllib.error, urllib.parse, argparse

BASE = "http://localhost:5000"

# (name, SMILES) pairs: SMILES feeds morgan/maccs/chemberta; name feeds nafm's q= resolver.
PROBES = [
    ("caffeine",    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ("curcumin",    "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O"),
    ("quinine",     "COc1ccc2nccc([C@@H](O)[C@H]3CC4CCN3C[C@@H]4C=C)c2c1"),
    ("morphine",    "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5"),
    ("resveratrol", "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1"),
    ("quercetin",   "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    ("berberine",   "COc1ccc2cc3[n+](cc2c1OC)CCc1cc2c(cc1-3)OCO2"),
    ("colchicine",  "COc1cc2c(c(OC)c1OC)-c1ccc(OC)c(=O)cc1[C@@H](NC(C)=O)CC2"),
    ("naringenin",  "O=C1C[C@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21"),
    ("kaempferol",  "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"),
]

def meminfo():
    d = {}
    for line in open("/proc/meminfo"):
        parts = line.split()
        d[parts[0].rstrip(":")] = int(parts[1]) // 1024   # kB -> MiB
    return d

def avail_mib():
    return meminfo().get("MemAvailable", 0)

def swap_counters():
    si = so = 0
    for line in open("/proc/vmstat"):
        if line.startswith("pswpin "):  si = int(line.split()[1])
        if line.startswith("pswpout "): so = int(line.split()[1])
    return si, so

def one_call(url):
    t0 = time.perf_counter()
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            r.read()
        return (time.perf_counter() - t0) * 1000.0, True
    except (urllib.error.URLError, TimeoutError, ConnectionError):
        return (time.perf_counter() - t0) * 1000.0, False

def bench(kind, url_fn, n, warm=20):
    lat, fails = [], 0
    for i in range(n + warm):
        _name, _smi = PROBES[i % len(PROBES)]
        ms, ok = one_call(url_fn(_name, _smi))
        if i < warm:
            continue
        if ok: lat.append(ms)
        else:  fails += 1
    if not lat:
        return {"kind": kind, "n": 0, "fails": fails, "mean": None,
                "median": None, "p95": None, "p99": None, "min": None, "max": None}
    lat.sort()
    pct = lambda p: lat[min(len(lat)-1, int(p/100*len(lat)))]
    return {"kind": kind, "n": len(lat), "fails": fails,
            "mean": round(statistics.mean(lat),1), "median": round(pct(50),1),
            "p95": round(pct(95),1), "p99": round(pct(99),1),
            "min": round(lat[0],1), "max": round(lat[-1],1)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10000)
    args = ap.parse_args()
    n = args.n

    avail0 = avail_mib()
    print(f"host memory before: MemAvailable={avail0} MiB")
    if avail0 < 500:
        print("WARNING: <500 MiB available. This script hits the RUNNING service so the")
        print("index is already resident, but if latencies look inflated, memory pressure")
        print("is the suspect — check the swap line at the end.")
    si0, so0 = swap_counters()
    t_start = time.time()

    results = []
    # embedding + fingerprint metrics that take a SMILES
    for m in ("morgan", "maccs", "chemberta"):
        results.append(bench(f"similarity:{m}",
            lambda name, smi, m=m: f"{BASE}/api/similarity?smiles={urllib.parse.quote(smi)}&metric={m}&limit=50", n))
        print(f"  done {m}: MemAvailable={avail_mib()} MiB")
    # nafm resolves by name via q=
    results.append(bench("similarity:nafm",
        lambda name, smi: f"{BASE}/api/similarity?q={urllib.parse.quote(name)}&metric=nafm&limit=50", n))
    print(f"  done nafm: MemAvailable={avail_mib()} MiB")
    # textual comparator (same workload, name search)
    results.append(bench("text:name",
        lambda name, smi: f"{BASE}/api/search?q={urllib.parse.quote(name)}&type=name&limit=50", n))

    si1, so1 = swap_counters()
    avail1 = avail_mib()
    elapsed = time.time() - t_start

    print("\n=== latency (ms) ===")
    for r in results:
        print(f"{r['kind']:20s} n={r['n']:6d} fails={r['fails']:4d} "
              f"mean={r['mean']} median={r['median']} p95={r['p95']} p99={r['p99']} "
              f"min={r['min']} max={r['max']}")
    print(f"\nwall={elapsed:.1f}s  MemAvailable after={avail1} MiB")
    pin, pout = si1 - si0, so1 - so0
    print(f"swap during run: pages_in={pin} pages_out={pout}")
    if pin > 0 or pout > 0:
        print("CAUTION: swap activity occurred — latencies may be inflated by paging.")
        print("Re-run when the host is quieter; do NOT quote these numbers as-is.")
    else:
        print("No swap activity — timings reflect in-memory retrieval, safe to quote.")

    json.dump({"results": results, "elapsed_s": round(elapsed,1),
               "mem_avail_before": avail0, "mem_avail_after": avail1,
               "swap_pages_in": pin, "swap_pages_out": pout, "n_requested": n},
              open("bench_latency_v135.json", "w"), indent=2)
    print("wrote bench_latency_v135.json")

if __name__ == "__main__":
    main()
