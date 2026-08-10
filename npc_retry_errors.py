#!/usr/bin/env python3
"""Retry only the __ERROR__ rows in npc_output.csv by re-querying the API,
and rewrite the file in place (atomic) with the recovered classifications.
Small set (~18k), so fast. Safe to re-run."""
import csv, os, sys, time, threading, queue
import requests

API = "https://npclassifier.gnps2.org/classify"
LOCK = threading.Lock()

def classify(smiles, session, retries=5):
    for attempt in range(retries):
        try:
            r = session.get(API, params={"smiles": smiles}, timeout=40)
            if r.status_code == 200:
                d = r.json()
                return (";".join(d.get("pathway_results") or []),
                        ";".join(d.get("superclass_results") or []),
                        ";".join(d.get("class_results") or []))
            if r.status_code in (429, 500, 502, 503):
                time.sleep(2 * (attempt + 1)); continue
            return ("", "", "")
        except Exception:
            time.sleep(2 * (attempt + 1))
    return ("__ERROR__", "", "")

def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else "/home/thorben.klamt/theobroma/npc_output.csv"
    in_path  = sys.argv[2] if len(sys.argv) > 2 else "/home/thorben.klamt/theobroma/npc_input.csv"
    workers  = int(sys.argv[3]) if len(sys.argv) > 3 else 24

    # smiles lookup
    smi = {}
    with open(in_path) as fh:
        r = csv.reader(fh); next(r, None)
        for row in r:
            if len(row) >= 2: smi[row[0]] = row[1]

    # load current results, find errors
    rows = {}
    with open(out_path) as fh:
        r = csv.reader(fh); next(r, None)
        for row in r:
            if row: rows[row[0]] = row[1:4] + [""]*(3-len(row[1:4]))
    errors = [cid for cid, v in rows.items() if v[0] == "__ERROR__"]
    print(f"error rows to retry: {len(errors)}", flush=True)
    if not errors:
        print("nothing to retry"); return

    q = queue.Queue()
    for cid in errors: q.put(cid)
    cnt = {"n": 0, "fixed": 0}
    t0 = time.time()

    def worker():
        sess = requests.Session()
        while True:
            try: cid = q.get_nowait()
            except queue.Empty: return
            pw, sc, cl = classify(smi.get(cid, ""), sess)
            with LOCK:
                rows[cid] = [pw, sc, cl]
                cnt["n"] += 1
                if pw not in ("__ERROR__", ""): cnt["fixed"] += 1
                if cnt["n"] % 500 == 0:
                    print(f"  {cnt['n']}/{len(errors)} fixed={cnt['fixed']} {cnt['n']/(time.time()-t0):.1f}/s", flush=True)
            q.task_done()

    ts = [threading.Thread(target=worker, daemon=True) for _ in range(workers)]
    for t in ts: t.start()
    for t in ts: t.join()

    # rewrite atomically
    tmp = out_path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["comp_id", "np_pathway", "np_superclass", "np_class"])
        for cid, v in rows.items():
            w.writerow([cid] + v)
    os.replace(tmp, out_path)
    still = sum(1 for v in rows.values() if v[0] == "__ERROR__")
    print(f"done: retried {cnt['n']}, newly fixed {cnt['fixed']}, still error {still} -> {out_path}", flush=True)

if __name__ == "__main__":
    main()
