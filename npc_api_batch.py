#!/usr/bin/env python3
"""Batch NPClassifier via the official GNPS2 API with concurrency, retries,
and resumability. Correct-by-construction (uses canonical NPClassifier).
Resumable: re-run after interruption and it skips already-done comp_ids.

Usage:
  python npc_api_batch.py --in npc_input.csv --out npc_output.csv --workers 32
"""
import argparse, csv, os, sys, time, threading, queue
import requests

API = "https://npclassifier.gnps2.org/classify"
LOCK = threading.Lock()

def classify(smiles, session, retries=4):
    for attempt in range(retries):
        try:
            r = session.get(API, params={"smiles": smiles}, timeout=30)
            if r.status_code == 200:
                d = r.json()
                return (";".join(d.get("pathway_results") or []),
                        ";".join(d.get("superclass_results") or []),
                        ";".join(d.get("class_results") or []))
            if r.status_code in (429, 500, 502, 503):
                time.sleep(1.5 * (attempt + 1)); continue
            return ("", "", "")  # 4xx other than rate-limit: unclassifiable
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return ("__ERROR__", "", "")

def load_done(out_path):
    done = set()
    if os.path.exists(out_path):
        with open(out_path) as fh:
            r = csv.reader(fh)
            next(r, None)
            for row in r:
                if row: done.add(row[0])
    return done

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default="/workspace/npc_input.csv")
    ap.add_argument("--out", dest="out", default="/workspace/npc_output.csv")
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    # read input
    rows = []
    with open(a.inp) as fh:
        r = csv.reader(fh); next(r, None)
        for row in r:
            if len(row) >= 2 and row[1].strip():
                rows.append((row[0], row[1]))
    if a.limit: rows = rows[:a.limit]

    done = load_done(a.out)
    todo = [(c, s) for c, s in rows if c not in done]
    print(f"total={len(rows)} done={len(done)} todo={len(todo)} workers={a.workers}", flush=True)

    # append mode; write header only if new file
    new = not os.path.exists(a.out)
    out_fh = open(a.out, "a", newline="")
    writer = csv.writer(out_fh)
    if new:
        writer.writerow(["comp_id", "np_pathway", "np_superclass", "np_class"])
        out_fh.flush()

    q = queue.Queue()
    for item in todo: q.put(item)
    counter = {"n": 0}
    t0 = time.time()

    def worker():
        sess = requests.Session()
        while True:
            try: cid, smi = q.get_nowait()
            except queue.Empty: return
            pw, sc, cl = classify(smi, sess)
            with LOCK:
                writer.writerow([cid, pw, sc, cl])
                counter["n"] += 1
                if counter["n"] % 500 == 0:
                    out_fh.flush()
                    rate = counter["n"] / (time.time() - t0)
                    eta = (len(todo) - counter["n"]) / rate / 3600 if rate else 0
                    print(f"  {counter['n']}/{len(todo)}  {rate:.1f}/s  ETA {eta:.1f}h", flush=True)
            q.task_done()

    threads = [threading.Thread(target=worker, daemon=True) for _ in range(a.workers)]
    for t in threads: t.start()
    for t in threads: t.join()
    out_fh.flush(); out_fh.close()
    print(f"done: {counter['n']} classified in {(time.time()-t0)/3600:.2f}h -> {a.out}", flush=True)

if __name__ == "__main__":
    main()
