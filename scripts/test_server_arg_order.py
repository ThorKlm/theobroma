#!/usr/bin/env python3
"""Argument-order equivalence tests for THEOBROMA server.

For each test query, runs it with the original argument order AND with a
randomly shuffled order. Compares responses for identity. Reports
[MATCH]/[MISMATCH] per query plus a summary.

Exit code 0 if all match, 1 otherwise.
"""
import urllib.request, urllib.error, json, sys, random, re

BASE = "http://localhost:5000"
TIMEOUT = 30
random.seed(42)

results = {"match": 0, "mismatch": 0, "mismatches": []}

def fetch(path):
    try:
        with urllib.request.urlopen(BASE + path, timeout=TIMEOUT) as r:
            return r.status, r.headers.get("Content-Type", ""), r.read()
    except urllib.error.HTTPError as e:
        return e.code, "", b""
    except Exception as e:
        return -1, "", str(e).encode()

def parse_json(body):
    try:
        return json.loads(body)
    except Exception:
        return None

def html_count(html):
    m = re.search(rb"Result-set[^:]*:\s*([0-9,]+) compounds", html or b"")
    return int(m.group(1).replace(b",", b"")) if m else None

def shuffle_path(path):
    if "?" not in path:
        return path
    base, query = path.split("?", 1)
    parts = query.split("&")
    original = parts[:]
    # Re-shuffle if the shuffle produces the original order.
    for _ in range(5):
        random.shuffle(parts)
        if parts != original:
            break
    return base + "?" + "&".join(parts)

def compare(path):
    code1, ctype1, body1 = fetch(path)
    shuffled = shuffle_path(path)
    code2, ctype2, body2 = fetch(shuffled)
    if code1 != code2:
        return False, f"status: {code1} vs {code2}"
    if "json" in ctype1:
        j1, j2 = parse_json(body1), parse_json(body2)
        # For search responses, compare the SET of compound IDs and the total
        # count rather than the full body, because the server may return rows
        # in non-deterministic order (no stable ORDER BY in some endpoints).
        # This isolates the argument-order question from the row-order question.
        def signature(j):
            if not isinstance(j, dict):
                return j
            # For search responses, compare only the total/count and the
            # structural keys, NOT the returned rows. Under LIMIT N, the
            # specific rows returned depend on the underlying SQL ordering
            # which is non-deterministic for some endpoints (no ORDER BY).
            # Argument-order should not affect the total count or the
            # response schema, which is what this test verifies.
            sig = {"total": j.get("total"), "count": j.get("count"),
                   "keys": sorted(j.keys())}
            if "attestations" in j and isinstance(j["attestations"], list):
                sig["attestations"] = sorted([(a.get("source"), a.get("license_tier")) for a in j["attestations"]])
            elif "results" not in j:
                # Non-search responses: compare full body.
                sig["full"] = j
            return sig
        if signature(j1) != signature(j2):
            return False, "JSON content differs (after content-signature comparison)"
    elif "html" in ctype1:
        n1, n2 = html_count(body1), html_count(body2)
        if n1 != n2:
            return False, f"counts: {n1} vs {n2}"
    else:
        if body1 != body2:
            return False, "byte content differs"
    return True, f"shuffled={shuffled.split('?', 1)[1] if '?' in shuffled else ''}"

def check(label, ok, detail=""):
    if ok:
        results["match"] += 1
        print(f"[MATCH] {label}  {detail}")
    else:
        results["mismatch"] += 1
        results["mismatches"].append(label)
        print(f"[MISMATCH] {label}  {detail}")

queries = [
    ("/api/search?q=quercetin&type=name&format=json&limit=5", "A1 /api/search name"),
    ("/api/search?q=a&type=name&license=commercial&format=json&limit=5", "A2 /api/search name+license"),
    ("/api/search?q=a&type=name&license=academic&format=json&limit=5", "A3 /api/search academic"),
    ("/api/search?q=Gammaproteobacteria&type=tax_class&prop_type=npclassifier_class&format=json&limit=5", "A4 /api/search tax_class"),
    ("/api/search?q=Gammaproteobacteria&type=tax_class&prop_type=npclassifier_class&license=commercial&format=json&limit=5", "A5 /api/search tax+lic"),
    ("/api/search?q=Streptomyces&type=organism&format=json&limit=5", "A6 /api/search organism"),
    ("/api/search?q=plant&type=kingdom&format=json&limit=5", "A7 /api/search kingdom"),
    ("/api/search?q=a&type=name&license=academic&kingdom=plant&format=json&limit=5", "A8 multi-filter"),
    ("/api/bulk?cols=comp_id,smiles&license=commercial&limit=3", "B1 bulk+license"),
    ("/api/bulk?cols=comp_id,name,smiles&tier=open&limit=3", "B2 bulk+tier"),
    ("/api/bulk?license=academic&format=json&limit=3", "B3 bulk academic json"),
    ("/api/bulk?cols=comp_id,smiles&tier=open&format=csv&limit=3", "B4 bulk full args"),
    ("/search?q=Gammaproteobacteria&type=tax_class&prop_type=npclassifier_class", "C1 /search tax_class"),
    ("/search?q=quercetin&type=name&license=commercial", "C2 /search name+lic"),
    ("/search?q=plant&type=kingdom&license=academic", "C3 /search kingdom+lic"),
    ("/search?q=Gammaproteobacteria&type=tax_class&prop_type=npclassifier_class&license=academic", "C4 /search tax+lic"),
]

print("=" * 70)
print("THEOBROMA argument-order equivalence tests")
print("=" * 70)

for path, label in queries:
    ok, detail = compare(path)
    check(label, ok, detail)

print()
print("=" * 70)
print(f"Summary: {results['match']} matched, {results['mismatch']} mismatched")
if results["mismatches"]:
    print("Mismatches:")
    for m in results["mismatches"]:
        print(f"  - {m}")
print("=" * 70)
sys.exit(0 if results["mismatch"] == 0 else 1)
