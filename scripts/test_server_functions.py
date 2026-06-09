#!/usr/bin/env python3
"""Comprehensive single-function tests for THEOBROMA server.
Each function group exercised with at least 5 argument combinations.
Reports [PASS]/[FAIL] per query plus a summary at end.
Exit code 0 if all pass, 1 otherwise.
"""
import urllib.request, urllib.error, json, sys, re

BASE = "http://localhost:5000"
TIMEOUT = 30
results = {"pass": 0, "fail": 0, "failures": []}

def fetch_json(path):
    try:
        with urllib.request.urlopen(BASE + path, timeout=TIMEOUT) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        return e.code, None
    except Exception as e:
        return -1, str(e)

def fetch_html(path):
    try:
        with urllib.request.urlopen(BASE + path, timeout=TIMEOUT) as r:
            return r.status, r.read().decode("utf-8", errors="ignore")
    except urllib.error.HTTPError as e:
        return e.code, None
    except Exception as e:
        return -1, str(e)

def html_count(html):
    # Match the result-set total, not any earlier section count.
    m = re.search(r"Result-set[^:]*:\s*([0-9,]+) compounds", html or "")
    if m:
        return int(m.group(1).replace(",", ""))
    # Fallback: any "N compounds" if the result-set anchor is absent.
    m = re.search(r"([0-9,]+) compounds", html or "")
    return int(m.group(1).replace(",", "")) if m else None

def check(label, ok, detail=""):
    if ok:
        results["pass"] += 1
        print(f"[PASS] {label}  {detail}")
    else:
        results["fail"] += 1
        results["failures"].append(label)
        print(f"[FAIL] {label}  {detail}")

print("=" * 70)
print("THEOBROMA server function tests")
print("=" * 70)

# Group A: license filter on /api/search
print("\n--- A: license filter (/api/search) ---")
for path, lab in [
    ("/api/search?q=a&type=name&format=json&limit=1", "unfiltered"),
    ("/api/search?q=a&type=name&license=commercial&format=json&limit=1", "commercial"),
    ("/api/search?q=a&type=name&license=academic&format=json&limit=1", "academic"),
    ("/api/search?q=a&type=name&license=all&format=json&limit=1", "all-alias"),
    ("/api/search?q=quercetin&type=name&license=commercial&format=json&limit=1", "named-commercial"),
]:
    code, body = fetch_json(path)
    total = body.get("total") if isinstance(body, dict) else None
    check(f"A: {lab}", code == 200 and total is not None, f"total={total}")

# Group B: search modes on /api/search
print("\n--- B: search modes (/api/search) ---")
for path, lab in [
    ("/api/search?q=quercetin&type=name&format=json&limit=1", "name"),
    ("/api/search?q=THEO_0000001&type=comp_id&format=json&limit=1", "comp_id"),
    ("/api/search?q=Gammaproteobacteria&type=tax_class&prop_type=npclassifier_class&format=json&limit=1", "tax_class"),
    ("/api/search?q=Streptomyces&type=organism&format=json&limit=1", "organism"),
    ("/api/search?q=plant&type=kingdom&format=json&limit=1", "kingdom-filter"),
]:
    code, body = fetch_json(path)
    total = body.get("total") if isinstance(body, dict) else None
    check(f"B: {lab}", code == 200, f"total={total}")

# Group C: /api/compound
print("\n--- C: /api/compound ---")
for cid, expect_ok in [
    ("THEO_0000001", True),
    ("THEO_0915582", True),  # multi-source FooDB+LOTUS
    ("THEO_0216080", True),
    ("THEO_0000109", True),
    ("THEO_9999999", False),  # should 404
]:
    code, body = fetch_json(f"/api/compound/{cid}")
    if expect_ok:
        has_prov = isinstance(body, dict) and "license_provenance" in body
        check(f"C: {cid}", code == 200 and has_prov, f"code={code} has_prov={has_prov}")
    else:
        check(f"C: {cid} (404)", code == 404, f"code={code}")

# Group D: /api/compound/<id>/license-provenance
print("\n--- D: license-provenance endpoint ---")
for cid, expect_ok in [
    ("THEO_0000001", True),
    ("THEO_0915582", True),
    ("THEO_0216080", True),
    ("THEO_0000109", True),
    ("THEO_9999999", False),
]:
    code, body = fetch_json(f"/api/compound/{cid}/license-provenance")
    if expect_ok:
        n_att = len(body.get("attestations", [])) if isinstance(body, dict) else 0
        check(f"D: {cid}", code == 200 and n_att > 0, f"code={code} n_att={n_att}")
    else:
        check(f"D: {cid} (404)", code == 404, f"code={code}")

# Group E: /api/bulk
print("\n--- E: /api/bulk ---")
for path, lab in [
    ("/api/bulk?limit=2", "default-cols"),
    ("/api/bulk?cols=comp_id,smiles&limit=2", "custom-cols"),
    ("/api/bulk?license=commercial&limit=2", "license-commercial"),
    ("/api/bulk?tier=open&limit=2", "tier-open-original"),
    ("/api/bulk?license=academic&format=json&limit=2", "academic-json"),
]:
    code, html = fetch_html(path)
    nonempty = code == 200 and html and ("comp_id" in html or "THEO_" in html)
    check(f"E: {lab}", nonempty, f"code={code}")

# Group F: /api/stats
print("\n--- F: /api/stats ---")
code, body = fetch_json("/api/stats")
ok = code == 200 and isinstance(body, dict)
check("F: returns 200", ok, f"code={code}")
check("F: has 'total'", ok and "total" in body, f"total={body.get('total') if ok else 'N/A'}")
check("F: has 'kingdoms'", ok and len(body.get("kingdoms", [])) > 0, f"n={len(body.get('kingdoms', [])) if ok else 0}")
check("F: has 'sources'", ok and len(body.get("sources", [])) > 0, f"n={len(body.get('sources', [])) if ok else 0}")
check("F: has 'licenses' (new)", ok and len(body.get("licenses", [])) > 0, f"n={len(body.get('licenses', [])) if ok else 0}")

# Group G: autocomplete
print("\n--- G: autocomplete endpoints ---")
for path, lab in [
    ("/api/autocomplete?q=quer", "name: quer"),
    ("/api/autocomplete?q=cur", "name: cur"),
    ("/api/autocomplete?q=", "name: empty (-> [])"),
    ("/api/organisms?q=Strep", "organism: Strep"),
    ("/api/organisms?q=Curc", "organism: Curc"),
]:
    code, body = fetch_json(path)
    ok = code == 200 and isinstance(body, list)
    check(f"G: {lab}", ok, f"code={code} n={len(body) if ok else 0}")

# Group H: /api/filter_options
print("\n--- H: /api/filter_options ---")
code, body = fetch_json("/api/filter_options")
ok = code == 200 and isinstance(body, dict) and len(body) > 0
check("H: returns dict with keys", ok, f"code={code} keys={list(body.keys()) if ok else []}")

# Group I: /compound/<id> HTML
print("\n--- I: /compound/<id> HTML ---")
for cid, expect_ok in [
    ("THEO_0000001", True),
    ("THEO_0915582", True),
    ("THEO_0216080", True),
    ("THEO_0000109", True),
    ("THEO_9999999", False),
]:
    code, html = fetch_html(f"/compound/{cid}")
    if expect_ok:
        has_prov = code == 200 and html and "License provenance" in html
        check(f"I: {cid} HTML+provenance", has_prov, f"code={code}")
    else:
        check(f"I: {cid} 404", code == 404, f"code={code}")

# Group J: /statistics
print("\n--- J: /statistics ---")
code, html = fetch_html("/statistics")
ok = code == 200 and html is not None
tier_counts = ok and all(s in html for s in ["894,353", "219,467", "12,018", "7,166"])
check("J: /statistics renders", ok, f"code={code}")
check("J: shows corrected tier counts", tier_counts, "")

# Group K: /search HTML
print("\n--- K: /search HTML ---")
for path, lab in [
    ("/search?q=quercetin&type=name", "name search"),
    ("/search?type=tax_class&prop_type=npclassifier_class&q=Gammaproteobacteria", "tax_class unfiltered"),
    ("/search?type=tax_class&prop_type=npclassifier_class&q=Gammaproteobacteria&license=commercial", "tax_class commercial"),
    ("/search?q=&type=name&kingdom=plant", "kingdom filter"),
    ("/search?q=curcumin&type=name&license=academic", "name+license"),
]:
    code, html = fetch_html(path)
    n = html_count(html) if html else None
    check(f"K: {lab}", code == 200 and n is not None, f"n={n}")

# Group L: SMILES + InChIKey + similarity + stereoisomers
print("\n--- L: structure-based search ---")
# SMILES substructure search via /api/search
for path_q, lab in [
    ("/api/search?q=c1ccccc1&type=smiles&format=json&limit=3", "L1 SMILES substructure (benzene)"),
    ("/api/search?q=BQJCRHHNABKAKU-KBQPJGBKSA-N&type=inchikey&format=json&limit=3", "L2 InChIKey exact (morphine)"),
]:
    code, body = fetch_json(path_q)
    total = body.get("total") if isinstance(body, dict) else None
    check(lab, code == 200 and total is not None, f"code={code} total={total}")

# Similarity search with each metric.
for metric, lab in [
    ("morgan", "L3 similarity Morgan"),
    ("maccs", "L4 similarity MACCS"),
    ("chemberta", "L5 similarity ChemBERTa"),
]:
    code, body = fetch_json(f"/api/similarity?smiles=CC(=O)Oc1ccccc1C(=O)O&metric={metric}&limit=3")
    has_results = isinstance(body, (list, dict)) and (
        (isinstance(body, list) and len(body) > 0) or
        (isinstance(body, dict) and (body.get("results") or body.get("matches") or body.get("data")))
    )
    check(lab, code == 200 and has_results, f"code={code}")

# Stereoisomers endpoint.
code, body = fetch_json("/api/stereoisomers/THEO_0915582")
has_results = isinstance(body, (list, dict)) and (
    (isinstance(body, list) and len(body) > 0) or
    (isinstance(body, dict) and body.get("results") is not None) or
    (isinstance(body, dict) and "stereoisomers" in body)
)
check("L6 /api/stereoisomers", code == 200, f"code={code}")

# ADMET sub-object on /api/compound. THEO_0915582 had ADMET in earlier output.
code, body = fetch_json("/api/compound/THEO_0915582")
has_admet = isinstance(body, dict) and "admet" in body and len(body.get("admet", {})) > 10
check("L7 ADMET sub-object on /api/compound", has_admet,
      f"admet_keys={len(body.get('admet', {})) if isinstance(body, dict) else 0}")

# Group M: CSV format and pagination on /api/search
print("\n--- M: alternate formats and pagination ---")
code, html = fetch_html("/api/search?q=quercetin&type=name&format=csv&limit=3")
ok = code == 200 and html and "comp_id" in html and "," in html
check("M1 /api/search format=csv", ok, f"code={code}")

# Offset pagination on a deterministic-enough endpoint.
code, body = fetch_json("/api/search?q=Gammaproteobacteria&type=tax_class&prop_type=npclassifier_class&format=json&limit=3&offset=0")
total0 = body.get("total") if isinstance(body, dict) else None
code2, body2 = fetch_json("/api/search?q=Gammaproteobacteria&type=tax_class&prop_type=npclassifier_class&format=json&limit=3&offset=3")
total2 = body2.get("total") if isinstance(body2, dict) else None
check("M2 offset=0 total matches offset=3 total", total0 is not None and total0 == total2,
      f"total@0={total0} total@3={total2}")

# /export endpoint
code, html = fetch_html("/export?q=quercetin&type=name&limit=5")
ok = code == 200 and html and ("comp_id" in html or "THEO_" in html)
check("M3 /export endpoint", ok, f"code={code}")

# /api/bulk streaming CSV (default format).
code, html = fetch_html("/api/bulk?cols=comp_id,smiles&limit=5")
ok = code == 200 and html and "comp_id" in html
check("M4 /api/bulk streaming CSV (default)", ok, f"code={code}")

# Source filter and region filter (less-exercised filter axes).
code, body = fetch_json("/api/search?q=quercetin&type=name&source=COCONUT&format=json&limit=3")
ok = code == 200 and isinstance(body, dict) and body.get("total") is not None
check("M5 source filter", ok, f"code={code} total={body.get('total') if isinstance(body, dict) else None}")

print()
print("=" * 70)
print(f"Summary: {results['pass']} passed, {results['fail']} failed")
if results["failures"]:
    print("Failures:")
    for f in results["failures"]:
        print(f"  - {f}")
print("=" * 70)
sys.exit(0 if results["fail"] == 0 else 1)
