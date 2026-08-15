#!/usr/bin/env python3
"""THEOBROMA deposit integrity suite, 25 checks.

Everything here concerns the artefacts that leave the machine: the dump, the manifest,
the classifier bundle and the checksums. The other five suites verify the live database
and the served application; none of them looks at what is actually deposited.

The organising principle is that a deposit must be self-describing and self-consistent
without access to the live system. A reader restores the dump, opens sources.yaml, and
must be able to reproduce every claim the manuscript makes from those two objects alone.

  PASS    invariant holds
  FAIL    the deposit would misrepresent the resource; fix before upload
  KNOWN   documented incompleteness, asserted so it cannot silently grow
  SKIP    artefact not present at the expected path
  ERROR   the check could not run

Usage:
  python3 deposit_integrity_suite.py [--deposit ~/theobroma_deposit] [--manifest sources.yaml]

Exit 1 on any FAIL or ERROR.
"""
import argparse, glob, hashlib, json, os, re, subprocess, sys, traceback, zipfile

CORPUS = 1132805
TIERS = {"CC BY 4.0": 891860, "CC0": 8243, "CC BY-NC 4.0": 225536, "Unspecified": 7166}
ATTESTATIONS = 1512594
FAMILIES = 486032
NONDEGENERATE = 137149
N_SOURCES = 29
EXPECTED_TABLES = [
    "compounds", "compound_taxonomy", "resolved_taxonomy",
    "per_source_license_attestation", "source_license_ref", "compound_synonyms",
    "compound_region_map", "admet", "scaffolds", "npc_class_parents",
    "npc_super_parents", "phylum_kingdom_map", "apg_clade_ref", "family_lineage_ref",
    "genus_lineage_ref", "genus_ref", "resolved_name_ref", "coconut_organisms",
]
CLASSIFIER_REQUIRED = [
    "xgb_v135.ubj", "pca_components.npy", "pca_mean.npy", "classes.json",
    "hpo_results_v135.json", "hpo_trial_log_v135.csv", "hpo_per_class_eval_v135.csv",
    "threshold_calibration_v135.csv",
]

TESTS, COUNT = [], {"PASS": 0, "FAIL": 0, "KNOWN": 0, "SKIP": 0, "ERROR": 0}


def test(tag, note):
    def deco(fn):
        TESTS.append((tag, note, fn))
        return fn
    return deco


def emit(tag, code, headline, note=None):
    COUNT[code] += 1
    print("[%-5s] %-4s %s" % (code, tag, headline))
    if note:
        print("             %s" % note)


class Ctx:
    def __init__(self, deposit, manifest):
        self.deposit = os.path.expanduser(deposit)
        self.manifest_path = os.path.expanduser(manifest)
        self._man = None
        self._toc = None

    def path(self, *parts):
        return os.path.join(self.deposit, *parts)

    def dump(self):
        c = sorted(glob.glob(self.path("*.dump")), key=os.path.getmtime)
        return c[-1] if c else None

    def manifest(self):
        if self._man is None:
            import yaml
            self._man = yaml.safe_load(open(self.manifest_path, encoding="utf-8"))
        return self._man

    def toc(self):
        """pg_restore -l output, cached."""
        if self._toc is None:
            d = self.dump()
            if not d:
                self._toc = ""
            else:
                self._toc = subprocess.run(["pg_restore", "-l", d], capture_output=True,
                                           text=True, timeout=300).stdout
        return self._toc


# ------------------------------------------------------- A. dump structure

@test("A01", "the dump exists and is a custom-format archive")
def a01(c):
    d = c.dump()
    if not d:
        return "SKIP", "no .dump in %s" % c.deposit
    head = open(d, "rb").read(5)
    mb = os.path.getsize(d) / 1e6
    return head.startswith(b"PGDMP"), "%s, %.0f MB, magic %s" % (
        os.path.basename(d), mb, head[:5].decode("ascii", "replace"))


@test("A02", "exactly the eighteen intended tables, no scratch tables")
def a02(c):
    toc = c.toc()
    if not toc:
        return "SKIP", "no dump to read"
    got = sorted({l.split()[-2] for l in toc.splitlines()
                  if "TABLE DATA" in l and len(l.split()) >= 2})
    missing = [t for t in EXPECTED_TABLES if t not in got]
    extra = [t for t in got if t not in EXPECTED_TABLES]
    return not missing and not extra, "%d tables" % len(got), \
        ("missing %s; extra %s" % (missing or "none", extra or "none")
         if (missing or extra) else "exactly the include list")


@test("A03", "no backup, archive or scratch table leaked in")
def a03(c):
    toc = c.toc()
    if not toc:
        return "SKIP", "no dump to read"
    pat = re.compile(r"(_bak|_pre_|preswap|_archive|_old|_buggy|_broken|^v135|_check$)")
    got = {l.split()[-2] for l in toc.splitlines() if "TABLE DATA" in l and len(l.split()) >= 2}
    bad = sorted(t for t in got if pat.search(t))
    return not bad, "suspicious table names %s" % (bad or "none"), \
        "resolved_taxonomy_buggy in a published dump would be actively misleading"


@test("A04", "the sha256 sidecar matches the dump")
def a04(c):
    d = c.dump()
    if not d:
        return "SKIP", "no dump"
    side = d + ".sha256"
    if not os.path.exists(side):
        return "FAIL", "no .sha256 sidecar beside the dump"
    declared = open(side).read().split()[0]
    h = hashlib.sha256()
    with open(d, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest() == declared, "declared %s… computed %s…" % (
        declared[:16], h.hexdigest()[:16])


@test("A05", "SHA256SUMS covers every deposited file")
def a05(c):
    s = c.path("SHA256SUMS")
    if not os.path.exists(s):
        return "FAIL", "SHA256SUMS absent"
    listed = {l.split()[-1].lstrip("*") for l in open(s) if l.strip()}
    present = {os.path.basename(p) for p in glob.glob(c.path("*"))
               if os.path.isfile(p) and os.path.basename(p) != "SHA256SUMS"}
    uncovered = sorted(present - listed)
    return not uncovered, "%d listed, %d present, %d uncovered" % (
        len(listed), len(present), len(uncovered)), ", ".join(uncovered[:4])


@test("A06", "no superseded dump left in the deposit directory")
def a06(c):
    d = sorted(glob.glob(c.path("*.dump")))
    stale = sorted(glob.glob(c.path("*.superseded")) + glob.glob(c.path("*.dump.bak")))
    return len(d) == 1 and not stale, "%d dump(s), %d superseded artefacts" % (
        len(d), len(stale)), "uploading two dumps invites the reader to pick the wrong one"


# ---------------------------------------------------- B. manifest integrity

@test("B01", "the manifest parses as YAML")
def b01(c):
    m = c.manifest()
    return isinstance(m, dict) and "sources" in m, "top-level keys %s" % sorted(m), \
        "the manuscript calls this machine-readable; it must load"


@test("B02", "source count matches the manuscript")
def b02(c):
    n = len(c.manifest()["sources"])
    return n == N_SOURCES, "sources %d" % n


@test("B03", "in_theobroma sums to the corpus")
def b03(c):
    t = sum(s.get("in_theobroma", 0) for s in c.manifest()["sources"])
    return t == CORPUS, "sum %s against corpus %s" % (format(t, ","), format(CORPUS, ",")), \
        "this makes the manifest a checkable decomposition in two lines"


@test("B04", "no source claims more integrated than it held")
def b04(c):
    bad = [(s["name"], s.get("raw_compounds"), s["in_theobroma"])
           for s in c.manifest()["sources"]
           if isinstance(s.get("raw_compounds"), int)
           and isinstance(s.get("in_theobroma"), int)
           and s["in_theobroma"] > s["raw_compounds"]]
    return not bad, "inversions %s" % (len(bad)), \
        "; ".join("%s %s<%s" % b for b in bad[:3]) or "none"


@test("B05", "every DOI is syntactically a DOI")
def b05(c):
    bad = []
    for s in c.manifest()["sources"]:
        cit = str(s.get("citation", ""))
        if cit.startswith("10.") and not re.match(r"^10\.\d{4,9}/[-._;()/:a-zA-Z0-9]+$", cit):
            bad.append((s["name"], cit))
    return not bad, "malformed DOIs %d" % len(bad), \
        "; ".join("%s %s" % b for b in bad[:3]) or "all well-formed"


@test("B06", "no DOI appears against two different sources")
def b06(c):
    seen = {}
    dup = []
    for s in c.manifest()["sources"]:
        cit = str(s.get("citation", ""))
        if cit.startswith("10."):
            if cit in seen:
                dup.append((seen[cit], s["name"], cit))
            seen[cit] = s["name"]
    return not dup, "duplicated DOIs %d" % len(dup), \
        "; ".join("%s/%s %s" % d for d in dup[:3]) or "none"


@test("B07", "sources without a DOI say so rather than guessing")
def b07(c):
    nod = [s["name"] for s in c.manifest()["sources"]
           if not str(s.get("citation", "")).startswith("10.")]
    return ("KNOWN" if nod else "PASS"), "%d without a DOI" % len(nod), \
        ", ".join(nod) + "  (each must have no publication, not an unverified one)"


@test("B08", "every source carries a licence and a download date")
def b08(c):
    bad = [s["name"] for s in c.manifest()["sources"]
           if not s.get("license_data") or not s.get("download_date")]
    return not bad, "incomplete provenance %s" % (bad or "none")


@test("B09", "conservative assignments record their basis")
def b09(c):
    flagged = [s["name"] for s in c.manifest()["sources"]
               if "conservative" in str(s.get("license_data", "")).lower()]
    return ("KNOWN" if flagged else "PASS"), "%d assignments flagged conservative" % len(flagged), \
        ", ".join(flagged) + "  (each departs from the source's own stated terms)"


@test("B10", "manifest licences agree with the deposited licence map")
def b10(c):
    import psycopg2
    con = psycopg2.connect("host=localhost dbname=theobroma user=theobroma")
    cur = con.cursor()
    cur.execute("SELECT lower(src), license_tier FROM source_license_ref")
    db = dict(cur.fetchall())
    con.close()
    bad = []
    for s in c.manifest()["sources"]:
        key = re.sub(r"[^a-z0-9_-]", "", s["name"].split()[0].lower())
        tier = db.get(key)
        m = re.match(r"(CC0|CC BY-NC-SA 4\.0|CC BY-NC 4\.0|CC BY 4\.0|Unspecified)",
                     str(s.get("license_data", "")))
        if tier and m and tier != m.group(1):
            bad.append((s["name"], m.group(1), tier))
    return not bad, "disagreements %d" % len(bad), \
        "; ".join("%s manifest %s db %s" % b for b in bad[:3]) or "all agree"


# ------------------------------------------- C. classifier reproducibility

@test("C01", "the classifier bundle is present")
def c01(c):
    z = c.path("classifier_reproducibility.zip")
    d = c.path("classifier_reproducibility")
    if os.path.exists(z):
        names = [os.path.basename(n) for n in zipfile.ZipFile(z).namelist()]
        return True, "zip with %d entries, %.0f MB" % (len(names), os.path.getsize(z) / 1e6)
    if os.path.isdir(d):
        return True, "directory with %d files" % len(os.listdir(d))
    return "SKIP", "no classifier bundle at either path"


@test("C02", "the model is unusable without the PCA basis, so both must ship")
def c02(c):
    z = c.path("classifier_reproducibility.zip")
    d = c.path("classifier_reproducibility")
    if os.path.exists(z):
        names = {os.path.basename(n) for n in zipfile.ZipFile(z).namelist()}
    elif os.path.isdir(d):
        names = set(os.listdir(d))
    else:
        return "SKIP", "no bundle"
    missing = [f for f in CLASSIFIER_REQUIRED if f not in names]
    return not missing, "%d of %d required files" % (
        len(CLASSIFIER_REQUIRED) - len(missing), len(CLASSIFIER_REQUIRED)), \
        "missing " + ", ".join(missing) if missing else "complete"


@test("C03", "the HPO record reports the deployed configuration")
def c03(c):
    z = c.path("classifier_reproducibility.zip")
    raw = None
    if os.path.exists(z):
        with zipfile.ZipFile(z) as zf:
            n = [x for x in zf.namelist() if x.endswith("hpo_results_v135.json")]
            if n:
                raw = zf.read(n[0])
    else:
        p = c.path("classifier_reproducibility", "hpo_results_v135.json")
        if os.path.exists(p):
            raw = open(p, "rb").read()
    if raw is None:
        return "SKIP", "hpo_results_v135.json not in the bundle"
    j = json.loads(raw)
    txt = json.dumps(j)
    return ("0.809" in txt or "0.80" in txt), "keys %s" % sorted(j)[:6], \
        "the main text quotes 0.809 CV and 0.825 held-out"


@test("C04", "no absolute paths from the build machine leak into the bundle")
def c04(c):
    z = c.path("classifier_reproducibility.zip")
    if not os.path.exists(z):
        return "SKIP", "no zip"
    bad = [n for n in zipfile.ZipFile(z).namelist()
           if n.startswith("/") or n.startswith("..") or "workspace" in n]
    return not bad, "leaking entries %d" % len(bad), ", ".join(bad[:3]) or "none"


# ------------------------------------------------- D. deposit self-description

@test("D01", "a README explains how to restore")
def d01(c):
    cand = glob.glob(c.path("*README*")) + glob.glob(c.path("*.md"))
    if not cand:
        return "FAIL", "no README in the deposit"
    txt = open(cand[0], encoding="utf-8", errors="replace").read().lower()
    has_pg = "pg_restore" in txt
    has_ext = "pg_trgm" in txt
    return has_pg and has_ext, "%s: pg_restore %s, pg_trgm %s" % (
        os.path.basename(cand[0]), has_pg, has_ext), \
        "restore emits one error without CREATE EXTENSION pg_trgm"


@test("D02", "the README states the corpus size and licence split")
def d02(c):
    cand = glob.glob(c.path("*README*")) + glob.glob(c.path("*.md"))
    if not cand:
        return "SKIP", "no README"
    txt = open(cand[0], encoding="utf-8", errors="replace").read()
    hits = [v for v in ("1,132,805", "1132805") if v in txt]
    tiers = sum(1 for v in ("891,860", "225,536", "8,243", "7,166") if v in txt)
    return bool(hits) and tiers >= 2, "corpus figure %s, tier figures %d of 4" % (
        bool(hits), tiers)


@test("D03", "the manifest is deposited alongside the dump")
def d03(c):
    p = c.path("sources.yaml")
    return os.path.exists(p), "sources.yaml in the deposit directory: %s" % os.path.exists(p), \
        "the manuscript cites it by name as the reproducibility manifest"


@test("D04", "no working files or editor residue in the deposit")
def d04(c):
    junk = []
    for p in glob.glob(c.path("*")):
        b = os.path.basename(p)
        if b.endswith(("~", ".swp", ".tmp", ".log", ".bak")) or b.startswith("."):
            junk.append(b)
    return not junk, "stray files %s" % (junk or "none")


@test("D05", "total deposit size is plausible for the stated contents")
def d05(c):
    tot = sum(os.path.getsize(p) for p in glob.glob(c.path("**/*"), recursive=True)
              if os.path.isfile(p))
    return 3e8 < tot < 3e9, "total %.0f MB" % (tot / 1e6), \
        "under 300 MB suggests a truncated dump; over 3 GB suggests the vectors crept back in"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deposit", default="~/theobroma_deposit")
    ap.add_argument("--manifest", default="sources.yaml")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    ctx = Ctx(a.deposit, a.manifest)
    print("THEOBROMA deposit integrity  %d checks  deposit=%s" % (len(TESTS), ctx.deposit))
    print("-" * 78)
    for tag, note, fn in TESTS:
        try:
            r = fn(ctx)
            code = r[0] if isinstance(r[0], str) else ("PASS" if r[0] else "FAIL")
            emit(tag, code, r[1], r[2] if len(r) > 2 else note)
        except Exception as e:
            emit(tag, "ERROR", "%s: %s" % (type(e).__name__, str(e).strip()[:70]), note)
            if a.verbose:
                traceback.print_exc()
    print("-" * 78)
    print("SUMMARY  PASS=%d FAIL=%d KNOWN=%d SKIP=%d ERROR=%d"
          % (COUNT["PASS"], COUNT["FAIL"], COUNT["KNOWN"], COUNT["SKIP"], COUNT["ERROR"]))
    return 1 if COUNT["FAIL"] or COUNT["ERROR"] else 0


if __name__ == "__main__":
    sys.exit(main())
