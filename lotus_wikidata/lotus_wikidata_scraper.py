#!/usr/bin/env python3
"""
Wikidata LOTUS SPARQL harvester v3 for THEOBROMA integration.

Avoids wdt:P171* (recursive taxonomy traversal) which times out on
Wikidata's 60s query deadline. Instead:
  Mode 1 (--all):   Download all NP compound-organism pairs, filter later
  Mode 2 (--genera): Query by explicit genus/species names (fast, no traversal)
  Mode 3 (--taxa-file): Query by a user-provided list of taxon QIDs

Usage:
    # Download ALL LOTUS NPs from Wikidata (~250k pairs, ~2-3 hours)
    python lotus_wd_harvest.py --all -o all_np.csv

    # Query specific fungal genera (fast, minutes)
    python lotus_wd_harvest.py --genera Ganoderma Hericium Inonotus Psilocybe -o mushroom_np.csv

    # Query from a file of genera (one per line)
    python lotus_wd_harvest.py --genera-file fungal_genera.txt -o fungi_np.csv

    # Query specific taxon QIDs (from previous harvest or manual curation)
    python lotus_wd_harvest.py --taxa-file taxon_qids.txt -o fungi_np.csv

    # With deduplication
    python lotus_wd_harvest.py --all -o novel.csv --dedup-file existing_inchikeys.txt

Dependencies: requests (pip install requests)
"""
from __future__ import annotations
import argparse, time, sys, csv, os

try:
    import requests
except ImportError:
    sys.exit("pip install requests")

SPARQL_URL = "https://query.wikidata.org/sparql"
UA = "THEOBROMA-NP-Harvester/3.0 (natural products database project)"
MAX_RETRIES = 6
RETRY_WAIT = 45
PAGE_SIZE = 5000

# Common basidiomycete genera for quick testing
BASIDIOMYCETE_GENERA = [
    "Ganoderma", "Hericium", "Inonotus", "Trametes", "Pleurotus",
    "Agaricus", "Psilocybe", "Omphalotus", "Antrodia", "Cordyceps",
    "Phellinus", "Fomes", "Fomitopsis", "Laetiporus", "Schizophyllum",
    "Coprinopsis", "Coprinus", "Lentinula", "Clitocybe", "Armillaria",
    "Boletus", "Russula", "Lactarius", "Cortinarius", "Amanita",
    "Lycoperdon", "Cyathus", "Strobilurus", "Oudemansiella", "Sarcodon",
    "Gymnopilus", "Panaeolus", "Conocybe", "Pluteus", "Inocybe",
    "Tricholoma", "Flammulina", "Grifola", "Auricularia", "Tremella",
    "Ustilago", "Puccinia", "Claviceps",
]

ASCOMYCETE_GENERA = [
    "Aspergillus", "Penicillium", "Fusarium", "Trichoderma", "Cladosporium",
    "Alternaria", "Chaetomium", "Xylaria", "Daldinia", "Hypoxylon",
    "Pestalotiopsis", "Colletotrichum", "Diaporthe", "Epicoccum",
    "Talaromyces", "Emericella", "Neosartorya", "Eurotium",
]


def sparql(query, attempt=0):
    """Run SPARQL with retry on 429/500/502/503/504/timeout."""
    headers = {"Accept": "application/json", "User-Agent": UA}
    try:
        r = requests.get(SPARQL_URL, params={"query": query},
                         headers=headers, timeout=90)
        if r.status_code in (429, 500, 502, 503, 504):
            if attempt < MAX_RETRIES:
                w = RETRY_WAIT * (attempt + 1)
                print(f"    HTTP {r.status_code}, retry in {w}s...", flush=True)
                time.sleep(w)
                return sparql(query, attempt + 1)
            raise RuntimeError(f"Max retries on HTTP {r.status_code}")
        r.raise_for_status()
        return r.json()["results"]["bindings"]
    except requests.exceptions.Timeout:
        if attempt < MAX_RETRIES:
            w = RETRY_WAIT * (attempt + 1)
            print(f"    Timeout, retry in {w}s...", flush=True)
            time.sleep(w)
            return sparql(query, attempt + 1)
        raise


def v(b, k):
    return b.get(k, {}).get("value", "")


def qid(uri):
    return uri.rsplit("/", 1)[-1] if uri else ""


# --- Query builders ----------------------------------------------------------

def query_all(offset, limit):
    """All compounds with SMILES + organism link. No taxonomy filter."""
    return f"""
SELECT ?compound ?compoundLabel ?smiles ?inchikey ?inchi ?taxon ?taxonLabel
WHERE {{
  {{
    SELECT ?compound ?smiles ?inchikey ?inchi ?taxon WHERE {{
      ?compound wdt:P233 ?smiles .
      ?compound wdt:P703 ?taxon .
      OPTIONAL {{ ?compound wdt:P235 ?inchikey . }}
      OPTIONAL {{ ?compound wdt:P234 ?inchi . }}
    }}
    LIMIT {limit} OFFSET {offset}
  }}
  SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
}}"""


def query_by_genus(genus_name, offset, limit):
    """Compounds from organisms whose label starts with a genus name.
    No taxonomy traversal needed -- pure string match on taxon label."""
    return f"""
SELECT ?compound ?compoundLabel ?smiles ?inchikey ?inchi ?taxon ?taxonLabel
WHERE {{
  {{
    SELECT ?compound ?smiles ?inchikey ?inchi ?taxon WHERE {{
      ?taxon rdfs:label ?taxonName .
      FILTER(LANG(?taxonName) = "en")
      FILTER(STRSTARTS(?taxonName, "{genus_name}"))
      ?compound wdt:P703 ?taxon .
      ?compound wdt:P233 ?smiles .
      OPTIONAL {{ ?compound wdt:P235 ?inchikey . }}
      OPTIONAL {{ ?compound wdt:P234 ?inchi . }}
    }}
    LIMIT {limit} OFFSET {offset}
  }}
  SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
}}"""


def query_by_taxa_qids(qids, offset, limit):
    """Compounds for an explicit set of taxon QIDs via VALUES."""
    values = " ".join(f"wd:{q}" for q in qids)
    return f"""
SELECT ?compound ?compoundLabel ?smiles ?inchikey ?inchi ?taxon ?taxonLabel
WHERE {{
  {{
    SELECT ?compound ?smiles ?inchikey ?inchi ?taxon WHERE {{
      VALUES ?taxon {{ {values} }}
      ?compound wdt:P703 ?taxon .
      ?compound wdt:P233 ?smiles .
      OPTIONAL {{ ?compound wdt:P235 ?inchikey . }}
      OPTIONAL {{ ?compound wdt:P234 ?inchi . }}
    }}
    LIMIT {limit} OFFSET {offset}
  }}
  SERVICE wikibase:label {{ bd:serviceParam wikibase:language "en". }}
}}"""


# --- Core harvest logic ------------------------------------------------------

FIELDS = ["compound_qid", "compound_name", "smiles", "inchikey",
          "inchi", "taxon_qid", "taxon_name"]


def load_dedup(path):
    keys = set()
    if path and os.path.exists(path):
        with open(path) as f:
            for line in f:
                k = line.strip()
                if len(k) == 27:
                    keys.add(k)
        print(f"Loaded {len(keys)} InChIKeys for dedup")
    return keys


def write_rows(writer, rows, seen, existing):
    """Process SPARQL result rows, write to CSV, return (written, skipped)."""
    written = skipped = 0
    for b in rows:
        cq, tq = qid(v(b, "compound")), qid(v(b, "taxon"))
        if (cq, tq) in seen:
            continue
        seen.add((cq, tq))
        ik = v(b, "inchikey")
        if ik and ik in existing:
            skipped += 1
            continue
        writer.writerow({
            "compound_qid": cq, "compound_name": v(b, "compoundLabel"),
            "smiles": v(b, "smiles"), "inchikey": ik,
            "inchi": v(b, "inchi"), "taxon_qid": tq,
            "taxon_name": v(b, "taxonLabel"),
        })
        written += 1
    return written, skipped


def paginate(query_fn, writer, seen, existing, label=""):
    """Generic pagination loop for a query builder function."""
    total_w = total_s = 0
    offset = 0
    while True:
        q = query_fn(offset, PAGE_SIZE)
        print(f"  {label}offset={offset}...", end=" ", flush=True)
        rows = sparql(q)
        print(f"{len(rows)} rows")
        if not rows:
            break
        w, s = write_rows(writer, rows, seen, existing)
        total_w += w
        total_s += s
        if len(rows) < PAGE_SIZE:
            break
        offset += PAGE_SIZE
        time.sleep(2)
    return total_w, total_s


def harvest_all(output, dedup_file):
    existing = load_dedup(dedup_file)
    seen = set()
    total_w = total_s = 0
    print("Harvesting ALL LOTUS/Wikidata NPs (no taxonomy filter)...")
    with open(output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        w, s = paginate(query_all, writer, seen, existing)
        total_w += w
        total_s += s
    print(f"\nDone. {total_w} pairs -> {output}")
    if total_s:
        print(f"  {total_s} skipped (dedup)")


def harvest_genera(genera, output, dedup_file):
    existing = load_dedup(dedup_file)
    seen = set()
    total_w = total_s = 0
    print(f"Harvesting compounds for {len(genera)} genera...")
    with open(output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for i, genus in enumerate(genera):
            label = f"[{i+1}/{len(genera)}] {genus}: "
            print(f"\n{label}")
            fn = lambda off, lim, g=genus: query_by_genus(g, off, lim)
            w, s = paginate(fn, writer, seen, existing, label="  ")
            total_w += w
            total_s += s
            time.sleep(1)
    print(f"\nDone. {total_w} pairs -> {output}")
    if total_s:
        print(f"  {total_s} skipped (dedup)")


def harvest_taxa_file(taxa_file, output, dedup_file):
    with open(taxa_file) as f:
        qids = [line.strip() for line in f if line.strip().startswith("Q")]
    if not qids:
        sys.exit(f"No QIDs found in {taxa_file}")
    existing = load_dedup(dedup_file)
    seen = set()
    total_w = total_s = 0
    batch_size = 40
    print(f"Harvesting compounds for {len(qids)} taxon QIDs...")
    with open(output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for i in range(0, len(qids), batch_size):
            batch = qids[i:i+batch_size]
            label = f"batch {i//batch_size+1}/{(len(qids)-1)//batch_size+1}: "
            fn = lambda off, lim, b=batch: query_by_taxa_qids(b, off, lim)
            w, s = paginate(fn, writer, seen, existing, label=label)
            total_w += w
            total_s += s
            time.sleep(1)
    print(f"\nDone. {total_w} pairs -> {output}")
    if total_s:
        print(f"  {total_s} skipped (dedup)")


def main():
    p = argparse.ArgumentParser(
        description="Harvest LOTUS NPs from Wikidata (v3, no taxonomy traversal)")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--all", action="store_true",
                   help="Download all NP-organism pairs (hours, ~250k+)")
    g.add_argument("--genera", nargs="+", metavar="GENUS",
                   help="Query by genus name(s), e.g. --genera Ganoderma Hericium")
    g.add_argument("--genera-file", type=str,
                   help="File with one genus name per line")
    g.add_argument("--preset", choices=["basidiomycetes", "ascomycetes", "all_fungi"],
                   help="Use built-in genus lists")
    g.add_argument("--taxa-file", type=str,
                   help="File with one Wikidata QID per line")
    p.add_argument("-o", "--output", default="lotus_wikidata_harvest.csv")
    p.add_argument("--dedup-file", type=str, default=None)
    args = p.parse_args()

    if args.all:
        harvest_all(args.output, args.dedup_file)
    elif args.genera:
        harvest_genera(args.genera, args.output, args.dedup_file)
    elif args.genera_file:
        with open(args.genera_file) as f:
            genera = [line.strip() for line in f if line.strip()]
        harvest_genera(genera, args.output, args.dedup_file)
    elif args.preset:
        if args.preset == "basidiomycetes":
            harvest_genera(BASIDIOMYCETE_GENERA, args.output, args.dedup_file)
        elif args.preset == "ascomycetes":
            harvest_genera(ASCOMYCETE_GENERA, args.output, args.dedup_file)
        elif args.preset == "all_fungi":
            harvest_genera(BASIDIOMYCETE_GENERA + ASCOMYCETE_GENERA,
                          args.output, args.dedup_file)
    elif args.taxa_file:
        harvest_taxa_file(args.taxa_file, args.output, args.dedup_file)


if __name__ == "__main__":
    main()