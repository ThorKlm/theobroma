#!/usr/bin/env python3
"""
FDP scraper (final) — all slugs confirmed working as of Apr 2026.
Covers LMDB (2145), TMDB (1276), SMDB (4084), CMDB (1069).

Usage:
    python fdp_scraper.py --db lmdb -o lmdb.csv
    python fdp_scraper.py --db smdb -o smdb.csv
    python fdp_scraper.py --db all -o fdp_all.csv
"""
from __future__ import annotations
import argparse, csv, json, time, sys

try:
    import requests
except ImportError:
    sys.exit("pip install requests")

API = "https://research-paper-27y2.onrender.com/api"
UA = "THEOBROMA-FDP-Scraper/5.0"
FIELDS = ["database", "slug_used", "species", "sn", "cid", "common_name",
          "iupac_name", "molecular_formula", "molecular_weight",
          "canonical_smiles", "reference", "super_class", "compound_class"]

# =============================================================================
# LMDB — dual approach: genus slugs + camelCase species slugs
# Path: /api/lichmdb/{slug}/
# =============================================================================
LMDB_SLUGS = [
    # Genus-level (catch-all for single-species genera)
    ("alectoria", "[genus]"), ("cetraria", "[genus]"), ("cladonia", "[genus]"),
    ("dimelaena", "[genus]"), ("evernia", "[genus]"), ("flavoparmelia", "[genus]"),
    ("heterodermia", "[genus]"), ("menegazzia", "[genus]"), ("ochrolechia", "[genus]"),
    ("oxneria", "[genus]"), ("platismatia", "[genus]"), ("ramalina", "[genus]"),
    ("ricasolia", "[genus]"), ("roccella", "[genus]"), ("sphaerophorus", "[genus]"),
    ("sulcaria", "[genus]"), ("teloschistes", "[genus]"), ("xanthoparmelia", "[genus]"),
    ("xalocoa", "[genus]"),
    # Species-level camelCase (catches species missed by genus slug)
    ("alectoriaSarmentosa", "Alectoria sarmentosa"),
    ("anziaHypoleucoides", "Anzia hypoleucoides"),
    ("candelariaConcolor", "Candelaria concolor"),
    ("cladoniaConvoluta", "Cladonia convoluta"),
    ("cladoniaCristatella", "Cladonia cristatella"),
    ("cladoniaCryptochlorophaea", "Cladonia cryptochlorophaea"),
    ("cladoniaFurcata", "Cladonia furcata"),
    ("cladoniaMacilenta", "Cladonia macilenta"),
    ("cladoniaMitis", "Cladonia mitis"),
    ("cladoniaRangiferina", "Cladonia rangiferina"),
    ("cladoniaVerticillata", "Cladonia verticillata"),
    ("coniocarponCinnabarinum", "Coniocarpon cinnabarinum"),
    ("dolichousnealongissima", "Dolichousnea longissima"),
    ("dolichousneaLongissima", "Dolichousnea longissima"),
    ("flavoparmeliaHaysomii", "Flavoparmelia haysomii"),
    ("hypogymniaEnteromorpha", "Hypogymnia enteromorpha"),
    ("hypotrachynaCirrhata", "Hypotrachyna cirrhata"),
    ("hypotrachynaLeiophylla", "Hypotrachyna leiophylla"),
    ("lecanoraFrustulosa", "Lecanora frustulosa"),
    ("lecanoraHybocarpa", "Lecanora hybocarpa"),
    ("lecanoraLividocinerea", "Lecanora lividocinerea"),
    ("leptogiumSaturninum", "Leptogium saturninum"),
    ("lethariellaCladonioides", "Lethariella cladonioides"),
    ("lethariellaSernanderi", "Lethariella sernanderi"),
    ("lobariaPulmonaria", "Lobaria pulmonaria"),
    ("lobariaScrobiculata", "Lobaria scrobiculata"),
    ("ophioparmaVentosa", "Ophioparma ventosa"),
    ("parmeliaOmphalodes", "Parmelia omphalodes"),
    ("parmeliaSaxatilis", "Parmelia saxatilis"),
    ("parmotremaNilgherrense", "Parmotrema nilgherrense"),
    ("parmotremaPraesorediosum", "Parmotrema praesorediosum"),
    ("parmotremaSanctiAngelii", "Parmotrema sancti-angelii"),
    ("parmotremaTinctorum", "Parmotrema tinctorum"),
    ("peltigeraCanina", "Peltigera canina"),
    ("peltigeraMalacea", "Peltigera malacea"),
    ("peltigeraVenosa", "Peltigera venosa"),
    ("pseudocyphellariaFaveolata", "Pseudocyphellaria faveolata"),
    ("pseudocyphellariaNorvegica", "Pseudocyphellaria norvegica"),
    ("pyrenulaPseudobufonia", "Pyrenula pseudobufonia"),
    ("ramalinaFarinacea", "Ramalina farinacea"),
    ("ramalinaHierrensis", "Ramalina hierrensis"),
    ("ramalinaPollinaria", "Ramalina pollinaria"),
    ("ramalinaSiliquosa", "Ramalina siliquosa"),
    ("ramalinaStenospora", "Ramalina stenospora"),
    ("roccellaFuciformis", "Roccella fuciformis"),
    ("stereocaulonAzoreum", "Stereocaulon azoreum"),
    ("stereocaulonCorticatulum", "Stereocaulon corticatulum"),
    ("stereocaulonJaponicum", "Stereocaulon japonicum"),
    ("stereocaulonPaschale", "Stereocaulon paschale"),
    ("stereocaulonRamulosum", "Stereocaulon ramulosum"),
    ("stereocaulonSasakii", "Stereocaulon sasakii"),
    ("stereocaulonVesuvianum", "Stereocaulon vesuvianum"),
    ("teloschistesFlavicans", "Teloschistes flavicans"),
    ("umbilicariaCinereorufescens", "Umbilicaria cinereorufescens"),
    ("umbilicariaCrustulosa", "Umbilicaria crustulosa"),
    ("umbilicariaEsculenta", "Umbilicaria esculenta"),
    ("umbilicariaHypococcinea", "Umbilicaria hypococcinea"),
    ("umbilicariaPolyphylla", "Umbilicaria polyphylla"),
    ("umbilicariaProboscidea", "Umbilicaria proboscidea"),
    ("umbilicariaSpodochroa", "Umbilicaria spodochroa"),
    ("usneaArticulata", "Usnea articulata"),
    ("usneaFlorida", "Usnea florida"),
    ("usneaUndulata", "Usnea undulata"),
    ("xanthoparmeliaCamtschadalis", "Xanthoparmelia camtschadalis"),
    ("xanthoparmeliaCompetita", "Xanthoparmelia competita"),
]

# =============================================================================
# SMDB — mixed Hindi + camelCase slugs (all 40 species confirmed)
# Path: /api/smdb/{slug}/
# =============================================================================
SMDB_SLUGS = [
    # Hindi slugs
    ("ratanjot", "Alkanna tinktoria"),
    ("lahsun", "Allium sativum"),
    ("pyaj", "Allium cepa"),
    ("dalchini", "Cinnamomum verum"),
    ("dhaniya", "Coriandrum sativum"),
    ("kesar", "Crocus sativus"),
    ("jeera", "Cuminum cyminum"),
    ("haldi", "Curcuma longa"),
    ("hing", "Ferula assafoetida"),
    ("pudina", "Mentha spicata"),
    ("kalonji", "Nigella sativa"),
    ("ajwain", "Trachyspermum ammi"),
    ("adrak", "Zingiber officinale"),
    ("kalimirch", "Piper nigrum"),
    ("elayachi", "Elettaria cardamomum"),
    ("jaifal", "Myristica fragrans"),
    ("loung", "Syzygium aromaticum"),
    ("tezpatta", "Cinnamomum tamala"),
    ("currypatta", "Murraya koengii"),
    # camelCase Latin slugs
    ("anthriscusCerefolium", "Anthriscus cerefolium"),
    ("apiumGraveolens", "Apium graveolens"),
    ("armoraciaRusticana", "Armoracia rusticana"),
    ("brassicaJuncea", "Brassica juncea"),
    ("capparisSpinosa", "Capparis spinosa"),
    ("capsicumAnnuum", "Capsicum annuum"),
    ("dracocephalumOfficinale", "Dracocephalum officinale"),
    ("foeniculumVulgare", "Foeniculum vulgare"),
    ("illiciumVerum", "Illicium verum"),
    ("laurusNobilis", "Laurus nobilis"),
    ("levisticumOfficinale", "Levisticum officinale"),
    ("mangiferaIndica", "Mangifera indica"),
    ("ocimumBasilicum", "Ocimum basilicum"),
    ("papaverSomniferum", "Papaver somniferum"),
    ("pimpinellaAnisum", "Pimpinella anisum"),
    ("punicaGranatum", "Punica granatum"),
    ("salviaRosmarinus", "Salvia rosmarinus"),
    ("saturejaHortensis", "Satureja hortensis"),
    ("sesamumIndicum", "Sesamum indicum"),
    ("tamarindusIndica", "Tamarindus indica"),
    ("thymusVulgaris", "Thymus vulgaris"),
]

# =============================================================================
# TMDB — /api/trichoderma/{epithet} (29/30 species confirmed)
# =============================================================================
TMDB_EPITHETS = [
    "album", "applanatum", "asperellum", "atroviride", "aureoviride",
    "avellaneum", "brevicompactum", "cerinum", "citrinoviride", "citrinum",
    "cornudamae", "crassum", "cremeum", "deliquescens", "gamsii",
    "hamatum", "harzianum", "hypocreamuroiana", "hypoxylon", "koningii",
    "koningiopsis", "longibrachiatum", "polysporum", "reesei",
    "saturnisporum", "spirale", "velutinum", "virens", "viride", "vinosum",
]

# =============================================================================
# CMDB — /api/{common_name}/ (all 7 confirmed)
# =============================================================================
CMDB_SLUGS = ["wheat", "rice", "sorghum", "ragi", "barley", "maize", "bajra"]

# =============================================================================
# Database configs
# =============================================================================
DB_CONFIG = {
    "lmdb": {"label": "LMDB (Lichen)",      "url_fn": lambda s: f"{API}/lichmdb/{s}/"},
    "smdb": {"label": "SMDB (Spice)",        "url_fn": lambda s: f"{API}/smdb/{s}/"},
    "tmdb": {"label": "TMDB (Trichoderma)",  "url_fn": lambda s: f"{API}/trichoderma/{s}"},
    "cmdb": {"label": "CMDB (Cereals)",      "url_fn": lambda s: f"{API}/{s}/"},
}

def get_slugs(db_key):
    """Return list of (slug, species_label) for a database."""
    if db_key == "lmdb":
        return LMDB_SLUGS
    elif db_key == "smdb":
        return SMDB_SLUGS
    elif db_key == "tmdb":
        return [(e, f"Trichoderma {e}") for e in TMDB_EPITHETS]
    elif db_key == "cmdb":
        return [(s, s.capitalize()) for s in CMDB_SLUGS]


def fetch(url, attempt=0):
    try:
        r = requests.get(url, headers={"User-Agent": UA, "Accept": "application/json"},
                         timeout=60)
        if r.status_code == 404:
            return None
        if r.status_code in (429, 500, 502, 503, 504):
            if attempt < 4:
                time.sleep(15 * (attempt + 1))
                return fetch(url, attempt + 1)
            return None
        r.raise_for_status()
        data = r.json()
        return data if isinstance(data, list) and len(data) > 0 else None
    except (requests.exceptions.RequestException, json.JSONDecodeError, ValueError):
        if attempt < 4:
            time.sleep(15 * (attempt + 1))
            return fetch(url, attempt + 1)
        return None


def scrape_database(db_key, writer):
    cfg = DB_CONFIG[db_key]
    slugs = get_slugs(db_key)
    seen_ids = set()
    total = 0
    active = 0
    print(f"\n{'='*60}")
    print(f"  {cfg['label']}  ({len(slugs)} slugs)")
    print(f"{'='*60}")
    for i, (slug, species) in enumerate(slugs):
        url = cfg["url_fn"](slug)
        print(f"  [{i+1}/{len(slugs)}] {slug}...", end=" ", flush=True)
        data = fetch(url)
        if not data:
            print("empty/404")
            continue
        active += 1
        new = 0
        for c in data:
            cid = c.get("_id", "") or c.get("cid", "")
            if cid in seen_ids:
                continue
            seen_ids.add(cid)
            sp_ret = (c.get("lichen") or c.get("species") or
                      c.get("plant") or c.get("cereal") or "")
            writer.writerow({
                "database": db_key.upper(),
                "slug_used": slug,
                "species": sp_ret or species,
                "sn": c.get("sn", ""),
                "cid": c.get("cid", ""),
                "common_name": c.get("commonName", ""),
                "iupac_name": c.get("IUPACName", ""),
                "molecular_formula": c.get("molecularFormula", ""),
                "molecular_weight": c.get("molecularWeight", ""),
                "canonical_smiles": c.get("canonicalSmiles", ""),
                "reference": c.get("Reference", ""),
                "super_class": c.get("superClass", ""),
                "compound_class": c.get("class", ""),
            })
            new += 1
        total += new
        print(f"{new} new ({len(data)} raw)")
        time.sleep(2)
    print(f"\n  >> {cfg['label']}: {total} unique compounds from {active}/{len(slugs)} slugs")
    return total


def main():
    p = argparse.ArgumentParser(description="Scrape FDP databases (fooddrugs.in)")
    p.add_argument("--db", required=True,
                   choices=list(DB_CONFIG.keys()) + ["all"])
    p.add_argument("-o", "--output", default="fdp_compounds.csv")
    args = p.parse_args()
    dbs = list(DB_CONFIG.keys()) if args.db == "all" else [args.db]
    grand = 0
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for db in dbs:
            grand += scrape_database(db, w)
    print(f"\n{'='*60}")
    print(f"  GRAND TOTAL: {grand} unique compounds -> {args.output}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()