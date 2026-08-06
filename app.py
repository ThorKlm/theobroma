"""
THEOBROMA v3 — An open multi-kingdom natural products database
spanning 27 sources across six kingdoms and six continents.
"""
from flask import (Flask, render_template, request, send_from_directory,
                   jsonify, abort, redirect, url_for, Response)
from config import Config
from config import VERSION_DISPLAY, VERSION_EXTERNAL, VERSION_INTERNAL
import psycopg2, psycopg2.extras, os, math, re, csv, io
import sys
import json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scripts.similarity import SimilarityEngine
from kingdom_thumbnail import kingdom_thumbnail_svg

sim_engine = SimilarityEngine(vectors_dir="data/vectors")
sim_engine.load()

# Warm the ChemBERTa embedding mmap so the first /api/similarity?metric=chemberta
# call does not pay the 30-60s cold-mmap cost. Reading a few scattered rows
# forces the OS to populate the page cache with the index pages.
try:
    import numpy as _np_warm
    _emb_path = "data/vectors/chemberta_embeddings.npy"
    if os.path.exists(_emb_path):
        _warm = _np_warm.load(_emb_path, mmap_mode="r")
        _ = _warm[0].copy(); _ = _warm[len(_warm)//2].copy(); _ = _warm[-1].copy()
        del _warm
except Exception as _e:
    print(f"[startup] WARN: could not warm chemberta mmap: {_e}")

app = Flask(__name__)

@app.context_processor
def inject_version():
    return {"theobroma_version": VERSION_DISPLAY,
            "theobroma_version_external": VERSION_EXTERNAL,
            "theobroma_version_internal": VERSION_INTERNAL,
            # Maintenance banner: on when env THEOBROMA_MAINTENANCE is set to
            # a truthy value ("1"/"true"/"on"). Flip off by unsetting it.
            "maintenance_notice": os.environ.get("THEOBROMA_MAINTENANCE", "").lower() in ("1", "true", "on", "yes")}
app.config.from_object(Config)

# Module-level cache of browse-page dropdown options. Loaded once at import
# time; refreshed only by service restart, which is the v33 corpus-update
# convention. Removes 3 redundant queries per /browse request (the dominant
# component of /browse latency in the v33 benchmark).
_BROWSE_CACHE = {"all_kingdoms": None, "all_sources_list": None, "all_regions": None}

def _load_browse_options():
    """Populate _BROWSE_CACHE. Called once at import time after sim_engine.load()."""
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""SELECT kingdom, COUNT(*) AS cnt FROM resolved_taxonomy
                           WHERE kingdom IS NOT NULL AND kingdom!=''
                           GROUP BY kingdom ORDER BY kingdom""")
            _BROWSE_CACHE["all_kingdoms"] = cur.fetchall()
            cur.execute("SELECT source_db, COUNT(*) AS cnt FROM compounds "
                        "GROUP BY source_db ORDER BY source_db")
            _BROWSE_CACHE["all_sources_list"] = cur.fetchall()
            cur.execute("SELECT DISTINCT macro_region AS reg FROM compound_region_map UNION SELECT 'global / unresolved' ORDER BY reg")
            _BROWSE_CACHE["all_regions"] = [r["reg"] for r in cur.fetchall()]

SORTABLE = {"comp_id","name","kingdom","source_db","region","source_organism",
            "mw","logp","tpsa","hba","hbd","n_rings","rotatable_bonds","license_tier"}
ALLOWED_PER_PAGE = {25, 50, 100, 1000}

def get_db():
    return psycopg2.connect(app.config["DB_URI"])

def _resolve_comp_id_for_nafm(raw_query, smiles=None):
    """Resolve a similarity query to a corpus comp_id for NaFM lookup.
    Order: explicit THEO_ id -> exact name match -> InChIKey of the query SMILES.
    Returns comp_id str or None."""
    if not raw_query and not smiles:
        return None
    q = (raw_query or "").strip()
    try:
        with get_db() as conn:
            with conn.cursor() as cur:
                if q.upper().startswith("THEO_"):
                    cur.execute("SELECT comp_id FROM compounds WHERE comp_id=%s LIMIT 1", (q,))
                    r = cur.fetchone()
                    if r: return r[0]
                if q:
                    cur.execute("SELECT comp_id FROM compounds WHERE LOWER(name)=%s LIMIT 1", (q.lower(),))
                    r = cur.fetchone()
                    if r: return r[0]
                if smiles:
                    try:
                        from rdkit import Chem
                        m = Chem.MolFromSmiles(smiles)
                        if m is not None:
                            ik = Chem.MolToInchiKey(m)
                            if ik:
                                cur.execute("SELECT comp_id FROM compounds WHERE inchikey=%s LIMIT 1", (ik,))
                                r = cur.fetchone()
                                if r: return r[0]
                    except Exception:
                        pass
    except Exception:
        return None
    return None

def resolve_taxon(rank, raw_query):
    """Normalize the query, find exact match in resolved_taxonomy.<col>, or fall back
    to the single closest distinct value by trigram similarity (>= 0.3).
    Returns the canonical taxon name to search for, or None if no acceptable match."""
    col_map = {"genus": "genus", "family": "family", "order": "taxorder",
               "tax_class": "taxclass", "clade": "taxclass", "phylum": "phylum", "kingdom": "kingdom"}
    col = col_map.get(rank)
    if not col or not raw_query:
        return raw_query
    import re as _re
    norm = _re.sub(r"[^a-zA-Z0-9]+", " ", raw_query).strip().lower()
    if not norm:
        return raw_query
    with get_db() as conn:
        with conn.cursor() as cur:
            # Exact match first (case-insensitive after normalization)
            cur.execute(
                f"SELECT DISTINCT {col} FROM resolved_taxonomy "
                f"WHERE LOWER(REGEXP_REPLACE({col}, '[^a-zA-Z0-9]+', ' ', 'g')) = %s LIMIT 1",
                (norm,)
            )
            row = cur.fetchone()
            if row:
                return row[0]
            # Fuzzy fallback via trigram similarity
            cur.execute(
                f"SELECT {col}, similarity(LOWER({col}), %s) AS sim "
                f"FROM resolved_taxonomy WHERE {col} IS NOT NULL "
                f"ORDER BY sim DESC LIMIT 1",
                (norm,)
            )
            row = cur.fetchone()
            if row and row[1] >= 0.3:
                return row[0]
    return None  # no acceptable match


def get_sort(default="comp_id"):
    s = request.args.get("sort", default)
    o = request.args.get("order", "asc").lower()
    if s not in SORTABLE: s = default
    if o not in ("asc","desc"): o = "asc"
    return s, o, f"ORDER BY {s} {o} NULLS LAST"

def get_per_page():
    try:
        pp = int(request.args.get("per_page", Config.PER_PAGE))
    except (ValueError, TypeError):
        pp = Config.PER_PAGE
    return pp if pp in ALLOWED_PER_PAGE else Config.PER_PAGE

def paginate(query, params, page, per_page, conn):
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        cur.execute(f"SELECT COUNT(*) FROM ({query}) AS subq", params)
        total = cur.fetchone()["count"]
        pages = max(1, math.ceil(total / per_page))
        cur.execute(query + " LIMIT %s OFFSET %s", params + (per_page, (page-1)*per_page))
        results = cur.fetchall()
    return results, total, pages

REGION_SQL = """CASE WHEN region IS NULL OR region='' OR region='global' THEN 'global / unresolved' ELSE region END"""

import re as _re
# Source-database identifier prefixes: names starting with these are NOT human
# chemical names and must never be the primary display name if any real name exists.
_ID_PREFIX = _re.compile(
    r'^(chebi:?|npc|npo|npa|sa_|mol_?|cid[0-9]|hmdb|zinc|pubchem|c\.i\.|e\s?\d|unii|ec\s?\d'
    r'|orb|kbiogr|kbioss|kbio|refchem|ncimech|schembl|mls\d|mcule|chembl|ambn|stk|bdbm'
    r'|nsc\d|cas[-_ ]?\d|acon|dtxsid|dtxcid|gtpl|lmfa|lmgp|lmpr|lmsp|lmst|q\d{4,})',
    _re.I)
_CAS = _re.compile(r'^\d{1,7}-\d{2}-\d$')                 # CAS registry number
_ALLCODE = _re.compile(r'^[^a-z]*$')                       # no lowercase letter -> code/acronym
_SYSTEMATIC = _re.compile(r'[0-9]+[,\-]|[\[\]\(\)]|(\b\d[EZRS]\b)|hepta|hexa|penta|tetra|dione|yl\b')

def is_id_like(s):
    """True if s is a source-database identifier / code rather than a human name."""
    if not s:
        return True
    t = s.strip()
    return bool(_ID_PREFIX.match(t) or _CAS.match(t) or _ALLCODE.match(t)
                or t in ("nan", "") or t.startswith("NPO"))

# Locant / stereo / structural descriptor prefixes whose case is meaningful and must
# NOT be capitalized (n-, o-, p-, m-, cis-, trans-, tert-, sec-, d-, l-, e-, z-, etc.).
_DESCRIPTOR_PREFIX = _re.compile(
    r'^(n|o|p|m|d|l|e|z|r|s|h|cis|trans|tert|sec|syn|anti|alpha|beta|gamma|delta|omega|ortho|meta|para)-',
    _re.I)

def display_name(s):
    """Presentation-only capitalization. Conservative: capitalize the first letter ONLY
    when the name is a plain trivial name (pure lowercase letters, no digits/punctuation,
    not a locant/descriptor prefix). Systematic names, descriptor-led names, and anything
    with internal case meaning are returned unchanged. Never lowercases anything."""
    if not s:
        return s
    t = str(s)
    # already has an uppercase letter, or does not start with a lowercase letter -> leave as-is
    if not t[:1].islower():
        return t
    # descriptor/locant-led (p-cymene, N-methyl, cis-..., L-dopa) -> leave as-is
    if _DESCRIPTOR_PREFIX.match(t):
        return t
    # pure lowercase letters only (curcumin, quercetin, geraniol) -> capitalize first letter
    if t.isalpha() and t.islower():
        return t[0].upper() + t[1:]
    # starts with a lowercase letter and the first token is an alphabetic word
    # (e.g. "curcumin dimer", "ginkgolide a"): capitalize first letter only.
    first_tok = t.split()[0] if t.split() else ""
    if first_tok.isalpha() and t[:1].islower():
        return t[0].upper() + t[1:]
    # otherwise (leading digit, punctuation, locant): leave unchanged
    return t

def synonym_tier(s):
    """Lower tier = better display name. 0 chemical common name, 1 systematic/IUPAC,
    2 vernacular/trade, 3 registry code/identifier. Used to order synonyms and pick a name."""
    if not s:
        return 4
    t = s.strip()
    if is_id_like(t):
        return 3
    if _SYSTEMATIC.search(t.lower()):
        return 1
    # short, single- or few-word, mostly-alpha names: treat as common/vernacular (tier 2),
    # promoted to 0 only when it is a plausible chemical common name (single token or ends in
    # a chemical suffix). This keeps "curcumin" above "Kacha haldi" without a curated lexicon.
    words = t.split()
    chem_suffix = _re.search(r'(ine|ol|one|ate|ide|acid|in|ene|ane|oside|genin)$', t.lower())
    if len(words) <= 1 and chem_suffix:
        return 0
    if chem_suffix and len(words) <= 2:
        return 0
    return 2

def order_synonyms(syns):
    """Stable-order a synonym list best-first by tier, then by length within tier."""
    seen, uniq = set(), []
    for s in syns:
        k = (s or "").strip().lower()
        if k and k not in seen:
            seen.add(k); uniq.append(s)
    return sorted(uniq, key=lambda s: (synonym_tier(s), len(s)))

def canonical_display_name(chebi_name, name, synonyms, chebi_iupac_name, query_norm=None):
    """Best display name: curated chebi_name > query-matched synonym > best-tier synonym
    > systematic IUPAC > existing name > empty. Returns '' if nothing usable."""
    if chebi_name and chebi_name.strip():
        return chebi_name.strip()
    nm = str(name or "")
    # Treat an ID-style primary name (Npc18984, SCHEMBL..., CHEBI:..., Mol_...) as a
    # placeholder so a real synonym or the systematic name is preferred over it.
    is_placeholder = (not nm or is_id_like(nm))
    ordered = order_synonyms(synonyms or [])
    if query_norm:
        for s in ordered:
            if normalize_query(s) == query_norm:
                return s
    if not is_placeholder:
        return nm
    for s in ordered:
        if synonym_tier(s) <= 2:
            return s
    if chebi_iupac_name and chebi_iupac_name.strip():
        return chebi_iupac_name.split(";")[0].strip()
    # Last resort: a systematic-looking synonym (tier 1) beats an ID or a dash.
    for s in ordered:
        if synonym_tier(s) <= 3 and not is_id_like(s):
            return s
    return nm if nm else (ordered[0] if ordered else "")

@app.template_filter("normalize_name")
def normalize_name(name):
    """Title-case all-uppercase compound names from source data; preserve intentional mixed case."""
    if not name:
        return name
    s_str = str(name)
    if any(c.islower() for c in s_str):
        return s_str
    return s_str[0].upper() + s_str[1:].lower()

@app.template_filter("region_label")
def region_label(v):
    return "global / unresolved" if not v or v in ("global","nan","") else v

@app.template_filter("kingdom_label")
def kingdom_label(v):
    return "unresolved" if not v or v in ("nan","") else v


# User-agents whose crawling traffic dominated v32 server time with no user
# benefit. GPTBot alone hit /admet 856k times in two months. The robots.txt
# Disallow is voluntary; this enforcement returns 429 so the bot eventually
# backs off. Substring match so version variants (GPTBot/1.3 etc.) are caught.
_BLOCKED_BOTS_HEAVY_ROUTES = ("GPTBot", "ClaudeBot", "anthropic-ai", "CCBot",
                              "Bytespider", "PerplexityBot")
_BOT_BLOCKED_PREFIXES = ("/admet", "/api/", "/export", "/similarity",
                         "/substructure", "/scaffolds")

@app.before_request
def block_heavy_bots():
    """Return 429 Too Many Requests for known scraper bots on heavy routes."""
    ua = request.headers.get("User-Agent", "")
    if any(bot in ua for bot in _BLOCKED_BOTS_HEAVY_ROUTES):
        if any(request.path.startswith(p) for p in _BOT_BLOCKED_PREFIXES):
            return Response("Bot traffic on heavy routes is rate-limited.\n",
                            status=429, mimetype="text/plain",
                            headers={"Retry-After": "86400"})


@app.after_request
def log_access(response):
    if request.path.startswith("/static/"):
        return response
    try:
        ip = request.headers.get("X-Forwarded-For", request.remote_addr)
        with get_db() as conn:
            with conn.cursor() as cur:
                cur.execute(
                    "INSERT INTO access_log (path, method, ip, user_agent) VALUES (%s,%s,%s,%s)",
                    (request.path, request.method, ip, request.headers.get("User-Agent","")[:200])
                )
            conn.commit()
    except:
        pass
    return response


def normalize_query(q):
    """Normalize search query to match the stored name_norm rule, which strips
    spaces and hyphens but preserves other punctuation (e.g. commas):
    'Withaferin A' -> 'withaferina', '2,3-Dihydro X' -> '2,3dihydrox'."""
    q = q.lower().strip()
    # Greek letter equivalents
    replacements = [
        ("α", "alpha"), ("β", "beta"), ("γ", "gamma"), ("δ", "delta"),
        ("ε", "epsilon"), ("ω", "omega"), ("μ", "mu"),
    ]
    for greek, latin in replacements:
        q = q.replace(greek, latin)
    # Remove hyphens AND spaces to match name_norm (commas etc. are preserved)
    q = q.replace("-", "").replace("–", "").replace("—", "")
    q = q.replace(" ", "")
    return q

_load_browse_options()

# --- Routes ---

@app.route("/")
def index():
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT COUNT(*) AS cnt FROM compounds")
            total = cur.fetchone()["cnt"]
            cur.execute("""SELECT kingdom, COUNT(*) AS cnt FROM resolved_taxonomy
                          WHERE kingdom IS NOT NULL AND kingdom!='' GROUP BY kingdom ORDER BY cnt DESC""")
            kingdoms = cur.fetchall()
            cur.execute("SELECT COUNT(DISTINCT source_db) AS cnt FROM compounds")
            n_sources = cur.fetchone()["cnt"]
            cur.execute("""SELECT reg, cnt FROM (
                   SELECT macro_region AS reg, COUNT(DISTINCT comp_id) AS cnt FROM compound_region_map GROUP BY 1
                   UNION ALL
                   SELECT 'global / unresolved' AS reg,
                          (SELECT COUNT(*) FROM compounds c WHERE NOT EXISTS
                             (SELECT 1 FROM compound_region_map m WHERE m.comp_id=c.comp_id)) AS cnt
                 ) t ORDER BY cnt DESC""")
            regions = cur.fetchall()
            n_regions = len([r for r in regions if r["reg"]!="unresolved"])
    home_thumb = kingdom_thumbnail_svg(kingdoms, total=total, size=180, title="Corpus kingdoms")
    home_region_counts = []
    try:
        import json as _json
        with open("/home/thorben.klamt/theobroma/static/compounds_by_country.json") as _f:
            _data = _json.load(_f)
            home_region_counts = _data.get("region_counts", [])
    except Exception:
        home_region_counts = []
    home_region_css, home_region_titles = build_region_color_css(home_region_counts, mode="light")
    linear_tree = compute_linear_tree()
    return render_template("index.html", total=total, kingdoms=kingdoms,
                           n_sources=n_sources, n_regions=n_regions, regions=regions,
                           home_thumb=home_thumb,
                           home_region_css=home_region_css,
                           home_region_counts=home_region_counts,
                           home_region_titles=home_region_titles,
                           linear_tree=linear_tree)

def _exact_flag():
    """Return True if the request asks for exact-match semantics."""
    return request.args.get("exact", "").lower() in ("true", "1", "yes")

@app.route("/search")
def search():
    q = request.args.get("q","").strip()
    exact = _exact_flag()
    st_preview = request.args.get("type","name")
    if st_preview == "smiles" and q:
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(q)
            if mol:
                ik = Chem.MolToInchiKey(mol)
                return redirect(url_for("search", q=ik, type="inchikey"))
        except:
            pass
    st = request.args.get("type","name")
    page = max(1, int(request.args.get("page",1)))
    per_page = get_per_page()
    sort, order, oc = get_sort()
    # Resolve "property" type to actual search type.
    # Accepts common aliases case-insensitively; unknown values fall back to
    # npclassifier_class rather than silently returning the unfiltered corpus.
    if st in ("property", "classification"):
        raw_prop_type = request.args.get("prop_type", "").strip().lower()
        prop_type_aliases = {
            "class": "npclassifier_class",
            "chemical_class": "npclassifier_class",
            "chem_class": "npclassifier_class",
            "npclassifier_class": "npclassifier_class",
            "npc_class": "npclassifier_class",
            "np_class": "npclassifier_class",
            "superclass": "npclassifier_superclass",
            "np_superclass": "npclassifier_superclass",
            "npclassifier_superclass": "npclassifier_superclass",
            "classyfire_class": "classyfire_class",
            "cf_class": "classyfire_class",
            "classyfire_superclass": "classyfire_class",
            "pathway": "pathway",
            "np_pathway": "pathway",
            "genus": "genus",
            "family": "family",
            "order": "order",
            "tax_class": "tax_class",
            "phylum": "phylum",
            "mw": "mw",
            "logp": "mw",
            "tpsa": "mw",
            "hba": "mw",
            "hbd": "mw",
            "n_rings": "mw",
            "rotatable_bonds": "mw",
        }
        st = prop_type_aliases.get(raw_prop_type, "npclassifier_class")
    has_extra = any(request.args.get(f"extra_type_{i}") and request.args.get(f"extra_q_{i}") for i in range(1, 11))
    has_range = any(request.args.get(f"{p}_min") or request.args.get(f"{p}_max")
                    for p in ["mw","logp","tpsa","hba","hbd","n_rings","rotatable_bonds",
                              "Lipinski","QED","stereo_centers","PAINS_alert","BRENK_alert","NIH_alert",
                              "AMES","BBB_Martins","Bioavailability_Ma",
                              "CYP1A2_Veith","CYP2C19_Veith","CYP2C9_Substrate_CarbonMangels","CYP2C9_Veith",
                              "CYP2D6_Substrate_CarbonMangels","CYP2D6_Veith","CYP3A4_Substrate_CarbonMangels","CYP3A4_Veith",
                              "Carcinogens_Lagunin","ClinTox","DILI","HIA_Hou",
                              "NR_AR_LBD","NR_AR","NR_AhR","NR_Aromatase","NR_ER_LBD","NR_ER","NR_PPAR_gamma",
                              "PAMPA_NCATS","Pgp_Broccatelli",
                              "SR_ARE","SR_ATAD5","SR_HSE","SR_MMP","SR_p53",
                              "Skin_Reaction","hERG","Caco2_Wang",
                              "Clearance_Hepatocyte_AZ","Clearance_Microsome_AZ","Half_Life_Obach",
                              "HydrationFreeEnergy_FreeSolv","LD50_Zhu","Lipophilicity_AstraZeneca","PPBR_AZ",
                              "Solubility_AqSolDB","VDss_Lombardo"])
    try:
        linear_tree = compute_linear_tree()
    except Exception:
        linear_tree = []
    if not q and not has_extra and not has_range and st not in ("mw",):
        # Whole-corpus widgets for the bare landing page (world map + circular tree
        # + linear cladogram), mirroring the home route and using cached data so no
        # live 1.13M aggregation happens on a query-less hit.
        landing_thumb = None
        landing_region_counts = []
        landing_region_css = ""
        landing_region_titles = {}
        try:
            with get_db() as conn_lp:
                with conn_lp.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur_lp:
                    cur_lp.execute("SELECT kingdom, COUNT(*) AS cnt FROM compounds WHERE kingdom IS NOT NULL AND kingdom!='' GROUP BY kingdom ORDER BY cnt DESC")
                    lp_kingdoms = cur_lp.fetchall()
                    cur_lp.execute("SELECT COUNT(*) AS cnt FROM compounds")
                    lp_total = cur_lp.fetchone()["cnt"]
            landing_thumb = kingdom_thumbnail_svg(lp_kingdoms, total=lp_total, size=150, title="Corpus kingdoms")
        except Exception:
            landing_thumb = None
        try:
            import json as _json_lp
            with open("/home/thorben.klamt/theobroma/static/compounds_by_country.json") as _f_lp:
                landing_region_counts = _json_lp.load(_f_lp).get("region_counts", [])
            landing_region_css, landing_region_titles = build_region_color_css(landing_region_counts)
        except Exception:
            landing_region_counts = []
            landing_region_css = ""
            landing_region_titles = {}
        return render_template("search.html", results=[], query="", search_type=st,
                               searched=False,
                               page=1, total=0, pages=0, sort=sort, order=order, per_page=per_page,
                               linear_tree=linear_tree,
                               thumb=landing_thumb,
                               region_counts_filtered=landing_region_counts,
                               region_css=landing_region_css,
                               region_titles=landing_region_titles)
    if not q and (has_extra or has_range):
        st = "_base"
    # THEO_id direct redirect
    if re.match(r"^THEO_\d{7}$", q):
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute("SELECT comp_id FROM compounds WHERE comp_id=%s LIMIT 1", (q,))
                if cur.fetchone():
                    return redirect(url_for("compound_detail", comp_id=q))

    # InChIKey direct redirect
    if st == "inchikey" and re.match(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$', q):
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute("SELECT comp_id FROM compounds WHERE inchikey=%s LIMIT 1", (q,))
                row = cur.fetchone()
                if row:
                    return redirect(url_for("compound_detail", comp_id=row["comp_id"]))
    # Fuzzy taxonomic resolution: must run BEFORE tq dict is built since params are baked in.
    if st in ("genus", "family", "order", "tax_class", "clade", "phylum"):
        resolved = resolve_taxon(st, q)
        if resolved is not None and resolved != q:
            q = resolved
        elif resolved is None:
            q = "__no_such_taxon__"
    tq = {
        "name":    ((f"""SELECT c.*, 0 AS rank_key FROM (
               SELECT DISTINCT ON (sn.comp_id) sn.comp_id, 0 AS relevance
               FROM search_names sn
               WHERE sn.name_norm = %s
             ) matched
             JOIN compounds c ON c.comp_id = matched.comp_id
             ORDER BY c.name""", (normalize_query(q),)) if exact else (f"""SELECT c.*, matched.rank_key,
               CASE WHEN (c.chebi_name IS NOT NULL AND c.chebi_name <> '')
                      OR (c.name IS NOT NULL AND c.name <> ''
                          AND c.name !~* '^(mol_|npo|sa_|npc|npo|schembl|orb|kbio|refchem|ncimech|mls[0-9]|cid[0-9]|chembl|zinc|hmdb|nan)')
                    THEN 0 ELSE 1 END AS has_name
             FROM (
               SELECT DISTINCT ON (sn.comp_id) sn.comp_id,
               levenshtein(sn.name_norm, %s) AS rank_key
               FROM search_names sn
               WHERE sn.name_norm LIKE %s
               ORDER BY sn.comp_id, levenshtein(sn.name_norm, %s)
             ) matched
             JOIN compounds c ON c.comp_id = matched.comp_id
             ORDER BY has_name, matched.rank_key, LENGTH(c.name), c.name""", (normalize_query(q), '%'+normalize_query(q)+'%', normalize_query(q)))),
        "smiles":  (f"""SELECT * FROM compounds WHERE inchikey = (
            SELECT inchikey FROM compounds WHERE smiles=%s LIMIT 1
          ) {oc}""", (q,)),
        "inchikey":(f"SELECT * FROM compounds WHERE inchikey=%s {oc}", (q,)),
        "source":  (f"SELECT * FROM compounds WHERE LOWER(source_db) = LOWER(%s) {oc}", (q,)),
        "organism":((f"SELECT * FROM compounds WHERE LOWER(source_organism) = LOWER(%s) {oc if request.args.get('sort') else ''}", (q,)) if exact else (f"SELECT * FROM compounds WHERE source_organism ILIKE %s {oc}", (f"%{q}%",))),
        "region":  (f"SELECT * FROM compounds WHERE EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = compounds.comp_id AND LOWER(crm.macro_region) = LOWER(%s)) {oc}", (q,)),
        "kingdom": (f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE (LOWER(rt.kingdom) = LOWER(%s) OR LOWER(%s) = ANY(rt.secondary_kingdoms)) {oc}", (q, q)),
        "genus":   ((f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE LOWER(rt.genus) = LOWER(%s) {oc}", (q,)) if exact else (f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE rt.genus ILIKE %s {oc}", (f"%{q}%",))),
        "family":  ((f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE LOWER(rt.family) = LOWER(%s) {oc}", (q,)) if exact else (f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE rt.family ILIKE %s {oc}", (f"%{q}%",))),
        "order":   ((f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE LOWER(rt.taxorder) = LOWER(%s) {oc}", (q,)) if exact else (f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE rt.taxorder ILIKE %s {oc}", (f"%{q}%",))),
        "tax_class":   ((f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE LOWER(rt.taxclass) = LOWER(%s) {oc}", (q,)) if exact else (f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE rt.taxclass ILIKE %s {oc}", (f"%{q}%",))),
        "clade":   ((f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE LOWER(rt.taxclass) = LOWER(%s) {oc}", (q,)) if exact else (f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE rt.taxclass ILIKE %s {oc}", (f"%{q}%",))),
        "phylum":  ((f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE LOWER(rt.phylum) = LOWER(%s) {oc}", (q,)) if exact else (f"SELECT c.* FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id WHERE rt.phylum ILIKE %s {oc}", (f"%{q}%",))),
        "class":   (f"SELECT * FROM compounds WHERE np_class ILIKE %s OR classyfire_superclass ILIKE %s OR inferred_class ILIKE %s OR np_superclass ILIKE %s OR np_pathway ILIKE %s {oc}", (f"%{q}%", f"%{q}%", f"%{q}%", f"%{q}%", f"%{q}%")),
        "npclassifier_class": (f"SELECT * FROM compounds WHERE np_class ILIKE %s OR inferred_class ILIKE %s {oc}", (f"%{q}%", f"%{q}%")),
        "npclassifier_superclass": (f"SELECT * FROM compounds WHERE np_superclass ILIKE %s {oc}", (f"%{q}%",)),
        "classyfire_class":   (f"SELECT * FROM compounds WHERE classyfire_superclass ILIKE %s {oc}", (f"%{q}%",)),
        "pathway": (f"SELECT * FROM compounds WHERE np_pathway ILIKE %s {oc}", (f"%{q}%",)),
        "mw":      (f"SELECT * FROM compounds WHERE 1=1 {oc}", ()),
        "_base":   (f"SELECT * FROM compounds WHERE 1=1 {oc}", ()),
    }
    query, params = tq.get(st, tq["name"])
    # Handle extra AND filters
    extra_clauses = []
    extra_params = []
    for i in range(1, 11):
        et = request.args.get(f"extra_type_{i}", "")
        eq = request.args.get(f"extra_q_{i}", "").strip()
        if not et or not eq:
            continue
        emap = {
            "name": "LOWER(name) = LOWER(%s)" if exact else "LOWER(name) LIKE %s",
            "organism": "LOWER(source_organism) = LOWER(%s)" if exact else "source_organism ILIKE %s",
            "kingdom": "EXISTS(SELECT 1 FROM resolved_taxonomy rt2 WHERE rt2.comp_id = base.comp_id AND (LOWER(rt2.kingdom) = LOWER(%s) OR LOWER(%s) = ANY(rt2.secondary_kingdoms)))",
            "region": "EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = base.comp_id AND LOWER(crm.macro_region) = LOWER(%s))",
            "source": "LOWER(source_db) = LOWER(%s)",
            "class": "(np_class ILIKE %s OR classyfire_superclass ILIKE %s)",
            "npclassifier_class": "(np_class ILIKE %s OR inferred_class ILIKE %s)",
            "classyfire_class": "classyfire_superclass ILIKE %s",
            "pathway": "np_pathway ILIKE %s",
            "genus":  "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = base.comp_id AND LOWER(rt.genus) = LOWER(%s))",
            "family": "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = base.comp_id AND LOWER(rt.family) = LOWER(%s))",
            "order":  "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = base.comp_id AND LOWER(rt.taxorder) = LOWER(%s))",
            "tax_class": "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = base.comp_id AND LOWER(rt.taxclass) = LOWER(%s))",
            "clade": "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = base.comp_id AND LOWER(rt.taxclass) = LOWER(%s))",
            "phylum": "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = base.comp_id AND LOWER(rt.phylum) = LOWER(%s))",
        }
        # Fuzzy resolution for taxonomic ranks in extras as well.
        if et in ("genus", "family", "order", "tax_class", "clade", "phylum"):
            resolved = resolve_taxon(et, eq)
            if resolved is not None and resolved != eq:
                eq = resolved
            elif resolved is None:
                # Force this extra to match nothing.
                eq = "__no_such_taxon__"
        esql = emap.get(et)
        if esql:
            if et in ("class", "npclassifier_class"):
                extra_clauses.append(esql)
                extra_params.extend([f"%{eq}%", f"%{eq}%"])
            elif et == "classyfire_class":
                extra_clauses.append(esql)
                extra_params.append(f"%{eq}%")
            elif exact and et in ("name", "organism"):
                extra_clauses.append(esql)
                extra_params.append(eq)
            elif et == "kingdom":
                # Kingdom uses primary OR secondary check; needs the value twice.
                extra_clauses.append(esql)
                extra_params.extend([eq, eq])
            elif et in ("region", "source"):
                extra_clauses.append(esql)
                extra_params.append(eq)
            elif et in ("genus", "family"):
                extra_clauses.append(esql)
                extra_params.append(eq)
            elif et in ("order", "tax_class", "clade", "phylum"):
                extra_clauses.append(esql)
                extra_params.append(eq)
            else:
                extra_clauses.append(esql)
                extra_params.append(f"%{eq}%")
    # Handle property range filters (physchem + ADMET)
    admet_cols = {"Lipinski","QED","stereo_centers","PAINS_alert","BRENK_alert","NIH_alert","AMES","BBB_Martins","Bioavailability_Ma","CYP1A2_Veith","CYP2C19_Veith","CYP2C9_Substrate_CarbonMangels","CYP2C9_Veith","CYP2D6_Substrate_CarbonMangels","CYP2D6_Veith","CYP3A4_Substrate_CarbonMangels","CYP3A4_Veith","Carcinogens_Lagunin","ClinTox","DILI","HIA_Hou","NR_AR_LBD","NR_AR","NR_AhR","NR_Aromatase","NR_ER_LBD","NR_ER","NR_PPAR_gamma","PAMPA_NCATS","Pgp_Broccatelli","SR_ARE","SR_ATAD5","SR_HSE","SR_MMP","SR_p53","Skin_Reaction","hERG","Caco2_Wang","Clearance_Hepatocyte_AZ","Clearance_Microsome_AZ","Half_Life_Obach","HydrationFreeEnergy_FreeSolv","LD50_Zhu","Lipophilicity_AstraZeneca","PPBR_AZ","Solubility_AqSolDB","VDss_Lombardo"}
    needs_admet_join = False
    for prop in ["mw", "logp", "tpsa", "hba", "hbd", "n_rings", "rotatable_bonds",
                  "Lipinski", "QED", "stereo_centers", "PAINS_alert", "BRENK_alert", "NIH_alert", "AMES", "BBB_Martins", "Bioavailability_Ma", "CYP1A2_Veith", "CYP2C19_Veith", "CYP2C9_Substrate_CarbonMangels", "CYP2C9_Veith", "CYP2D6_Substrate_CarbonMangels", "CYP2D6_Veith", "CYP3A4_Substrate_CarbonMangels", "CYP3A4_Veith", "Carcinogens_Lagunin", "ClinTox", "DILI", "HIA_Hou", "NR_AR_LBD", "NR_AR", "NR_AhR", "NR_Aromatase", "NR_ER_LBD", "NR_ER", "NR_PPAR_gamma", "PAMPA_NCATS", "Pgp_Broccatelli", "SR_ARE", "SR_ATAD5", "SR_HSE", "SR_MMP", "SR_p53", "Skin_Reaction", "hERG", "Caco2_Wang", "Clearance_Hepatocyte_AZ", "Clearance_Microsome_AZ", "Half_Life_Obach", "HydrationFreeEnergy_FreeSolv", "LD50_Zhu", "Lipophilicity_AstraZeneca", "PPBR_AZ", "Solubility_AqSolDB", "VDss_Lombardo"]:
        pmin = request.args.get(f"{prop}_min", "")
        pmax = request.args.get(f"{prop}_max", "")
        if pmin or pmax:
            if prop in admet_cols:
                needs_admet_join = True
                col = f'admet."{prop}"'
            else:
                col = prop
            if pmin:
                extra_clauses.append(f"{col} >= %s")
                extra_params.append(float(pmin))
            if pmax:
                extra_clauses.append(f"{col} <= %s")
                extra_params.append(float(pmax))
    # License-tier filter applied alongside typed-search extras. Mirrors the
    # basic-browse license logic at the /search route's lower branch (line ~764).
    license_filter = request.args.get("license", "all")
    if license_filter == "commercial":
        extra_clauses.append("tier_rank <= 1")
    elif license_filter == "academic":
        extra_clauses.append("tier_rank <= 4")
    # Named-only toggle: filter to compounds with non-empty name
    named_only = request.args.get("named", "")
    if named_only:
        extra_clauses.append("name IS NOT NULL AND name != ''")
    if extra_clauses:
        extra_where = " AND ".join(extra_clauses)
        rank_order = " ORDER BY base.has_name, base.rank_key, LENGTH(base.name), base.name" if st == "name" and not exact else ""
        if needs_admet_join:
            query = f"SELECT base.* FROM ({query}) AS base JOIN admet ON base.comp_id=admet.comp_id WHERE {extra_where}{rank_order}"
        else:
            query = f"SELECT * FROM ({query}) AS base WHERE {extra_where}{rank_order}"
        params = params + tuple(extra_params)
    with get_db() as conn:
        results, total, pages = paginate(query, params, page, per_page, conn)
    if st == "name" and results:
        seen = {}
        deduped = []
        for r in results:
            ik14 = (r.get("inchikey") or "")[:14]
            name_key = (r.get("name") or "").lower().strip()
            dedup_key = (ik14, name_key) if ik14 else (name_key,)
            if dedup_key in seen:
                continue
            seen[dedup_key] = True
            nm = r.get("name") or ""
            if nm and nm[0].islower():
                r["name"] = nm[0].upper() + nm[1:]
            deduped.append(r)
        results = deduped
    # Display name for every row: canonical helper prefers curated chebi_name, then a
    # synonym matching the query, then best-tier synonym, then systematic IUPAC.
    q_norm = normalize_query(q) if st == "name" else None
    # Rows that may need a synonym lookup: no chebi_name and (placeholder OR name-search).
    def _is_placeholder(nm):
        return is_id_like(str(nm or ""))
    need_syn = [r for r in results
                if not (r.get("chebi_name") or "").strip()
                and (_is_placeholder(r.get("name")) or st == "name")]
    syn_map = {}
    iks = [r.get("inchikey") for r in need_syn if r.get("inchikey")]
    if iks:
        with get_db() as conn2:
            with conn2.cursor() as cur2:
                ph = ",".join(["%s"]*len(iks))
                cur2.execute(f"SELECT inchikey, synonym FROM compound_synonyms WHERE inchikey IN ({ph})", tuple(iks))
                for ik, syn in cur2.fetchall():
                    syn_map.setdefault(ik, []).append(syn)
    for r in results:
        disp = canonical_display_name(r.get("chebi_name"), r.get("name"),
                                      syn_map.get(r.get("inchikey"), []),
                                      r.get("chebi_iupac_name"), q_norm)
        if disp:
            r["name"] = display_name(disp)
    # Kingdom donut for the full search result set.
    # Special case A: explicit kingdom filter (type=kingdom OR extra_*=kingdom).
    # Special case B: taxonomic filter (genus/family/order/tax_class/phylum) where
    # the result set has a clear majority kingdom (>=95%). In both cases show 100%
    # of that kingdom to reflect that the query implies a single biological kingdom.
    thumb = None
    if total > 0:
        try:
            extras_kingdom_q = None
            for i in range(1, 11):
                et = request.args.get(f"extra_type_{i}", "")
                eq_v = request.args.get(f"extra_q_{i}", "").strip()
                if et == "kingdom" and eq_v:
                    extras_kingdom_q = eq_v
                    break
            kingdom_label_override = None
            if st == "kingdom" and q:
                kingdom_label_override = q.lower()
            elif extras_kingdom_q:
                kingdom_label_override = extras_kingdom_q.lower()
            # Primary-kingdom breakdown SQL (used for inspection + fallback).
            # When a taxonomic filter is active (genus/family/order/tax_class/phylum
            # as the search type OR in any extra slot), drop the secondary-kingdom
            # union: the user filtered at a taxonomic level and expects the donut
            # to reflect primary kingdoms of result compounds only. Cross-kingdom
            # secondary attestations get filtered out.
            tax_filter_active = (
                st in ('genus','family','order','tax_class','clade','phylum') or
                any(request.args.get(f'extra_type_{i}','') in ('genus','family','order','tax_class','clade','phylum')
                    for i in range(1, 11))
            )
            # Result-set kingdom donut: count PRIMARY kingdom only, matching the
            # expanded taxonomy tree's primary-lineage representation. A name query
            # for a plant compound (e.g. curcumin) shows plant-only, not incidental
            # secondary-kingdom attestations. (The browse/search kingdom FILTER still
            # matches primary+secondary for findability; this is the visualization.)
            kingdom_sql = f"""
                WITH result_set AS ({query})
                SELECT rt.kingdom AS kingdom, COUNT(DISTINCT rs.comp_id) AS cnt
                FROM result_set rs JOIN resolved_taxonomy rt ON rt.comp_id = rs.comp_id
                WHERE rt.kingdom IS NOT NULL AND rt.kingdom != ''
                GROUP BY 1 ORDER BY cnt DESC"""
            with get_db() as conn:
                with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                    cur.execute(kingdom_sql, params)
                    thumb_kingdoms = cur.fetchall()
            # If no explicit kingdom override but the query is taxonomic and one
            # kingdom dominates, treat that kingdom as implied.
            if kingdom_label_override is None and st in ("genus", "family", "order", "tax_class", "clade", "phylum") and thumb_kingdoms:
                top = thumb_kingdoms[0]
                if top["cnt"] / total >= 0.95:
                    kingdom_label_override = top["kingdom"]
            if kingdom_label_override:
                thumb_kingdoms = [{"kingdom": kingdom_label_override, "cnt": total}]
            thumb = kingdom_thumbnail_svg(thumb_kingdoms, total=total, size=150, title="Result-set kingdoms")
        except Exception:
            thumb = None
    # Region distribution for the world-map (same subquery pattern as kingdom thumb)
    region_counts_filtered = []
    region_css = ""
    region_titles = {}
    if total > 0:
        try:
            region_sql = f"SELECT crm.macro_region AS region, COUNT(DISTINCT sub.comp_id) AS count FROM ({query}) AS sub JOIN compound_region_map crm ON crm.comp_id = sub.comp_id GROUP BY crm.macro_region ORDER BY count DESC"
            with get_db() as conn_rs:
                with conn_rs.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur_rs:
                    cur_rs.execute(region_sql, params)
                    region_counts_filtered = cur_rs.fetchall()
            region_css, region_titles = build_region_color_css(region_counts_filtered)
        except Exception:
            region_counts_filtered = []
            region_css = ""
            region_titles = {}
    # Recompute linear_tree scoped to the search result-set
    try:
        lt_sql = f"""
            WITH result_set AS ({query})
            SELECT rt.taxorder AS name, rt.kingdom AS kingdom,
                   COUNT(DISTINCT rs.comp_id) AS count
            FROM result_set rs JOIN resolved_taxonomy rt ON rt.comp_id = rs.comp_id
            WHERE rt.taxorder IS NOT NULL
            GROUP BY 1, 2 ORDER BY count DESC LIMIT 15"""
        with get_db() as conn_lt:
            with conn_lt.cursor() as cur_lt:
                cur_lt.execute(lt_sql, params)
                linear_tree = [{"name": r[0], "kingdom": r[1], "count": r[2]} for r in cur_lt.fetchall()]
    except Exception:
        pass
    return render_template("search.html", results=results, query=q, search_type=st,
                           searched=bool(q or has_extra or has_range),
                           page=page, total=total, pages=pages, sort=sort, order=order, per_page=per_page,
                           thumb=thumb, linear_tree=linear_tree,
                           region_css=region_css, region_counts_filtered=region_counts_filtered, region_titles=region_titles)


# Region-to-ISO2 mapping (mirrors scripts/build_compound_region_map.py)
_REGION_TO_ISO2 = {
    "East Asia":      ["cn","jp","kr","kp","tw","mn","hk","mo"],
    "South Asia":     ["in","pk","bd","np","lk","bt","mv","af"],
    "Southeast Asia": ["id","th","vn","my","ph","sg","kh","la","mm","bn","tl"],
    "Europe":         ["de","fr","gb","it","es","nl","be","ch","at","se","no","dk","fi","is","pt","gr","pl","cz","sk","hu","ro","bg","hr","si","rs","ba","me","mk","al","xk","ua","by","md","lt","lv","ee","ie","lu","mt","cy","gl"],
    "North America":  ["us","ca","mx"],
    "Latin America":  ["br","ar","co","cl","pe","ve","ec","bo","py","uy","gy","sr","gf","cu","do","ht","jm","tt","pr","bs","bb","bz","cr","gt","hn","ni","pa","sv"],
    "Africa":         ["ng","et","eg","cd","tz","ke","za","ug","dz","sd","ma","ao","mz","gh","mg","cm","ci","ne","bf","ml","mw","zm","sn","so","td","zw","gn","rw","bj","tn","bi","ss","tg","sl","ly","cf","cg","lr","mr","er","na","bw","gm","ls","gw","gq","mu","dj","km","cv","st","sc"],
    "Middle East":    ["sa","ir","iq","il","jo","ae","kw","lb","om","qa","sy","ye","bh","ps","tr"],
    "Russia/CIS":     ["ru","kz","uz","tm","kg","tj","am","az","ge"],
    "Central Asia":   ["kz","uz","tm","kg","tj"],
    "Australia":      ["au"],
    "New Zealand":    ["nz"],
    "Oceania":        ["pg","fj","sb","vu","nc","pf","ws","to","ki","mh","fm","pw","nr","tv","ck","as","gu","mp"],
}

def build_region_color_css(region_counts, mode="dark"):
    """Given a list of {region, count} dicts, return CSS style rules for the world-map SVG.
    Sqrt scale of share-of-max. Bacteria-blue (#0077b6) ramp.
    mode="dark": low=black; mode="light": low=light blue (legacy index-page look).
    Hover lightens 25% toward white."""
    counts_by_iso2 = {}
    region_by_iso2 = {}
    for entry in region_counts:
        region = entry.get("region") or entry.get("kingdom")
        count = entry.get("count") or entry.get("cnt") or entry.get("n") or 0
        iso2_list = _REGION_TO_ISO2.get(region, [])
        for iso2 in iso2_list:
            # Dominant-region colouring: a country shared by >1 macro-region
            # (e.g. Kazakhstan in both Russia/CIS and Central Asia) takes its
            # strongest region's count, not the sum, to avoid colour inflation.
            if count > counts_by_iso2.get(iso2, 0):
                counts_by_iso2[iso2] = count
                region_by_iso2[iso2] = region
    if not counts_by_iso2:
        return "", {}
    max_count = max(counts_by_iso2.values())
    import math
    end_color = (0, 119, 182)  # bacteria-blue
    if mode == "light":
        start_color = (232, 241, 248)
        zero_color = "#f0f0f0"
    else:
        start_color = (180, 188, 196)
        zero_color = "#b0b8c0"
    def color(c):
        if c <= 0: return zero_color
        t = c / max_count  # linear share-of-max
        return f"#{int(start_color[0]+(end_color[0]-start_color[0])*t):02x}{int(start_color[1]+(end_color[1]-start_color[1])*t):02x}{int(start_color[2]+(end_color[2]-start_color[2])*t):02x}"
    def hover(c):
        if c <= 0: return "#bfbfbf" if mode == "light" else "#c8d0d8"
        t = c / max_count
        r0 = start_color[0]+(end_color[0]-start_color[0])*t
        g0 = start_color[1]+(end_color[1]-start_color[1])*t
        b0 = start_color[2]+(end_color[2]-start_color[2])*t
        r = int(0.75*r0 + 0.25*255); g = int(0.75*g0 + 0.25*255); b = int(0.75*b0 + 0.25*255)
        return f"#{r:02x}{g:02x}{b:02x}"
    rules = []
    titles = {}
    for iso2, count in counts_by_iso2.items():
        rules.append(f"a#a_{iso2} path, a#a_{iso2} g path {{ fill: {color(count)}; }}")
        rules.append(f"a#a_{iso2}:hover path, a#a_{iso2}:hover g path {{ fill: {hover(count)}; }}")
        titles[iso2] = f"{region_by_iso2[iso2]}: {count:,} compounds"
    rules.append("path { stroke: #fff; stroke-width: 0.3; transition: fill 0.15s; }")
    rules.append("a:hover { cursor: pointer; }")
    return "\n".join(rules), titles


def compute_linear_tree(where_sql=None, where_params=None, limit=15):
    """Return top-N orders by compound count for given WHERE context.
    where_sql should reference 'compounds c'. If None, returns global top-N from static cache."""
    if where_sql is None:
        # Static / unfiltered: read precomputed cache if available
        cache_path = os.path.join(app.static_folder, "linear_tree_global.json")
        if os.path.exists(cache_path):
            try:
                with open(cache_path) as f:
                    return json.load(f)[:limit]
            except Exception:
                pass
    sql = """
        SELECT rt.taxorder AS name,
               rt.kingdom AS kingdom,
               COUNT(DISTINCT c.comp_id) AS count
        FROM compounds c
        JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id
        WHERE rt.taxorder IS NOT NULL
    """
    params = list(where_params or [])
    if where_sql:
        # Rewrite unqualified `compounds.` references (from /browse and /search clause
        # builders that assume table name) to the local alias `c.`
        where_sql_local = where_sql.replace('compounds.', 'c.')
        sql += " AND " + where_sql_local
    sql += """ GROUP BY rt.taxorder, rt.kingdom ORDER BY count DESC LIMIT %s"""
    params.append(limit)
    with get_db() as conn:
        with conn.cursor() as cur:
            cur.execute(sql, tuple(params))
            return [{"name": r[0], "kingdom": r[1], "count": r[2]} for r in cur.fetchall()]


def compute_chem_tree():
    """Read the precomputed NPClassifier chemistry-tree cache (static/chem_tree.json).
    Never builds in-request: the cache is generated offline by
    scripts/build_chem_tree_cache.py at deploy time. Returns empty if absent."""
    cache_path = os.path.join(app.static_folder, "chem_tree.json")
    try:
        with open(cache_path) as f:
            return json.load(f)
    except Exception:
        return {"pathways": [], "top_superclass": None}

@app.route("/tree")
def tree_page():
    """Dedicated cladogram page with all the search/browse filter context."""
    linear_tree = compute_linear_tree()
    chem_tree = compute_chem_tree()
    return render_template("tree.html", linear_tree=linear_tree, chem_tree=chem_tree)


@app.route("/browse")
def browse():
    # If the URL uses /search-style params (type=, q=, extra_*), redirect there
    # so the user gets the result they expected. Preserves all current params.
    has_search_params = (
        request.args.get("type") or request.args.get("q") or
        any(request.args.get(f"extra_type_{i}") or request.args.get(f"extra_q_{i}") for i in range(1, 11))
    )
    if has_search_params:
        # Translate region=X into extra_type_N=region&extra_q_N=X for /search
        args = dict(request.args)
        region_val = args.pop("region", None)
        if region_val:
            for i in range(1, 11):
                if not args.get(f"extra_type_{i}"):
                    args[f"extra_type_{i}"] = "region"
                    args[f"extra_q_{i}"] = region_val
                    break
        return redirect(url_for("search", **args))
    page = max(1, int(request.args.get("page",1)))
    per_page = get_per_page()
    kingdom = request.args.get("kingdom","")
    source = request.args.get("source","")
    region = request.args.get("region","")
    named = request.args.get("named","")
    chem_class = request.args.get("class","")
    sort, order, oc = get_sort()
    clauses, params = [], ()
    if kingdom:
        clauses.append("EXISTS(SELECT 1 FROM resolved_taxonomy rt2 WHERE rt2.comp_id = compounds.comp_id AND (LOWER(rt2.kingdom) = LOWER(%s) OR LOWER(%s) = ANY(rt2.secondary_kingdoms)))"); params += (kingdom, kingdom)
    if source:
        clauses.append("(LOWER(source_db) = LOWER(%s) OR all_sources LIKE %s)"); params += (source, f"%{source}%")
    if region:
        if region in ("unresolved", "global / unresolved"):
            clauses.append("(region IS NULL OR region='' OR region='global')")
        else:
            clauses.append("EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = compounds.comp_id AND LOWER(crm.macro_region) = LOWER(%s))"); params += (region,)
    license_filter = request.args.get("license", "all")
    if license_filter == "commercial":
        clauses.append("tier_rank <= 1")
    elif license_filter == "academic":
        clauses.append("tier_rank <= 4")
    if named:
        clauses.append("name IS NOT NULL AND name != ''")
    if chem_class:
        clauses.append("(np_class ILIKE %s OR classyfire_superclass ILIKE %s)")
        params += (f"%{chem_class}%", f"%{chem_class}%")
    where = "WHERE "+" AND ".join(clauses) if clauses else ""
    query = f"SELECT * FROM compounds {where} {oc}"
    with get_db() as conn:
        results, total, pages = paginate(query, params, page, per_page, conn)
        if kingdom:
            # User explicitly filtered by kingdom: every returned compound IS
            # that kingdom (primary or secondary). Show 100% queried kingdom.
            thumb_kingdoms = [{"kingdom": kingdom.lower(), "cnt": total}]
        elif source or region or named or chem_class or license_filter != "all":
            # Other filter active: show primary-kingdom breakdown of result set.
            count_query = f"""
                WITH result_set AS (SELECT compounds.comp_id FROM compounds {where}),
                expanded AS (
                    SELECT rs.comp_id, rt.kingdom AS kingdom
                    FROM result_set rs JOIN resolved_taxonomy rt ON rt.comp_id = rs.comp_id
                    WHERE rt.kingdom IS NOT NULL AND rt.kingdom != ''
                    UNION ALL
                    SELECT rs.comp_id, unnest(rt.secondary_kingdoms) AS kingdom
                    FROM result_set rs JOIN resolved_taxonomy rt ON rt.comp_id = rs.comp_id
                    WHERE rt.secondary_kingdoms IS NOT NULL AND array_length(rt.secondary_kingdoms, 1) > 0
                )
                SELECT kingdom, COUNT(DISTINCT comp_id) AS cnt FROM expanded
                GROUP BY 1 ORDER BY cnt DESC"""
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute(count_query, params)
                thumb_kingdoms = cur.fetchall()
        else:
            thumb_kingdoms = _BROWSE_CACHE["all_kingdoms"]
    thumb = kingdom_thumbnail_svg(thumb_kingdoms, total=total, size=150, title="Result-set kingdoms")
    compounds_by_country = []
    region_counts_filtered = []
    try:
        any_filter = bool(kingdom or source or region or named or chem_class or license_filter != "all")
        if any_filter:
            with get_db() as conn_rm:
                with conn_rm.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur_rm:
                    region_q = f"SELECT crm.macro_region AS region, COUNT(DISTINCT compounds.comp_id) AS count FROM compounds JOIN compound_region_map crm ON crm.comp_id = compounds.comp_id {where} GROUP BY crm.macro_region ORDER BY count DESC" if where else "SELECT crm.macro_region AS region, COUNT(DISTINCT compounds.comp_id) AS count FROM compounds JOIN compound_region_map crm ON crm.comp_id = compounds.comp_id GROUP BY crm.macro_region ORDER BY count DESC"
                    cur_rm.execute(region_q, params)
                    region_counts_filtered = cur_rm.fetchall()
        else:
            with open("/home/thorben.klamt/theobroma/static/compounds_by_country.json") as _f:
                _data = json.load(_f)
                region_counts_filtered = _data.get("region_counts", [])
    except Exception:
        region_counts_filtered = []
    region_css, region_titles = build_region_color_css(region_counts_filtered)
    linear_tree = compute_linear_tree(
        where_sql=(" AND ".join(clauses) if clauses else None),
        where_params=params if clauses else None)
    return render_template("browse.html", results=results, page=page, total=total,
                           pages=pages, kingdom=kingdom, source=source, region=region,
                           all_kingdoms=_BROWSE_CACHE["all_kingdoms"],
                           all_sources_list=_BROWSE_CACHE["all_sources_list"],
                           thumb=thumb,
                           all_regions=_BROWSE_CACHE["all_regions"],
                           sort=sort, order=order,
                           named=named, per_page=per_page, chem_class=chem_class,
                           compounds_by_country=compounds_by_country,
                           region_css=region_css, region_counts_filtered=region_counts_filtered, region_titles=region_titles,
                           linear_tree=linear_tree)

@app.route("/compound/<comp_id>")
def compound_detail(comp_id):
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT * FROM compounds WHERE comp_id=%s", (comp_id,))
            c = cur.fetchone()
    if not c: abort(404)
    src_list = [s.strip() for s in (c.get("all_sources") or "").split("|") if s.strip()]
    synonyms = []
    admet_data = None
    if c.get("inchikey"):
        with get_db() as conn2:
            with conn2.cursor() as cur2:
                cur2.execute("SELECT synonym FROM compound_synonyms WHERE inchikey=%s LIMIT 20", (c["inchikey"],))
                synonyms = [r[0] for r in cur2.fetchall()]
    # Canonical display name (chebi_name > synonym > IUPAC), then presentation-capitalize.
    disp = canonical_display_name(c.get("chebi_name"), c.get("name"), synonyms, c.get("chebi_iupac_name"))
    if disp:
        c["name"] = display_name(disp)
    # Order synonyms, presentation-capitalize each, and drop any that duplicate the
    # primary name case-insensitively (so "curcumin" name + "Curcumin" synonym collapse).
    name_key = (c.get("name") or "").strip().casefold()
    seen_syn = set()
    ordered_syns = []
    for s in order_synonyms(synonyms):
        k = (s or "").strip().casefold()
        if not k or k == name_key or k in seen_syn:
            continue
        seen_syn.add(k)
        ordered_syns.append(display_name(s))
    synonyms = ordered_syns
    with get_db() as conn3:
        with conn3.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur3:
            try:
                cur3.execute("SELECT * FROM admet WHERE comp_id=%s", (c["comp_id"],))
                admet_data = cur3.fetchone()
            except:
                conn3.rollback()
    taxonomy = []
    with get_db() as conn4:
        with conn4.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur4:
            try:
                cur4.execute("""SELECT token, accepted_taxon, genus, family, native_distribution
                                FROM compound_taxonomy
                                WHERE comp_id = %s
                                ORDER BY (genus IS NULL), token""", (c["comp_id"],))
                taxonomy = cur4.fetchall()
            except Exception:
                conn4.rollback()
    # Resolved taxonomic lineage from resolved_taxonomy: kingdom-restricted
    # majority-voted lineage, distinct from the full multi-attestation list
    # in `taxonomy`. Rendered as a "Resolved lineage" block on the compound
    # page so the canonical resolver output is visible alongside the
    # all-attestations listing.
    resolved_lineage = None
    with get_db() as conn_rt:
        with conn_rt.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur_rt:
            try:
                cur_rt.execute("""SELECT kingdom, phylum, taxclass, taxorder,
                                         family, genus, secondary_kingdoms
                                  FROM resolved_taxonomy WHERE comp_id = %s""",
                               (c["comp_id"],))
                resolved_lineage = cur_rt.fetchone()
            except Exception:
                conn_rt.rollback()
    class_hierarchy = []
    if c.get("np_class"):
        atomic_classes = [a.strip() for a in c["np_class"].split(" $ ") if a.strip()]
        if atomic_classes:
            with get_db() as conn5:
                with conn5.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur5:
                    try:
                        placeholders = ",".join(["%s"] * len(atomic_classes))
                        cur5.execute(
                            f"SELECT DISTINCT class_name, superclass_name, pathway_name FROM npc_class_parents WHERE class_name IN ({placeholders}) ORDER BY class_name",
                            tuple(atomic_classes)
                        )
                        class_hierarchy = cur5.fetchall()
                    except Exception:
                        conn5.rollback()
    stereoisomers = []
    if c.get("inchikey"):
        ik_prefix = c["inchikey"][:14]
        with get_db() as conn_si:
            with conn_si.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur_si:
                try:
                    cur_si.execute("""SELECT comp_id, name, smiles, inchikey, source_db, license_tier
                                      FROM compounds
                                      WHERE SUBSTRING(inchikey,1,14)=%s
                                      ORDER BY (name IS NULL), comp_id""", (ik_prefix,))
                    stereoisomers = cur_si.fetchall()
                except Exception:
                    conn_si.rollback()
    stereoisomer_groups = {'reference': [], 'stereo': [], 'unspecified': [], 'protonation': []}
    if c.get("inchikey") and len(c["inchikey"]) >= 27:
        ref_stereo = c["inchikey"][15:25]
        ref_proton = c["inchikey"][26]
        for si in stereoisomers:
            ik = si.get("inchikey") or ""
            if len(ik) < 27:
                stereoisomer_groups['stereo'].append(si)
                continue
            sib_stereo = ik[15:25]
            sib_proton = ik[26]
            if si['comp_id'] == c['comp_id']:
                stereoisomer_groups['reference'].append(si)
            elif sib_stereo == 'UHFFFAOYSA' and ref_stereo != 'UHFFFAOYSA':
                stereoisomer_groups['unspecified'].append(si)
            elif sib_stereo == ref_stereo and sib_proton != ref_proton:
                stereoisomer_groups['protonation'].append(si)
            else:
                stereoisomer_groups['stereo'].append(si)
    # License provenance: per-source attestation chain to display alongside
    # the resolved license tier. Lets users audit the resolution decision
    # without leaving the compound detail page.
    license_attestations = []
    with get_db() as conn_lic:
        with conn_lic.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur_lic:
            try:
                cur_lic.execute("""SELECT source, license_tier, attested_at
                                   FROM per_source_license_attestation
                                   WHERE comp_id = %s
                                   ORDER BY source""", (c["comp_id"],))
                license_attestations = cur_lic.fetchall()
            except Exception:
                conn_lic.rollback()
    return render_template("compound.html", c=c, all_sources_list=src_list, synonyms=synonyms, admet=admet_data, taxonomy=taxonomy, resolved_lineage=resolved_lineage, class_hierarchy=class_hierarchy, stereoisomers=stereoisomers, stereoisomer_groups=stereoisomer_groups, license_attestations=license_attestations)

@app.route("/statistics")
def statistics():
    # Prefer the offline-built cache JSON if present and recent.
    import os, json, time
    cache_path = "/home/thorben.klamt/theobroma/static/statistics_cache.json"
    cache_max_age_sec = 7 * 24 * 3600  # 7 days; cron should refresh daily but we tolerate gaps
    cached = None
    if os.path.exists(cache_path):
        try:
            age = time.time() - os.path.getmtime(cache_path)
            if age < cache_max_age_sec:
                with open(cache_path) as f:
                    cached = json.load(f)
        except Exception:
            cached = None
    if cached is not None:
        visitors_by_country = []
        visitors_meta = {}
        try:
            with open("/home/thorben.klamt/theobroma/static/visitors_by_country.json") as f:
                visitors_by_country = json.load(f)
            with open("/home/thorben.klamt/theobroma/static/visitors_meta.json") as f:
                visitors_meta = json.load(f)
        except Exception:
            pass
        return render_template("statistics.html",
                               total=cached.get("total", 0),
                               kingdoms=cached.get("kingdoms", []),
                               sources=cached.get("sources", []),
                               regions=cached.get("regions", []),
                               prop_stats=cached.get("prop_stats", {}),
                               licenses=cached.get("licenses", []),
                               multi_source=cached.get("multi_source", 0),
                               admet_stats=cached.get("admet_stats", {}),
                               visitors_by_country=visitors_by_country,
                               visitors_meta=visitors_meta)
    # Live-query fallback below (used if cache missing/stale)
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT COUNT(*) AS cnt FROM compounds")
            total = cur.fetchone()["cnt"]
            cur.execute("SELECT kingdom, COUNT(*) AS cnt FROM resolved_taxonomy WHERE kingdom IS NOT NULL AND kingdom!='' GROUP BY kingdom ORDER BY cnt DESC")
            kingdoms = cur.fetchall()
            cur.execute("SELECT source_db, COUNT(*) AS cnt FROM compounds GROUP BY source_db ORDER BY cnt DESC")
            sources = cur.fetchall()
            cur.execute("""SELECT region, cnt FROM (
                   SELECT macro_region AS region, COUNT(DISTINCT comp_id) AS cnt FROM compound_region_map GROUP BY 1
                   UNION ALL
                   SELECT 'global / unresolved' AS region,
                          (SELECT COUNT(*) FROM compounds c WHERE NOT EXISTS
                             (SELECT 1 FROM compound_region_map m WHERE m.comp_id=c.comp_id)) AS cnt
                 ) t ORDER BY cnt DESC LIMIT 30""")
            regions = cur.fetchall()
            cur.execute("SELECT AVG(mw) AS avg_mw, AVG(logp) AS avg_logp, AVG(tpsa) AS avg_tpsa, AVG(hba) AS avg_hba, AVG(hbd) AS avg_hbd FROM compounds WHERE mw IS NOT NULL")
            prop_stats = cur.fetchone()
            cur.execute("SELECT license_tier, COUNT(*) AS cnt FROM compounds GROUP BY license_tier ORDER BY cnt DESC")
            licenses = cur.fetchall()
            cur.execute("SELECT COUNT(*) AS cnt FROM compounds WHERE all_sources LIKE '%%|%%'")
            multi_source = cur.fetchone()["cnt"]
    visitors_by_country = []
    visitors_meta = {}
    try:
        import json as _json
        with open("/home/thorben.klamt/theobroma/static/visitors_by_country.json") as _f:
            visitors_by_country = _json.load(_f)
        with open("/home/thorben.klamt/theobroma/static/visitors_meta.json") as _f:
            visitors_meta = _json.load(_f)
    except Exception:
        pass
    # ADMET summary stats
    admet_stats = {}
    try:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            key_admet = [
                ("hERG", "hERG risk"),
                ("AMES", "AMES mutagenicity"),
                ("BBB_Martins", "BBB penetration"),
                ("HIA_Hou", "Intestinal absorption"),
                ("Caco2_Wang", "Caco-2 permeability"),
                ("Solubility_AqSolDB", "Solubility"),
                ("Lipophilicity_AstraZeneca", "Lipophilicity"),
                ("CYP2D6_Veith", "CYP2D6 inhib."),
                ("CYP3A4_Veith", "CYP3A4 inhib."),
                ("CYP2C9_Veith", "CYP2C9 inhib."),
                ("ClinTox", "Clinical tox."),
                ("DILI", "DILI risk"),
                ("Bioavailability_Ma", "Bioavailability"),
                ("Pgp_Broccatelli", "P-gp inhibition"),
                ("VDss_Lombardo", "Volume of distribution"),
            ]
            for col, label in key_admet:
                try:
                    cur.execute(f'SELECT AVG("{col}"), MIN("{col}"), MAX("{col}") FROM admet WHERE "{col}" IS NOT NULL')
                    row = cur.fetchone()
                    if row and row["avg"] is not None:
                        admet_stats[label] = {"avg": round(float(row["avg"]), 3), "min": round(float(row["min"]), 3), "max": round(float(row["max"]), 3)}
                except:
                    conn.rollback()
    except:
        pass
    return render_template("statistics.html", total=total, kingdoms=kingdoms,
                           sources=sources, regions=regions, prop_stats=prop_stats,
                           licenses=licenses, multi_source=multi_source, admet_stats=admet_stats, visitors_by_country=visitors_by_country, visitors_meta=visitors_meta)

def _build_filename(args, ext):
    """Construct a descriptive filename from active query args.
    Format: theobroma_<key>_<val>_..._<ext>. Falls back to theobroma_export."""
    keys = ("q", "type", "kingdom", "region", "source", "class",
            "license", "tier", "named")
    parts = []
    for k in keys:
        v = args.get(k, "").strip()
        if v and v not in ("all", ""):
            safe = "".join(c if c.isalnum() else "_" for c in v)[:24]
            parts.append(f"{k}_{safe}")
    base = "theobroma_" + "_".join(parts) if parts else "theobroma_export"
    return f"{base}.{ext}"

@app.route("/api/cladogram")
def api_cladogram():
    """Return per-order compound counts with optional family breakdown."""
    where_clauses = []
    where_params = []
    kingdom = request.args.get("kingdom")
    source = request.args.get("source")
    region = request.args.get("region")
    license_f = request.args.get("license", "all")
    named = request.args.get("named")
    if kingdom:
        where_clauses.append("EXISTS(SELECT 1 FROM resolved_taxonomy rt2 WHERE rt2.comp_id=c.comp_id AND (LOWER(rt2.kingdom)=LOWER(%s) OR LOWER(%s)=ANY(rt2.secondary_kingdoms)))")
        where_params.extend([kingdom, kingdom])
    if source:
        where_clauses.append("(LOWER(c.source_db)=LOWER(%s) OR c.all_sources ILIKE %s)")
        where_params.extend([source, "%"+source+"%"])
    if region:
        where_clauses.append("EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = c.comp_id AND LOWER(crm.macro_region) = LOWER(%s))")
        where_params.append(region)
    if license_f == "commercial":
        where_clauses.append("c.tier_rank <= 1")
    elif license_f == "academic":
        where_clauses.append("c.tier_rank <= 4")
    if named:
        where_clauses.append("c.name IS NOT NULL AND c.name != ''")
    # Search-page filters: type/q + extra_type_N/extra_q_N
    search_type = request.args.get("type", "").strip()
    search_q = request.args.get("q", "").strip()
    exact = request.args.get("exact", "false").lower() == "true"
    # Resolve property/classification meta-type to the concrete search type
    # (mirrors the /search route) so type=property&prop_type=X scopes the tree.
    if search_type in ("property", "classification"):
        _prop_aliases = {
            "class": "npclassifier_class", "chemical_class": "npclassifier_class",
            "chem_class": "npclassifier_class", "npclassifier_class": "npclassifier_class",
            "npc_class": "npclassifier_class", "np_class": "npclassifier_class",
            "superclass": "npclassifier_superclass", "np_superclass": "npclassifier_superclass",
            "npclassifier_superclass": "npclassifier_superclass",
            "classyfire_class": "classyfire_class", "cf_class": "classyfire_class",
            "classyfire_superclass": "classyfire_class",
            "pathway": "pathway", "np_pathway": "pathway",
            "genus": "genus", "family": "family", "order": "order",
            "tax_class": "tax_class", "phylum": "phylum",
        }
        search_type = _prop_aliases.get(request.args.get("prop_type", "").strip().lower(), "npclassifier_class")


    def _add_search_clause(t, qval):
        """Append a WHERE clause for a search-type/value pair using same patterns as /search."""
        if not t or not qval:
            return
        # Fuzzy taxonomic resolution
        local_q = qval
        if t in ("genus", "family", "order", "tax_class", "clade", "phylum"):
            try:
                resolved = resolve_taxon(t, qval)
                if resolved is not None and resolved != qval:
                    local_q = resolved
                elif resolved is None:
                    local_q = "__no_such_taxon__"
            except Exception:
                pass
        clauses = {
            "name": (("LOWER(c.name) = LOWER(%s)", [local_q]) if exact else
                     ("c.name ILIKE %s", ["%"+local_q+"%"])),
            "organism": (("LOWER(c.source_organism) = LOWER(%s)", [local_q]) if exact else
                         ("c.source_organism ILIKE %s", ["%"+local_q+"%"])),
            "kingdom": ("(LOWER(rt.kingdom)=LOWER(%s) OR LOWER(%s)=ANY(rt.secondary_kingdoms))", [local_q, local_q]),
            "source": ("(LOWER(c.source_db)=LOWER(%s) OR c.all_sources ILIKE %s)", [local_q, "%"+local_q+"%"]),
            "region": ("EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = c.comp_id AND LOWER(crm.macro_region) = LOWER(%s))", [local_q]),
            "genus": (("LOWER(rt.genus)=LOWER(%s)", [local_q]) if exact else
                      ("rt.genus ILIKE %s", ["%"+local_q+"%"])),
            "family": (("LOWER(rt.family)=LOWER(%s)", [local_q]) if exact else
                       ("rt.family ILIKE %s", ["%"+local_q+"%"])),
            "order": (("LOWER(rt.taxorder)=LOWER(%s)", [local_q]) if exact else
                      ("rt.taxorder ILIKE %s", ["%"+local_q+"%"])),
            "clade": (("LOWER(rt.taxclass)=LOWER(%s)", [local_q]) if exact else ("rt.taxclass ILIKE %s", [f"%{local_q}%"])),
            "tax_class": (("LOWER(rt.taxclass)=LOWER(%s)", [local_q]) if exact else
                          ("rt.taxclass ILIKE %s", ["%"+local_q+"%"])),
            "phylum": (("LOWER(rt.phylum)=LOWER(%s)", [local_q]) if exact else
                       ("rt.phylum ILIKE %s", ["%"+local_q+"%"])),
            "inchikey": ("c.inchikey = %s", [local_q]),
            "smiles": ("c.smiles = %s", [local_q]),
        }
        c = clauses.get(t)
        if c:
            where_clauses.append(c[0])
            where_params.extend(c[1])

    _add_search_clause(search_type, search_q)
    for i in range(1, 11):
        _add_search_clause(request.args.get(f"extra_type_{i}", "").strip(),
                           request.args.get(f"extra_q_{i}", "").strip())

    where_sql = ("WHERE " + " AND ".join(where_clauses)) if where_clauses else ""
    conn_join = ("AND" if where_sql else "WHERE")
    order_sql = (
        "SELECT rt.taxorder, rt.kingdom, COUNT(DISTINCT c.comp_id) AS n "
        "FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id "
        + where_sql + " " + conn_join + " rt.taxorder IS NOT NULL "
        "GROUP BY rt.taxorder, rt.kingdom ORDER BY n DESC LIMIT 50"
    )
    family_sql = (
        "SELECT rt.taxorder, rt.family, COUNT(DISTINCT c.comp_id) AS n "
        "FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id "
        + where_sql + " " + conn_join + " rt.taxorder IS NOT NULL AND rt.family IS NOT NULL "
        "GROUP BY rt.taxorder, rt.family ORDER BY rt.taxorder, n DESC"
    )
    orders, families_by_order, total = [], {}, 0
    with get_db() as conn:
        with conn.cursor() as cur:
            cur.execute(order_sql, tuple(where_params))
            for row in cur.fetchall():
                orders.append({"name": row[0], "kingdom": row[1], "count": row[2], "has_families": False, "families": []})
                total += row[2]
            cur.execute(family_sql, tuple(where_params))
            for row in cur.fetchall():
                families_by_order.setdefault(row[0], []).append({"name": row[1], "count": row[2]})
    for o in orders:
        if o["name"] in families_by_order:
            o["families"] = families_by_order[o["name"]][:25]
            o["has_families"] = len(o["families"]) > 0
    return jsonify({"orders": orders, "total": total, "n_orders": len(orders)})



@app.route("/api/cladogram/export_pdf", methods=["POST"])
def api_cladogram_export_pdf():
    """Receive serialized SVG markup, return rendered PDF via cairosvg."""
    import cairosvg
    from io import BytesIO
    svg_markup = request.get_data(as_text=True)
    if not svg_markup or "<svg" not in svg_markup:
        return jsonify({"error": "invalid or empty SVG payload"}), 400
    buf = BytesIO()
    try:
        cairosvg.svg2pdf(bytestring=svg_markup.encode("utf-8"), write_to=buf)
    except Exception as e:
        return jsonify({"error": "PDF conversion failed: " + str(e)}), 500
    pdf_bytes = buf.getvalue()
    resp = Response(pdf_bytes, mimetype="application/pdf")
    resp.headers["Content-Disposition"] = 'attachment; filename="theobroma_cladogram.pdf"'
    return resp


@app.route("/api/taxonomy_tree")
def api_taxonomy_tree():
    """Return taxonomic tree. If no filter params, serve precomputed cache.
    Otherwise compute filtered tree on the fly from compounds + compound_taxonomy."""
    filter_keys = ("kingdom", "source", "region", "license", "named", "type", "q",
                   "extra_type_1", "extra_type_2", "extra_type_3", "extra_type_4", "extra_type_5",
                   "comp_ids")
    has_filter = any(request.args.get(k) for k in filter_keys)
    cache_path = os.path.join(app.static_folder, "taxonomy_cache.json")
    if not has_filter and os.path.exists(cache_path):
        return send_from_directory(app.static_folder, "taxonomy_cache.json", mimetype="application/json")
    where = []
    params = []
    # comp_ids scoping for similarity-on-tree (#4): restrict the tree to an
    # explicit compound set passed as comma-separated THEO_ ids in ?comp_ids=.
    _comp_ids_raw = request.args.get("comp_ids", "").strip()
    if _comp_ids_raw:
        _cid_list = [c for c in _comp_ids_raw.split(",") if c.startswith("THEO_")]
        if _cid_list:
            where.append("c.comp_id = ANY(%s)")
            params.append(_cid_list)
    kingdom = request.args.get("kingdom")
    source = request.args.get("source")
    region = request.args.get("region")
    license_f = request.args.get("license", "all")
    named = request.args.get("named")
    # Search-page filter params for /search context
    search_type = request.args.get("type", "")
    search_q = request.args.get("q", "").strip()
    exact = request.args.get("exact", "false").lower() == "true"
    # Resolve property/classification meta-type to the concrete search type
    # (mirrors the /search route) so type=property&prop_type=X scopes the tree.
    if search_type in ("property", "classification"):
        _prop_aliases = {
            "class": "npclassifier_class", "chemical_class": "npclassifier_class",
            "chem_class": "npclassifier_class", "npclassifier_class": "npclassifier_class",
            "npc_class": "npclassifier_class", "np_class": "npclassifier_class",
            "superclass": "npclassifier_superclass", "np_superclass": "npclassifier_superclass",
            "npclassifier_superclass": "npclassifier_superclass",
            "classyfire_class": "classyfire_class", "cf_class": "classyfire_class",
            "classyfire_superclass": "classyfire_class",
            "pathway": "pathway", "np_pathway": "pathway",
            "genus": "genus", "family": "family", "order": "order",
            "tax_class": "tax_class", "phylum": "phylum",
        }
        search_type = _prop_aliases.get(request.args.get("prop_type", "").strip().lower(), "npclassifier_class")
    # Apply search-page "extra_type_N / extra_q_N" filters (1..5) by reusing
    # the same per-rank SQL patterns we use for the primary search clause.
    for i in range(1, 11):
        et = request.args.get(f"extra_type_{i}", "").strip()
        eq = request.args.get(f"extra_q_{i}", "").strip()
        if not et or not eq:
            continue
        extra_clauses = {
            "name": (("LOWER(c.name) = LOWER(%s)", (eq,)) if exact else ("c.name ILIKE %s", (f"%{eq}%",))),
            "organism": (("LOWER(c.source_organism) = LOWER(%s)", (eq,)) if exact else ("c.source_organism ILIKE %s", (f"%{eq}%",))),
            "kingdom": ("(LOWER(rt.kingdom) = LOWER(%s) OR LOWER(%s) = ANY(rt.secondary_kingdoms))", (eq, eq)),
            "source": ("(LOWER(c.source_db) = LOWER(%s) OR c.all_sources LIKE %s)", (eq, f"%{eq}%")),
            "region": ("EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = c.comp_id AND LOWER(crm.macro_region) = LOWER(%s))", (eq,)),
            "genus": (("LOWER(rt.genus) = LOWER(%s)", (eq,)) if exact else ("rt.genus ILIKE %s", (f"%{eq}%",))),
            "family": (("LOWER(rt.family) = LOWER(%s)", (eq,)) if exact else ("rt.family ILIKE %s", (f"%{eq}%",))),
            "order": (("LOWER(rt.taxorder) = LOWER(%s)", (eq,)) if exact else ("rt.taxorder ILIKE %s", (f"%{eq}%",))),
            "tax_class": (("LOWER(rt.taxclass) = LOWER(%s)", (eq,)) if exact else ("rt.taxclass ILIKE %s", (f"%{eq}%",))),
            "clade": (("LOWER(rt.taxclass) = LOWER(%s)", (eq,)) if exact else ("rt.taxclass ILIKE %s", (f"%{eq}%",))),
            "phylum": (("LOWER(rt.phylum) = LOWER(%s)", (eq,)) if exact else ("rt.phylum ILIKE %s", (f"%{eq}%",))),
        }
        # Fuzzy resolution for taxonomic ranks
        if et in ("genus", "family", "order", "tax_class", "clade", "phylum"):
            resolved = resolve_taxon(et, eq)
            if resolved is not None and resolved != eq:
                eq = resolved
            elif resolved is None:
                eq = "__no_such_taxon__"
        ec = extra_clauses.get(et)
        if ec:
            where.append(ec[0])
            params.extend(ec[1])
    if search_type and search_q:
        # Fuzzy taxonomic resolution
        if search_type in ("genus", "family", "order", "tax_class", "clade", "phylum"):
            resolved = resolve_taxon(search_type, search_q)
            if resolved is not None and resolved != search_q:
                search_q = resolved
            elif resolved is None:
                search_q = "__no_such_taxon__"
        type_clauses = {
            "name": (("LOWER(c.name) = LOWER(%s)", (search_q,)) if exact else ("c.name ILIKE %s", (f"%{search_q}%",))),
            "organism": (("LOWER(c.source_organism) = LOWER(%s)", (search_q,)) if exact else ("c.source_organism ILIKE %s", (f"%{search_q}%",))),
            "kingdom": ("(LOWER(rt.kingdom) = LOWER(%s) OR LOWER(%s) = ANY(rt.secondary_kingdoms))", (search_q, search_q)),
            "source": ("(LOWER(c.source_db) = LOWER(%s) OR c.all_sources LIKE %s)", (search_q, f"%{search_q}%")),
            "region": ("EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = c.comp_id AND LOWER(crm.macro_region) = LOWER(%s))", (search_q,)),
            "smiles": ("c.smiles = %s", (search_q,)),
            "inchikey": ("c.inchikey = %s", (search_q,)),
            "genus": (("LOWER(rt.genus) = LOWER(%s)", (search_q,)) if exact else ("rt.genus ILIKE %s", (f"%{search_q}%",))),
            "family": (("LOWER(rt.family) = LOWER(%s)", (search_q,)) if exact else ("rt.family ILIKE %s", (f"%{search_q}%",))),
            "order": (("LOWER(rt.taxorder) = LOWER(%s)", (search_q,)) if exact else ("rt.taxorder ILIKE %s", (f"%{search_q}%",))),
            "tax_class": (("LOWER(rt.taxclass) = LOWER(%s)", (search_q,)) if exact else ("rt.taxclass ILIKE %s", (f"%{search_q}%",))),
            "clade": (("LOWER(rt.taxclass) = LOWER(%s)", (search_q,)) if exact else ("rt.taxclass ILIKE %s", (f"%{search_q}%",))),
            "phylum": (("LOWER(rt.phylum) = LOWER(%s)", (search_q,)) if exact else ("rt.phylum ILIKE %s", (f"%{search_q}%",))),
            "class": ("(c.np_class ILIKE %s OR c.classyfire_superclass ILIKE %s OR c.inferred_class ILIKE %s)", (f"%{search_q}%", f"%{search_q}%", f"%{search_q}%")),
            "npclassifier_class": ("(c.np_class ILIKE %s OR c.inferred_class ILIKE %s)", (f"%{search_q}%", f"%{search_q}%")),
            "classyfire_class": ("c.classyfire_superclass ILIKE %s", (f"%{search_q}%",)),
            "pathway": ("c.np_pathway ILIKE %s", (f"%{search_q}%",)),
        }
        clause_pair = type_clauses.get(search_type)
        if clause_pair:
            clause_sql, clause_params = clause_pair
            where.append(clause_sql)
            params.extend(clause_params)
    # Detect if extras already provide kingdom/source/region; if so, skip the
    # legacy single-param to avoid AND-combining the same dimension twice
    # (which leads to empty results when values differ, e.g. region=East+Asia
    # alongside extra_q_1=Australia).
    extras_have = set()
    for ii in range(1, 11):
        et_chk = request.args.get(f"extra_type_{ii}", "").strip()
        eq_chk = request.args.get(f"extra_q_{ii}", "").strip()
        if et_chk and eq_chk:
            extras_have.add(et_chk)
    if kingdom and "kingdom" not in extras_have:
        where.append("rt.kingdom = %s")
        params.append(kingdom)
    if source and "source" not in extras_have:
        where.append("(LOWER(c.source_db) = LOWER(%s) OR c.all_sources LIKE %s)")
        params.append(source)
        params.append(f"%{source}%")
    if region and "region" not in extras_have:
        where.append("EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = c.comp_id AND LOWER(crm.macro_region) = LOWER(%s))")
        params.append(region)
    if named:
        where.append("c.name IS NOT NULL AND c.name != ''")
    if license_f == "commercial":
        where.append("c.tier_rank <= 1")
    elif license_f == "academic":
        where.append("c.tier_rank <= 4")
    # Range filters: physchem (in compounds) + ADMET (in admet table)
    admet_cols_taxapi = {"Lipinski","QED","stereo_centers","PAINS_alert","BRENK_alert","NIH_alert","AMES","BBB_Martins","Bioavailability_Ma","CYP1A2_Veith","CYP2C19_Veith","CYP2C9_Substrate_CarbonMangels","CYP2C9_Veith","CYP2D6_Substrate_CarbonMangels","CYP2D6_Veith","CYP3A4_Substrate_CarbonMangels","CYP3A4_Veith","Carcinogens_Lagunin","ClinTox","DILI","HIA_Hou","NR_AR_LBD","NR_AR","NR_AhR","NR_Aromatase","NR_ER_LBD","NR_ER","NR_PPAR_gamma","PAMPA_NCATS","Pgp_Broccatelli","SR_ARE","SR_ATAD5","SR_HSE","SR_MMP","SR_p53","Skin_Reaction","hERG","Caco2_Wang","Clearance_Hepatocyte_AZ","Clearance_Microsome_AZ","Half_Life_Obach","HydrationFreeEnergy_FreeSolv","LD50_Zhu","Lipophilicity_AstraZeneca","PPBR_AZ","Solubility_AqSolDB","VDss_Lombardo"}
    needs_admet_join = False
    for prop in ["mw","logp","tpsa","hba","hbd","n_rings","rotatable_bonds"] + list(admet_cols_taxapi):
        pmin = request.args.get(f"{prop}_min", "")
        pmax = request.args.get(f"{prop}_max", "")
        if not pmin and not pmax:
            continue
        if prop in admet_cols_taxapi:
            needs_admet_join = True
            col = f'a."{prop}"'
        else:
            col = f"c.{prop}"
        if pmin:
            where.append(f"{col} >= %s")
            params.append(float(pmin))
        if pmax:
            where.append(f"{col} <= %s")
            params.append(float(pmax))
    admet_join_sql = " LEFT JOIN admet a ON a.comp_id = c.comp_id" if needs_admet_join else ""
    where_sql = ("WHERE " + " AND ".join(where)) if where else ""
    # When the user filtered by kingdom (directly via type=kingdom, an extra,
    # or implicitly because a taxonomic filter dominates a single kingdom),
    # the tree should only show that kingdom's slice. Detect:
    #   (a) search_type == "kingdom"
    #   (b) any extra_type_i == "kingdom"
    #   (c) search_type is taxonomic AND base set is >=95% dominated by one kingdom
    kingdom_override = None
    if search_type == "kingdom" and search_q:
        kingdom_override = search_q
    elif kingdom:
        # /browse-style URL param ?kingdom=X (no type=kingdom involved)
        kingdom_override = kingdom
    else:
        for i in range(1, 11):
            et = request.args.get(f"extra_type_{i}", "")
            eq_v = request.args.get(f"extra_q_{i}", "").strip()
            if et == "kingdom" and eq_v:
                kingdom_override = eq_v
                break
    if kingdom_override is None and search_type in ("genus", "family", "order", "tax_class", "clade", "phylum") and search_q:
        # Quick dominance check using the where_sql already built
        try:
            with get_db() as conn:
                with conn.cursor() as cur:
                    dom_sql = f"""SELECT rt.kingdom, COUNT(*) AS n
                                   FROM compounds c LEFT JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id{admet_join_sql}
                                   {where_sql}
                                   GROUP BY 1 ORDER BY 2 DESC LIMIT 1"""
                    cur.execute(dom_sql, params)
                    row = cur.fetchone()
                    if row and row[0]:
                        cur.execute(f"SELECT COUNT(*) FROM compounds c LEFT JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id{admet_join_sql} {where_sql}", params)
                        total_for_dom = cur.fetchone()[0]
                        if total_for_dom > 0 and row[1] / total_for_dom >= 0.95:
                            kingdom_override = row[0]
        except Exception:
            kingdom_override = None
    # Snapshot params BEFORE appending kingdom_override; sql_compounds (below)
    # uses only the base where_sql and must NOT see the extra placeholder.
    params_base = list(params)
    if kingdom_override:
        kingdom_filter = "WHERE LOWER(theobroma_kingdom) = LOWER(%s)"
        params.append(kingdom_override)
    else:
        kingdom_filter = ""

    # Detect taxonomic filter (incl. extras)
    api_tax_filter_active = (
        search_type in ('genus','family','order','tax_class','clade','phylum') or
        any(request.args.get(f'extra_type_{i}','') in ('genus','family','order','tax_class','clade','phylum')
            for i in range(1, 11))
        or bool(_comp_ids_raw)  # comp_ids set -> primary-lineage only (match thumbnail)
    )
    secondary_union = '' if api_tax_filter_active else (
        "UNION ALL SELECT comp_id, unnest(secondary_kingdoms), NULL, NULL, NULL, NULL, NULL "
        "FROM base WHERE secondary_kingdoms IS NOT NULL AND secondary_kingdoms <> '{{}}'"
    )
    sql = f"""
        WITH base AS (
            SELECT c.comp_id, rt.kingdom AS primary_k, rt.secondary_kingdoms,
                   rt.phylum, rt.taxclass, rt.taxorder, rt.family, rt.genus
            FROM compounds c
            LEFT JOIN resolved_taxonomy rt ON c.comp_id = rt.comp_id{admet_join_sql}
            {where_sql}
        ),
        expanded AS (
            -- When a taxonomic filter is active, only primary kingdoms count.
            -- Otherwise expand to include secondary attestations.
            SELECT comp_id, primary_k AS theobroma_kingdom, phylum, taxclass, taxorder, family, genus FROM base
            {secondary_union}
        ),
        narrowed AS (
            SELECT * FROM expanded {kingdom_filter}
        )
        SELECT
            theobroma_kingdom,
            NULL::text AS lineage_kingdom,
            phylum AS phylum,
            taxclass AS class,
            taxorder AS taxon_order,
            family AS family,
            genus AS genus,
            COUNT(DISTINCT comp_id) AS n
        FROM narrowed
        GROUP BY 1,2,3,4,5,6,7
    """
    with get_db() as conn:
        with conn.cursor() as cur:
            cur.execute(sql, params)
            rows = cur.fetchall()
    def insert(tree, path, count):
        node = tree
        for level in path:
            if level is None or level == '':
                level = '(unresolved)'
            if 'children' not in node:
                node['children'] = {}
            if level not in node['children']:
                node['children'][level] = {'name': level, 'value': 0}
            node = node['children'][level]
            node['value'] += count
    root = {'name': 'THEOBROMA', 'value': 0}
    total = 0
    for r in rows:
        theobroma_k, lineage_k, phylum, klass, order_, family, genus, n = r
        n = int(n)
        total += n
        insert(root, [theobroma_k, phylum, klass, order_, family, genus], n)
    root['value'] = total
    # When the result set is small, fetch compound names so the sunburst can
    # render compound nodes as a 7th ring under each genus.
    compounds_by_path = {}
    show_compounds = total > 0 and total <= 250
    if show_compounds:
        sql_compounds = f"""
            SELECT
                c.comp_id,
                COALESCE(c.name, '') AS name,
                rt.kingdom AS theobroma_kingdom,
                rt.phylum AS phylum,
                rt.taxclass AS class,
                rt.taxorder AS taxon_order,
                rt.family AS family,
                rt.genus AS genus
            FROM compounds c
            LEFT JOIN resolved_taxonomy rt ON c.comp_id = rt.comp_id{admet_join_sql}
            {where_sql}
            ORDER BY c.comp_id
        """
        with get_db() as conn2:
            with conn2.cursor() as cur2:
                cur2.execute(sql_compounds, params_base)
                seen_per_path = set()
                for row in cur2.fetchall():
                    cid, nm, tk, ph, cl, od, fa, gn = row
                    if (cid, tk, ph, cl, od, fa, gn) in seen_per_path: continue
                    seen_per_path.add((cid, tk, ph, cl, od, fa, gn))
                    path = tuple(x or '(unresolved)' for x in (tk, ph, cl, od, fa, gn))
                    compounds_by_path.setdefault(path, []).append({'comp_id': cid, 'name': nm})
    def to_d3(node, path=()):
        out = {'name': node['name'], 'value': node['value']}
        if 'children' in node:
            ch = [to_d3(c, path + (node['name'],)) for c in node['children'].values()]
            ch.sort(key=lambda x: -x['value'])
            out['children'] = ch
        elif show_compounds:
            # Leaf node (genus level) -- attach compound leaves if available.
            # path is (root, kingdom, phylum, class, order, family); we are at genus.
            leaf_path = path[1:] + (node['name'],)  # drop 'THEOBROMA'
            comps = compounds_by_path.get(leaf_path, [])
            if comps:
                out['children'] = [
                    {'name': (c['name'] or c['comp_id']),
                     'value': 1,
                     'comp_id': c['comp_id'],
                     'is_compound': True}
                    for c in comps
                ]
        return out
    return jsonify({"tree": to_d3(root), "total_compounds": total, "distinct_paths": len(rows), "show_compounds": show_compounds})

@app.route("/download")
def download_page():
    d = app.config["DATA_DIR"]
    files = []
    for fn in ["theobroma_final.csv","theobroma_all.sdf"]:
        p = os.path.join(d, fn)
        if os.path.exists(p):
            files.append({"name":fn, "fmt":fn.split(".")[-1].upper(), "size":f"{os.path.getsize(p)/1024/1024:.1f} MB", "url":f"/download/{fn}"})
    files.append({"name":"theobroma_all.json", "fmt":"JSON", "size":"streamed", "url":"/api/bulk?format=json"})
    return render_template("download.html", files=files)

@app.route("/download/<filename>")
def download_file(filename):
    d = app.config["DATA_DIR"]
    if not os.path.exists(os.path.join(d, filename)): abort(404)
    return send_from_directory(d, filename, as_attachment=True)

@app.route("/help")
def help_page():
    return render_template("help.html")

# --- API routes ---

@app.route("/api/search")
def api_search():
    """JSON API for programmatic access. Supports name, smiles, inchikey, kingdom, organism, region, source searches."""
    q = request.args.get("q","").strip()
    st = request.args.get("type","name")
    limit = min(10000, max(1, int(request.args.get("limit",50))))
    offset = max(0, int(request.args.get("offset",0)))
    fmt = request.args.get("format","json")
    if not q: return jsonify({"error":"query parameter 'q' required", "usage":{
        "endpoint": "/api/search",
        "params": {"q":"search query (required)", "type":"name|smiles|inchikey|kingdom|organism|region|source",
                    "limit":"max results (1-10000, default 50)", "offset":"pagination offset (default 0)",
                    "format":"json|csv (default json)"},
        "examples": ["/api/search?q=curcumin&type=name", "/api/search?q=fungi&type=kingdom&limit=100",
                     "/api/search?q=Streptomyces&type=organism&format=csv"]
    }}), 400
    cols = "comp_id,name,smiles,inchikey,kingdom,source_db,all_sources,source_organism,region,mw,logp,license_tier"
    exact = _exact_flag()
    if st == "name":
        nq = normalize_query(q)
        cols = "c.comp_id,c.name,c.smiles,c.inchikey,c.kingdom,c.source_db,c.all_sources,c.source_organism,c.region,c.mw,c.logp,c.license_tier"
        if exact:
            base = f"""SELECT {cols} FROM (
              SELECT DISTINCT ON (sn.comp_id) sn.comp_id, 0 AS relevance
              FROM search_names sn
              WHERE sn.name_norm = %s
            ) matched
            JOIN compounds c ON c.comp_id = matched.comp_id
            ORDER BY c.name"""
            pm_tuple = (nq,)
        else:
            base = f"""SELECT {cols} FROM (
              SELECT DISTINCT ON (sn.comp_id) sn.comp_id,
              CASE WHEN sn.name_norm = %s THEN 0
                   WHEN sn.name_norm LIKE %s THEN 1
                   ELSE 2 END AS relevance
              FROM search_names sn
              WHERE sn.name_norm LIKE %s
            ) matched
            JOIN compounds c ON c.comp_id = matched.comp_id
            ORDER BY matched.relevance, LENGTH(c.name), c.name"""
            pm_tuple = (nq, nq+'%', '%'+nq+'%')
        license_filter = request.args.get("license", "all")
        if license_filter == "commercial":
            base = base.replace("ORDER BY", "WHERE c.tier_rank <= 1 ORDER BY")
        elif license_filter == "academic":
            base = base.replace("ORDER BY", "WHERE c.tier_rank <= 4 ORDER BY")
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute(f"SELECT COUNT(*) FROM ({base}) sq", pm_tuple)
                total = cur.fetchone()["count"]
                cur.execute(base + " LIMIT %s OFFSET %s", pm_tuple + (limit, offset))
                results = cur.fetchall()
    else:
        tc = {"smiles":"smiles=%s","inchikey":"inchikey=%s",
              "kingdom":"EXISTS(SELECT 1 FROM resolved_taxonomy rt2 WHERE rt2.comp_id = compounds.comp_id AND (LOWER(rt2.kingdom) = LOWER(%s) OR LOWER(%s) = ANY(rt2.secondary_kingdoms)))","organism": ("LOWER(source_organism) = LOWER(%s)" if exact else "source_organism ILIKE %s"),
              "genus": ("EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.genus) = LOWER(%s))" if exact else "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND rt.genus ILIKE %s)"),
              "family": ("EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.family) = LOWER(%s))" if exact else "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND rt.family ILIKE %s)"),
              "order": ("EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.taxorder) = LOWER(%s))" if exact else "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND rt.taxorder ILIKE %s)"),
              "tax_class": ("EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.taxclass) = LOWER(%s))" if exact else "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND rt.taxclass ILIKE %s)"),
              "clade": ("EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.taxclass) = LOWER(%s))" if exact else "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND rt.taxclass ILIKE %s)"),
              "phylum": ("EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.phylum) = LOWER(%s))" if exact else "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND rt.phylum ILIKE %s)"),
              "region":"EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = compounds.comp_id AND LOWER(crm.macro_region) = LOWER(%s))","source":"LOWER(source_db) = LOWER(%s)",
              "class":"np_class ILIKE %s OR classyfire_superclass ILIKE %s OR inferred_class ILIKE %s OR np_superclass ILIKE %s OR np_pathway ILIKE %s",
              "npclassifier_class":"np_class ILIKE %s OR inferred_class ILIKE %s",
              "npclassifier_superclass":"np_superclass ILIKE %s",
              "classyfire_class":"classyfire_superclass ILIKE %s",
              "pathway":"np_pathway ILIKE %s"}
        cl = tc.get(st, "LOWER(name) LIKE %s")
        if st == "pathway":
            pm = (f"%{q}%",)
        elif st == "class":
            pm = (f"%{q}%", f"%{q}%", f"%{q}%", f"%{q}%", f"%{q}%")
        elif st == "npclassifier_class":
            pm = (f"%{q}%", f"%{q}%")
        elif st in ("classyfire_class", "npclassifier_superclass"):
            pm = (f"%{q}%",)
        elif st in ("tax_class", "clade", "order", "phylum"):
            pm = (q,) if exact else (f"%{q}%",)
        elif st == "kingdom":
            pm = (q, q)
        else:
            pm = (q if exact else f"%{q.lower()}%") if st == "organism" else (q if st in ("region","source") else q)
        # Extra AND-filters (extra_type_N / extra_q_N): mirror the /search page so
        # the API honours multi-filter intersection instead of silently ignoring
        # it. Conditions reference compounds.comp_id (the API queries FROM
        # compounds directly). Values are always bound parameters.
        extra_cl = []
        extra_pm = []
        emap = {
            "name": ("LOWER(name) = LOWER(%s)" if exact else "LOWER(name) LIKE %s"),
            "organism": ("LOWER(source_organism) = LOWER(%s)" if exact else "source_organism ILIKE %s"),
            "kingdom": "EXISTS(SELECT 1 FROM resolved_taxonomy rt2 WHERE rt2.comp_id = compounds.comp_id AND (LOWER(rt2.kingdom) = LOWER(%s) OR LOWER(%s) = ANY(rt2.secondary_kingdoms)))",
            "region": "EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = compounds.comp_id AND LOWER(crm.macro_region) = LOWER(%s))",
            "source": "LOWER(source_db) = LOWER(%s)",
            "class": "(np_class ILIKE %s OR classyfire_superclass ILIKE %s)",
            "npclassifier_class": "(np_class ILIKE %s OR inferred_class ILIKE %s)",
            "classyfire_class": "classyfire_superclass ILIKE %s",
            "npclassifier_superclass": "np_superclass ILIKE %s",
            "pathway": "np_pathway ILIKE %s",
            "genus":  "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.genus) = LOWER(%s))",
            "family": "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.family) = LOWER(%s))",
            "order":  "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.taxorder) = LOWER(%s))",
            "tax_class": "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.taxclass) = LOWER(%s))",
            "clade": "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.taxclass) = LOWER(%s))",
            "phylum": "EXISTS(SELECT 1 FROM resolved_taxonomy rt WHERE rt.comp_id = compounds.comp_id AND LOWER(rt.phylum) = LOWER(%s))",
        }
        for _i in range(1, 11):
            _et = request.args.get(f"extra_type_{_i}", "")
            _eq = request.args.get(f"extra_q_{_i}", "").strip()
            if not _et or not _eq:
                continue
            if _et in ("genus", "family", "order", "tax_class", "clade", "phylum"):
                _res = resolve_taxon(_et, _eq)
                if _res is not None and _res != _eq:
                    _eq = _res
                elif _res is None:
                    _eq = "__no_such_taxon__"
            _sql = emap.get(_et)
            if not _sql:
                continue
            extra_cl.append(_sql)
            if _et in ("class", "npclassifier_class"):
                extra_pm.extend([f"%{_eq}%", f"%{_eq}%"])
            elif _et in ("classyfire_class", "pathway", "npclassifier_superclass"):
                extra_pm.append(f"%{_eq}%")
            elif _et == "kingdom":
                extra_pm.extend([_eq, _eq])
            elif exact and _et in ("name", "organism"):
                extra_pm.append(_eq)
            elif _et in ("region", "source", "genus", "family", "order", "tax_class", "clade", "phylum"):
                extra_pm.append(_eq)
            else:
                extra_pm.append(f"%{_eq}%")
        extra_sql = ("".join(" AND (" + c + ")" for c in extra_cl))
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                params = pm if isinstance(pm, tuple) else (pm,)
                params = params + tuple(extra_pm)
                license_filter = request.args.get("license", "all")
                lic_clause = ""
                if license_filter == "commercial":
                    lic_clause = " AND tier_rank <= 1"
                elif license_filter == "academic":
                    lic_clause = " AND tier_rank <= 4"
                cur.execute(f"SELECT COUNT(*) FROM compounds WHERE ({cl}){extra_sql}{lic_clause}", params)
                total = cur.fetchone()["count"]
                cur.execute(f"SELECT {cols} FROM compounds WHERE ({cl}){extra_sql}{lic_clause} LIMIT %s OFFSET %s", params + (limit, offset))
                results = cur.fetchall()
    if fmt == "csv":
        si = io.StringIO()
        if results:
            w = csv.DictWriter(si, fieldnames=results[0].keys())
            w.writeheader()
            w.writerows(results)
        return Response(si.getvalue(), mimetype="text/csv",
                       headers={"Content-Disposition": f"attachment; filename=theobroma_{st}_{q[:20]}.csv"})
    if fmt == "json" and request.args.get("download", "").lower() in ("true", "1", "yes"):
        payload = json.dumps({"count":len(results), "total":total, "offset":offset, "limit":limit, "results":[dict(r) for r in results]}, default=str)
        return Response(payload, mimetype="application/json",
                       headers={"Content-Disposition": f"attachment; filename=theobroma_{st}_{q[:20]}.json"})
    return jsonify({"count":len(results), "total":total, "offset":offset, "limit":limit, "results":results})

@app.route("/api/compound/<comp_id>")
def api_compound(comp_id):
    """Get full compound details by comp_id, including ADMET predictions."""
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT * FROM compounds WHERE comp_id=%s", (comp_id,))
            c = cur.fetchone()
            if not c: return jsonify({"error":"compound not found"}), 404
            # Resolved taxonomic lineage merged into the response so API
            # consumers do not need a separate /resolved_taxonomy lookup.
            cur.execute("""SELECT phylum, taxclass, taxorder, family, genus,
                                  secondary_kingdoms
                           FROM resolved_taxonomy WHERE comp_id = %s""",
                        (comp_id,))
            rt = cur.fetchone()
            if rt:
                c["phylum"]             = rt.get("phylum")
                c["taxclass"]           = rt.get("taxclass")
                c["taxorder"]           = rt.get("taxorder")
                c["family"]             = rt.get("family")
                c["genus"]              = rt.get("genus")
                c["secondary_kingdoms"] = rt.get("secondary_kingdoms")
            cur.execute("SELECT * FROM admet WHERE comp_id=%s", (comp_id,))
            a = cur.fetchone()
            if a:
                a.pop("comp_id", None)
                c["admet"] = dict(a)
            # License provenance: per-source attestation list plus resolution
            # metadata. Returned alongside the main compound payload so a
            # downstream user can audit the resolved tier against the
            # underlying source attestations without a separate API call.
            cur.execute("""
                SELECT source, license_tier, attested_at
                FROM per_source_license_attestation
                WHERE comp_id = %s
                ORDER BY source
            """, (comp_id,))
            attestations = [
                {"source": r["source"],
                 "license_tier": r["license_tier"],
                 "attested_at": r["attested_at"].isoformat() if r["attested_at"] else None}
                for r in cur.fetchall()
            ]
            c["license_provenance"] = {
                "resolved_tier": c.get("license_tier"),
                "resolution_rule": "most-permissive-wins (structure-as-fact)",
                "resolution_note": "A structure is a fact obtainable from whichever source is least restrictive, so the resolved tier is the most permissive license under which the structure is available across its sources. Unspecified is treated as unknown (no signal): it never overrides a known license, and a compound is Unspecified only if all of its sources are.",
                "precedence_order": ["CC0", "CC BY 4.0", "CC BY-NC 4.0", "CC BY-NC-SA 4.0", "CC BY-NC-ND 4.0", "Unspecified"],
                "attestations": attestations,
            }
    return jsonify(c)

@app.route("/api/compound/<comp_id>/license-provenance")
def api_license_provenance(comp_id):
    """Per-compound license provenance for compliance auditing.
    Returns only the per-source attestation list and resolution metadata
    without the full compound payload, suitable for bulk audit workflows
    that do not need ADMET, taxonomy, or structural data per compound.
    """
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT comp_id, license_tier FROM compounds WHERE comp_id=%s", (comp_id,))
            c = cur.fetchone()
            if not c:
                return jsonify({"error": "compound not found"}), 404
            cur.execute("""
                SELECT source, license_tier, attested_at
                FROM per_source_license_attestation
                WHERE comp_id = %s
                ORDER BY source
            """, (comp_id,))
            attestations = [
                {"source": r["source"],
                 "license_tier": r["license_tier"],
                 "attested_at": r["attested_at"].isoformat() if r["attested_at"] else None}
                for r in cur.fetchall()
            ]
    return jsonify({
        "comp_id": comp_id,
        "resolved_tier": c["license_tier"],
        "resolution_rule": "most-permissive-wins (structure-as-fact)",
        "resolution_note": "A structure is a fact obtainable from whichever source is least restrictive, so the resolved tier is the most permissive license under which the structure is available across its sources. Unspecified is treated as unknown (no signal): it never overrides a known license, and a compound is Unspecified only if all of its sources are.",
        "precedence_order": ["CC0", "CC BY 4.0", "CC BY-NC 4.0", "CC BY-NC-SA 4.0", "CC BY-NC-ND 4.0", "Unspecified"],
        "attestations": attestations,
    })

@app.route("/api/autocomplete")
def api_autocomplete():
    """Name autocomplete for search box."""
    q = request.args.get("q","").strip()
    if len(q) < 2: return jsonify([])
    with get_db() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT DISTINCT LOWER(name) AS lname, MIN(name) AS sample FROM compounds WHERE LOWER(name) LIKE %s AND name IS NOT NULL AND name != '' GROUP BY LOWER(name) ORDER BY lname LIMIT 12",
                        (f"{q.lower()}%",))
            rows = cur.fetchall()
    # Normalize first letter to uppercase for display
    results = []
    for _, sample in rows:
        if sample and sample[0].islower():
            sample = sample[0].upper() + sample[1:]
        results.append(sample)
    return jsonify(results)

@app.route("/api/organisms")
def api_organisms():
    """Organism autocomplete."""
    q = request.args.get("q","").strip()
    if len(q) < 2: return jsonify([])
    with get_db() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT DISTINCT source_organism FROM compounds WHERE source_organism ILIKE %s AND source_organism IS NOT NULL AND source_organism != '' ORDER BY source_organism LIMIT 12",
                        (f"%{q}%",))
            results = [r[0] for r in cur.fetchall()]
    return jsonify(results)

@app.route("/api/stats")
def api_stats():
    """Database statistics as JSON."""
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT COUNT(*) AS total FROM compounds")
            total = cur.fetchone()["total"]
            cur.execute("SELECT kingdom, COUNT(*) AS cnt FROM compounds WHERE kingdom IS NOT NULL GROUP BY kingdom ORDER BY cnt DESC")
            kingdoms = cur.fetchall()
            cur.execute("SELECT source_db, COUNT(*) AS cnt FROM compounds GROUP BY source_db ORDER BY cnt DESC")
            sources = cur.fetchall()
            cur.execute("SELECT COUNT(DISTINCT source_db) AS n_sources FROM compounds")
            n_sources = cur.fetchone()["n_sources"]
            # License-tier breakdown from the resolver-applied compounds table.
            cur.execute("SELECT license_tier, COUNT(*) AS cnt FROM compounds GROUP BY license_tier ORDER BY cnt DESC")
            licenses = cur.fetchall()
    return jsonify({"total":total, "n_sources":n_sources, "kingdoms":kingdoms, "sources":sources, "licenses":licenses})

@app.route("/export")
def export_results():
    """Export current search/browse results as CSV (up to 10k rows)."""
    q = request.args.get("q","")
    st = request.args.get("type","name")
    kingdom = request.args.get("kingdom","")
    source = request.args.get("source","")
    region = request.args.get("region","")
    named = request.args.get("named","")
    limit = min(10000, max(1, int(request.args.get("limit",10000))))
    # Build query based on context
    clauses, params = [], ()
    if q:
        tmap = {"name":"LOWER(name) LIKE %s","smiles":"smiles=%s","inchikey":"inchikey=%s",
                "source":"LOWER(source_db) = LOWER(%s)","organism":"source_organism ILIKE %s",
                "region":"LOWER(region) = LOWER(%s)","kingdom":"EXISTS(SELECT 1 FROM resolved_taxonomy rt2 WHERE rt2.comp_id = compounds.comp_id AND (LOWER(rt2.kingdom) = LOWER(%s) OR LOWER(%s) = ANY(rt2.secondary_kingdoms)))"}
        cl = tmap.get(st, tmap["name"])
        pm = f"%{q.lower()}%" if st in ("name","organism") else (q if st in ("kingdom","region","source") else q)
        clauses.append(cl); params += (pm, pm) if st == "kingdom" else (pm,)
    if kingdom:
        clauses.append("EXISTS(SELECT 1 FROM resolved_taxonomy rt2 WHERE rt2.comp_id = compounds.comp_id AND (LOWER(rt2.kingdom) = LOWER(%s) OR LOWER(%s) = ANY(rt2.secondary_kingdoms)))"); params += (kingdom, kingdom)
    if source:
        clauses.append("(LOWER(source_db) = LOWER(%s) OR all_sources LIKE %s)"); params += (source, f"%{source}%")
    if region and region != "unresolved":
        clauses.append("LOWER(region) = LOWER(%s)"); params += (region,)
    elif region == "unresolved":
        clauses.append("(region IS NULL OR region='' OR region='global')")
    license_filter = request.args.get("license", "all")
    if license_filter == "commercial":
        clauses.append("tier_rank <= 1")
    elif license_filter == "academic":
        clauses.append("tier_rank <= 4")
    if named:
        clauses.append("name IS NOT NULL AND name != ''")
    where = "WHERE "+" AND ".join(clauses) if clauses else ""
    cols = "comp_id,name,smiles,inchikey,kingdom,source_db,all_sources,source_organism,region,mw,logp,tpsa,hba,hbd,license_tier"
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(f"SELECT {cols} FROM compounds {where} LIMIT %s", params + (limit,))
            results = cur.fetchall()
    fmt = request.args.get("format", "csv").lower()
    if fmt == "json":
        payload = json.dumps({"count": len(results), "results": [dict(r) for r in results]}, default=str)
        return Response(payload, mimetype="application/json",
                       headers={"Content-Disposition": f"attachment; filename={_build_filename(request.args, 'json')}"})
    si = io.StringIO()
    if results:
        w = csv.DictWriter(si, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)
    return Response(si.getvalue(), mimetype="text/csv",
                   headers={"Content-Disposition": f"attachment; filename={_build_filename(request.args, 'csv')}"})

@app.route("/similarity")
def similarity():
    query_input = (request.args.get("smiles") or request.args.get("q") or "").strip()
    query_smiles = query_input
    query_comp_id = None  # Set when input resolves to a corpus compound; enables click-through on Query structure card.
    # Resolve name or comp_id to SMILES, and capture the matching comp_id.
    if query_input and not any(c in query_input for c in "()=#@[]"):
        with get_db() as conn:
            with conn.cursor() as cur:
                if query_input.startswith("THEO_"):
                    cur.execute("SELECT comp_id, smiles FROM compounds WHERE comp_id=%s", (query_input,))
                    row = cur.fetchone()
                    if row:
                        query_comp_id = row[0]
                        query_smiles  = row[1]
                else:
                    cur.execute("SELECT comp_id, smiles FROM compounds WHERE LOWER(name)=%s LIMIT 1", (query_input.lower(),))
                    row = cur.fetchone()
                    if row:
                        query_comp_id = row[0]
                        query_smiles  = row[1]
    # Raw SMILES case: try to resolve to a corpus compound via InChIKey lookup.
    if query_comp_id is None and query_smiles and any(c in query_smiles for c in "()=#@[]"):
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(query_smiles)
            if mol is not None:
                ik = Chem.MolToInchiKey(mol)
                if ik:
                    with get_db() as conn:
                        with conn.cursor() as cur:
                            cur.execute("SELECT comp_id FROM compounds WHERE inchikey=%s LIMIT 1", (ik,))
                            row = cur.fetchone()
                            if row:
                                query_comp_id = row[0]
        except Exception:
            pass
    top_n = min(200, max(1, int(request.args.get("top_n", 50))))
    threshold = max(0.0, min(1.0, float(request.args.get("threshold", "0.3"))))
    results = []
    error = None
    if query_smiles:
        if not sim_engine.loaded:
            error = "Similarity search not available (vectors not loaded)."
        else:
            metric = request.args.get("metric", "morgan")
            if metric == "maccs":
                hits = sim_engine.maccs_search(query_smiles, top_n=top_n, threshold=threshold)
            elif metric == "chemberta":
                hits = sim_engine.chemberta_search(query_smiles, top_n=top_n)
            elif metric == "nafm":
                _cid = _resolve_comp_id_for_nafm(request.args.get("q","").strip(), query_smiles)
                hits = sim_engine.nafm_search_by_comp_id(_cid, top_n=top_n) if _cid else []
                if not hits and not _cid:
                    error = "NaFM similarity needs a compound in the database (search by name or THEO_ id, or a SMILES that matches a corpus compound). For novel structures, use Morgan or ChemBERTa."
            else:
                hits = sim_engine.tanimoto_search(query_smiles, top_n=top_n, threshold=threshold)
            if not hits:
                error = "Invalid SMILES or no similar compounds found."
            else:
                comp_ids = [h["comp_id"] for h in hits]
                scores = {h["comp_id"]: h["tanimoto"] for h in hits}
                ph = ",".join(["%s"] * len(comp_ids))
                with get_db() as conn:
                    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                        cur.execute(f"SELECT * FROM compounds WHERE comp_id IN ({ph})", tuple(comp_ids))
                        db_rows = {r["comp_id"]: r for r in cur.fetchall()}
                for cid in comp_ids:
                    if cid in db_rows:
                        row = db_rows[cid]
                        row["tanimoto"] = scores[cid]
                        results.append(row)
    # Geographic analog: filter by region (multi-valued via compound_region_map).
    region_filter = request.args.get("region", "")
    if region_filter and results:
        _rids = [r.get("comp_id") for r in results if r.get("comp_id")]
        _in_region = set()
        if _rids:
            _ph = ",".join(["%s"]*len(_rids))
            with get_db() as _conn:
                with _conn.cursor() as _cur:
                    _cur.execute(
                        f"SELECT DISTINCT comp_id FROM compound_region_map "
                        f"WHERE comp_id IN ({_ph}) AND LOWER(macro_region)=LOWER(%s)",
                        tuple(_rids)+(region_filter,))
                    _in_region = {row[0] for row in _cur.fetchall()}
        results = [r for r in results if r.get("comp_id") in _in_region]
    # Apply license filter
    license_filter = request.args.get("license", "all")
    if license_filter == "commercial" and results:
        results = [r for r in results if r.get("license_tier") in ("CC BY 4.0", "CC0")]
    elif license_filter == "academic" and results:
        results = [r for r in results if r.get("license_tier") in ("CC BY 4.0", "CC0", "CC BY-NC 4.0")]
    # Apply post-filters
    region_filter = request.args.get("region", "")
    kingdom_filter = request.args.get("kingdom", "")
    class_filter = request.args.get("class_filter", "")
    if region_filter and results:
        results = [r for r in results if r.get("region") == region_filter]
    if kingdom_filter and results:
        results = [r for r in results if r.get("kingdom") == kingdom_filter]
    if class_filter and results:
        results = [r for r in results if class_filter.lower() in (r.get("np_class") or "").lower() or class_filter.lower() in (r.get("classyfire_superclass") or "").lower()]
    # Deduplicate stereo variants (keep highest similarity)
    seen_connectivity = {}
    deduped_results = []
    for r in results:
        ik = (r.get("inchikey") or "")[:14]
        if ik and ik in seen_connectivity:
            continue
        if ik:
            seen_connectivity[ik] = True
        deduped_results.append(r)
    results = deduped_results
    # Deduplicate stereo variants (keep highest similarity)
    seen_connectivity = {}
    deduped_results = []
    for r in results:
        ik = (r.get("inchikey") or "")[:14]
        if ik and ik in seen_connectivity:
            continue
        if ik:
            seen_connectivity[ik] = True
        deduped_results.append(r)
    results = deduped_results
    metric = request.args.get("metric", "morgan")
    # --- three-widget data (map + linear tree + circular tree), scoped to similarity results ---
    sim_comp_ids = [r["comp_id"] for r in results if r.get("comp_id")]
    sim_thumb = None
    sim_linear_tree = []
    sim_region_counts = []
    sim_region_css = ""
    sim_region_titles = {}
    if sim_comp_ids:
        try:
            with get_db() as _sc:
                with _sc.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as _scur:
                    _scur.execute("""WITH result_set AS (SELECT unnest(%s::text[]) AS comp_id)
                        SELECT rt.kingdom AS kingdom, COUNT(DISTINCT rs.comp_id) AS cnt
                        FROM result_set rs JOIN resolved_taxonomy rt ON rt.comp_id = rs.comp_id
                        WHERE rt.kingdom IS NOT NULL AND rt.kingdom != ''
                        GROUP BY 1 ORDER BY cnt DESC""", (sim_comp_ids,))
                    _sk = _scur.fetchall()
                    _scur.execute("""WITH result_set AS (SELECT unnest(%s::text[]) AS comp_id)
                        SELECT c.region AS region, COUNT(*) AS count
                        FROM result_set rs JOIN compounds c ON c.comp_id = rs.comp_id
                        WHERE c.region IS NOT NULL AND c.region != ''
                        GROUP BY 1 ORDER BY count DESC""", (sim_comp_ids,))
                    sim_region_counts = _scur.fetchall()
            sim_thumb = kingdom_thumbnail_svg(_sk, total=len(sim_comp_ids), size=150, title="Result-set kingdoms")
            sim_linear_tree = compute_linear_tree(where_sql="c.comp_id = ANY(%s)", where_params=[sim_comp_ids])
            sim_region_css, sim_region_titles = build_region_color_css(sim_region_counts)
        except Exception:
            sim_thumb = None
    return render_template("similarity.html", query_smiles=query_smiles, query_input=query_input,
                           query_comp_id=query_comp_id,
                           results=results, metric=metric,
                           region_filter=region_filter, kingdom_filter=kingdom_filter, class_filter=class_filter,
                           top_n=top_n, threshold=threshold, error=error,
                           engine_loaded=sim_engine.loaded,
                           thumb=sim_thumb, linear_tree=sim_linear_tree,
                           region_css=sim_region_css, region_counts_filtered=sim_region_counts,
                           region_titles=sim_region_titles,
                           sim_comp_ids_str=",".join(sim_comp_ids))

@app.route("/api/similarity")
def api_similarity():
    raw = (request.args.get("smiles") or request.args.get("q") or "").strip()
    top_n = min(200, max(1, int(request.args.get("top_n", 50))))
    threshold = max(0.0, min(1.0, float(request.args.get("threshold", "0.3"))))
    if not raw:
        return jsonify({"error": "smiles or q parameter required"}), 400
    if not sim_engine.loaded:
        return jsonify({"error": "similarity search not available"}), 503
    smiles = raw
    if not any(c in raw for c in "()=#@[]"):
        with get_db() as conn:
            with conn.cursor() as cur:
                if raw.startswith("THEO_"):
                    cur.execute("SELECT smiles FROM compounds WHERE comp_id=%s", (raw,))
                else:
                    cur.execute("SELECT smiles FROM compounds WHERE LOWER(name)=%s LIMIT 1", (raw.lower(),))
                row = cur.fetchone()
                if row: smiles = row[0]
    metric = request.args.get("metric", "morgan")
    if metric == "chemberta":
        hits = sim_engine.chemberta_search(smiles, top_n=top_n)
    elif metric == "maccs":
        hits = sim_engine.maccs_search(smiles, top_n=top_n, threshold=threshold)
    elif metric == "nafm":
        _cid = _resolve_comp_id_for_nafm(request.args.get("q","").strip(), smiles)
        hits = sim_engine.nafm_search_by_comp_id(_cid, top_n=top_n) if _cid else []
    else:
        hits = sim_engine.tanimoto_search(smiles, top_n=top_n, threshold=threshold)
    comp_ids = [h["comp_id"] for h in hits]
    scores = {h["comp_id"]: h["tanimoto"] for h in hits}
    if not comp_ids:
        return jsonify({"count": 0, "results": []})
    ph = ",".join(["%s"] * len(comp_ids))
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(f"SELECT comp_id,name,smiles,inchikey,kingdom,source_db,source_organism,region,mw FROM compounds WHERE comp_id IN ({ph})", tuple(comp_ids))
            db_rows = {r["comp_id"]: r for r in cur.fetchall()}
    results = []
    seen = {}
    for cid in comp_ids:
        if cid in db_rows:
            r = db_rows[cid]
            ik14 = (r.get("inchikey") or "")[:14]
            if ik14 and ik14 in seen:
                continue
            if ik14:
                seen[ik14] = True
            r["tanimoto"] = scores[cid]
            results.append(r)
    return jsonify({"count": len(results), "results": results})


@app.route("/substructure")
def substructure():
    query = request.args.get("smarts", "").strip()
    max_results = min(500, max(1, int(request.args.get("max_results", 100))))
    results = []
    error = None
    # If the query happens to be a valid SMILES of a corpus compound, capture
    # the matching comp_id so the Query structure card can link to its page.
    query_comp_id = None
    if query:
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(query)
            if mol is not None:
                ik = Chem.MolToInchiKey(mol)
                if ik:
                    with get_db() as conn:
                        with conn.cursor() as cur:
                            cur.execute("SELECT comp_id FROM compounds WHERE inchikey=%s LIMIT 1", (ik,))
                            row = cur.fetchone()
                            if row:
                                query_comp_id = row[0]
        except Exception:
            pass
    if query:
        if not sim_engine.loaded:
            error = "Substructure search not available."
        else:
            from rdkit import Chem
            hits = sim_engine.substructure_search(query, max_results=max_results)
            if not hits:
                error = "Invalid SMARTS/SMILES or no matches found."
            else:
                comp_ids = [h["comp_id"] for h in hits]
                ph = ",".join(["%s"] * len(comp_ids))
                with get_db() as conn:
                    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                        cur.execute(f"SELECT * FROM compounds WHERE comp_id IN ({ph})", tuple(comp_ids))
                        db_rows = {r["comp_id"]: r for r in cur.fetchall()}
                # Verify substructure match with RDKit
                query_mol = Chem.MolFromSmarts(query)
                if query_mol is None:
                    query_mol = Chem.MolFromSmiles(query)
                for cid in comp_ids:
                    if cid in db_rows:
                        row = db_rows[cid]
                        mol = Chem.MolFromSmiles(row["smiles"])
                        if mol and mol.HasSubstructMatch(query_mol):
                            results.append(row)
                            if len(results) >= max_results:
                                break
    return render_template("substructure.html", query=query,
                           query_comp_id=query_comp_id,
                           results=results,
                           max_results=max_results, error=error,
                           engine_loaded=sim_engine.loaded)


@app.route("/scaffolds")
def scaffold_browser():
    page = max(1, int(request.args.get("page", 1)))
    scaffold = request.args.get("scaffold", "").strip()
    per_page = get_per_page()
    if scaffold:
        query = """SELECT c.* FROM compounds c JOIN scaffolds s ON c.comp_id=s.comp_id
                   WHERE s.scaffold=%s ORDER BY c.comp_id"""
        with get_db() as conn:
            results, total, pages = paginate(query, (scaffold,), page, per_page, conn)
        return render_template("scaffold_detail.html", scaffold=scaffold, results=results,
                               page=page, total=total, pages=pages, per_page=per_page)
    else:
        query = """SELECT scaffold, COUNT(*) AS cnt FROM scaffolds
                   WHERE scaffold IS NOT NULL AND scaffold != ''
                   GROUP BY scaffold ORDER BY cnt DESC"""
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute("SELECT COUNT(DISTINCT scaffold) FROM scaffolds WHERE scaffold != ''")
                total_scaffolds = cur.fetchone()["count"]
                cur.execute(query + " LIMIT %s OFFSET %s", (per_page, (page-1)*per_page))
                scaffolds = cur.fetchall()
                pages = max(1, -(-total_scaffolds // per_page))
        return render_template("scaffolds.html", scaffolds=scaffolds, total=total_scaffolds,
                               page=page, pages=pages, per_page=per_page)


@app.route("/admet")
def admet_browser():
    comp_id = request.args.get("comp_id", "").strip()
    if comp_id:
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute("SELECT * FROM admet WHERE comp_id=%s", (comp_id,))
                admet_data = cur.fetchone()
                cur.execute("SELECT * FROM compounds WHERE comp_id=%s", (comp_id,))
                compound = cur.fetchone()
        return render_template("admet_detail.html", compound=compound, admet=admet_data)
    # Filter mode
    filters = {}
    clauses = []
    params = []
    filter_defs = [
        ("hERG", "hERG cardiotoxicity", 0, 1),
        ("AMES", "AMES mutagenicity", 0, 1),
        ("BBB_Martins", "BBB penetration", 0, 1),
        ("HIA_Hou", "Human intestinal absorption", 0, 1),
        ("Caco2_Wang", "Caco-2 permeability", -8, -4),
        ("Solubility_AqSolDB", "Aqueous solubility", -10, 2),
        ("DILI", "DILI risk", 0, 1),
        ("Bioavailability_Ma", "Bioavailability", 0, 1),
        ("Pgp_Broccatelli", "P-gp inhibition", 0, 1),
        ("ClinTox", "Clinical toxicity", 0, 1),
        ("Lipophilicity_AstraZeneca", "Lipophilicity", -5, 5),
        ("VDss_Lombardo", "Volume of distribution", 0, 10),
    ]
    for col, label, default_min, default_max in filter_defs:
        lo = request.args.get(f"{col}_min", "")
        hi = request.args.get(f"{col}_max", "")
        if lo:
            clauses.append(f'a."{col}" >= %s')
            params.append(float(lo))
            filters[f"{col}_min"] = lo
        if hi:
            clauses.append(f'a."{col}" <= %s')
            params.append(float(hi))
            filters[f"{col}_max"] = hi
    page = max(1, int(request.args.get("page", 1)))
    per_page = get_per_page()
    results = []
    total = 0
    pages = 0
    if clauses:
        where = "WHERE " + " AND ".join(clauses)
        count_sql = f"SELECT COUNT(*) FROM admet a {where}"
        query_sql = f"""SELECT c.*, a.* FROM compounds c JOIN admet a ON c.comp_id=a.comp_id
                        {where} ORDER BY c.comp_id LIMIT %s OFFSET %s"""
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute(count_sql, tuple(params))
                total = cur.fetchone()["count"]
                pages = max(1, -(-total // per_page))
                cur.execute(query_sql, tuple(params) + (per_page, (page-1)*per_page))
                results = cur.fetchall()
    return render_template("admet.html", results=results, filters=filters,
                           filter_defs=filter_defs, total=total, page=page,
                           pages=pages, per_page=per_page)



@app.route("/api/usage")
def api_usage():
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT COUNT(*) AS total FROM access_log")
            total = cur.fetchone()["total"]
            cur.execute("SELECT COUNT(*) AS today FROM access_log WHERE ts::date = CURRENT_DATE")
            today = cur.fetchone()["today"]
            cur.execute("SELECT COUNT(DISTINCT ip) AS unique_ips FROM access_log")
            unique = cur.fetchone()["unique_ips"]
            cur.execute("SELECT path, COUNT(*) AS cnt FROM access_log WHERE path NOT LIKE '/static/%%' GROUP BY path ORDER BY cnt DESC LIMIT 10")
            top_pages = cur.fetchall()
    return jsonify({"total_requests": total, "today": today, "unique_visitors": unique, "top_pages": top_pages})



@app.route("/advanced")
def advanced_search():
    """Advanced search with AND logic across multiple filters + property ranges."""
    clauses, params = [], []
    filters = {}
    # Text filters (up to 5 AND conditions)
    for i in range(5):
        field = request.args.get(f"field_{i}", "")
        value = request.args.get(f"value_{i}", "").strip()
        if not field or not value:
            continue
        filters[f"field_{i}"] = field
        filters[f"value_{i}"] = value
        field_map = {
            "name": "LOWER(name) LIKE %s",
            "organism": "source_organism ILIKE %s",
            "kingdom": "EXISTS(SELECT 1 FROM resolved_taxonomy rt2 WHERE rt2.comp_id = compounds.comp_id AND (LOWER(rt2.kingdom) = LOWER(%s) OR LOWER(%s) = ANY(rt2.secondary_kingdoms)))",
            "region": "EXISTS(SELECT 1 FROM compound_region_map crm WHERE crm.comp_id = compounds.comp_id AND LOWER(crm.macro_region) = LOWER(%s))",
            "source": "LOWER(source_db) = LOWER(%s)",
            "class": "(np_class ILIKE %s OR classyfire_superclass ILIKE %s)",
            "pathway": "np_pathway ILIKE %s",
        }
        sql = field_map.get(field)
        if sql:
            if field == "pathway":
                clauses.append(sql)
                params.append(f"%{value}%")
            elif field == "class":
                clauses.append(sql)
                params.extend([f"%{value}%", f"%{value}%"])
            elif field == "kingdom":
                clauses.append(sql)
                params.extend([value, value])
            elif field in ("region", "source"):
                clauses.append(sql)
                params.append(value)
            else:
                clauses.append(sql)
                params.append(f"%{value}%")
    # License filter
    license_filter = request.args.get("license", "all")
    if license_filter == "commercial":
        clauses.append("tier_rank <= 1")
    elif license_filter == "academic":
        clauses.append("tier_rank <= 4")
    # Property range filters
    for prop in ["mw", "logp", "tpsa", "hba", "hbd", "n_rings", "rotatable_bonds"]:
        lo = request.args.get(f"{prop}_min", "")
        hi = request.args.get(f"{prop}_max", "")
        if lo:
            clauses.append(f"{prop} >= %s")
            params.append(float(lo))
            filters[f"{prop}_min"] = lo
        if hi:
            clauses.append(f"{prop} <= %s")
            params.append(float(hi))
            filters[f"{prop}_max"] = hi
    page = max(1, int(request.args.get("page", 1)))
    per_page = get_per_page()
    sort, order, oc = get_sort()
    results = []
    total = 0
    pages = 0
    if clauses:
        where = "WHERE " + " AND ".join(clauses)
        query = f"SELECT * FROM compounds {where} {oc}"
        with get_db() as conn:
            results, total, pages = paginate(query, tuple(params), page, per_page, conn)
    import os, json
    hist_path = os.path.join("static", "histograms.json")
    histograms = {}
    if os.path.exists(hist_path):
        with open(hist_path) as hf:
            histograms = json.load(hf)
    return render_template("advanced.html", results=results, filters=filters,
                           total=total, page=page, pages=pages, per_page=per_page,
                           sort=sort, order=order, histograms=histograms)



_filter_options_cache = {}

@app.route("/api/filter_options")
def api_filter_options():
    """Return available values for each filter type. Class values are split
    by ontology: npclassifier_class (NPClassifier-derived, atomic class names
    from np_class and inferred_class), classyfire_class (ClassyFire
    superclass), plus a unified 'class' for backward compatibility.
    Cached in module-level dict; restart service to refresh."""
    if _filter_options_cache:
        return jsonify(_filter_options_cache)
    with get_db() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT DISTINCT kingdom FROM compounds WHERE kingdom IS NOT NULL AND kingdom != '' ORDER BY kingdom")
            kingdoms = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT region FROM compounds WHERE region IS NOT NULL AND region != '' AND region != 'global' ORDER BY region")
            regions = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT source_db FROM compounds WHERE source_db IS NOT NULL ORDER BY source_db")
            sources = [r[0] for r in cur.fetchall()]
            # NPClassifier classes: atomic (split on ' $ ') from np_class.
            cur.execute("SELECT DISTINCT TRIM(c) AS cls FROM compounds, regexp_split_to_table(np_class, ' *[|$] *') AS c WHERE np_class IS NOT NULL AND np_class != '' AND TRIM(c) != '' ORDER BY cls")
            npc_classes = [r[0] for r in cur.fetchall()]
            # ClassyFire superclasses: stored as plain single values, not dollar-split.
            cur.execute("SELECT DISTINCT classyfire_superclass FROM compounds WHERE classyfire_superclass IS NOT NULL AND classyfire_superclass != '' ORDER BY classyfire_superclass")
            cf_classes = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT genus FROM compound_taxonomy WHERE genus IS NOT NULL ORDER BY genus")
            genera = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT family FROM compound_taxonomy WHERE family IS NOT NULL ORDER BY family")
            families = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT taxorder FROM resolved_taxonomy WHERE taxorder IS NOT NULL AND taxorder != '' ORDER BY taxorder")
            orders = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT taxclass FROM resolved_taxonomy WHERE taxclass IS NOT NULL AND taxclass != '' ORDER BY taxclass")
            tax_classes = [r[0] for r in cur.fetchall()]
            # APG clades now live in taxclass (APG fix used overwrite-in-place),
            # so 'clade' mirrors tax_class rather than a separate taxclass column.
            clades = tax_classes
            cur.execute("SELECT DISTINCT phylum FROM resolved_taxonomy WHERE phylum IS NOT NULL AND phylum != '' ORDER BY phylum")
            phyla = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT TRIM(sc) AS scls FROM compounds, regexp_split_to_table(np_superclass, ' *[|$] *') AS sc WHERE np_superclass IS NOT NULL AND np_superclass != '' AND TRIM(sc) != '' ORDER BY scls")
            npc_superclasses = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT TRIM(pw) AS pwy FROM compounds, regexp_split_to_table(np_pathway, ' *[|$] *') AS pw WHERE np_pathway IS NOT NULL AND np_pathway != '' AND TRIM(pw) != '' ORDER BY pwy")
            np_pathways = [r[0] for r in cur.fetchall()]
    # Union of both for backward-compatible 'class' key (deduplicated).
    classes = sorted(set(npc_classes) | set(cf_classes))
    result = {"kingdom": kingdoms, "region": regions, "source": sources,
              "class": classes,
              "npclassifier_class": npc_classes,
              "npclassifier_superclass": npc_superclasses,
              "classyfire_class": cf_classes,
              "pathway": np_pathways,
              "genus": genera, "family": families,
              "order": orders, "tax_class": tax_classes, "clade": clades, "phylum": phyla}
    _filter_options_cache.update(result)
    return jsonify(result)



@app.route("/api/depict")
def api_depict():
    """Render 2D structure as SVG via RDKit."""
    smiles = request.args.get("smiles", "")
    w = int(request.args.get("w", 300))
    h = int(request.args.get("h", 200))
    if not smiles:
        return "", 404
    from rdkit import Chem
    from rdkit.Chem import Draw
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "", 404
    from rdkit.Chem.Draw import rdMolDraw2D
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 0.09 if w < 150 else 0.2
    opts.minFontSize = 5 if w < 150 else 8
    opts.maxFontSize = 8 if w < 150 else 12
    opts.padding = 0.05
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    import re
    svg = re.sub(r"<rect[^>]*fill:#FFFFFF[^<]*</rect>", "", svg, flags=re.DOTALL)
    svg = re.sub(r"<rect[^>]*fill:#FFFFFF[^/]*/>", "", svg)
    return svg, 200, {"Content-Type": "image/svg+xml"}



@app.route("/api/family_figure/<comp_id>")
def api_family_figure(comp_id):
    """Compose all stereoisomer family members for a compound into one SVG/PDF/PNG.

    Query params:
      scale  : float, 0.5-3.0 (default 1.0), scales all dimensions and font sizes
      format : svg | pdf | png (default svg)
      download : if truthy, set Content-Disposition: attachment
    """
    import math, re
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    def _get_scale(name, default=1.0, lo=0.5, hi=3.0):
        try:
            v = float(request.args.get(name, str(default)))
        except (TypeError, ValueError):
            v = default
        return max(lo, min(hi, v))
    legacy = request.args.get("scale")
    if legacy is not None:
        try:
            legacy_v = max(0.5, min(3.0, float(legacy)))
        except (TypeError, ValueError):
            legacy_v = 1.0
        font_scale = legacy_v
        cell_scale = legacy_v
    else:
        font_scale = _get_scale("font_scale")
        cell_scale = _get_scale("cell_scale")
    scale = 1.0
    fmt = (request.args.get("format") or "svg").lower()
    if fmt not in ("svg", "pdf", "png"):
        fmt = "svg"
    embed = request.args.get("embed", "").lower() in ("1", "true", "yes")
    if embed:
        # Embed defaults: larger fonts, wider side margin so labels don't truncate
        if request.args.get("font_scale") is None and legacy is None:
            font_scale = 1.6
        if request.args.get("cell_scale") is None and legacy is None:
            cell_scale = 1.0
    try:
        name_chars = int(request.args.get("name_chars", "30"))
    except (TypeError, ValueError):
        name_chars = 30
    name_chars = max(8, min(80, name_chars))
    title_case = request.args.get("title_case", "").lower() in ("1", "true", "yes")
    try:
        mol_font_scale = float(request.args.get("mol_font_scale", str(font_scale)))
    except (TypeError, ValueError):
        mol_font_scale = font_scale
    mol_font_scale = max(0.3, min(3.0, mol_font_scale))

    with get_db() as conn_ff:
        with conn_ff.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur_ff:
            cur_ff.execute("SELECT inchikey FROM compounds WHERE comp_id=%s", (comp_id,))
            row = cur_ff.fetchone()
            if not row or not row.get("inchikey"):
                return "", 404
            ik_prefix = row["inchikey"][:14]
            cur_ff.execute("""SELECT comp_id, name, smiles, inchikey
                              FROM compounds WHERE SUBSTRING(inchikey,1,14)=%s
                              ORDER BY (name IS NULL), comp_id""", (ik_prefix,))
            members = cur_ff.fetchall()
    if len(members) < 2:
        return "", 404
    members.sort(key=lambda m: 0 if m["comp_id"] == comp_id else 1)
    ref = members[0]
    sats = members[1:]
    N = len(sats)

    # apply scale uniformly to every geometric/typographic constant
    W = 800
    H = 880
    title_h = 80
    cx = W / 2
    cy = title_h + (W / 2)
    sat_w = max(100 * cell_scale, min(160 * cell_scale, (200 - N * 6) * cell_scale))
    sat_h = max(80 * cell_scale, min(130 * cell_scale, (165 - N * 5) * cell_scale))
    ref_w, ref_h = sat_w, sat_h
    radius = (W / 2) - max(160, 130) / 2 - (90 if embed else 30)

    def render_mol(smiles, mw, mh):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        drawer = rdMolDraw2D.MolDraw2DSVG(int(mw), int(mh))
        opts = drawer.drawOptions()
        opts.bondLineWidth = 1.6
        opts.minFontSize = int(10 * mol_font_scale)
        opts.maxFontSize = int(16 * mol_font_scale)
        opts.padding = 0.04
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg = re.sub(r"<rect[^>]*fill:#FFFFFF[^<]*</rect>", "", svg, flags=re.DOTALL)
        svg = re.sub(r"<rect[^>]*fill:#FFFFFF[^/]*/>", "", svg)
        m_inner = re.search(r"<svg[^>]*>(.*)</svg>", svg, flags=re.DOTALL)
        return m_inner.group(1) if m_inner else None

    ref_name = ref["name"] or ref["comp_id"]
    if len(ref_name) > 60:
        ref_name = ref_name[:58] + "..."

    fs_title = 22 * font_scale
    fs_subtitle = 18 * font_scale
    fs_footer = 9 * font_scale
    fs_label = 11 * font_scale

    if embed:
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.1f} {H:.1f}" width="{W:.1f}" height="{H:.1f}" font-family="sans-serif">',
        ]
    else:
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.1f} {H:.1f}" width="{W:.1f}" height="{H:.1f}" font-family="sans-serif">',
            f'<text x="32" y="36" text-anchor="start" font-size="{fs_title:.1f}" font-weight="600" fill="#222">Stereoisomer family</text>',
            f'<text x="32" y="62" text-anchor="start" font-size="{fs_subtitle:.1f}" fill="#444">{ref_name}</text>',
            f'<text x="{W-12:.0f}" y="{H-10:.0f}" text-anchor="end" font-size="{fs_footer:.1f}" fill="#aaa">THEOBROMA {VERSION_DISPLAY}</text>'
        ]
    line_start = radius / 3
    line_end = 2 * radius / 3
    for i, sat in enumerate(sats):
        ang = (2 * math.pi * i / N) - math.pi / 2
        x1 = cx + line_start * math.cos(ang)
        y1 = cy + line_start * math.sin(ang)
        x2 = cx + line_end * math.cos(ang)
        y2 = cy + line_end * math.sin(ang)
        parts.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#888" stroke-width="1" opacity="0.5"/>')
    def emit_cell(m, x, y, w, h):
        inner = render_mol(m["smiles"], w, h - 20 * font_scale) or ""
        parts.append(f'<g data-comp-id="{m["comp_id"]}" transform="translate({x:.1f},{y:.1f})">')
        parts.append(f'<svg viewBox="0 0 {int(w)} {int(h-20*font_scale)}" width="{int(w)}" height="{int(h-20*font_scale)}">{inner}</svg>')
        name = m["name"] or m["comp_id"]
        if name and not any(c.islower() for c in name):
            name = name[0].upper() + name[1:].lower()
        if title_case and name and name[0].islower():
            name = name[0].upper() + name[1:]
        if len(name) > name_chars:
            name = name[:name_chars-2] + "..."
        parts.append(f'<text x="{w/2:.1f}" y="{h-5*font_scale:.1f}" text-anchor="middle" font-size="{fs_label:.1f}" fill="#222">{name}</text>')
        parts.append('</g>')
    emit_cell(ref, cx - ref_w / 2, cy - ref_h / 2, ref_w, ref_h)
    for i, sat in enumerate(sats):
        ang = (2 * math.pi * i / N) - math.pi / 2
        x = cx + radius * math.cos(ang) - sat_w / 2
        y = cy + radius * math.sin(ang) - sat_h / 2
        emit_cell(sat, x, y, sat_w, sat_h)
    parts.append("</svg>")
    svg_text = "".join(parts)

    download = request.args.get("download")

    if fmt == "svg":
        headers = {"Content-Type": "image/svg+xml"}
        if download:
            headers["Content-Disposition"] = f'attachment; filename="theobroma_family_{comp_id}.svg"'
        return svg_text, 200, headers

    if fmt == "pdf":
        import cairosvg
        from io import BytesIO
        buf = BytesIO()
        try:
            cairosvg.svg2pdf(bytestring=svg_text.encode("utf-8"), write_to=buf)
        except Exception as e:
            return jsonify({"error": "PDF conversion failed: " + str(e)}), 500
        resp = Response(buf.getvalue(), mimetype="application/pdf")
        resp.headers["Content-Disposition"] = f'attachment; filename="theobroma_family_{comp_id}.pdf"'
        return resp

    # png
    import cairosvg
    from io import BytesIO
    buf = BytesIO()
    try:
        cairosvg.svg2png(bytestring=svg_text.encode("utf-8"),
                         write_to=buf, output_width=int(W * 3))
    except Exception as e:
        return jsonify({"error": "PNG conversion failed: " + str(e)}), 500
    resp = Response(buf.getvalue(), mimetype="image/png")
    resp.headers["Content-Disposition"] = f'attachment; filename="theobroma_family_{comp_id}.png"'
    return resp


@app.route("/api/bulk")
def api_bulk():
    """Bulk export of compound subsets with selectable columns.
    Streams as CSV for downloads of any size.
    Params: cols (comma-separated), tier (open/nc/all), limit (default unlimited)
    Example: /api/bulk?cols=comp_id,smiles
    """
    # Default to a richer field set when cols is unspecified, matching the
    # manuscript Section 5 description of "full annotation payloads."
    default_cols = "comp_id,name,smiles,inchikey,kingdom,source_db,all_sources,source_organism,region,license_tier,np_class,classyfire_superclass,mw,logp"
    cols_param = request.args.get("cols", default_cols)
    # Accept both 'tier' and 'license' as parameter names. 'license' matches
    # the /search and /api/search convention; 'tier' is the original /api/bulk
    # parameter. Map the /search-style values (commercial, academic) to the
    # /api/bulk-style values (open, nc) for consistency.
    license_param = request.args.get("license", "")
    license_to_tier = {"commercial": "open", "academic": "nc", "all": "all"}
    tier = license_to_tier.get(license_param, request.args.get("tier", "all"))
    limit = request.args.get("limit", "")

    allowed_cols = {
        "comp_id", "name", "smiles", "inchi", "inchikey",
        "source_db", "kingdom", "region", "source_organism",
        "mw", "logp", "tpsa", "hba", "hbd", "n_rings", "rotatable_bonds",
        "license_tier", "all_sources", "np_class", "classyfire_superclass",
        "inferred_class", "reference_doi", "trad_medicine", "trust_score"
    }
    cols = [c.strip() for c in cols_param.split(",") if c.strip() in allowed_cols]
    if not cols:
        return "No valid columns requested.", 400

    tier_map = {
        "open": ["CC BY 4.0", "CC0"],
        "nc":   ["CC BY 4.0", "CC0", "CC BY-NC 4.0"],
        "all":  None,
    }
    allowed_tier = tier_map.get(tier)

    col_sql = ",".join(cols)
    where = ""
    params = ()
    if allowed_tier:
        ph = ",".join(["%s"] * len(allowed_tier))
        where = f"WHERE license_tier IN ({ph})"
        params = tuple(allowed_tier)
    limit_sql = f"LIMIT {int(limit)}" if limit.isdigit() else ""

    fmt = request.args.get("format", "csv").lower()
    def stream_csv():
        import csv, io
        buf = io.StringIO()
        w = csv.writer(buf)
        w.writerow(cols)
        yield buf.getvalue()
        buf.seek(0); buf.truncate()
        with get_db() as conn:
            with conn.cursor(name="bulk_cursor") as cur:
                cur.itersize = 10000
                cur.execute(f"SELECT {col_sql} FROM compounds {where} {limit_sql}", params)
                for row in cur:
                    w.writerow(["" if v is None else str(v) for v in row])
                    if buf.tell() > 65536:
                        yield buf.getvalue()
                        buf.seek(0); buf.truncate()
                yield buf.getvalue()
    def stream_json():
        yield "["
        first = True
        with get_db() as conn:
            with conn.cursor(name="bulk_cursor_json") as cur:
                cur.itersize = 10000
                cur.execute(f"SELECT {col_sql} FROM compounds {where} {limit_sql}", params)
                for row in cur:
                    obj = {k: (None if v is None else str(v)) for k, v in zip(cols, row)}
                    prefix = "" if first else ","
                    first = False
                    yield prefix + json.dumps(obj)
        yield "]"
    if fmt == "json":
        return Response(stream_json(), mimetype="application/json",
                        headers={"Content-Disposition": f"attachment; filename=theobroma_bulk_{tier}.json"})
    return Response(stream_csv(), mimetype="text/csv",
                    headers={"Content-Disposition": f"attachment; filename=theobroma_bulk_{tier}.csv"})

@app.errorhandler(404)
def not_found(e):
    return render_template("404.html"), 404


@app.route("/api/docs")
def api_docs():
    return send_from_directory("static", "openapi.yaml", mimetype="text/yaml")

@app.route("/api")
def api_index():
    return """<!DOCTYPE html><html><head><title>THEOBROMA API</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/swagger-ui-dist/swagger-ui.css">
    </head><body><div id="swagger-ui"></div>
    <script src="https://cdn.jsdelivr.net/npm/swagger-ui-dist/swagger-ui-bundle.js"></script>
    <script>SwaggerUIBundle({url:"/api/docs",dom_id:"#swagger-ui"})</script>
    </body></html>"""


@app.route("/api/stereoisomers/<comp_id>")
def api_stereoisomers(comp_id):
    conn = get_db()
    cur = conn.cursor()
    cur.execute("SELECT inchikey FROM compounds WHERE comp_id=%s", (comp_id,))
    row = cur.fetchone()
    if not row:
        return jsonify({"error": "not found"}), 404
    ik_prefix = row[0][:14]
    cur.execute("""SELECT comp_id, name, smiles, inchikey, source_db, np_class, trust_score
                   FROM compounds WHERE SUBSTRING(inchikey,1,14)=%s ORDER BY comp_id""", (ik_prefix,))
    cols = [d[0] for d in cur.description]
    results = [dict(zip(cols, r)) for r in cur.fetchall()]
    return jsonify({"inchikey_prefix": ik_prefix, "count": len(results), "stereoisomers": results})


@app.route("/annotate")
def annotate_page():
    return render_template("annotate.html")

@app.route("/api/annotate", methods=["POST"])
def api_annotate():
    """Batch lookup of SMILES/InChIKey inputs against the corpus. Returns matched
    rows with full annotation and an unmatched list. Cap 1000 inputs per call;
    larger batches are chunked client-side."""
    data = request.get_json(force=True, silent=True) or {}
    inputs = data.get("inputs", [])
    if not isinstance(inputs, list) or not inputs:
        return jsonify({"error": "inputs must be a non-empty list"}), 400
    if len(inputs) > 1000:
        return jsonify({"error": "max 1000 inputs per call"}), 400
    from rdkit import Chem
    IK_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
    resolved = []
    for i, item in enumerate(inputs):
        rid = item.get("id", i)
        ik = (item.get("inchikey") or "").strip().upper()
        smi = (item.get("smiles") or "").strip()
        if ik and IK_RE.match(ik):
            resolved.append({"id": rid, "in_smi": smi, "in_ik": ik, "prefix": ik[:14], "err": None})
        elif smi:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                resolved.append({"id": rid, "in_smi": smi, "in_ik": "", "prefix": None, "err": "invalid SMILES"})
                continue
            ik_c = Chem.MolToInchiKey(mol)
            resolved.append({"id": rid, "in_smi": smi, "in_ik": ik_c, "prefix": ik_c[:14], "err": None})
        else:
            resolved.append({"id": rid, "in_smi": "", "in_ik": "", "prefix": None, "err": "no smiles or inchikey"})
    prefixes = sorted({r["prefix"] for r in resolved if r["prefix"]})
    by_prefix, tax_map = {}, {}
    if prefixes:
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute("""SELECT comp_id, name, smiles, inchikey, source_db, all_sources,
                                      kingdom, region, source_organism, mw, logp, tpsa, hba, hbd,
                                      n_rings, rotatable_bonds, license_tier, np_class, np_superclass,
                                      np_pathway, classyfire_superclass, inferred_class,
                                      inferred_confidence, trad_medicine, trust_score,
                                      novelty_morgan, novelty_maccs, sa_score
                               FROM compounds
                               WHERE SUBSTRING(inchikey,1,14) = ANY(%s::text[])""", (prefixes,))
                for row in cur:
                    by_prefix.setdefault(row["inchikey"][:14], []).append(dict(row))
                comp_ids = [c["comp_id"] for fam in by_prefix.values() for c in fam]
                if comp_ids:
                    cur.execute("""SELECT ct.comp_id,
                                          string_agg(DISTINCT rt.phylum, ', ') AS phyla,
                                          string_agg(DISTINCT rt.taxclass, ', ') AS tax_classes,
                                          string_agg(DISTINCT rt.taxorder, ', ') AS orders,
                                          string_agg(DISTINCT ct.family, ', ') AS families,
                                          string_agg(DISTINCT ct.genus, ', ') AS genera
                                   FROM compound_taxonomy ct
                                   LEFT JOIN resolved_taxonomy rt ON rt.comp_id = ct.comp_id
                                   WHERE ct.comp_id = ANY(%s::text[])
                                   GROUP BY ct.comp_id""", (comp_ids,))
                    for row in cur:
                        tax_map[row["comp_id"]] = {k: row[k] for k in
                                                    ("phyla","tax_classes","orders","families","genera")}
    matched, unmatched = [], []
    for r in resolved:
        if not r["prefix"]:
            unmatched.append({"id": r["id"], "input_smiles": r["in_smi"],
                              "input_inchikey": r["in_ik"], "reason": r["err"]})
            continue
        fam = by_prefix.get(r["prefix"], [])
        if not fam:
            unmatched.append({"id": r["id"], "input_smiles": r["in_smi"],
                              "input_inchikey": r["in_ik"], "reason": "no match in corpus"})
            continue
        exact = next((c for c in fam if c["inchikey"] == r["in_ik"]), None)
        chosen = exact if exact else fam[0]
        out = {"id": r["id"], "input_smiles": r["in_smi"], "input_inchikey": r["in_ik"],
               "match_type": "exact" if exact else "family",
               "family_size": len(fam),
               "family_members": [{"comp_id": c["comp_id"], "name": c["name"], "inchikey": c["inchikey"]}
                                  for c in fam if c["comp_id"] != chosen["comp_id"]]}
        out.update(chosen)
        out.update(tax_map.get(chosen["comp_id"], {}))
        matched.append(out)
    return jsonify({"matched": matched, "unmatched": unmatched})


@app.route("/sources")
def sources():
    return redirect(url_for("statistics"))


@app.route("/robots.txt")
def robots_txt():
    return send_from_directory("static", "robots.txt", mimetype="text/plain")

@app.route("/admin/refresh_cache")
def admin_refresh_cache():
    """Reload the browse-options cache. Localhost-only so SSH-tunneled admin
    requests work but external traffic cannot trigger cache rebuilds. Use after
    any bulk corpus mutation (INSERT, UPDATE, DELETE on compounds) that was
    applied without a service restart."""
    if request.remote_addr not in ("127.0.0.1", "::1"):
        abort(403)
    _load_browse_options()
    return jsonify({"ok": True,
                    "kingdoms": len(_BROWSE_CACHE["all_kingdoms"]),
                    "sources": len(_BROWSE_CACHE["all_sources_list"]),
                    "regions": len(_BROWSE_CACHE["all_regions"])})

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
