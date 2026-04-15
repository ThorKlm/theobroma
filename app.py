"""
THEOBROMA v3 — An open multi-kingdom natural products database
spanning 27 sources across six kingdoms and six continents.
"""
from flask import (Flask, render_template, request, send_from_directory,
                   jsonify, abort, redirect, url_for, Response)
from config import Config
import psycopg2, psycopg2.extras, os, math, re, csv, io
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scripts.similarity import SimilarityEngine

sim_engine = SimilarityEngine(vectors_dir="data/vectors")
sim_engine.load()

app = Flask(__name__)
app.config.from_object(Config)

SORTABLE = {"comp_id","name","kingdom","source_db","region","source_organism",
            "mw","logp","tpsa","hba","hbd","n_rings","rotatable_bonds","license_tier"}
ALLOWED_PER_PAGE = {25, 50, 100, 1000}

def get_db():
    return psycopg2.connect(app.config["DB_URI"])

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

@app.template_filter("region_label")
def region_label(v):
    return "global / unresolved" if not v or v in ("global","nan","") else v

@app.template_filter("kingdom_label")
def kingdom_label(v):
    return "unresolved" if not v or v in ("nan","") else v


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
    """Normalize search query: strip hyphens, handle alpha/greek equivalents."""
    q = q.lower().strip()
    # Greek letter equivalents
    replacements = [
        ("α", "alpha"), ("β", "beta"), ("γ", "gamma"), ("δ", "delta"),
        ("ε", "epsilon"), ("ω", "omega"), ("μ", "mu"),
    ]
    for greek, latin in replacements:
        q = q.replace(greek, latin)
    # Remove hyphens for fuzzy matching
    q = q.replace("-", "").replace("–", "").replace("—", "")
    return q

# --- Routes ---

@app.route("/")
def index():
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT COUNT(*) AS cnt FROM compounds")
            total = cur.fetchone()["cnt"]
            cur.execute("""SELECT kingdom, COUNT(*) AS cnt FROM compounds
                          WHERE kingdom IS NOT NULL AND kingdom!='' GROUP BY kingdom ORDER BY cnt DESC""")
            kingdoms = cur.fetchall()
            cur.execute("SELECT COUNT(DISTINCT source_db) AS cnt FROM compounds")
            n_sources = cur.fetchone()["cnt"]
            cur.execute(f"SELECT {REGION_SQL} AS reg, COUNT(*) AS cnt FROM compounds GROUP BY 1 ORDER BY cnt DESC")
            regions = cur.fetchall()
            n_regions = len([r for r in regions if r["reg"]!="unresolved"])
    return render_template("index.html", total=total, kingdoms=kingdoms,
                           n_sources=n_sources, n_regions=n_regions, regions=regions)

@app.route("/search")
def search():
    q = request.args.get("q","").strip()
    st = request.args.get("type","name")
    page = max(1, int(request.args.get("page",1)))
    per_page = get_per_page()
    sort, order, oc = get_sort()
    # Resolve "property" type to actual search type
    if st == "property":
        prop_type = request.args.get("prop_type", "class")
        if prop_type == "class":
            st = "class"
        else:
            st = "mw"  # handled by range filters
    has_extra = any(request.args.get(f"extra_type_{i}") and request.args.get(f"extra_q_{i}") for i in range(1,6))
    has_range = any(request.args.get(f"{p}_min") or request.args.get(f"{p}_max")
                    for p in ["mw","logp","tpsa","hba","hbd","n_rings","rotatable_bonds"])
    if not q and not has_extra and not has_range and st not in ("mw",):
        return render_template("search.html", results=[], query="", search_type=st,
                               page=1, total=0, pages=0, sort=sort, order=order, per_page=per_page)
    # InChIKey direct redirect
    if st == "inchikey" and re.match(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$', q):
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute("SELECT comp_id FROM compounds WHERE inchikey=%s LIMIT 1", (q,))
                row = cur.fetchone()
                if row:
                    return redirect(url_for("compound_detail", comp_id=row["comp_id"]))
    tq = {
        "name":    (f"""SELECT * FROM (
               SELECT DISTINCT ON (c.comp_id) c.*,
               CASE WHEN LOWER(c.name) = %s THEN 0
                    WHEN LOWER(c.name) LIKE %s THEN 1
                    WHEN LOWER(c.name) LIKE %s THEN 2
                    ELSE 3 END AS relevance
               FROM compounds c
               LEFT JOIN compound_synonyms s ON c.inchikey = s.inchikey
               WHERE LOWER(c.name) LIKE %s OR LOWER(s.synonym) LIKE %s
               OR REPLACE(REPLACE(LOWER(c.name),'-',''),' ','') LIKE %s
               OR REPLACE(REPLACE(LOWER(s.synonym),'-',''),' ','') LIKE %s
             ) sub ORDER BY relevance, LENGTH(name), name""", (q.lower(), f"{q.lower()}%", f"% {q.lower()}%", f"%{q.lower()}%", f"%{q.lower()}%", f"%{normalize_query(q)}%", f"%{normalize_query(q)}%")),
        "smiles":  (f"SELECT * FROM compounds WHERE smiles=%s {oc}", (q,)),
        "inchikey":(f"SELECT * FROM compounds WHERE inchikey=%s {oc}", (q,)),
        "source":  (f"SELECT * FROM compounds WHERE source_db ILIKE %s {oc}", (f"%{q}%",)),
        "organism":(f"SELECT * FROM compounds WHERE source_organism ILIKE %s {oc}", (f"%{q}%",)),
        "region":  (f"SELECT * FROM compounds WHERE region ILIKE %s {oc}", (f"%{q}%",)),
        "kingdom": (f"SELECT * FROM compounds WHERE kingdom ILIKE %s {oc}", (f"%{q}%",)),
        "class":   (f"SELECT * FROM compounds WHERE np_class ILIKE %s OR classyfire_superclass ILIKE %s {oc}", (f"%{q}%", f"%{q}%")),
        "mw":      (f"SELECT * FROM compounds WHERE 1=1 {oc}", ()),
    }
    query, params = tq.get(st, tq["name"])
    # Handle extra AND filters
    extra_clauses = []
    extra_params = []
    for i in range(1, 6):
        et = request.args.get(f"extra_type_{i}", "")
        eq = request.args.get(f"extra_q_{i}", "").strip()
        if not et or not eq:
            continue
        emap = {
            "name": "LOWER(name) LIKE %s",
            "organism": "source_organism ILIKE %s",
            "kingdom": "kingdom ILIKE %s",
            "region": "region ILIKE %s",
            "source": "source_db ILIKE %s",
            "class": "(np_class ILIKE %s OR classyfire_superclass ILIKE %s)",
        }
        esql = emap.get(et)
        if esql:
            if et == "class":
                extra_clauses.append(esql)
                extra_params.extend([f"%{eq}%", f"%{eq}%"])
            else:
                extra_clauses.append(esql)
                extra_params.append(f"%{eq}%")
    # Handle property range filters
    for prop in ["mw", "logp", "tpsa", "hba", "hbd", "n_rings", "rotatable_bonds"]:
        pmin = request.args.get(f"{prop}_min", "")
        pmax = request.args.get(f"{prop}_max", "")
        if pmin:
            extra_clauses.append(f"{prop} >= %s")
            extra_params.append(float(pmin))
        if pmax:
            extra_clauses.append(f"{prop} <= %s")
            extra_params.append(float(pmax))
    if extra_clauses:
        extra_where = " AND ".join(extra_clauses)
        query = f"SELECT * FROM ({query}) AS base WHERE {extra_where}"
        params = params + tuple(extra_params)
    with get_db() as conn:
        results, total, pages = paginate(query, params, page, per_page, conn)
    return render_template("search.html", results=results, query=q, search_type=st,
                           page=page, total=total, pages=pages, sort=sort, order=order, per_page=per_page)

@app.route("/browse")
def browse():
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
        clauses.append("kingdom=%s"); params += (kingdom,)
    if source:
        clauses.append("(source_db=%s OR all_sources LIKE %s)"); params += (source, f"%{source}%")
    if region:
        if region in ("unresolved", "global / unresolved"):
            clauses.append("(region IS NULL OR region='' OR region='global')")
        else:
            clauses.append("region=%s"); params += (region,)
    if named:
        clauses.append("name IS NOT NULL AND name != ''")
    if chem_class:
        clauses.append("(np_class ILIKE %s OR classyfire_superclass ILIKE %s)")
        params += (f"%{chem_class}%", f"%{chem_class}%")
    where = "WHERE "+" AND ".join(clauses) if clauses else ""
    query = f"SELECT * FROM compounds {where} {oc}"
    with get_db() as conn:
        results, total, pages = paginate(query, params, page, per_page, conn)
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""SELECT kingdom, COUNT(*) AS cnt FROM compounds
                          WHERE kingdom IS NOT NULL AND kingdom!='' GROUP BY kingdom ORDER BY kingdom""")
            all_kingdoms = cur.fetchall()
            # Count each source by appearances in all_sources (full credit to smaller DBs)
            cur.execute("""
                        SELECT TRIM(src) AS source_db, COUNT(*) AS cnt
                        FROM (SELECT unnest(string_to_array(COALESCE(NULLIF(all_sources, ''), source_db), '|')) AS src
                              FROM compounds) t
                        WHERE TRIM(src) != '' AND LENGTH(TRIM(src)) > 0
                        GROUP BY TRIM (src)
                        ORDER BY cnt DESC
                        """)
            raw_sources = cur.fetchall()

            # For COCONUT, show only exclusive compounds (not in any other source)
            cur.execute("""
                        SELECT COUNT(*) AS cnt
                        FROM compounds
                        WHERE (all_sources IS NULL OR all_sources = '' OR all_sources = 'COCONUT'
                            OR all_sources NOT LIKE '%%|%%')
                          AND (source_db = 'COCONUT' OR source_db IS NULL)
                        """)
            coconut_exclusive = cur.fetchone()["cnt"]

            sources = []
            for s in raw_sources:
                if s["source_db"] == "COCONUT":
                    sources.append({"source_db": "COCONUT (exclusive)", "cnt": coconut_exclusive})
                else:
                    sources.append(s)
            sources.sort(key=lambda x: -x["cnt"])
            cur.execute("SELECT source_db, COUNT(*) AS cnt FROM compounds GROUP BY source_db ORDER BY source_db")
            all_sources = cur.fetchall()
            cur.execute(f"SELECT DISTINCT {REGION_SQL} AS reg FROM compounds ORDER BY reg")
            all_regions = [r["reg"] for r in cur.fetchall()]
    return render_template("browse.html", results=results, page=page, total=total,
                           pages=pages, kingdom=kingdom, source=source, region=region,
                           all_kingdoms=all_kingdoms,
                           all_sources_list=all_sources,
                           all_regions=all_regions, sort=sort, order=order,
                           named=named, per_page=per_page, chem_class=chem_class)

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
    with get_db() as conn3:
        with conn3.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur3:
            try:
                cur3.execute("SELECT * FROM admet WHERE comp_id=%s", (c["comp_id"],))
                admet_data = cur3.fetchone()
            except:
                conn3.rollback()
    return render_template("compound.html", c=c, all_sources_list=src_list, synonyms=synonyms, admet=admet_data)

@app.route("/statistics")
def statistics():
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT COUNT(*) AS cnt FROM compounds")
            total = cur.fetchone()["cnt"]
            cur.execute("SELECT kingdom, COUNT(*) AS cnt FROM compounds WHERE kingdom IS NOT NULL AND kingdom!='' GROUP BY kingdom ORDER BY cnt DESC")
            kingdoms = cur.fetchall()
            cur.execute("SELECT source_db, COUNT(*) AS cnt FROM compounds GROUP BY source_db ORDER BY cnt DESC")
            sources = cur.fetchall()
            cur.execute(f"SELECT {REGION_SQL} AS region, COUNT(*) AS cnt FROM compounds GROUP BY 1 ORDER BY cnt DESC LIMIT 30")
            regions = cur.fetchall()
            cur.execute("SELECT AVG(mw) AS avg_mw, AVG(logp) AS avg_logp, AVG(tpsa) AS avg_tpsa, AVG(hba) AS avg_hba, AVG(hbd) AS avg_hbd FROM compounds WHERE mw IS NOT NULL")
            prop_stats = cur.fetchone()
            cur.execute("SELECT license_tier, COUNT(*) AS cnt FROM compounds GROUP BY license_tier ORDER BY cnt DESC")
            licenses = cur.fetchall()
            cur.execute("SELECT COUNT(*) AS cnt FROM compounds WHERE all_sources LIKE '%%|%%'")
            multi_source = cur.fetchone()["cnt"]
    return render_template("statistics.html", total=total, kingdoms=kingdoms,
                           sources=sources, regions=regions, prop_stats=prop_stats,
                           licenses=licenses, multi_source=multi_source)

@app.route("/download")
def download_page():
    d = app.config["DATA_DIR"]
    files = []
    for fn in ["theobroma_final.csv","theobroma_all.sdf"]:
        p = os.path.join(d, fn)
        if os.path.exists(p):
            files.append({"name":fn, "fmt":fn.split(".")[-1].upper(), "size":f"{os.path.getsize(p)/1024/1024:.1f} MB"})
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
    if st == "name":
        nq = normalize_query(q)
        base = f"""SELECT * FROM (
              SELECT DISTINCT ON (c.comp_id) c.comp_id,c.name,c.smiles,c.inchikey,c.kingdom,c.source_db,c.all_sources,c.source_organism,c.region,c.mw,c.logp,c.license_tier,
              CASE WHEN LOWER(c.name) = %s THEN 0
                   WHEN LOWER(c.name) LIKE %s THEN 1
                   WHEN LOWER(c.name) LIKE %s THEN 2
                   ELSE 3 END AS relevance
              FROM compounds c LEFT JOIN compound_synonyms s ON c.inchikey=s.inchikey
              WHERE LOWER(c.name) LIKE %s OR LOWER(s.synonym) LIKE %s
              OR REPLACE(REPLACE(LOWER(c.name),'-',''),' ','') LIKE %s
              OR REPLACE(REPLACE(LOWER(s.synonym),'-',''),' ','') LIKE %s
            ) sub ORDER BY relevance, LENGTH(name), name"""
        pm_tuple = (q.lower(), f"{q.lower()}%", f"% {q.lower()}%", f"%{q.lower()}%", f"%{q.lower()}%", f"%{nq}%", f"%{nq}%")
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                cur.execute(f"SELECT COUNT(*) FROM ({base}) sq", pm_tuple)
                total = cur.fetchone()["count"]
                cur.execute(base + " LIMIT %s OFFSET %s", pm_tuple + (limit, offset))
                results = cur.fetchall()
    else:
        tc = {"smiles":"smiles=%s","inchikey":"inchikey=%s",
              "kingdom":"kingdom ILIKE %s","organism":"source_organism ILIKE %s",
              "region":"region ILIKE %s","source":"source_db ILIKE %s",
              "class":"np_class ILIKE %s OR classyfire_superclass ILIKE %s"}
        cl = tc.get(st, "LOWER(name) LIKE %s")
        if st == "class":
            pm = (f"%{q}%", f"%{q}%")
        else:
            pm = f"%{q.lower()}%" if st == "organism" else (f"%{q}%" if st in ("kingdom","region","source") else q)
        with get_db() as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
                params = pm if isinstance(pm, tuple) else (pm,)
                cur.execute(f"SELECT COUNT(*) FROM compounds WHERE {cl}", params)
                total = cur.fetchone()["count"]
                cur.execute(f"SELECT {cols} FROM compounds WHERE {cl} LIMIT %s OFFSET %s", params + (limit, offset))
                results = cur.fetchall()
    if fmt == "csv":
        si = io.StringIO()
        if results:
            w = csv.DictWriter(si, fieldnames=results[0].keys())
            w.writeheader()
            w.writerows(results)
        return Response(si.getvalue(), mimetype="text/csv",
                       headers={"Content-Disposition": f"attachment; filename=theobroma_{st}_{q[:20]}.csv"})
    return jsonify({"count":len(results), "total":total, "offset":offset, "limit":limit, "results":results})

@app.route("/api/compound/<comp_id>")
def api_compound(comp_id):
    """Get full compound details by comp_id."""
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT * FROM compounds WHERE comp_id=%s", (comp_id,))
            c = cur.fetchone()
    if not c: return jsonify({"error":"compound not found"}), 404
    return jsonify(c)

@app.route("/api/autocomplete")
def api_autocomplete():
    """Name autocomplete for search box."""
    q = request.args.get("q","").strip()
    if len(q) < 2: return jsonify([])
    with get_db() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT DISTINCT name FROM compounds WHERE LOWER(name) LIKE %s AND name IS NOT NULL AND name != '' ORDER BY name LIMIT 12",
                        (f"{q.lower()}%",))
            results = [r[0] for r in cur.fetchall()]
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
    return jsonify({"total":total, "n_sources":n_sources, "kingdoms":kingdoms, "sources":sources})

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
                "source":"source_db ILIKE %s","organism":"source_organism ILIKE %s",
                "region":"region ILIKE %s","kingdom":"kingdom ILIKE %s"}
        cl = tmap.get(st, tmap["name"])
        pm = f"%{q.lower()}%" if st in ("name","organism") else (f"%{q}%" if st in ("kingdom","region","source") else q)
        clauses.append(cl); params += (pm,)
    if kingdom:
        clauses.append("kingdom=%s"); params += (kingdom,)
    if source:
        clauses.append("(source_db=%s OR all_sources LIKE %s)"); params += (source, f"%{source}%")
    if region and region != "unresolved":
        clauses.append("region=%s"); params += (region,)
    elif region == "unresolved":
        clauses.append("(region IS NULL OR region='' OR region='global')")
    if named:
        clauses.append("name IS NOT NULL AND name != ''")
    where = "WHERE "+" AND ".join(clauses) if clauses else ""
    cols = "comp_id,name,smiles,inchikey,kingdom,source_db,all_sources,source_organism,region,mw,logp,tpsa,hba,hbd,license_tier"
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(f"SELECT {cols} FROM compounds {where} LIMIT %s", params + (limit,))
            results = cur.fetchall()
    si = io.StringIO()
    if results:
        w = csv.DictWriter(si, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)
    return Response(si.getvalue(), mimetype="text/csv",
                   headers={"Content-Disposition": "attachment; filename=theobroma_export.csv"})

@app.route("/similarity")
def similarity():
    query_input = request.args.get("smiles", "").strip()
    query_smiles = query_input
    # Resolve name or comp_id to SMILES
    if query_input and not any(c in query_input for c in "()=#@[]"):
        with get_db() as conn:
            with conn.cursor() as cur:
                # Try comp_id first
                if query_input.startswith("THEO_"):
                    cur.execute("SELECT smiles FROM compounds WHERE comp_id=%s", (query_input,))
                    row = cur.fetchone()
                    if row: query_smiles = row[0]
                else:
                    # Try name search
                    cur.execute("SELECT smiles FROM compounds WHERE LOWER(name)=%s LIMIT 1", (query_input.lower(),))
                    row = cur.fetchone()
                    if row: query_smiles = row[0]
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
    # Geographic analog: filter by region
    region_filter = request.args.get("region", "")
    if region_filter and results:
        results = [r for r in results if r.get("region") == region_filter]
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
    return render_template("similarity.html", query_smiles=query_smiles, query_input=query_input, results=results, metric=metric,
                           region_filter=region_filter, kingdom_filter=kingdom_filter, class_filter=class_filter,
                           top_n=top_n, threshold=threshold, error=error,
                           engine_loaded=sim_engine.loaded)

@app.route("/api/similarity")
def api_similarity():
    raw = request.args.get("smiles", "").strip()
    top_n = min(200, max(1, int(request.args.get("top_n", 50))))
    threshold = max(0.0, min(1.0, float(request.args.get("threshold", "0.3"))))
    if not raw:
        return jsonify({"error": "smiles parameter required"}), 400
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
    return render_template("substructure.html", query=query, results=results,
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
        ("hERG_Karim-et-al", "hERG risk", 0, 1),
        ("AMES_Li-et-al", "AMES mutagenicity", 0, 1),
        ("BBB_Martins-et-al", "BBB penetration", 0, 1),
        ("HIA_Hou-et-al", "Human intestinal absorption", 0, 1),
        ("Caco2_Wang-et-al", "Caco-2 permeability", -8, -4),
        ("Solubility_AqSolDB", "Aqueous solubility", -10, 2),
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
            "kingdom": "kingdom ILIKE %s",
            "region": "region ILIKE %s",
            "source": "source_db ILIKE %s",
            "class": "(np_class ILIKE %s OR classyfire_superclass ILIKE %s)",
        }
        sql = field_map.get(field)
        if sql:
            if field == "class":
                clauses.append(sql)
                params.extend([f"%{value}%", f"%{value}%"])
            else:
                clauses.append(sql)
                params.append(f"%{value}%" if field != "kingdom" else f"%{value}%")
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



@app.route("/api/filter_options")
def api_filter_options():
    """Return available values for each filter type."""
    with get_db() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT DISTINCT kingdom FROM compounds WHERE kingdom IS NOT NULL AND kingdom != '' ORDER BY kingdom")
            kingdoms = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT region FROM compounds WHERE region IS NOT NULL AND region != '' AND region != 'global' ORDER BY region")
            regions = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT source_db FROM compounds WHERE source_db IS NOT NULL ORDER BY source_db")
            sources = [r[0] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT np_class FROM compounds WHERE np_class IS NOT NULL AND np_class != '' ORDER BY np_class LIMIT 100")
            classes = [r[0] for r in cur.fetchall()]
    return jsonify({"kingdom": kingdoms, "region": regions, "source": sources, "class": classes})


@app.errorhandler(404)
def not_found(e):
    return render_template("404.html"), 404

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
