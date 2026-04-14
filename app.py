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
    if not q:
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
        "name":    (f"SELECT * FROM compounds WHERE LOWER(name) LIKE %s {oc}", (f"%{q.lower()}%",)),
        "smiles":  (f"SELECT * FROM compounds WHERE smiles=%s {oc}", (q,)),
        "inchikey":(f"SELECT * FROM compounds WHERE inchikey=%s {oc}", (q,)),
        "source":  (f"SELECT * FROM compounds WHERE source_db ILIKE %s {oc}", (f"%{q}%",)),
        "organism":(f"SELECT * FROM compounds WHERE source_organism ILIKE %s {oc}", (f"%{q}%",)),
        "region":  (f"SELECT * FROM compounds WHERE region ILIKE %s {oc}", (f"%{q}%",)),
        "kingdom": (f"SELECT * FROM compounds WHERE kingdom ILIKE %s {oc}", (f"%{q}%",)),
    }
    query, params = tq.get(st, tq["name"])
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
    sort, order, oc = get_sort()
    clauses, params = [], ()
    if kingdom:
        clauses.append("kingdom=%s"); params += (kingdom,)
    if source:
        clauses.append("source_db=%s"); params += (source,)
    if region:
        if region in ("unresolved", "global / unresolved"):
            clauses.append("(region IS NULL OR region='' OR region='global')")
        else:
            clauses.append("region=%s"); params += (region,)
    if named:
        clauses.append("name IS NOT NULL AND name != ''")
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
            cur.execute(f"SELECT DISTINCT {REGION_SQL} AS reg FROM compounds ORDER BY reg")
            all_regions = [r["reg"] for r in cur.fetchall()]
    return render_template("browse.html", results=results, page=page, total=total,
                           pages=pages, kingdom=kingdom, source=source, region=region,
                           all_kingdoms=all_kingdoms,
                           all_sources_list=all_sources_list,
                           all_regions=all_regions, sort=sort, order=order,
                           named=named, per_page=per_page)

@app.route("/compound/<comp_id>")
def compound_detail(comp_id):
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT * FROM compounds WHERE comp_id=%s", (comp_id,))
            c = cur.fetchone()
    if not c: abort(404)
    src_list = [s.strip() for s in (c.get("all_sources") or "").split("|") if s.strip()]
    return render_template("compound.html", c=c, all_sources_list=src_list)

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
    tc = {"name":"LOWER(name) LIKE %s","smiles":"smiles=%s","inchikey":"inchikey=%s",
          "kingdom":"kingdom ILIKE %s","organism":"source_organism ILIKE %s",
          "region":"region ILIKE %s","source":"source_db ILIKE %s"}
    cl = tc.get(st, tc["name"])
    pm = f"%{q.lower()}%" if st in ("name","organism") else (f"%{q}%" if st in ("kingdom","region","source") else q)
    cols = "comp_id,name,smiles,inchikey,kingdom,source_db,all_sources,source_organism,region,mw,logp,license_tier"
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(f"SELECT COUNT(*) FROM compounds WHERE {cl}", (pm,))
            total = cur.fetchone()["count"]
            cur.execute(f"SELECT {cols} FROM compounds WHERE {cl} LIMIT %s OFFSET %s", (pm, limit, offset))
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
        clauses.append("source_db=%s"); params += (source,)
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
    query_smiles = request.args.get("smiles", "").strip()
    top_n = min(200, max(1, int(request.args.get("top_n", 50))))
    threshold = max(0.0, min(1.0, float(request.args.get("threshold", "0.3"))))
    results = []
    error = None
    if query_smiles:
        if not sim_engine.loaded:
            error = "Similarity search not available (vectors not loaded)."
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
    return render_template("similarity.html", query_smiles=query_smiles, results=results,
                           top_n=top_n, threshold=threshold, error=error,
                           engine_loaded=sim_engine.loaded)

@app.route("/api/similarity")
def api_similarity():
    smiles = request.args.get("smiles", "").strip()
    top_n = min(200, max(1, int(request.args.get("top_n", 50))))
    threshold = max(0.0, min(1.0, float(request.args.get("threshold", "0.3"))))
    if not smiles:
        return jsonify({"error": "smiles parameter required"}), 400
    if not sim_engine.loaded:
        return jsonify({"error": "similarity search not available"}), 503
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
    for cid in comp_ids:
        if cid in db_rows:
            r = db_rows[cid]
            r["tanimoto"] = scores[cid]
            results.append(r)
    return jsonify({"count": len(results), "results": results})

@app.errorhandler(404)
def not_found(e):
    return render_template("404.html"), 404

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
