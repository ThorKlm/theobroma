"""
THEOBROMA — An open multi-kingdom natural products database
spanning 21+ sources across five kingdoms and six continents
for virtual screening and drug discovery.
"""
from flask import (Flask, render_template, request, send_from_directory,
                   jsonify, abort)
from config import Config
import psycopg2, psycopg2.extras, os, math

app = Flask(__name__)
app.config.from_object(Config)

SORTABLE = {"comp_id","name","kingdom","source_db","region","source_organism",
            "mw","logp","tpsa","hba","hbd","n_rings","rotatable_bonds","license_tier"}

def get_db():
    return psycopg2.connect(app.config["DB_URI"])

def get_sort(default="comp_id"):
    s = request.args.get("sort", default)
    o = request.args.get("order", "asc").lower()
    if s not in SORTABLE: s = default
    if o not in ("asc","desc"): o = "asc"
    return s, o, f"ORDER BY {s} {o} NULLS LAST"

def paginate(query, params, page, per_page, conn):
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        cur.execute(f"SELECT COUNT(*) FROM ({query}) AS subq", params)
        total = cur.fetchone()["count"]
        pages = max(1, math.ceil(total / per_page))
        cur.execute(query + " LIMIT %s OFFSET %s", params + (per_page, (page-1)*per_page))
        results = cur.fetchall()
    return results, total, pages

REGION_SQL = """CASE WHEN region IS NULL OR region='' OR region='global' THEN 'unresolved' ELSE region END"""

@app.template_filter("region_label")
def region_label(v):
    return "unresolved" if not v or v in ("global","nan","") else v

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
    sort, order, oc = get_sort()
    if not q:
        return render_template("search.html", results=[], query="", search_type=st,
                               page=1, total=0, pages=0, sort=sort, order=order)
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
        results, total, pages = paginate(query, params, page, Config.PER_PAGE, conn)
    return render_template("search.html", results=results, query=q, search_type=st,
                           page=page, total=total, pages=pages, sort=sort, order=order)

@app.route("/browse")
def browse():
    page = max(1, int(request.args.get("page",1)))
    kingdom = request.args.get("kingdom","")
    source = request.args.get("source","")
    region = request.args.get("region","")
    sort, order, oc = get_sort()
    clauses, params = [], ()
    if kingdom:
        clauses.append("kingdom=%s"); params += (kingdom,)
    if source:
        clauses.append("source_db=%s"); params += (source,)
    if region:
        if region == "unresolved":
            clauses.append("(region IS NULL OR region='' OR region='global')")
        else:
            clauses.append("region=%s"); params += (region,)
    where = "WHERE "+" AND ".join(clauses) if clauses else ""
    query = f"SELECT * FROM compounds {where} {oc}"
    with get_db() as conn:
        results, total, pages = paginate(query, params, page, Config.PER_PAGE, conn)
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT DISTINCT kingdom FROM compounds WHERE kingdom IS NOT NULL AND kingdom!='' ORDER BY kingdom")
            all_kingdoms = [r["kingdom"] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT source_db FROM compounds ORDER BY source_db")
            all_sources = [r["source_db"] for r in cur.fetchall()]
            cur.execute(f"SELECT DISTINCT {REGION_SQL} AS reg FROM compounds ORDER BY reg")
            all_regions = [r["reg"] for r in cur.fetchall()]
    return render_template("browse.html", results=results, page=page, total=total,
                           pages=pages, kingdom=kingdom, source=source, region=region,
                           all_kingdoms=all_kingdoms, all_sources=all_sources,
                           all_regions=all_regions, sort=sort, order=order)

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

@app.route("/api/search")
def api_search():
    q = request.args.get("q","").strip()
    st = request.args.get("type","name")
    limit = min(1000, max(1, int(request.args.get("limit",50))))
    if not q: return jsonify({"error":"query parameter 'q' required"}), 400
    tc = {"name":"LOWER(name) LIKE %s","smiles":"smiles=%s","inchikey":"inchikey=%s",
          "kingdom":"kingdom ILIKE %s","organism":"source_organism ILIKE %s"}
    cl = tc.get(st, tc["name"])
    pm = f"%{q.lower()}%" if st in ("name","organism") else (f"%{q}%" if st=="kingdom" else q)
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(f"SELECT comp_id,name,smiles,inchikey,kingdom,source_db,all_sources,source_organism,region,license_tier FROM compounds WHERE {cl} LIMIT %s", (pm, limit))
            results = cur.fetchall()
    return jsonify({"count":len(results), "results":results})

@app.errorhandler(404)
def not_found(e):
    return render_template("404.html"), 404

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
