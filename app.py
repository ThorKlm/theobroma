"""
THEOBROMA — An open global multi-kingdom natural products database
for virtual screening and drug discovery.

Flask web application serving ~960K compounds from 21+ sources
across five kingdoms and six continents.
"""
from flask import (Flask, render_template, request, send_from_directory,
                   jsonify, abort, redirect, url_for)
from config import Config
import psycopg2
import psycopg2.extras
import os
import math

app = Flask(__name__)
app.config.from_object(Config)


def get_db():
    return psycopg2.connect(app.config["DB_URI"])


def paginate(query, params, page, per_page, conn):
    """Execute a paginated query, returning (results, total, pages)."""
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        count_q = f"SELECT COUNT(*) FROM ({query}) AS subq"
        cur.execute(count_q, params)
        total = cur.fetchone()["count"]
        pages = max(1, math.ceil(total / per_page))
        offset = (page - 1) * per_page
        cur.execute(query + " LIMIT %s OFFSET %s", params + (per_page, offset))
        results = cur.fetchall()
    return results, total, pages


# --- Routes ---

@app.route("/")
def index():
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT COUNT(*) AS cnt FROM compounds")
            total = cur.fetchone()["cnt"]
            cur.execute("""SELECT kingdom, COUNT(*) AS cnt FROM compounds
                          WHERE kingdom IS NOT NULL AND kingdom != ''
                          GROUP BY kingdom ORDER BY cnt DESC""")
            kingdoms = cur.fetchall()
            cur.execute("SELECT COUNT(DISTINCT source_db) AS cnt FROM compounds")
            n_sources = cur.fetchone()["cnt"]
            cur.execute("SELECT COUNT(DISTINCT region) AS cnt FROM compounds WHERE region IS NOT NULL AND region != ''")
            n_regions = cur.fetchone()["cnt"]
    return render_template("index.html", total=total, kingdoms=kingdoms,
                           n_sources=n_sources, n_regions=n_regions)


@app.route("/search")
def search():
    q = request.args.get("q", "").strip()
    search_type = request.args.get("type", "name")
    page = max(1, int(request.args.get("page", 1)))
    if not q:
        return render_template("search.html", results=[], query="",
                               search_type=search_type, page=1, total=0, pages=0)
    type_to_query = {
        "name": ("SELECT * FROM compounds WHERE LOWER(name) LIKE %s ORDER BY name", (f"%{q.lower()}%",)),
        "smiles": ("SELECT * FROM compounds WHERE smiles = %s ORDER BY comp_id", (q,)),
        "inchikey": ("SELECT * FROM compounds WHERE inchikey = %s ORDER BY comp_id", (q,)),
        "source": ("SELECT * FROM compounds WHERE source_db ILIKE %s ORDER BY comp_id", (f"%{q}%",)),
        "organism": ("SELECT * FROM compounds WHERE source_organism ILIKE %s ORDER BY name", (f"%{q}%",)),
        "region": ("SELECT * FROM compounds WHERE region ILIKE %s ORDER BY name", (f"%{q}%",)),
        "kingdom": ("SELECT * FROM compounds WHERE kingdom ILIKE %s ORDER BY name", (f"%{q}%",)),
    }
    query, params = type_to_query.get(search_type, type_to_query["name"])
    with get_db() as conn:
        results, total, pages = paginate(query, params, page, Config.PER_PAGE, conn)
    return render_template("search.html", results=results, query=q,
                           search_type=search_type, page=page, total=total, pages=pages)


@app.route("/browse")
def browse():
    page = max(1, int(request.args.get("page", 1)))
    kingdom = request.args.get("kingdom", "")
    source = request.args.get("source", "")
    clauses, params = [], ()
    if kingdom:
        clauses.append("kingdom = %s")
        params += (kingdom,)
    if source:
        clauses.append("source_db = %s")
        params += (source,)
    where = "WHERE " + " AND ".join(clauses) if clauses else ""
    query = f"SELECT * FROM compounds {where} ORDER BY comp_id"
    with get_db() as conn:
        results, total, pages = paginate(query, params, page, Config.PER_PAGE, conn)
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT DISTINCT kingdom FROM compounds WHERE kingdom IS NOT NULL AND kingdom != '' ORDER BY kingdom")
            all_kingdoms = [r["kingdom"] for r in cur.fetchall()]
            cur.execute("SELECT DISTINCT source_db FROM compounds ORDER BY source_db")
            all_sources = [r["source_db"] for r in cur.fetchall()]
    return render_template("browse.html", results=results, page=page, total=total,
                           pages=pages, kingdom=kingdom, source=source,
                           all_kingdoms=all_kingdoms, all_sources=all_sources)


@app.route("/compound/<comp_id>")
def compound_detail(comp_id):
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT * FROM compounds WHERE comp_id = %s", (comp_id,))
            compound = cur.fetchone()
    if not compound:
        abort(404)
    return render_template("compound.html", c=compound)


@app.route("/statistics")
def statistics():
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT COUNT(*) AS cnt FROM compounds")
            total = cur.fetchone()["cnt"]
            cur.execute("""SELECT kingdom, COUNT(*) AS cnt FROM compounds
                          WHERE kingdom IS NOT NULL AND kingdom != ''
                          GROUP BY kingdom ORDER BY cnt DESC""")
            kingdoms = cur.fetchall()
            cur.execute("""SELECT source_db, COUNT(*) AS cnt FROM compounds
                          GROUP BY source_db ORDER BY cnt DESC""")
            sources = cur.fetchall()
            cur.execute("""SELECT region, COUNT(*) AS cnt FROM compounds
                          WHERE region IS NOT NULL AND region != ''
                          GROUP BY region ORDER BY cnt DESC LIMIT 30""")
            regions = cur.fetchall()
            cur.execute("""SELECT
                AVG(mw) AS avg_mw, AVG(logp) AS avg_logp, AVG(tpsa) AS avg_tpsa,
                AVG(hba) AS avg_hba, AVG(hbd) AS avg_hbd
                FROM compounds WHERE mw IS NOT NULL""")
            prop_stats = cur.fetchone()
    return render_template("statistics.html", total=total, kingdoms=kingdoms,
                           sources=sources, regions=regions, prop_stats=prop_stats)


@app.route("/download")
def download_page():
    data_dir = app.config["DATA_DIR"]
    files = []
    for fmt in ["csv", "sdf"]:
        fname = f"theobroma_all.{fmt}"
        path = os.path.join(data_dir, fname)
        if os.path.exists(path):
            size_mb = os.path.getsize(path) / (1024 * 1024)
            files.append({"name": fname, "fmt": fmt.upper(), "size": f"{size_mb:.1f} MB"})
    return render_template("download.html", files=files)


@app.route("/download/<filename>")
def download_file(filename):
    data_dir = app.config["DATA_DIR"]
    path = os.path.join(data_dir, filename)
    if not os.path.exists(path):
        abort(404)
    return send_from_directory(data_dir, filename, as_attachment=True)


@app.route("/help")
def help_page():
    return render_template("help.html")


@app.route("/api/search")
def api_search():
    """JSON API for programmatic access."""
    q = request.args.get("q", "").strip()
    search_type = request.args.get("type", "name")
    limit = min(1000, max(1, int(request.args.get("limit", 50))))
    if not q:
        return jsonify({"error": "query parameter 'q' required"}), 400
    type_to_col = {
        "name": "LOWER(name) LIKE %s",
        "smiles": "smiles = %s",
        "inchikey": "inchikey = %s",
        "kingdom": "kingdom ILIKE %s",
        "organism": "source_organism ILIKE %s",
    }
    clause = type_to_col.get(search_type, type_to_col["name"])
    param = f"%{q.lower()}%" if search_type in ("name", "organism") else (f"%{q}%" if search_type == "kingdom" else q)
    with get_db() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(f"SELECT comp_id, name, smiles, inchikey, kingdom, source_db, source_organism, region FROM compounds WHERE {clause} LIMIT %s",
                        (param, limit))
            results = cur.fetchall()
    return jsonify({"count": len(results), "results": results})


@app.errorhandler(404)
def not_found(e):
    return render_template("404.html"), 404


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
