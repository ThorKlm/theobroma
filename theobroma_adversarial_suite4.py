#!/usr/bin/env python3
"""THEOBROMA adversarial suite IV, 50 checks.

Suites I to III cover integrity, taxonomy repair, licensing, query semantics,
adversarial input, similarity, caches, exports, cross-view agreement and edges.
This one covers what the compound-page and kingdom-badge work opened up: the
agreement between what a filter matches and what a page displays, taxonomy
invariants at every rank, the attested-lineage layer, and the interface
behaviour that only shows up on specific compounds.

  PASS   invariant holds
  FAIL   regression, investigate before the deposit
  KNOWN  documented limitation or expected divergence, asserted so it cannot drift
  ERROR  the check itself could not run

Usage:  python3 theobroma_adversarial_suite4.py [--base http://localhost:5000]
Exit 1 if any FAIL or ERROR.
"""
import argparse, json, re, sys, time, traceback
import urllib.error, urllib.parse, urllib.request

import psycopg2

DSN = "host=localhost dbname=theobroma user=theobroma"
CORPUS = 1132805
KINGDOMS = ("plant", "animal", "fungi", "bacteria", "unresolved")

TESTS = []
COUNT = {"PASS": 0, "FAIL": 0, "KNOWN": 0, "ERROR": 0}


def test(tag, note):
    def deco(fn):
        TESTS.append((tag, note, fn))
        return fn
    return deco


def emit(tag, code, headline, detail=None):
    COUNT[code] += 1
    print("[%-5s] %-4s %s" % (code, tag, headline))
    for line in (detail or [])[:2]:
        if line:
            print("             %s" % line)


class Ctx:
    def __init__(self, con, base):
        self.con, self.base = con, base
        self.cur = con.cursor()

    def q(self, sql, args=None):
        self.cur.execute(sql, args or ())
        return self.cur.fetchall()

    def q1(self, sql, args=None):
        r = self.q(sql, args)
        return r[0][0] if r else None

    def get(self, path, timeout=60):
        try:
            t0 = time.time()
            with urllib.request.urlopen(self.base.rstrip("/") + path, timeout=timeout) as r:
                return r.status, r.read().decode("utf-8", "replace"), time.time() - t0
        except urllib.error.HTTPError as e:
            return e.code, "", 0.0
        except Exception as e:
            return 0, str(e), 0.0

    def js(self, path, timeout=60):
        st, body, dt = self.get(path, timeout)
        if st != 200:
            return None, st, dt
        try:
            return json.loads(body), st, dt
        except Exception:
            return None, st, dt

    def total(self, query):
        d, _, _ = self.js("/api/search?%s&limit=1" % query)
        return (d or {}).get("total")

    def page(self, cid):
        return self.get("/compound/%s" % cid)[1]


# ------------------------------------ R. filter matches vs page displays

@test("R01", "kingdom badge appears only where a secondary is attested")
def r01(c):
    with_sec = c.q1("""SELECT comp_id FROM resolved_taxonomy
                       WHERE array_length(secondary_kingdoms,1) >= 1 LIMIT 1""")
    without = c.q1("""SELECT comp_id FROM resolved_taxonomy
                      WHERE secondary_kingdoms IS NULL
                         OR array_length(secondary_kingdoms,1) IS NULL LIMIT 1""")
    a = c.page(with_sec).count("badge-secondary-att")
    b = c.page(without).count("badge-secondary-att")
    return a >= 1 and b == 0, "with secondary %s badges, without %s" % (a, b), \
        "the title chip must not badge compounds that have no secondary attestation"


@test("R02", "every attested secondary kingdom is rendered")
def r02(c):
    cid, sec = c.q("""SELECT comp_id, secondary_kingdoms FROM resolved_taxonomy
                      WHERE array_length(secondary_kingdoms,1) >= 2 LIMIT 1""")[0]
    html = c.page(cid)
    shown = [k for k in sec if ('>%s<' % k) in html or ('from %s"' % k) in html]
    return len(shown) == len(sec), "%s: db %s, rendered %s" % (cid, len(sec), len(shown)), \
        "a partially rendered array is worse than none, the page then contradicts the filter"


@test("R03", "a non-kingdom search still shows secondary kingdoms on the compound page")
def r03(c):
    cid = c.q1("""SELECT rt.comp_id FROM resolved_taxonomy rt
                  WHERE array_length(rt.secondary_kingdoms,1) >= 1 LIMIT 1""")
    html = c.page(cid)
    return "badge-secondary-att" in html, \
        "compound page badges secondaries independently of the search type", \
        "the compound page has no search context, so the badge must not depend on one"


@test("R04", "kingdom results mark only rows matched on a secondary")
def r04(c):
    st, html, _ = c.get("/search?type=kingdom&q=fungi")
    ids = list(dict.fromkeys(re.findall(r"/compound/(THEO_\d+)", html)))
    if not ids:
        return "ERROR", "no result rows parsed", ""
    expect = c.q1("""SELECT count(*) FROM resolved_taxonomy
                     WHERE comp_id = ANY(%s) AND lower(kingdom) <> 'fungi'
                       AND 'fungi' = ANY(secondary_kingdoms)""", (ids,))
    got = html.count("badge-secondary-att")
    return got == expect, "%s rows, %s badges, %s expected" % (len(ids), got, expect), \
        "a badge on a primary-matched row would misreport how the hit was made"


@test("R05", "a non-kingdom search shows no secondary-attestation marks in results")
def r05(c):
    st, html, _ = c.get("/search?type=npclassifier_class&q=Flavonols")
    return st == 200 and html.count("badge-secondary-att") == 0, \
        "class search secondary marks: %s" % html.count("badge-secondary-att"), \
        "the mark answers 'why did this row match', which only applies to a kingdom filter"


@test("R06", "the region row shows every mapped macro-region")
def r06(c):
    cid, n = c.q("""SELECT comp_id, count(*) FROM compound_region_map
                    GROUP BY 1 ORDER BY 2 DESC LIMIT 1""")[0]
    regs = [r[0] for r in c.q("""SELECT DISTINCT macro_region FROM compound_region_map
                                 WHERE comp_id = %s""", (cid,))]
    html = c.page(cid)
    missing = [r for r in regs if r not in html]
    return not missing, "%s: %s regions, %s absent from the page" % (cid, n, len(missing)), \
        "compounds.region is a single legacy value while the filter uses the map"


@test("R07", "attested-lineage block renders where more than one kingdom is attested")
def r07(c):
    cid = c.q1("""SELECT comp_id FROM resolved_taxonomy
                  WHERE array_length(secondary_kingdoms,1) >= 1 LIMIT 1""")
    return "Attested lineages" in c.page(cid), "%s shows the lineage block" % cid, \
        "without it the page shows a voted kingdom above a contradicting lineage"


@test("R08", "single-kingdom compounds do not carry the lineage block")
def r08(c):
    cid = c.q1("""SELECT comp_id FROM resolved_taxonomy
                  WHERE secondary_kingdoms IS NULL
                     OR array_length(secondary_kingdoms,1) IS NULL LIMIT 1""")
    return "Attested lineages" not in c.page(cid), "%s omits the lineage block" % cid, \
        "the block is only informative where the lineages differ"


# ------------------------------------------- S. taxonomy rank invariants

@test("S01", "one organism token resolves to one genus and family")
def s01(c):
    n = c.q1("""SELECT count(*) FROM (SELECT token FROM compound_taxonomy
                WHERE coalesce(token,'') <> '' GROUP BY token
                HAVING count(DISTINCT coalesce(genus,'')) > 1
                    OR count(DISTINCT coalesce(family,'')) > 1) x""")
    return n == 0, "tokens with inconsistent lineage %s" % n, \
        "organism to lineage is a fact about taxonomy, not about the compound"


@test("S02", "genera spanning more than one family are nomenclatural homonyms")
def s02(c):
    n = c.q1("""SELECT count(*) FROM (SELECT genus FROM compound_taxonomy
                WHERE coalesce(genus,'') <> '' AND coalesce(family,'') <> ''
                GROUP BY genus HAVING count(DISTINCT family) > 1) x""")
    return ("KNOWN" if n is not None and n <= 20 else "FAIL"), \
        "genera under several families %s" % n, \
        "ICZN and ICN are separate codes, so a plant and an animal may share a genus name"


@test("S03", "families spanning more than one order")
def s03(c):
    n = c.q1("""SELECT count(*) FROM (SELECT family FROM resolved_taxonomy
                WHERE coalesce(family,'') <> '' AND coalesce(taxorder,'') <> ''
                GROUP BY family HAVING count(DISTINCT taxorder) > 1) x""")
    return ("KNOWN" if n is not None and n <= 10 else "FAIL"), \
        "families under several orders %s" % n, \
        "downstream of the genus homonyms rather than independent"


@test("S04", "orders spanning more than one class")
def s04(c):
    n = c.q1("""SELECT count(*) FROM (SELECT taxorder FROM resolved_taxonomy
                WHERE coalesce(taxorder,'') <> '' AND coalesce(taxclass,'') <> ''
                GROUP BY taxorder HAVING count(DISTINCT taxclass) > 1) x""")
    return ("KNOWN" if n is not None and n <= 26 else "FAIL"), \
        "orders under several classes %s" % n, \
        "the eudicots and magnoliids rank-vocabulary question, deliberately not decided"


@test("S07", "the compound page renders no unnormalized kingdom name")
def s07(c):
    cid = c.q1("""SELECT ct.comp_id FROM compound_taxonomy ct
                  WHERE lower(ct.kingdom_any) IN ('bacillati','animalia','viridiplantae',
                                                  'archaea','chromista') LIMIT 1""")
    html = c.page(cid) if cid else ""
    bad = [w for w in ("bacillati", "animalia", "viridiplantae", "chromista") if w in html]
    return not bad, "%s: unnormalized names on the page %s" % (cid, bad or "none"), \
        "an NCBI superkingdom name has no matching badge class and reads as a different kingdom"


@test("S08", "monocot orders stay under liliopsida")
def s08(c):
    n = c.q1("""SELECT count(*) FROM resolved_taxonomy
        WHERE taxorder IN ('acorales','alismatales','arecales','asparagales','commelinales',
                           'dioscoreales','liliales','pandanales','poales','zingiberales',
                           'petrosaviales') AND taxclass IS DISTINCT FROM 'liliopsida'""")
    return n == 0, "monocot rows not under liliopsida %s" % n


@test("S09", "one phylum name for land plants")
def s09(c):
    n = c.q1("SELECT count(*) FROM resolved_taxonomy WHERE phylum='tracheophyta'")
    return n == 0, "tracheophyta rows %s" % n


@test("S10", "every phylum resolves to a display kingdom")
def s10(c):
    r = c.q("""SELECT DISTINCT rt.phylum FROM resolved_taxonomy rt
               LEFT JOIN phylum_kingdom_map m ON m.phylum = rt.phylum
               WHERE coalesce(rt.phylum,'') <> '' AND m.phylum IS NULL""")
    return not r, "unmapped phyla %s" % len(r), ", ".join(x[0] for x in r[:5]) or "none"


@test("S11", "taxonomy case consistency")
def s11(c):
    r = c.q("""SELECT count(*) FILTER (WHERE taxclass <> lower(taxclass)),
                      count(*) FILTER (WHERE taxorder <> lower(taxorder)),
                      count(*) FILTER (WHERE phylum   <> lower(phylum))
               FROM resolved_taxonomy""")[0]
    return sum(r) == 0, "mixed-case class/order/phylum %s" % (tuple(r),)


@test("S12", "lineage kingdom against voted kingdom, phylum level")
def s12(c):
    n = c.q1("""SELECT count(*) FROM resolved_taxonomy rt
                JOIN phylum_kingdom_map m ON m.phylum = rt.phylum
                WHERE m.lineage_kingdom IS DISTINCT FROM rt.kingdom""")
    return ("KNOWN" if n is not None and 4000 < n < 9000 else "FAIL"), \
        "compounds where the two disagree %s" % n, \
        "the tree groups by lineage, Table S1 reports the vote"


@test("S13", "lineage kingdom against voted kingdom, organism level")
def s13(c):
    n = c.q1("""WITH m AS (
        SELECT ct.comp_id,
               CASE WHEN ct.ncbi_lineage::text ILIKE '%%Metazoa%%'       THEN 'animal'
                    WHEN ct.ncbi_lineage::text ILIKE '%%Fungi%%'         THEN 'fungi'
                    WHEN ct.ncbi_lineage::text ILIKE '%%Bacteria%%'      THEN 'bacteria'
                    WHEN ct.ncbi_lineage::text ILIKE '%%Viridiplantae%%' THEN 'plant' END AS lk,
               rt.kingdom AS vote
        FROM compound_taxonomy ct JOIN resolved_taxonomy rt ON rt.comp_id = ct.comp_id
        WHERE ct.ncbi_lineage IS NOT NULL)
        SELECT count(DISTINCT comp_id) FROM m WHERE lk IS NOT NULL AND lk IS DISTINCT FROM vote""")
    return ("KNOWN" if n is not None and n < 60000 else "FAIL"), \
        "compounds with an attesting organism contradicting the vote %s" % n, \
        "six times the phylum-level figure, since it counts any attesting organism"


@test("S14", "genuine cross-kingdom attestation by authority lineage")
def s14(c):
    n = c.q1("""WITH k AS (
        SELECT ct.comp_id, count(DISTINCT
          CASE WHEN ct.ncbi_lineage::text ILIKE '%%Metazoa%%'       THEN 'animal'
               WHEN ct.ncbi_lineage::text ILIKE '%%Fungi%%'         THEN 'fungi'
               WHEN ct.ncbi_lineage::text ILIKE '%%Bacteria%%'      THEN 'bacteria'
               WHEN ct.ncbi_lineage::text ILIKE '%%Viridiplantae%%' THEN 'plant' END) AS n
        FROM compound_taxonomy ct WHERE ct.ncbi_lineage IS NOT NULL GROUP BY 1)
        SELECT count(*) FROM k WHERE n > 1""")
    return ("KNOWN" if n is not None else "FAIL"), \
        "compounds attested across kingdoms by lineage %s" % n, \
        "independent of the vote, and the number S2 could report alongside it"


@test("S15", "secondary kingdoms are drawn from the closed vocabulary")
def s15(c):
    r = c.q("""SELECT DISTINCT unnest(secondary_kingdoms) FROM resolved_taxonomy
               WHERE secondary_kingdoms IS NOT NULL""")
    bad = [x[0] for x in r if x[0] not in KINGDOMS]
    return not bad, "secondary kingdom values outside the vocabulary %s" % (bad or "none")


@test("S16", "a compound never lists its primary kingdom as a secondary")
def s16(c):
    n = c.q1("""SELECT count(*) FROM resolved_taxonomy
                WHERE kingdom = ANY(secondary_kingdoms)""")
    return n == 0, "compounds with their primary repeated as secondary %s" % n


@test("S17", "at most three secondary kingdoms")
def s17(c):
    n = c.q1("SELECT count(*) FROM resolved_taxonomy WHERE array_length(secondary_kingdoms,1) > 3")
    return n == 0, "compounds with more than three secondaries %s" % n, \
        "four resolvable kingdoms imply a ceiling of three"


# ------------------------------------------- T. attested-lineage layer

@test("T01", "every compound with a family has a resolvable lineage row")
def t01(c):
    n = c.q1("""SELECT count(*) FROM resolved_taxonomy rt
                WHERE coalesce(rt.family,'') <> ''
                  AND NOT EXISTS (SELECT 1 FROM compound_taxonomy ct
                                  WHERE ct.comp_id = rt.comp_id AND coalesce(ct.family,'') <> '')""")
    return ("KNOWN" if n is not None and n < 1000 else "FAIL"), \
        "resolved families with no attested counterpart %s" % n


@test("T02", "the resolved family is among the attested families")
def t02(c):
    n = c.q1("""SELECT count(*) FROM resolved_taxonomy rt
                WHERE coalesce(rt.family,'') <> ''
                  AND EXISTS (SELECT 1 FROM compound_taxonomy ct
                              WHERE ct.comp_id = rt.comp_id AND coalesce(ct.family,'') <> '')
                  AND NOT EXISTS (SELECT 1 FROM compound_taxonomy ct
                                  WHERE ct.comp_id = rt.comp_id
                                    AND lower(ct.family) = lower(rt.family))""")
    return n == 0, "resolved families absent from the attested set %s" % n, \
        "the resolver must select from the evidence, not invent a rank"


@test("T03", "the resolved genus is among the attested genera")
def t03(c):
    n = c.q1("""SELECT count(*) FROM resolved_taxonomy rt
                WHERE coalesce(rt.genus,'') <> ''
                  AND EXISTS (SELECT 1 FROM compound_taxonomy ct
                              WHERE ct.comp_id = rt.comp_id AND coalesce(ct.genus,'') <> '')
                  AND NOT EXISTS (SELECT 1 FROM compound_taxonomy ct
                                  WHERE ct.comp_id = rt.comp_id
                                    AND lower(ct.genus) = lower(rt.genus))""")
    return n == 0, "resolved genera absent from the attested set %s" % n



# ---------------------------------------------- U. interface behaviour

@test("U01", "every kingdom filter value returns rows")
def u01(c):
    out = {k: c.total("type=kingdom&q=" + k) for k in KINGDOMS}
    return all(v for v in out.values()), \
        "; ".join("%s %s" % (k, v) for k, v in out.items())


@test("U02", "kingdom filter matches primary plus secondary")
def u02(c):
    api = c.total("type=kingdom&q=fungi")
    db = c.q1("""SELECT count(*) FROM resolved_taxonomy
                 WHERE lower(kingdom)='fungi' OR 'fungi' = ANY(secondary_kingdoms)""")
    return api == db, "api %s, db %s" % (api, db), \
        "the filter is attestation-wide while Table S1 reports the primary vote"


@test("U03", "the kingdom column header says primary")
def u03(c):
    st, html, _ = c.get("/search?type=kingdom&q=fungi")
    return "Kingdom (primary)" in html, "results header names the primary reading", \
        "the column shows one value where the filter matched any attestation"


@test("U04", "the kingdom filter control says any attestation")
def u04(c):
    st, html, _ = c.get("/search?type=kingdom&q=fungi")
    return "Kingdom (any attestation)" in html, "filter control names the matching rule"


@test("U05", "the compound page renders for every kingdom")
def u05(c):
    out = {}
    for k in KINGDOMS:
        cid = c.q1("SELECT comp_id FROM resolved_taxonomy WHERE lower(kingdom)=%s LIMIT 1", (k,))
        out[k] = c.get("/compound/%s" % cid)[0] if cid else 0
    return all(v == 200 for v in out.values()), "status per kingdom %s" % out


@test("U06", "the tree renders for a secondary-only kingdom query")
def u06(c):
    d, st, dt = c.js("/api/taxonomy_tree?type=kingdom&q=bacteria")
    n = (d or {}).get("total_compounds")
    return st == 200 and n, "tree returns %s compounds in %.1fs" % (n, dt)


@test("U08", "long attestation lists are capped with an expander")
def u08(c):
    cid = c.q1("""SELECT comp_id FROM compounds
                  WHERE array_length(string_to_array(source_organism,'; '),1) > 40 LIMIT 1""")
    html = c.page(cid)
    return "more, show all" in html, "%s renders the inline expander" % cid, \
        "197 organisms unpaginated pushed the page to five columns-worth of list"


@test("U09", "the expander is absent where lists are short")
def u09(c):
    cid = c.q1("""SELECT comp_id FROM compounds
                  WHERE coalesce(source_organism,'') <> ''
                    AND array_length(string_to_array(source_organism,'; '),1) <= 3 LIMIT 1""")
    html = c.page(cid)
    return "more, show all" not in html, "%s renders no expander" % cid


@test("U10", "a compound with no organism renders")
def u10(c):
    cid = c.q1("SELECT comp_id FROM compounds WHERE coalesce(source_organism,'')='' LIMIT 1")
    st, html, _ = c.get("/compound/%s" % cid)
    return st == 200 and len(html) > 2000, "%s renders, %s bytes" % (cid, len(html))


@test("U11", "licence provenance echoes the resolved tier")
def u11(c):
    cid = c.q1("""SELECT comp_id FROM per_source_license_attestation
                  GROUP BY 1 HAVING count(*) > 5 LIMIT 1""")
    tier = c.q1("SELECT license_tier FROM compounds WHERE comp_id=%s", (cid,))
    d, st, _ = c.js("/api/compound/%s/license-provenance" % cid)
    return st == 200 and tier and tier in json.dumps(d or {}), \
        "%s tier %s echoed by the endpoint" % (cid, tier)


@test("U12", "stereoisomer family members share the connectivity prefix")
def u12(c):
    cid = c.q1("""SELECT comp_id FROM compounds WHERE left(inchikey,14) IN
                  (SELECT left(inchikey,14) FROM compounds GROUP BY 1 HAVING count(*) > 2)
                  LIMIT 1""")
    d, st, _ = c.js("/api/stereoisomers/%s" % cid)
    mem = (d or {}).get("stereoisomers") or []
    pref = {(m.get("inchikey") or "")[:14] for m in mem}
    echoed = (d or {}).get("inchikey_prefix") or ""
    return st == 200 and len(mem) >= 2 and pref == {echoed}, \
        "%s members, %s distinct prefixes" % (len(mem), len(pref))


# --------------------------------------------- V. classification layer

@test("V01", "no PLACEHOLDER survives in the materialized class")
def v01(c):
    n = c.q1("SELECT count(*) FROM compounds WHERE effective_class='PLACEHOLDER'")
    return n == 0, "PLACEHOLDER rows %s" % n


@test("V02", "the three tiers partition the corpus")
def v02(c):
    r = dict(c.q("SELECT coalesce(class_source,'(null)'), count(*) FROM compounds GROUP BY 1"))
    return sum(r.values()) == CORPUS, "tier sum %s" % format(sum(r.values()), ","), \
        ", ".join("%s %s" % (k, format(v, ",")) for k, v in sorted(r.items(), key=lambda x: -x[1]))


@test("V03", "model-tier classes stay inside the curated vocabulary")
def v03(c):
    n = c.q1("""SELECT count(DISTINCT c.effective_class) FROM compounds c
                WHERE c.class_source='inferred_xgb' AND coalesce(c.effective_class,'')<>''
                  AND NOT EXISTS (SELECT 1 FROM compounds t WHERE t.class_source='curated'
                                  AND t.effective_class = c.effective_class)""")
    return n == 0, "model classes absent from the curated set %s" % n


@test("V04", "separator conventions hold per rank")
def v04(c):
    r = c.q("""SELECT
        count(*) FILTER (WHERE effective_class LIKE '%%  %%' OR effective_class LIKE '%% | %%'),
        count(*) FILTER (WHERE np_pathway LIKE '%%  %%' OR np_pathway LIKE '%% | %%'),
        count(*) FILTER (WHERE inferred_class LIKE '%%  %%' OR inferred_class LIKE '%% | %%')
        FROM compounds""")[0]
    return sum(r) == 0, "tight-separator columns with artifacts %s" % (tuple(r),), \
        "effective_superclass and effective_pathway use ' | ' by convention"


@test("V05", "ontology parents are ranked and unique")
def v05(c):
    a = c.q1("SELECT count(*) FROM npc_class_parents WHERE parent_rank IS NULL")
    b = c.q1("""SELECT count(*) FROM (SELECT class_name FROM npc_class_parents
                WHERE parent_rank=1 GROUP BY 1 HAVING count(*)>1) x""")
    return a == 0 and b == 0, "unranked %s, duplicate rank 1 %s" % (a, b)


@test("V06", "the two pathway denominators")
def v06(c):
    a = c.q1("SELECT count(*) FROM compounds WHERE coalesce(np_pathway,'')<>''")
    b = c.q1("SELECT count(*) FROM compounds WHERE coalesce(effective_pathway,'')<>''")
    return ("KNOWN" if a == 1101638 else "FAIL"), \
        "np_pathway %s, effective_pathway %s" % (format(a, ","), format(b, ",")), \
        "Figure S3 quotes the source-supplied column, guard H53 the derived one"


# ------------------------------------------------ W. licensing invariants

@test("W01", "most-restrictive-wins holds")
def w01(c):
    n = c.q1("""SELECT count(*) FROM compounds c
        JOIN (SELECT a.comp_id, max(r.tier_rank) AS worst
              FROM per_source_license_attestation a
              JOIN source_license_ref r ON lower(r.src)=lower(a.source) GROUP BY 1) s
          ON s.comp_id = c.comp_id
        WHERE c.tier_rank IS DISTINCT FROM s.worst""")
    return n == 0, "compounds not at the most restrictive attested tier %s" % n


@test("W02", "least-restrictive bound holds")
def w02(c):
    n = c.q1("""SELECT count(*) FROM compounds c
        JOIN (SELECT a.comp_id, min(r.tier_rank) AS best
              FROM per_source_license_attestation a
              JOIN source_license_ref r ON lower(r.src)=lower(a.source) GROUP BY 1) s
          ON s.comp_id = c.comp_id
        WHERE c.tier_rank_min IS DISTINCT FROM s.best""")
    return n == 0, "compounds not at the least restrictive attested tier %s" % n, \
        "the second bound must track the attestation minimum or it drifts silently"


@test("W03", "the interval never inverts")
def w03(c):
    n = c.q1("SELECT count(*) FROM compounds WHERE tier_rank_min > tier_rank")
    return n == 0, "compounds where the least exceeds the most restrictive %s" % n


@test("W04", "the interval is degenerate for single-tier compounds")
def w04(c):
    n = c.q1("""SELECT count(*) FROM compounds c
        WHERE c.tier_rank_min <> c.tier_rank
          AND (SELECT count(DISTINCT r.tier_rank) FROM per_source_license_attestation a
               JOIN source_license_ref r ON lower(r.src)=lower(a.source)
               WHERE a.comp_id = c.comp_id) < 2""")
    return n == 0, "single-tier compounds with a non-degenerate interval %s" % n


@test("W05", "the public-domain tier rests on one source")
def w05(c):
    r = c.q("""SELECT r.src, count(*) FROM per_source_license_attestation a
               JOIN source_license_ref r ON lower(r.src)=lower(a.source)
               WHERE r.tier_rank = 0 GROUP BY 1""")
    return ("KNOWN" if len(r) == 1 else "FAIL"), \
        "CC0 sources with attestations %s" % len(r), \
        "%s alone determines the public-domain tier at both bounds" % (r[0][0] if r else "-")


@test("W06", "every attestation joins the licence map")
def w06(c):
    n = c.q1("""SELECT count(*) FROM per_source_license_attestation a
                LEFT JOIN source_license_ref r ON lower(r.src)=lower(a.source)
                WHERE r.src IS NULL""")
    return n == 0, "attestations with no licence-map match %s" % n


@test("W07", "conflicting licence-map variants are unreachable")
def w07(c):
    n = c.q1("""SELECT count(*) FROM (
        SELECT r.src FROM source_license_ref r
        WHERE EXISTS (SELECT 1 FROM per_source_license_attestation a
                      WHERE lower(a.source) = lower(r.src))
        GROUP BY r.src) x""")
    tot = c.q1("SELECT count(*) FROM source_license_ref")
    return ("KNOWN" if n and n < tot else "FAIL"), \
        "%s of %s licence-map rows are reachable" % (n, tot), \
        "orphan variants at conflicting tiers are a hazard for any re-run"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="http://localhost:5000")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    con = psycopg2.connect(DSN)
    con.set_session(readonly=True, autocommit=True)
    ctx = Ctx(con, a.base)
    print("THEOBROMA adversarial suite IV  %d checks  base=%s" % (len(TESTS), a.base))
    print("-" * 74)
    t0 = time.time()
    for tag, note, fn in TESTS:
        try:
            r = fn(ctx)
            code = r[0] if isinstance(r[0], str) else ("PASS" if r[0] else "FAIL")
            emit(tag, code, r[1], list(r[2:]))
        except Exception as e:
            emit(tag, "ERROR", "%s: %s" % (type(e).__name__, str(e).strip()[:70]), [note])
            if a.verbose:
                traceback.print_exc()
            try:
                ctx.cur.close()
                ctx.cur = con.cursor()
            except Exception:
                pass
    print("-" * 74)
    print("SUMMARY  PASS=%d FAIL=%d KNOWN=%d ERROR=%d  (%.0fs)" % (
        COUNT["PASS"], COUNT["FAIL"], COUNT["KNOWN"], COUNT["ERROR"], time.time() - t0))
    con.close()
    return 1 if COUNT["FAIL"] or COUNT["ERROR"] else 0


if __name__ == "__main__":
    sys.exit(main())
