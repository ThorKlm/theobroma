#!/usr/bin/env python3
"""THEOBROMA v1.35 regression guards.

Standalone so the main suite is untouched this close to submission. Merge these
into theobroma_test_suite.py afterwards.

Covers the failure modes found on 2026-08-09 that the existing suite could not
catch, plus one existing test that passes for the wrong reason:

  H52  effective_class never holds PLACEHOLDER. The original gate checked
       np_class while the corruption lived in effective_class, which is why
       300,162 rows reached production.
  H53  class and pathway coverage above 99.9 percent.
  H54  effective_superclass and effective_pathway agree with the primary-parent
       rule. Generalises H52: catches any divergence, not one token.
  H55  npc_class_parents.parent_rank populated and rank 1 unique per class.
       The precedence convention must not depend on physical row order.
  H56  effective_pathway tokens are the seven canonical NPClassifier pathways.
  H57  no double-space separator artifacts in effective_class.
  H58  tier 3 classes are a subset of tier 1, since the model can only predict
       classes present in its training labels.
  G45b license resolver, with the case-insensitive join. The suite's G45 joins
       source_license_ref.src (lowercase) to per_source_license_attestation
       .source (proper case), matches zero rows, and passes vacuously. The
       licensing rule is the paper's headline contribution and has therefore
       never actually been tested.

Usage:
  cd ~/theobroma && source venv/bin/activate && python3 regression_guards.py
Exit code 1 if any check fails.
"""
import os, sys
import psycopg2

DSN = os.environ.get("THEO_DB_DSN", "postgresql://theobroma:theobroma@localhost:5432/theobroma")
CANONICAL = ("Alkaloids", "Terpenoids", "Shikimates and Phenylpropanoids", "Polyketides",
             "Fatty acids", "Amino acids and Peptides", "Carbohydrates")
results = []

conn = psycopg2.connect(DSN)
conn.set_session(readonly=True, autocommit=True)

def q1(sql, args=None):
    cur = conn.cursor()
    cur.execute(sql, args) if args else cur.execute(sql)
    v = cur.fetchone(); cur.close()
    return v[0] if v else None

def check(cid, ok, msg):
    results.append((cid, "PASS" if ok else "FAIL", msg))
    print(f"[{'PASS' if ok else 'FAIL'}] {cid:6s} {msg}", flush=True)

# --- H52 placeholder ---------------------------------------------------------
n = q1("SELECT count(*) FROM compounds WHERE effective_class = 'PLACEHOLDER'")
check("H52", n == 0, f"effective_class PLACEHOLDER rows: {n:,}")

# --- H53 coverage ------------------------------------------------------------
tot  = q1("SELECT count(*) FROM compounds")
cls  = q1("SELECT count(*) FROM compounds WHERE effective_class IS NOT NULL AND effective_class <> ''")
path = q1("SELECT count(*) FROM compounds WHERE effective_pathway IS NOT NULL AND effective_pathway <> ''")
check("H53", cls / tot > 0.999 and path / tot > 0.999,
      f"class {cls:,} ({100*cls/tot:.3f}%) pathway {path:,} ({100*path/tot:.3f}%) of {tot:,}")

# --- H54 ontology agreement --------------------------------------------------
n = q1("""SELECT count(*) FROM compounds c
          JOIN (SELECT lower(trim(class_name)) k, superclass_name, pathway_name
                FROM npc_class_parents WHERE parent_rank = 1) o
            ON o.k = lower(trim(c.effective_class))
          WHERE c.effective_class NOT LIKE '%|%'
            AND c.effective_superclass IS NOT NULL AND c.effective_pathway IS NOT NULL
            AND (c.effective_superclass <> o.superclass_name
                 OR c.effective_pathway <> o.pathway_name)""")
agree = q1("""SELECT count(*) FROM compounds c
              JOIN (SELECT lower(trim(class_name)) k, superclass_name, pathway_name
                    FROM npc_class_parents WHERE parent_rank = 1) o
                ON o.k = lower(trim(c.effective_class))
              WHERE c.effective_class NOT LIKE '%|%'
                AND c.effective_superclass = o.superclass_name
                AND c.effective_pathway = o.pathway_name""")
check("H54", n == 0, f"rows diverging from primary-parent rule: {n:,} (agreeing: {agree:,})")

# --- H55 precedence column ---------------------------------------------------
null_rank = q1("SELECT count(*) FROM npc_class_parents WHERE parent_rank IS NULL")
dup_first = q1("""SELECT count(*) FROM (SELECT lower(trim(class_name)) k FROM npc_class_parents
                  WHERE parent_rank = 1 GROUP BY 1 HAVING count(*) > 1) t""")
classes   = q1("SELECT count(DISTINCT lower(trim(class_name))) FROM npc_class_parents")
rank1     = q1("SELECT count(*) FROM npc_class_parents WHERE parent_rank = 1")
check("H55", null_rank == 0 and dup_first == 0 and rank1 == classes,
      f"parent_rank: {null_rank} unranked, {dup_first} classes with duplicate rank 1, "
      f"{rank1} rank-1 rows for {classes} classes")

# --- H56 canonical pathway vocabulary ----------------------------------------
cur = conn.cursor()
cur.execute("""SELECT DISTINCT trim(t) FROM compounds, unnest(string_to_array(effective_pathway,'|')) t
               WHERE effective_pathway IS NOT NULL AND effective_pathway <> ''""")
toks = {r[0] for r in cur.fetchall()}; cur.close()
bad = toks - set(CANONICAL)
check("H56", not bad, f"{len(toks)} distinct pathway tokens, non-canonical: {sorted(bad) or 'none'}")

# --- H57 separator hygiene ---------------------------------------------------
n = q1("SELECT count(*) FROM compounds WHERE effective_class LIKE '%  %' "
       "OR effective_pathway LIKE '%  %' OR effective_superclass LIKE '%  %'")
check("H57", n == 0, f"double-space separator artifacts: {n:,}")

# --- H58 tier 3 vocabulary is a subset of tier 1 -----------------------------
n = q1("""WITH t1 AS (SELECT DISTINCT trim(t) c FROM compounds,
                      unnest(string_to_array(effective_class,'|')) t WHERE classification_tier = 1),
               t3 AS (SELECT DISTINCT trim(t) c FROM compounds,
                      unnest(string_to_array(effective_class,'|')) t WHERE classification_tier = 3)
          SELECT count(*) FROM t3 WHERE c NOT IN (SELECT c FROM t1)""")
check("H58", n == 0, f"tier-3 classes absent from tier-1 training vocabulary: {n}")

# --- G45b license resolver, case-insensitive join ----------------------------
# The stored per-compound tier must equal the most restrictive attestation.
joined = q1("""SELECT count(*) FROM per_source_license_attestation a
               JOIN source_license_ref r ON lower(trim(r.src)) = lower(trim(a.source))""")
total  = q1("SELECT count(*) FROM per_source_license_attestation")
viol   = q1("""WITH best AS (
                 SELECT a.comp_id, max(r.tier_rank) b
                 FROM per_source_license_attestation a
                 JOIN source_license_ref r ON lower(trim(r.src)) = lower(trim(a.source))
                 GROUP BY 1)
               SELECT count(*) FROM compounds c JOIN best ON best.comp_id = c.comp_id
               WHERE c.tier_rank <> best.b""")
check("G45b", joined > 0 and viol == 0,
      f"attestations joined {joined:,}/{total:,}, resolver violations {viol:,}")

# --- summary -----------------------------------------------------------------
fails = [r for r in results if r[1] == "FAIL"]
print("-" * 72)
print(f"SUMMARY  PASS={len(results)-len(fails)}  FAIL={len(fails)}")
sys.exit(1 if fails else 0)
