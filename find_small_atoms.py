#!/usr/bin/env python3
"""Identify compounds with <=3 TOTAL atoms (heavy + implicit H) using RDKit,
which correctly infers implicit hydrogens (SQL cannot). Read-only: writes a
CSV of the comp_ids to remove/flag; does not modify the DB.

Run on the server (venv with rdkit + psycopg2):
    python3 find_small_atoms.py
Outputs: small_atoms.csv  (comp_id,name,mw,smiles,total_atoms)
"""
import csv, psycopg2
from rdkit import Chem

MAXTOTAL = 3  # remove compounds with total atom count <= this

conn = psycopg2.connect(host="localhost", dbname="theobroma", user="theobroma", password="theobroma")
cur = conn.cursor(name="stream")  # server-side cursor to stream 1.1M rows
cur.itersize = 20000
cur.execute("SELECT comp_id, name, mw, smiles FROM compounds WHERE smiles IS NOT NULL AND smiles<>''")

hits = []
checked = 0
for comp_id, name, mw, smiles in cur:
    checked += 1
    # quick prefilter: only bother parsing very short SMILES (long ones can't be <=3 atoms)
    if len(smiles) > 12:
        continue
    try:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            continue
        m = Chem.AddHs(m)          # make implicit H explicit
        total = m.GetNumAtoms()    # heavy + H
        if total <= MAXTOTAL:
            hits.append((comp_id, name or "", round(mw,2) if mw else "", smiles, total))
    except Exception:
        continue

cur.close(); conn.close()

with open("small_atoms.csv","w",newline="") as f:
    w = csv.writer(f); w.writerow(["comp_id","name","mw","smiles","total_atoms"])
    w.writerows(hits)

# summary
from collections import Counter
byn = Counter(h[4] for h in hits)
print(f"scanned {checked} compounds; matched {len(hits)} with <= {MAXTOTAL} total atoms")
for n in sorted(byn):
    print(f"  {n} atom(s): {byn[n]}")
print("wrote small_atoms.csv")
# show a sample
for h in sorted(hits, key=lambda x: (x[4], x[2] if isinstance(x[2],float) else 0))[:30]:
    print(h)
