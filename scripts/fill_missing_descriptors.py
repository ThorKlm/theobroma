import psycopg2
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
RDLogger.DisableLog("rdApp.*")

con = psycopg2.connect("host=localhost dbname=theobroma user=theobroma")
cur = con.cursor()
cur.execute("""SELECT comp_id, smiles FROM compounds
               WHERE logp IS NULL AND coalesce(smiles,'') <> ''""")
rows = cur.fetchall()
print("candidates:", len(rows))

ok, bad, upd = 0, 0, []
for cid, smi in rows:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        bad += 1
        continue
    upd.append((
        Descriptors.ExactMolWt(m), Crippen.MolLogP(m), rdMolDescriptors.CalcTPSA(m),
        float(rdMolDescriptors.CalcNumHBA(m)), float(rdMolDescriptors.CalcNumHBD(m)),
        float(rdMolDescriptors.CalcNumRings(m)),
        float(rdMolDescriptors.CalcNumRotatableBonds(m)), cid))
    ok += 1
print("parsed:", ok, "failed:", bad)

cur.execute("CREATE TABLE IF NOT EXISTS compounds_desc_bak2_20260812 AS "
            "SELECT comp_id, mw, logp, tpsa, hba, hbd, n_rings, rotatable_bonds "
            "FROM compounds WHERE logp IS NULL")
cur.executemany("""UPDATE compounds SET mw=%s, logp=%s, tpsa=%s, hba=%s, hbd=%s,
                   n_rings=%s, rotatable_bonds=%s WHERE comp_id=%s""", upd)
con.commit()
cur.execute("SELECT count(*) FROM compounds WHERE logp IS NULL")
print("remaining NULL mw:", cur.fetchone()[0])
con.close()
