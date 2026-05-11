"""Regenerate Morgan and MACCS fingerprints from a SMILES column.

Reads compound SMILES from the compounds table, generates 2048-bit Morgan
(radius 2) and 167-bit MACCS keys via RDKit, and writes them to two
scipy.sparse CSR matrices on disk. The output indexing matches Space F
(see docs/vectors.md): only compounds whose SMILES parse with RDKit
contribute rows, indexed in compound-table order.

Outputs:
    data/vectors/morgan_fps.npz   (CSR, n_f x 2048)
    data/vectors/maccs_fps.npz    (CSR, n_f x 167)
    data/vectors/fingerprint_indices.npy   (int64, length n_f; positions
        in compounds.comp_ids that produced valid fingerprints)
"""
import os, sys
import numpy as np
import psycopg2
from scipy import sparse
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys

VEC = os.path.expanduser("~/theobroma/data/vectors")
DSN = "postgresql://theobroma:theobroma@localhost:5432/theobroma"

def main():
    morgan_rows, maccs_rows, valid_idx = [], [], []
    with psycopg2.connect(DSN) as conn, conn.cursor() as cur:
        cur.execute("SELECT comp_id, smiles FROM compounds ORDER BY comp_id")
        for i, (comp_id, smi) in enumerate(cur):
            mol = Chem.MolFromSmiles(smi) if smi else None
            if mol is None:
                continue
            mfp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            kfp = MACCSkeys.GenMACCSKeys(mol)
            morgan_rows.append(np.array(mfp, dtype=np.uint8))
            maccs_rows.append(np.array(kfp, dtype=np.uint8))
            valid_idx.append(i)
            if i % 50000 == 0:
                print(f"  {i:,} rows scanned, {len(valid_idx):,} valid")
    morgan = sparse.csr_matrix(np.vstack(morgan_rows))
    maccs = sparse.csr_matrix(np.vstack(maccs_rows))
    sparse.save_npz(os.path.join(VEC, "morgan_fps.npz"), morgan)
    sparse.save_npz(os.path.join(VEC, "maccs_fps.npz"), maccs)
    np.save(os.path.join(VEC, "fingerprint_indices.npy"),
            np.array(valid_idx, dtype=np.int64))
    print(f"saved {morgan.shape[0]:,} morgan + maccs rows, {len(valid_idx):,} valid_idx")

if __name__ == "__main__":
    main()
