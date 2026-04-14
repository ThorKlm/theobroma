"""Similarity search engine using precomputed Morgan fingerprints."""
import numpy as np
from scipy import sparse
import os

class SimilarityEngine:
    def __init__(self, vectors_dir="data/vectors"):
        self.loaded = False
        self.vectors_dir = vectors_dir

    def load(self):
        d = self.vectors_dir
        if not os.path.exists(os.path.join(d, "morgan_fps.npz")):
            print("[SimilarityEngine] No vectors found, disabled")
            return
        print("[SimilarityEngine] Loading Morgan fingerprints...")
        self.comp_ids = np.load(os.path.join(d, "comp_ids.npy"), allow_pickle=True)
        self.morgan = sparse.load_npz(os.path.join(d, "morgan_fps.npz"))
        self.valid_idx = np.load(os.path.join(d, "valid_indices.npy"))
        self.loaded = True
        print(f"[SimilarityEngine] {len(self.valid_idx):,} valid compounds")

    def smiles_to_morgan(self, smiles):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros(2048, dtype=np.uint8)
        for bit in fp.GetOnBits():
            arr[bit] = 1
        return sparse.csr_matrix(arr.reshape(1, -1))

    def tanimoto_search(self, query_smiles, top_n=50, threshold=0.3):
        if not self.loaded:
            return []
        qfp = self.smiles_to_morgan(query_smiles)
        if qfp is None:
            return []
        valid_morgan = self.morgan[self.valid_idx]
        intersection = valid_morgan.dot(qfp.T).toarray().flatten()
        query_bits = qfp.sum()
        target_bits = np.array(valid_morgan.sum(axis=1)).flatten()
        union = query_bits + target_bits - intersection
        union[union == 0] = 1
        tanimoto = intersection / union
        top_idx = np.argsort(tanimoto)[-top_n:][::-1]
        top_idx = top_idx[tanimoto[top_idx] >= threshold] if threshold > 0 else top_idx
        if len(top_idx) == 0:
            top_idx = np.argsort(tanimoto)[-top_n:][::-1]
        results = []
        for i in top_idx:
            results.append({
                "comp_id": str(self.comp_ids[self.valid_idx[i]]),
                "tanimoto": round(float(tanimoto[i]), 4)
            })
        return results
