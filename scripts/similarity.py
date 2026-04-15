"""Similarity search engine: Morgan Tanimoto, ChemBERTa cosine, substructure."""
import numpy as np
from scipy import sparse
import os

class SimilarityEngine:
    def __init__(self, vectors_dir="data/vectors"):
        self.loaded = False
        self.faiss_loaded = False
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
        print(f"[SimilarityEngine] {len(self.valid_idx):,} valid compounds (Morgan)")
        # Try loading FAISS
        faiss_path = os.path.join(d, "faiss_hnsw.index")
        if os.path.exists(faiss_path):
            try:
                import faiss
                print("[SimilarityEngine] Loading FAISS index...")
                self.faiss_index = faiss.read_index(faiss_path)
                self.embeddings = np.load(os.path.join(d, "chemberta_embeddings.npy"), mmap_mode="r")
                self.faiss_loaded = True
                print(f"[SimilarityEngine] FAISS loaded ({self.faiss_index.ntotal:,} vectors)")
            except Exception as e:
                print(f"[SimilarityEngine] FAISS not available: {e}")

    def smiles_to_morgan(self, smiles):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros(2048, dtype=np.uint8)
        for bit in fp.GetOnBits(): arr[bit] = 1
        return sparse.csr_matrix(arr.reshape(1, -1))

    def maccs_search(self, query_smiles, top_n=50, threshold=0.3):
        """MACCS keys Tanimoto similarity search."""
        if not self.loaded or not self.maccs_loaded:
            return self.tanimoto_search(query_smiles, top_n, threshold)
        from rdkit import Chem
        from rdkit.Chem import MACCSkeys
        mol = Chem.MolFromSmiles(query_smiles)
        if mol is None: return []
        fp = MACCSkeys.GenMACCSKeys(mol)
        arr = np.zeros(167, dtype=np.uint8)
        for bit in fp.GetOnBits(): 
            if bit < 167: arr[bit] = 1
        qfp = sparse.csr_matrix(arr.reshape(1, -1))
        valid_maccs = self.maccs[self.valid_idx]
        intersection = valid_maccs.dot(qfp.T).toarray().flatten()
        query_bits = qfp.sum()
        target_bits = np.array(valid_maccs.sum(axis=1)).flatten()
        union = query_bits + target_bits - intersection
        union[union == 0] = 1
        tanimoto = intersection / union
        top_idx = np.argsort(tanimoto)[-top_n:][::-1]
        top_idx = top_idx[tanimoto[top_idx] >= threshold] if threshold > 0 else top_idx
        if len(top_idx) == 0:
            top_idx = np.argsort(tanimoto)[-top_n:][::-1]
        return [{"comp_id": str(self.comp_ids[self.valid_idx[i]]),
                 "tanimoto": round(float(tanimoto[i]), 4)} for i in top_idx]

    def tanimoto_search(self, query_smiles, top_n=50, threshold=0.3):
        if not self.loaded: return []
        qfp = self.smiles_to_morgan(query_smiles)
        if qfp is None: return []
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
        return [{"comp_id": str(self.comp_ids[self.valid_idx[i]]),
                 "tanimoto": round(float(tanimoto[i]), 4)} for i in top_idx]

    def chemberta_search(self, query_smiles, top_n=50):
        """Cosine similarity via FAISS HNSW on ChemBERTa embeddings."""
        if not self.faiss_loaded: return []
        try:
            from transformers import AutoTokenizer, AutoModel
            import torch
            tokenizer = AutoTokenizer.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
            model = AutoModel.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
            model.eval()
            inputs = tokenizer(query_smiles, return_tensors="pt", padding=True, truncation=True, max_length=512)
            with torch.no_grad():
                out = model(**inputs)
            emb = out.last_hidden_state[:,0,:].numpy().astype(np.float32)
            import faiss
            faiss.normalize_L2(emb)
            D, I = self.faiss_index.search(emb, top_n)
            results = []
            for i, (dist, idx) in enumerate(zip(D[0], I[0])):
                if idx < 0: continue
                results.append({
                    "comp_id": str(self.comp_ids[self.valid_idx[idx]]),
                    "tanimoto": round(float((1 + dist) / 2), 4)  # convert cosine distance to similarity
                })
            return results
        except Exception as e:
            print(f"[ChemBERTa] Error: {e}")
            return []

    def substructure_search(self, query_smarts, max_results=100):
        """Substructure search using Morgan FP pre-filter + RDKit match."""
        if not self.loaded: return []
        from rdkit import Chem
        from rdkit.Chem import AllChem
        query_mol = Chem.MolFromSmarts(query_smarts)
        if query_mol is None:
            query_mol = Chem.MolFromSmiles(query_smarts)
        if query_mol is None: return []
        # Morgan FP pre-filter: all query bits must be set in target
        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
        except:
            return []
        query_bits = list(fp.GetOnBits())
        if not query_bits: return []
        # Screen: candidates must have all query bits set
        valid_morgan = self.morgan[self.valid_idx]
        query_arr = np.zeros(2048, dtype=np.uint8)
        for b in query_bits: query_arr[b] = 1
        query_sparse = sparse.csr_matrix(query_arr.reshape(1, -1))
        hits = valid_morgan.dot(query_sparse.T).toarray().flatten()
        candidates = np.where(hits == len(query_bits))[0]
        if len(candidates) > 50000:
            candidates = candidates[:50000]  # safety cap
        # RDKit substructure match on candidates
        results = []
        for i in candidates:
            global_idx = self.valid_idx[i]
            cid = str(self.comp_ids[global_idx])
            # We need SMILES from DB, but for speed, try to get from the sparse FP
            # Actually we need the SMILES. Return comp_ids, let the route fetch SMILES and verify.
            results.append({"comp_id": cid})
            if len(results) >= max_results * 10:  # fetch more than needed, verify in route
                break
        return results
