"""Subset feature blocks to labelled rows, resolving comp_id to feature row via valid_indices."""
import numpy as np, pandas as pd, os
from scipy import sparse
V, OUT = "data/vectors", "v135_bundle"
os.makedirs(OUT, exist_ok=True)
cid = np.load(f"{V}/comp_ids.npy", allow_pickle=True).astype(str)
vi = np.load(f"{V}/valid_indices.npy", allow_pickle=True)
mo = sparse.load_npz(f"{V}/morgan_fps.npz").tocsr()
ma = sparse.load_npz(f"{V}/maccs_fps.npz").tocsr()
emb = np.load(f"{V}/chemberta_embeddings.npy", mmap_mode="r")
n_feat = mo.shape[0]
assert ma.shape[0] == n_feat and emb.shape[0] == n_feat, (mo.shape, ma.shape, emb.shape)
if vi.dtype == bool:
    feat_cid = cid[vi]
elif len(vi) == n_feat:
    feat_cid = cid[vi.astype(np.int64)]
else:
    raise SystemExit(f"cannot map: valid_indices len {len(vi)} vs feature rows {n_feat}")
row_of = {c: i for i, c in enumerate(feat_cid)}
lab = pd.read_csv(f"{V}/labels_v135.tsv", sep="\t")
lab["feat_row"] = lab.comp_id.map(row_of)
missing = lab.feat_row.isna().sum()
print(f"labels {len(lab):,}  unmapped {missing:,}")
lab = lab.dropna(subset=["feat_row"]).copy()
idx = np.sort(lab.feat_row.astype(np.int64).unique())
remap = {int(o): i for i, o in enumerate(idx)}
lab["row_idx"] = lab.feat_row.astype(np.int64).map(remap)
lab.drop(columns=["feat_row"]).to_csv(f"{OUT}/labels_v135.tsv", sep="\t", index=False)
sparse.save_npz(f"{OUT}/morgan_fps.npz", mo[idx]); print("morgan", mo[idx].shape)
sparse.save_npz(f"{OUT}/maccs_fps.npz", ma[idx]); print("maccs", ma[idx].shape)
out = np.empty((len(idx), emb.shape[1]), dtype=np.float32)
for s in range(0, len(idx), 20000):
    out[s:s+20000] = emb[idx[s:s+20000]]
np.save(f"{OUT}/chemberta_embeddings.npy", out)
print("emb", out.shape, f"{out.nbytes/1e9:.2f} GB")
