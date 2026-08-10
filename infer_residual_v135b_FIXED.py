"""Corpus-wide inference for the INFERRED tier (Tier 3), matched to the ACTUAL HPO
artifacts: xgb_v135.ubj (booster), pca_components.npy + pca_mean.npy (PCA), classes.json.

Scores only the residual (comp_ids with neither a curated nor an NPClassifier class),
assembling the same 2,471-d vector (Morgan 2048 + MACCS 167 + PCA-256 of ChemBERTa),
and writes comp_id, inferred_class, inferred_confidence.
"""
import argparse, os, json
import numpy as np, pandas as pd
from scipy import sparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", default="out_v135b")
    ap.add_argument("--vectors", default="v135b_bundle_fullcorpus")  # see note below
    ap.add_argument("--residual", default="residual_comp_ids.tsv")
    ap.add_argument("--out", default="inferred_v135b.tsv")
    ap.add_argument("--batch", type=int, default=50000)
    args = ap.parse_args()
    import xgboost as xgb

    md = args.model_dir
    booster = xgb.Booster(); booster.load_model(os.path.join(md, "xgb_v135.ubj"))
    classes = np.asarray(json.load(open(os.path.join(md, "classes.json"))))
    comp = np.asarray(pca_components := np.load(os.path.join(md, "pca_components.npy")))  # (256,768)
    pca_mean = np.load(os.path.join(md, "pca_mean.npy"))                                  # (768,)
    def project(emb):  # replicate sklearn PCA.transform: (emb-mean) @ components_.T
        return (emb - pca_mean) @ comp.T

    V = args.vectors
    cid = np.load(f"{V}/comp_ids.npy", allow_pickle=True).astype(str)
    vi  = np.load(f"{V}/valid_indices.npy").astype(np.int64)
    fi  = np.load(f"{V}/fingerprint_indices.npy").astype(np.int64)
    emb = np.load(f"{V}/chemberta_embeddings.npy", mmap_mode="r")
    mo  = sparse.load_npz(f"{V}/morgan_fps.npz").tocsr()
    ma  = sparse.load_npz(f"{V}/maccs_fps.npz").tocsr()
    e_row = {c: i for i, c in enumerate(cid[vi])}
    f_row = {c: i for i, c in enumerate(cid[fi])}

    resid = pd.read_csv(args.residual, sep="\t").comp_id.astype(str).tolist()
    resid = [c for c in resid if c in e_row and c in f_row]
    print(f"[in] residual scored (with features): {len(resid):,}")
    rows = []
    for s in range(0, len(resid), args.batch):
        chunk = resid[s:s+args.batch]
        ei = np.array([e_row[c] for c in chunk]); fj = np.array([f_row[c] for c in chunk])
        proj = project(np.asarray(emb[ei], dtype=np.float32)).astype(np.float32)
        X = np.hstack([mo[fj].toarray().astype(np.float32),
                       ma[fj].toarray().astype(np.float32), proj])
        p = booster.predict(xgb.DMatrix(X))
        top = p.argmax(1); conf = p.max(1)
        for c, t, cf in zip(chunk, top, conf):
            rows.append((c, classes[t], float(cf)))
        print(f"  {min(s+args.batch,len(resid)):,}/{len(resid):,}")
    out = pd.DataFrame(rows, columns=["comp_id", "inferred_class", "inferred_confidence"])
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[out] {len(out):,} inferred -> {args.out}")

if __name__ == "__main__":
    main()
