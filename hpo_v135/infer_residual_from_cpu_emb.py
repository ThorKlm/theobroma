"""Phase-3 inference using residual-only CPU embeddings + full-corpus fingerprints.

Consumes the output of embed_residual_cpu.py (residual_embeddings.npy +
residual_comp_ids.npy, keyed by comp_id) and the full-corpus Morgan/MACCS already on
the box (data/vectors), assembles the 2,471-d vector (Morgan 2048 + MACCS 167 +
PCA-256 of ChemBERTa) per residual compound, and predicts class + confidence.

CPU or GPU for the XGBoost predict (booster device set from availability); the heavy
embedding was already done on CPU by embed_residual_cpu.py.
"""
import argparse, os, json
import numpy as np, pandas as pd
from scipy import sparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", default="out_v135b")
    ap.add_argument("--emb-dir", default="residual_vectors", help="embed_residual_cpu.py output")
    ap.add_argument("--vectors", default="data/vectors", help="full-corpus morgan/maccs + comp_ids")
    ap.add_argument("--out", default="inferred_v135b.tsv")
    ap.add_argument("--batch", type=int, default=50000)
    args = ap.parse_args()
    import xgboost as xgb

    md = args.model_dir
    booster = xgb.Booster(); booster.load_model(os.path.join(md, "xgb_v135.ubj"))
    classes = np.asarray(json.load(open(os.path.join(md, "classes.json"))))
    comp = np.load(os.path.join(md, "pca_components.npy"))   # (256,768)
    mean = np.load(os.path.join(md, "pca_mean.npy"))         # (768,)
    print(f"[model] {len(classes)} classes")

    # residual embeddings keyed by comp_id
    r_emb = np.load(os.path.join(args.emb_dir, "residual_embeddings.npy"))
    r_cid = np.load(os.path.join(args.emb_dir, "residual_comp_ids.npy"), allow_pickle=True).astype(str)
    emb_of = {c: i for i, c in enumerate(r_cid)}
    print(f"[emb] residual embeddings {r_emb.shape}")

    # full-corpus fingerprints, mapped by comp_id via fingerprint_indices
    V = args.vectors
    cid = np.load(f"{V}/comp_ids.npy", allow_pickle=True).astype(str)
    fi  = np.load(f"{V}/fingerprint_indices.npy").astype(np.int64)
    mo  = sparse.load_npz(f"{V}/morgan_fps.npz").tocsr()
    ma  = sparse.load_npz(f"{V}/maccs_fps.npz").tocsr()
    f_row = {c: i for i, c in enumerate(cid[fi])}

    # score compounds that have BOTH an embedding and a fingerprint
    ids = [c for c in r_cid if c in f_row]
    print(f"[in] residual with emb {len(r_cid):,}; also with fingerprint {len(ids):,}")

    rows = []
    for s in range(0, len(ids), args.batch):
        chunk = ids[s:s+args.batch]
        ei = np.fromiter((emb_of[c] for c in chunk), dtype=np.int64, count=len(chunk))
        fj = np.fromiter((f_row[c] for c in chunk), dtype=np.int64, count=len(chunk))
        proj = ((r_emb[ei] - mean) @ comp.T).astype(np.float32)
        X = np.hstack([mo[fj].toarray().astype(np.float32),
                       ma[fj].toarray().astype(np.float32), proj])
        p = booster.predict(xgb.DMatrix(X))
        top = p.argmax(1); conf = p.max(1)
        for c, t, cf in zip(chunk, top, conf):
            rows.append((c, str(classes[t]), float(cf)))
        print(f"  scored {min(s+args.batch,len(ids)):,}/{len(ids):,}")
    out = pd.DataFrame(rows, columns=["comp_id", "inferred_class", "inferred_confidence"])
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[out] {len(out):,} rows -> {args.out}")
    print("[dist] top predicted classes:")
    print(out.inferred_class.value_counts().head(10).to_string())
    print(f"[conf] mean {out.inferred_confidence.mean():.3f} "
          f"min {out.inferred_confidence.min():.3f} max {out.inferred_confidence.max():.3f}")

if __name__ == "__main__":
    main()
