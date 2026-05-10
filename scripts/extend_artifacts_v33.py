"""Phase B-2: extend production artifacts with the TIPdb Phase B bundle.

Inputs (on L3S, after Phase B bundle has landed):
    data/vectors/{chemberta_embeddings,comp_ids,valid_indices,morgan_fps,
                  maccs_fps,faiss_hnsw}.* (current production baseline)
    data/tipdb_phase_b/{chemberta_embeddings,comp_ids,morgan_fps,maccs_fps,
                        valid_mask}.* (Phase B bundle from vast.ai)
    PostgreSQL compounds table (for fetching the production-side THEO_ ids
        of the new TIPdb rows; the bundle's comp_ids carry locally-assigned
        ids that do not match what the production-merge SQL inserted).

Outputs (with .post-cleanup backups for the modified files):
    data/vectors/comp_ids.npy           (extended by 472 new entries)
    data/vectors/valid_indices.npy      (extended by 472 new entries)
    data/vectors/chemberta_embeddings.npy  (extended by 472 rows)
    data/vectors/morgan_fps.npz         (extended by 472 rows)
    data/vectors/maccs_fps.npz          (extended by 472 rows)
    data/vectors/faiss_hnsw.index       (rebuilt from augmented embeddings)

The FAISS rebuild on 1.07M+ vectors is the slowest step and may push the
L3S box close to its memory ceiling. To reduce the chance of OOM (which
killed an earlier dry-run when production gunicorn was holding 9 GB) the
script can be invoked with --skip-faiss to defer the FAISS rebuild to a
quieter window. In that case, similarity search via ChemBERTa returns
results from the pre-extension index until FAISS is rebuilt.

Run from the repo root with the production venv:
    cd ~/theobroma
    venv/bin/python scripts/extend_artifacts_v33.py --dry-run
    venv/bin/python scripts/extend_artifacts_v33.py
    venv/bin/python scripts/extend_artifacts_v33.py --skip-faiss
"""
import argparse, os, shutil, sys, time
import numpy as np
import psycopg2
import psycopg2.extras
from scipy import sparse

REPO = os.path.expanduser("~/theobroma")
VEC = os.path.join(REPO, "data", "vectors")
BUNDLE = os.path.join(REPO, "data", "tipdb_phase_b")
DB_URI = os.environ.get("DATABASE_URL",
    "postgresql://theobroma:theobroma@localhost:5432/theobroma")

EMB = os.path.join(VEC, "chemberta_embeddings.npy")
COMP_IDS = os.path.join(VEC, "comp_ids.npy")
VALID = os.path.join(VEC, "valid_indices.npy")
MORGAN = os.path.join(VEC, "morgan_fps.npz")
MACCS = os.path.join(VEC, "maccs_fps.npz")
FAISS_IDX = os.path.join(VEC, "faiss_hnsw.index")

BUNDLE_EMB = os.path.join(BUNDLE, "chemberta_embeddings.npy")
BUNDLE_CID = os.path.join(BUNDLE, "comp_ids.npy")
BUNDLE_MORGAN = os.path.join(BUNDLE, "morgan_fps.npz")
BUNDLE_MACCS = os.path.join(BUNDLE, "maccs_fps.npz")
BUNDLE_VALID = os.path.join(BUNDLE, "valid_mask.npy")

MODIFIED = [EMB, COMP_IDS, VALID, MORGAN, MACCS, FAISS_IDX]
BACKUP_SUFFIX = ".post-cleanup"

def fetch_production_ids():
    """Get the 472 newly-inserted TIPdb rows and their production-side
    THEO_ ids and inchikeys. Order matches the bundle's CSV row order
    via the staging-table sort by inchikey that the production-merge
    SQL applied."""
    print("[db] fetching production THEO_ ids for TIPdb-novel rows")
    conn = psycopg2.connect(DB_URI, cursor_factory=psycopg2.extras.RealDictCursor)
    with conn.cursor() as cur:
        cur.execute("""SELECT comp_id, inchikey
                       FROM compounds
                       WHERE source_db = 'TIPdb'
                       ORDER BY inchikey""")
        rows = cur.fetchall()
    conn.close()
    if len(rows) != 472:
        sys.exit(f"expected 472 TIPdb rows, got {len(rows)}")
    print(f"  {len(rows)} rows; first comp_id: {rows[0]['comp_id']}, "
          f"last: {rows[-1]['comp_id']}")
    ik_to_compid = {r["inchikey"]: r["comp_id"] for r in rows}
    return ik_to_compid

def backup_originals(dry_run):
    for src in MODIFIED:
        if not os.path.exists(src):
            print(f"  WARNING: {src} not found")
            continue
        dst = src + BACKUP_SUFFIX
        if os.path.exists(dst):
            print(f"  backup exists: {os.path.basename(dst)}")
            continue
        size_gb = os.path.getsize(src) / 1024**3
        print(f"  backup: {os.path.basename(src)} ({size_gb:.2f} GB) -> {BACKUP_SUFFIX}")
        if not dry_run:
            shutil.copy2(src, dst)

def load_bundle_with_compids(ik_to_compid):
    """Load the bundle and remap its comp_ids to production-side ids by
    joining on inchikey. The bundle's comp_ids carry locally-assigned ids
    (THEO_1078682..) that will not match the database (THEO_1080994..)."""
    print(f"[bundle] loading from {BUNDLE}")
    bundle_emb = np.load(BUNDLE_EMB)
    bundle_cid = np.load(BUNDLE_CID, allow_pickle=True).astype(str)
    bundle_morgan = sparse.load_npz(BUNDLE_MORGAN)
    bundle_maccs = sparse.load_npz(BUNDLE_MACCS)
    bundle_valid_mask = np.load(BUNDLE_VALID)
    print(f"  embeddings: {bundle_emb.shape}")
    print(f"  comp_ids: {len(bundle_cid)}")
    print(f"  morgan: {bundle_morgan.shape}")
    print(f"  maccs: {bundle_maccs.shape}")
    print(f"  valid_mask: {bundle_valid_mask.shape}, "
          f"true count: {int(bundle_valid_mask.sum())}")

    # Read the input CSV that produced the bundle; the row order there
    # is the same as the bundle arrays. From it we get inchikey per row.
    import pandas as pd
    smiles_csv = os.path.join(REPO, "data", "tipdb_novel_smiles.csv")
    df = pd.read_csv(smiles_csv)
    print(f"  input CSV: {len(df)} rows")
    if len(df) != len(bundle_cid):
        sys.exit(f"row-count mismatch: CSV has {len(df)}, bundle has "
                 f"{len(bundle_cid)}")

    # Map each bundle row to its production-side comp_id via inchikey.
    new_cids = []
    for ik in df["inchikey"]:
        if ik not in ik_to_compid:
            sys.exit(f"inchikey {ik} from bundle not found in production "
                     f"compounds table")
        new_cids.append(ik_to_compid[ik])
    new_cids = np.array(new_cids, dtype=object)
    print(f"  remapped comp_ids: first={new_cids[0]}, last={new_cids[-1]}")

    return (bundle_emb, new_cids, bundle_morgan, bundle_maccs,
            bundle_valid_mask)

def extend_artifacts(bundle_emb, new_cids, bundle_morgan, bundle_maccs,
                     bundle_valid_mask, dry_run):
    """Append bundle arrays to baseline arrays. valid_indices grows by
    the number of bundle entries with valid_mask=True; values are the
    new positions in the augmented comp_ids array (where the new
    compounds will sit)."""
    print("[load] baseline arrays")
    base_cid = np.load(COMP_IDS, allow_pickle=True).astype(str)
    base_valid = np.load(VALID)
    print(f"  base comp_ids: {len(base_cid)}")
    print(f"  base valid_indices: {len(base_valid)}")

    print("[extend] comp_ids")
    new_full_cid = np.concatenate([base_cid, new_cids])
    print(f"  new length: {len(new_full_cid)} ({len(base_cid)} + {len(new_cids)})")
    if not dry_run:
        np.save(COMP_IDS, new_full_cid)

    print("[extend] valid_indices")
    # New embedding entries sit at positions len(base_cid) onwards; only
    # those with valid_mask=True actually got embeddings (all 472 in our
    # case, but the script is robust to partial validity).
    n_base_cid = len(base_cid)
    new_valid_entries = np.where(bundle_valid_mask)[0] + n_base_cid
    new_full_valid = np.concatenate([base_valid, new_valid_entries.astype(base_valid.dtype)])
    print(f"  new length: {len(new_full_valid)} "
          f"({len(base_valid)} + {len(new_valid_entries)})")
    if not dry_run:
        np.save(VALID, new_full_valid)

    print("[extend] embeddings")
    base_emb = np.load(EMB, mmap_mode="r")
    valid_bundle_emb = bundle_emb[bundle_valid_mask]
    print(f"  base: {base_emb.shape}, bundle valid: {valid_bundle_emb.shape}")
    new_emb = np.concatenate([np.array(base_emb), valid_bundle_emb], axis=0)
    print(f"  new: {new_emb.shape}")
    if not dry_run:
        np.save(EMB, new_emb)

    print("[extend] morgan_fps")
    base_morgan = sparse.load_npz(MORGAN)
    print(f"  base: {base_morgan.shape}, bundle: {bundle_morgan.shape}")
    new_morgan = sparse.vstack([base_morgan, bundle_morgan]).tocsr()
    print(f"  new: {new_morgan.shape}")
    if not dry_run:
        sparse.save_npz(MORGAN, new_morgan)

    print("[extend] maccs_fps")
    base_maccs = sparse.load_npz(MACCS)
    print(f"  base: {base_maccs.shape}, bundle: {bundle_maccs.shape}")
    new_maccs = sparse.vstack([base_maccs, bundle_maccs]).tocsr()
    print(f"  new: {new_maccs.shape}")
    if not dry_run:
        sparse.save_npz(MACCS, new_maccs)

    return new_emb

def rebuild_faiss(emb, dry_run):
    import faiss
    print("[faiss] L2-normalize embeddings")
    valid_emb = emb.astype(np.float32, copy=False)
    norms = np.linalg.norm(valid_emb, axis=1, keepdims=True)
    norms[norms == 0] = 1
    valid_emb = valid_emb / norms
    d = valid_emb.shape[1]
    print(f"[faiss] building HNSW: shape={valid_emb.shape}, M=32, efC=200")
    print("        warning: holds ~5-10 GB RAM during construction; if "
          "production gunicorn is using significant RAM concurrently this "
          "may OOM-kill. Consider stopping gunicorn first.")
    t0 = time.time()
    index = faiss.IndexHNSWFlat(d, 32)
    index.hnsw.efConstruction = 200
    index.hnsw.efSearch = 128
    index.add(valid_emb)
    print(f"[faiss] built {index.ntotal:,} vectors in {time.time()-t0:.1f}s")
    if not dry_run:
        faiss.write_index(index, FAISS_IDX)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-faiss", action="store_true",
                    help="skip the FAISS rebuild (run later in a quieter window)")
    args = ap.parse_args()
    print(f"=== Phase B-2 ({'DRY RUN' if args.dry_run else 'WRITE'}) ===\n")

    if not os.path.exists(BUNDLE):
        sys.exit(f"bundle not found: {BUNDLE}")
    for path in [BUNDLE_EMB, BUNDLE_CID, BUNDLE_MORGAN, BUNDLE_MACCS,
                 BUNDLE_VALID]:
        if not os.path.exists(path):
            sys.exit(f"bundle missing file: {path}")

    ik_to_compid = fetch_production_ids()
    print()
    backup_originals(args.dry_run)
    print()
    bundle_data = load_bundle_with_compids(ik_to_compid)
    print()
    new_emb = extend_artifacts(*bundle_data, dry_run=args.dry_run)
    print()
    if args.skip_faiss:
        print("[faiss] skipped (--skip-faiss); rebuild later via separate run")
    else:
        rebuild_faiss(new_emb, args.dry_run)
    print(f"\n=== Phase B-2 complete ===")

if __name__ == "__main__":
    main()
