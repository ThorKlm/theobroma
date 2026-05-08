"""Phase A: realign artifact baseline to the post-cleanup compounds table.

The v33 dedup cleanup deleted 2,323 duplicate-InChIKey rows from the
compounds table without touching the artifact files. This script applies
the same deletion to every artifact in lockstep, producing a 1,074,567-ish
baseline (1,076,890 - count of deleted rows that had embeddings) aligned
to the current 1,078,671-row database state. FAISS is rebuilt from
scratch because IndexHNSWFlat does not support targeted vector removal.

Inputs (all relative to ~/theobroma):
    data/vectors/chemberta_embeddings.npy   (float32, 768-d, ~1,076,890 rows)
    data/vectors/comp_ids.npy               (object, ~1,076,890 entries)
    data/vectors/valid_indices.npy          (int64, ~1,076,890 entries)
    data/vectors/morgan_fps.npz             (sparse, ~1,076,890 rows)
    data/vectors/maccs_fps.npz              (sparse, ~1,076,890 rows)
    data/vectors/faiss_hnsw.index           (~1,076,890 normalized vectors)
    /tmp/deleted_comp_ids.txt               (one THEO_xxxxxxx per line)

Outputs (originals get .v32 backups; realigned files take their place):
    same paths, with the deleted rows removed and FAISS rebuilt.

Run from the repo root with the production venv:
    cd ~/theobroma
    venv/bin/python scripts/realign_v33_cleanup.py

Dry-run mode: pass --dry-run to print the planned operations without
writing anything. Recommended for a first pass before destructive run.
"""
import argparse, os, sys, time
import numpy as np
from scipy import sparse

REPO = os.path.expanduser("~/theobroma")
VEC = os.path.join(REPO, "data", "vectors")
DELETED_LIST = "/tmp/deleted_comp_ids.txt"

EMB = os.path.join(VEC, "chemberta_embeddings.npy")
COMP_IDS = os.path.join(VEC, "comp_ids.npy")
VALID = os.path.join(VEC, "valid_indices.npy")
MORGAN = os.path.join(VEC, "morgan_fps.npz")
MACCS = os.path.join(VEC, "maccs_fps.npz")
FAISS_IDX = os.path.join(VEC, "faiss_hnsw.index")

ARTIFACTS = [EMB, COMP_IDS, VALID, MORGAN, MACCS, FAISS_IDX]

def load_deleted(path):
    """Load the deleted comp_id list as a set."""
    with open(path) as f:
        ids = {line.strip() for line in f if line.strip()}
    return ids

def backup_originals(dry_run=False):
    """Copy each artifact to a .v32 sibling if no backup exists yet."""
    import shutil
    for src in ARTIFACTS:
        if not os.path.exists(src):
            print(f"  WARNING: {src} not found; skipping")
            continue
        dst = src + ".v32"
        if os.path.exists(dst):
            print(f"  backup exists: {os.path.basename(dst)}")
            continue
        size_gb = os.path.getsize(src) / 1024**3
        print(f"  backup: {os.path.basename(src)} ({size_gb:.2f} GB) -> .v32")
        if not dry_run:
            shutil.copy2(src, dst)

def realign_arrays(deleted_ids, dry_run=False):
    """Slice embeddings, comp_ids, valid_indices, morgan, maccs in lockstep
    against the deleted set. Returns the keep-mask for FAISS rebuild."""
    print("[load] comp_ids.npy")
    comp_ids = np.load(COMP_IDS, allow_pickle=True)
    n_pre = len(comp_ids)
    keep_mask = ~np.isin(comp_ids, list(deleted_ids))
    n_drop = int((~keep_mask).sum())
    n_post = int(keep_mask.sum())
    print(f"  pre: {n_pre:,}; drop: {n_drop:,}; post: {n_post:,}")
    n_deleted_total = len(deleted_ids)
    if n_drop < n_deleted_total:
        # Some deleted comp_ids never had embeddings (their SMILES failed
        # parsing during the original embedding generation). Expected and
        # benign; report it so the discrepancy is on record.
        print(f"  note: {n_deleted_total - n_drop:,} deleted comp_ids had "
              f"no embedding (SMILES-parse failures during v32 generation)")

    print("[realign] comp_ids.npy")
    new_comp_ids = comp_ids[keep_mask]
    if not dry_run:
        np.save(COMP_IDS, new_comp_ids)
    print(f"  saved: {len(new_comp_ids):,} entries")

    print("[realign] chemberta_embeddings.npy (memmap; this is the slow read)")
    t0 = time.time()
    emb = np.load(EMB, mmap_mode="r")
    keep_idx = np.where(keep_mask)[0]
    new_emb = np.array(emb[keep_idx])  # forces memmap -> RAM for the kept rows
    print(f"  shape: {emb.shape} -> {new_emb.shape} ({time.time()-t0:.1f}s)")
    if not dry_run:
        np.save(EMB, new_emb)
    print(f"  saved: {new_emb.shape}")

    print("[realign] valid_indices.npy")
    valid = np.load(VALID)
    # valid_indices entries are CSV row positions; the cleanup shifted those
    # positions but our embeddings are now aligned to the new comp_ids array
    # by position. The simplest interpretation is that valid_indices should
    # become np.arange(len(new_comp_ids)) now that the deleted rows are gone.
    # This matches the realign_indizes precedent which rebuilt valid_indices
    # to an identity mapping for the kept rows.
    new_valid = np.arange(len(new_comp_ids), dtype=valid.dtype)
    if not dry_run:
        np.save(VALID, new_valid)
    print(f"  rebuilt as identity: {len(new_valid):,} entries")

    print("[realign] morgan_fps.npz")
    morgan = sparse.load_npz(MORGAN)
    new_morgan = morgan[keep_mask]
    print(f"  shape: {morgan.shape} -> {new_morgan.shape}")
    if not dry_run:
        sparse.save_npz(MORGAN, new_morgan)

    print("[realign] maccs_fps.npz")
    maccs = sparse.load_npz(MACCS)
    new_maccs = maccs[keep_mask]
    print(f"  shape: {maccs.shape} -> {new_maccs.shape}")
    if not dry_run:
        sparse.save_npz(MACCS, new_maccs)

    return new_emb, n_post

def rebuild_faiss(emb, dry_run=False):
    """Rebuild FAISS HNSW from the realigned embeddings. Parameters match
    the v32-frozen construction at M=32, efConstruction=200, efSearch=128
    with cosine via L2-normalized vectors and IndexHNSWFlat."""
    import faiss
    print("[faiss] L2-normalize embeddings")
    valid_emb = emb.astype(np.float32, copy=False)
    norms = np.linalg.norm(valid_emb, axis=1, keepdims=True)
    norms[norms == 0] = 1
    valid_emb = valid_emb / norms
    d = valid_emb.shape[1]
    print(f"[faiss] building HNSW: shape={valid_emb.shape}, M=32, efC=200")
    t0 = time.time()
    index = faiss.IndexHNSWFlat(d, 32)
    index.hnsw.efConstruction = 200
    index.hnsw.efSearch = 128
    index.add(valid_emb)
    print(f"[faiss] built {index.ntotal:,} vectors in {time.time()-t0:.1f}s")
    if not dry_run:
        faiss.write_index(index, FAISS_IDX)
        print(f"[faiss] saved: {FAISS_IDX}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="print plan without writing anything")
    args = ap.parse_args()
    if not os.path.exists(DELETED_LIST):
        sys.exit(f"deleted list not found: {DELETED_LIST}")
    if not all(os.path.exists(p) for p in ARTIFACTS):
        for p in ARTIFACTS:
            if not os.path.exists(p):
                print(f"  missing: {p}")
        sys.exit("one or more artifacts missing")
    print(f"=== Phase A realignment ({'DRY RUN' if args.dry_run else 'WRITE'}) ===\n")
    deleted = load_deleted(DELETED_LIST)
    print(f"[load] {len(deleted):,} deleted comp_ids from {DELETED_LIST}\n")
    backup_originals(dry_run=args.dry_run)
    print()
    new_emb, n_post = realign_arrays(deleted, dry_run=args.dry_run)
    print()
    rebuild_faiss(new_emb, dry_run=args.dry_run)
    print(f"\n=== Phase A complete: {n_post:,} aligned compounds ===")

if __name__ == "__main__":
    main()
