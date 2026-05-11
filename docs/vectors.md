# Vector Artifacts (Space C, Space F, Space E)

The three artifact arrays in this directory have related but non-identical
row counts, reflecting the cascade of filters during corpus preparation:

- **Space C (corpus)**: `comp_ids.npy` indexes every row in the compounds
  table that has artifact-bearing data. As of v33 augmented, 1,077,361 entries
  (1,074,567 from the realigned post-cleanup baseline plus 472 TIPdb novel,
  with internal cleanup adjustments).

- **Space F (fingerprints)**: `morgan_fps.npz` and `maccs_fps.npz` index the
  subset of Space C that produced valid RDKit fingerprints. SMILES that
  failed RDKit parsing are excluded.

- **Space E (embeddings)**: `chemberta_embeddings.npy` and `faiss_hnsw.index`
  index the subset of Space F that also produced ChemBERTa embeddings
  (some SMILES exceed the tokenizer's max length or contain unsupported
  atoms; these are dropped).

The `valid_indices.npy` array maps Space E back to Space C positions.
At similarity-query time, `scripts/similarity.py` looks up
`comp_ids[valid_indices[hit_index]]` to recover the production compound
identifier for a FAISS or fingerprint hit.

## Files

- `chemberta_embeddings.npy`: float32 (n_e, 768) ChemBERTa mean-pooled CLS embeddings
- `comp_ids.npy`: object array of THEO_xxxxxxx strings, length n_c
- `valid_indices.npy`: int64 array of Space C positions for which embeddings exist, length n_e
- `morgan_fps.npz`: scipy.sparse CSR (n_f, 2048) Morgan-r2 fingerprints
- `maccs_fps.npz`: scipy.sparse CSR (n_f, 167) MACCS keys
- `faiss_hnsw.index`: FAISS IndexHNSWFlat, M=32, efConstruction=200, efSearch=128, n_e vectors
- `admet_predictions.csv`: full 104-endpoint ADMET-AI predictions (only 4 surfaced in the admet DB table)

## Regeneration

Realigning Space C after a compounds-table cleanup: `scripts/realign_v33_cleanup.py`.
Extending all three spaces with new compounds: `scripts/extend_artifacts_v33.py`.
