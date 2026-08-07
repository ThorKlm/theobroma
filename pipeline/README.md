# THEOBROMA build pipeline (consolidated)

Ordered reproduction sequence for the THEOBROMA corpus. This directory holds the
consolidated SQL stage scripts (assembled verbatim from the applied one-off scripts);
the Python stages that were already single canonical scripts are referenced in place
under `scripts/`. See ../BUILD.md for the full provenance manifest and rationale.

The deposited database is the authoritative artifact; these scripts document the method.
Each SQL stage was originally staged with an operator COMMIT after inspecting verification
output; the trailing `COMMIT;` is left commented exactly as in the source. Run each stage,
review its verification output, then COMMIT.

## Order

| # | Stage | Script | Notes |
|---|-------|--------|-------|
| 01 | Source assembly | `convert.py` + `dedup_and_merge.py` + `merge_with_coconut.py` (repo root) | Per-source CSV -> InChIKey dedup (regional priority, most-permissive license) -> `data/theobroma_final.csv`. RDKit standardization at convert time. |
| 02 | Load | `scripts/load_data.py` (+ `scripts/load_sdf.py`) | DROP+CREATE `compounds`; load merged CSV. Destructive; runs first. |
| 03 | Compound-taxonomy | `scripts/populate_compound_taxonomy.py` | Download WCVP (Kew, May 2026 release) -> `wcvp_names`; tokenize source_organism; resolve genus/family; carry NCBI/WoRMS/WCVP lineage + native_distribution. LOTUS occurrences enter here (NCBI+WoRMS taxon_source). |
| 04 | Resolved taxonomy | `scripts/build_resolved_taxonomy_v2.py` | Majority-vote kingdom + three-tier joint lineage (one consistent tuple per compound). Canonical resolver (v2 supersedes v1). |
| 04b | Genus/family joint fix | `04b_resolved_taxonomy_joint_fix.sql` | Correct residual non-co-attested genus/family pairs from compound_taxonomy. |
| 05a | Plant-lineage backfill | `scripts/backfill_plant_lineage.py` | Fill phylum/class/order for plant families (APG IV / PPG I embedded map). |
| 05b | Coherence pass 1 | `05b_coherence_genus_family.sql` | BUILD `genus_lineage_ref`; snap incoherent genus->family lineages. |
| 05c | Coherence pass 2 | `05c_coherence_family_order.sql` | BUILD `family_lineage_ref`; snap family->order/class lineages. |
| 05d | APG taxclass relabel | `05d_apg_taxclass_relabel.sql` (verbatim `apg_taxclass_fix.sql`) | Deprecated magnoliopsida/liliopsida -> modern APG clades by taxorder. |
| 06 | Classification | `06_classification.sql` (from `npc_load_apply.sql`) | Load full-corpus NPClassifier re-run (GNPS2 API via `npc_api_batch.py`); overwrite np_pathway/superclass/class. Fixes the v134 scramble. |
| 06b | Pathway straggler fix | `06b_classification_pathway_fix.sql` | 27 Triterpenoid-superclass compounds: pathway Alkaloids -> Terpenoids. |
| 07 | Organism curation | `07_organism_curation.sql` | Build `source_organism_curated` (LOTUS-enriched, normalized binomials). Canonical; supersedes two earlier variants. |
| 08 | Region | `08_region.sql` | Multi-valued `compound_region_map` via WCVP 52->13 crosswalk + regional sources + superset. |
| 09 | License | `09_license.sql` (from `license_apply_v134.sql`) | Per-source most-restrictive; per-compound most-permissive across sources; writes license_tier + tier_rank. |
| 09b | Reconcile attestations | `09b_reconcile_attestations.sql` | Align attestation tiers to source_license_ref (resolved tiers invariant). |
| 10 | Enrichment (specify) | `scripts/generate_admet.py`, `scripts/compute_trust_score.py`, `scripts/regen_fingerprints.py`, `scripts/extend_artifacts_v33.py`/`build_nafm_index.py`; scaffold/SA/target from the Phase-B bundle (no repo script) | ADMET-AI, trust, Morgan/MACCS, ChemBERTa (seyonec/ChemBERTa-zinc-base-v1) + FAISS HNSW (M=32, efC=200, efS=128). GPU on compute host. Documented, not re-run for validation. |
| 11 | Late-source merge | `scripts/assemble_production_merge.py` + `scripts/tipdb_production_merge.sql` | 472 novel TipDB rows + 7,495 overlap enrichments; server-assigned THEO_ ids. |
| 12 | Serving caches | `scripts/update_taxonomy_cache.py`, `scripts/build_chem_tree_cache.py`, `scripts/update_statistics_cache.py`, `scripts/build_compound_region_map.py` | Regenerate `static/*.json` after data is final. |

## Source manifest
`sources.yaml` (repo root) is the authoritative 29-source manifest (license, region scope,
converter mapping); it is the citation-grade provenance file and ships with the deposit.

## Reference tables built along the way
- `genus_lineage_ref` (05b), `family_lineage_ref` (05c), `apg_clade_ref` (from apg_clade_ref.sql / supplement).
- `npc_class_parents`, `npc_super_parents` (NPClassifier ontology; no standalone builder script -- loaded from the ontology alongside stage 06 / from the enrichment bundle).

## Environment (v1.35)
Python 3.13.5; RDKit 2026.3.1; transformers 5.5.4 (seyonec/ChemBERTa-zinc-base-v1);
torch 2.11.0; scikit-learn 1.8.0; scipy 1.17.1; faiss-cpu 1.13.2; PostgreSQL 17.

## Notes
- These stages reflect the corrected v1.35 corpus. Numbers (esp. license distribution)
  supersede the v1.34 Zenodo deposit; update Zenodo/HF/manuscript to v1.35 independently.
- The deposited DB is the ground truth; where a consolidated script approximates an applied
  correction, the DB carries the exact result. See ../BUILD.md section 13 (corrections ledger).
