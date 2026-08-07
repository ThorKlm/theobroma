# Archive index — superseded / one-off scripts

These scripts were applied to produce the live corpus and are preserved for provenance.
Their results are baked into the deposited database. They are NOT part of the runnable
consolidated pipeline (see ../pipeline/README.md); each is either superseded by a canonical
stage script or a standalone correction captured in the artifact.

## Superseded by a canonical pipeline stage
| Script | Superseded by | What it did |
|--------|---------------|-------------|
| build_source_organism_curated_apply.sql | pipeline/07 | initial organism curation (no LOTUS/normalization) |
| rebuild_source_organism_curated_apply.sql | pipeline/07 | curation rebuilt on corrected taxonomy (pre-LOTUS) |
| build_source_organism_curated_dryrun.sql | pipeline/07 | dry-run of the above |
| rebuild_curated_lotus_merge_dryrun.sql | pipeline/07 | dry-run of the LOTUS merge |
| fix_genus_family_joint_dryrun.sql | pipeline/04b | dry-run of the joint fix |
| fix_pathway_stragglers_dryrun.sql | pipeline/06b | dry-run of the straggler fix |
| region_derive_stage2_backfill.sql | pipeline/08 | earlier region build variant |
| region_derive_stage1_sample.py | pipeline/08 | region derivation sampling study (input) |
| build_resolved_taxonomy.py (v1) | pipeline/04 (v2) | independent-per-rank resolver (buggy pairs) |

## Standalone corrections (captured in the artifact; not folded)
| Script | What it did |
|--------|-------------|
| fix_secondary_kingdoms_apply.sql | remove primary kingdom redundantly present in secondary_kingdoms[] |
| cleanup_taxonomy_majority.py | majority resolution of residual lineage inconsistencies |
| remove_inorganic.sql | remove 199 non-NP entries (<=3 atoms: single elements, tiny inorganics) |
| null_placeholder_names.sql | null literal placeholder strings in compounds.name |
| rebuild_apg_archive.sql | rebuild pre-APG archive from current state + mapping |
| anne_fix_class_organism.sql | early classification/organism fixes (superseded by 06 + 07) |
| anne_fix_precise.sql | organism truncation diagnostics/fix (superseded by 07) |
| npc_batch.py / npclassifier_baseline.py | exploratory NPClassifier variants (canonical: npc_api_batch.py) |
| license_permissive_v134.sql | discarded most-permissive license option (not applied) |
| restore_v134.sql | full restore-to-clean-v1.34 utility |
| v135b_merge_relicense.sql | v1.34+v1.35 merge / relicense staging |

## App/runtime patchers (already in committed app code)
The many fix_help_*.py, fix_api_*.py, fix_region_*.py, fix_similarity_*.py, patch_*.py,
and revert_*.py scripts modified app.py/templates and are reflected in the committed
application code (commit 29c55f8 and predecessors). They are not data-pipeline stages.

## Evaluation / exploration (not build)
eval_*.sh, explore_*.sh/.sql, investigate_*.sh/.sql, verify_*.sh/.sql, validation_probes*.sql,
version_scan.sh, survey_rebuild_inputs.sh — diagnostics and audits, not corpus-producing.
