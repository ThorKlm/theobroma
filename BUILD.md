# THEOBROMA — Build and Provenance Manifest (BUILD.md)

**What this document is.** A provenance record of how the THEOBROMA corpus was produced,
from raw source databases to the deployed production database, together with a plan for
consolidating the scattered build and correction scripts into a small number of coherent
stage scripts. It is a capture-and-document manifest, not a claim of one-command,
bit-verified reproducibility.

**Authoritative artifact.** The deposited final database is the ground truth. Every
correction described here is already reflected in that database. The scripts document how
the state was produced; where a one-off correction is not folded into a consolidated stage
script, it is archived and listed in the corrections ledger (Section 13), so the method is
recoverable from the deposited database plus the versioned scripts at that state. Nothing
material is lost even where the scripts are not wired into a single runnable chain.

**Corpus at time of writing.** Internal version 1.35 (public release 1.0). ~1,132,805
compounds, 29 source databases, four resolved kingdoms plus an unresolved category, 13
macro-regions (multi-valued). Deduplicated by full 27-character InChIKey, stereoisomer
families exposed by the 14-character connectivity prefix. Live: theobroma.l3s.uni-hannover.de.
Repository: github.com/ThorKl/theobroma.

**Environment (pinned, v1.35).** Python 3.13.5; RDKit 2026.3.1; transformers 5.5.4 with
checkpoint seyonec/ChemBERTa-zinc-base-v1 (768-d, mean-pooled); torch 2.11.0 (CPU on serving,
GPU on the embedding-compute host); scikit-learn 1.8.0; scipy 1.17.1; pandas; xgboost and
Optuna on the training host; faiss-cpu 1.13.2 (HNSW). PostgreSQL 17, database `theobroma`,
role `theobroma`. Serving venv at ~/theobroma/venv. Offline SQL-driving scripts shell out to
psql (pattern: update_taxonomy_cache.py) so they run without the Python DB driver.

---

## 0. Conventions

- **Stage scripts** live in `scripts/` (the maintained pipeline) or were run as dated
  one-offs in the repo root (the v135 production build and subsequent corrections).
- **Consolidation target** notes, per stage, which scattered scripts should be merged into
  one canonical stage script and which are archived as-is.
- **Ground-truth principle.** Because the deposited DB reflects every applied correction, a
  consolidated script that approximates a correction is acceptable; the corrections ledger
  (Section 13) is the completeness checklist ensuring nothing is silently dropped.

---

## 1. Source assembly (pre-database, on the build host)

Raw source databases are converted to a uniform CSV schema, then deduplicated and merged
before loading into PostgreSQL.

- **Converters.** `convert.py` (unified converter for most sources: foodb, fdp, cmaup,
  npass, gbif enrichment, region normalization, and more), `simple_converter.py` (FDP
  scraper output split by sub-database), and source-specific converters
  (`convert_foodb.py`, `convert_microbial.py`, `convert_tipdb.py`, `convert_local.py`).
  Output: per-source CSVs in `data/converted/` with the uniform column set
  (comp_id, name, smiles, inchi, inchikey, source_db, all_sources, kingdom, region,
  source_organism, physicochemical descriptors, license_tier).
- **Dedup + merge.** `dedup_and_merge.py` deduplicates by InChIKey (computing missing keys
  from SMILES via RDKit MolToInchi -> InchiToInchiKey), applying regional-source priority
  (a SOURCE_PRIORITY table: IMPPAT/Phyto4Health/VIETHERB/ANPDB/... over COCONUT) and
  most-permissive-license-wins per InChIKey. `merge_with_coconut.py` performs the final
  cross-dedup against COCONUT and writes `data/theobroma_final.csv`.
- **Standardization.** RDKit canonical SMILES, neutralization, salt stripping (per the
  Zenodo methods summary). Full-stereochemistry InChIKey is the dedup key; the 14-char
  connectivity prefix defines stereoisomer families.
- **Raw inputs on disk** (`data/`): COCONUT SDF (coconut_prepped.sdf.gz, 636 MB), plus
  bacteria/fungi/marine_algae/plant_supplementary unique SDFs, compound_synonyms.csv, and
  the assembled theobroma_final.csv.

**Consolidation target:** already coherent. Document the three-step chain
(convert -> dedup_and_merge -> merge_with_coconut) as the source-assembly stage; no merging
needed. `sources.yaml` (repo root; verified) is the authoritative per-source manifest -- the
  29 sources with license, region scope, and converter mapping -- the citation-grade provenance
  file, and ships with the deposit.

---

## 2. Load into PostgreSQL

- **`scripts/load_data.py`** creates the `compounds` schema (DROP + CREATE, so this runs
  first and is destructive) and loads the converted/merged CSVs; carries a SOURCE_META table
  mapping each of the ~29 sources to (scope, region, default license). `scripts/load_sdf.py`
  handles SDF-format sources.
- Physicochemical descriptors (mw, logp, tpsa, hba, hbd, n_rings, rotatable_bonds) are
  loaded from the converter output (computed via RDKit at convert time).

**Outputs:** `compounds` (base columns). **Consolidation:** single script already; canonical.

---

## 3. Compound-taxonomy population (organism resolution)

- **`scripts/populate_compound_taxonomy.py`** downloads the WCVP names dump from Kew
  (sftp.kew.org; May 2026 WCVP release used for v33), loads it into `wcvp_names`
  (~1.4M rows), tokenizes each compound's source_organism on '; ', and resolves each token
  to genus/family via WCVP (accepted directly; synonyms followed to the accepted row).
  NPO/NPC identifier pollution (NPAtlas internal IDs) is filtered at populate time.
- The table also carries NCBI and WoRMS lineage JSON and provenance (`taxon_source`), which
  is how LOTUS occurrences enter: LOTUS structure-organism pairs resolved through NCBI+WoRMS
  are tagged `NCBI+WoRMS` (226,558 rows, the LOTUS set).

**Outputs:** `wcvp_names`; `compound_taxonomy` (1,504,425 rows / 355,970 compounds; 22 cols:
token, accepted_taxon, genus, family, ncbi_/worms_/wcvp_lineage JSON, native_distribution,
phylum_any/class_any/order_any/kingdom_any, taxon_source).
**taxon_source breakdown:** WCVP 706,405; blank 402,576; NCBI+WoRMS 226,558 (LOTUS);
NCBI 152,301; WoRMS 14,032; WCVP+WoRMS 2,553.
**Consolidation:** single script; canonical. (WCVP is an external download; pin the release.)

---

## 4. Resolved taxonomy (kingdom + joint lineage)

- **Canonical resolver: `scripts/build_resolved_taxonomy_v2.py`** (v2, written after v1;
  v1 `build_resolved_taxonomy.py` is superseded, keep for history only).
- **Method:** majority-vote kingdom across attestations (kingdom-mapping table for post-2022
  NCBI clade names), with a three-tier lineage fallback. Tier 1 (mappable NCBI/WoRMS kingdom)
  and Tier 2 (WCVP-only plants) both select, per compound, a single most-specific lineage row
  via `ROW_NUMBER() ... ORDER BY specificity DESC, family ASC` and take row 1, so genus and
  family come from ONE consistent organism lineage (this is the joint resolution). Tier 3
  (no attestations) uses compounds.kingdom with NULL lineage. Secondary kingdoms aggregated.
  Swaps `resolved_taxonomy` in place, rebuilds lowercase indexes, GRANTs to the role.
- **Correction reflected:** the independent-per-rank voting bug (impossible genus/family
  pairs, e.g. Gramine Dioscorea + Poaceae) is absent from v2's tuple selection. A later
  one-off (`fix_genus_family_joint_apply.sql`, this session) corrected residual live
  inconsistencies; the live DB now shows zero non-co-attested genus/family pairs.

**Outputs:** `resolved_taxonomy` (comp_id, kingdom, secondary_kingdoms[], phylum, taxclass,
taxorder, family, genus). **Checkpoints:** non-co-attested genus/family pairs = 0.
**Consolidation:** canonical resolver = v2. Append the joint-fix + the v135 coherence fixes
(below) as a documented "taxonomy corrections" companion, OR archive them beside v2 with a
note that the DB reflects them. Ground-truth principle applies.

---

## 5. Lineage backbone, backfill, and coherence

- **Reference backbones.** `genus_lineage_ref` (genus_lc -> family/taxorder/taxclass/phylum/
  apg_clade + support), `family_lineage_ref` (family_lc -> ...), `apg_clade_ref`
  (family -> clade). Built from the corpus + APG IV / PPG I references; ref_build.log records
  a resolved_name_ref/genus_ref build (30,283 / 6,810 rows; WCVP name-resolution helper tables).
  genus_lineage_ref and family_lineage_ref are built by the coherence scripts (stages 05b/05c);
  the resolved_name_ref/genus_ref helper is part of the populate/ref build (ref_build.log).
  Non-blocking, recoverable from the tables if ever needed.
- **Backfill.** `scripts/backfill_plant_lineage.py` fills phylum/class/order for plant
  families from a large embedded APG IV / PPG I family map (Streptophyta phylum, Magnoliopsida
  /Liliopsida/Polypodiopsida classes), via an UPDATE FROM a temp `_plant_lineage` table.
- **Coherence fixes.** `v135_coherence_fix.sql` and `v135_coherence_fix2.sql`
  (multi-organism lineage stitching: family <-> order/class), `apg_taxclass_fix.sql` +
  `rebuild_apg_archive.sql` (replace deprecated magnoliopsida/liliopsida class labels),
  `cleanup_taxonomy_majority.py` (majority resolution of residual inconsistencies).
- **Consistency guarantee (audited):** backfillable gaps = 0; genus->family consistency =
  6,134 single-family genera, 17 multi-family (14 cross-kingdom homonyms + 3 accepted
  revisions), zero errors. Residual non-backfillable stragglers: genus-no-family 529,
  family-no-order 714, order-no-class 536, class-no-phylum 246.

**Consolidation:** fold backfill + coherence fixes + apg taxclass fix into one canonical
"lineage backbone + coherence" script, OR archive the v135_coherence/apg one-offs beside
backfill_plant_lineage.py with a note. Ground-truth principle applies.

---

## 6. Chemical classification

- **NPClassifier (canonical).** `npc_load_apply.sql` loads `npc_output.csv` (full-corpus
  re-run via the canonical GNPS2 API, `npc_api_batch.py`) into `npc_reclassified`
  (comp_id, np_pathway, np_superclass, np_class), then applies to `compounds`. This is the
  classification of record and explicitly fixes the v134 residual classification scramble
  (Anne's point). Staged load -> coverage -> archive -> update -> verify.
- **ClassyFire.** `classyfire_superclass` carried from source where present.
- **Inferred-class fallback.** `inferred_class` / `inferred_class_source` / `inferred_confidence`
  from an XGBoost ensemble over a 2,471-d feature vector (Morgan + MACCS + RDKit descriptors +
  ChemBERTa pooled embedding), trained on 432 curated classes; Optuna 50-trial HPO, seed 1,
  best macro-F1 0.7401; predictions split into high-confidence and exploratory tiers.
- **Ontology tables.** `npc_class_parents` (class -> superclass -> pathway),
  `npc_super_parents` (superclass -> pathway). 7 pathways / 76 superclasses / 687 classes.
  No standalone builder script on the server: these were populated from the NPClassifier
  ontology alongside the stage-06 load (inline / from the enrichment bundle). Treated as
  provided reference data; a rebuild reloads the ontology, it is not regenerated by a repo script.
  Six classes map to two superclasses (biosynthetic hybrids); 17 superclasses to two pathways.
- **Pathway-straggler fix (this session).** `fix_pathway_stragglers_apply.sql`: 27 compounds
  (21 Fusidane, 6 Quassinoids) with impossible pathway/superclass combinations corrected to
  follow the ontology; residual = 0.

**Consolidation:** canonical = npc_load_apply.sql + the ontology-table builder. Append the
27-row straggler fix as an ontology-consistency step. Archive the exploratory npc_batch
variants.

---

## 7. Source-organism curation (includes the LOTUS merge)

- **Canonical (final) script:** `rebuild_curated_lotus_merge_apply.sql` (this session).
  Builds `source_organism_curated` from the union of compound_taxonomy tokens and raw
  source_organism tokens, normalized to clean `Genus species` binomials via leading-binomial
  extraction (regexp_match on `^([A-Za-z]+)\s+([A-Za-z]+)`; exclude epithets sp/spp/var/
  subsp/cf/aff/f/x/l/ssp/nov/gen; require epithet length > 1; initcap genus, lower epithet),
  deduplicated, ordered resolved-genus-first; NULL where no binomial (template falls back to
  raw source_organism). This surfaces the LOTUS occurrences into the display.
- **Superseded earlier variants** (same tiering idea, no LOTUS/normalization):
  `build_source_organism_curated_apply.sql` (initial), `rebuild_source_organism_curated_apply.sql`
  (rebuilt on corrected taxonomy). `v135_stage6_organisms.sql` builds the COCONUT-only
  `coconut_organisms` occurrence table (a separate input, not the curated field).
- **Effect (audited):** curated non-null 163,827 -> 335,530; Gramine 36 -> 49; Morphine leads
  Papaver; 10,480 nulled-to-fallback with identical raw, no data loss.

**Consolidation:** single canonical curation script = the LOTUS-merge version; supersedes the
two earlier variants (archive them). This is the one stage whose canonical logic lived only in
root one-offs; it is the primary consolidation win.

---

## 8. Geographic region derivation (multi-valued, DONE)

- **Canonical: `region_derive_stage2_persist.sql`.** DROPs + CREATEs `compound_region_map`
  (comp_id, macro_region), multi-valued, then dedups to distinct (comp_id, macro_region).
- **Method:** per compound, the SET of { crosswalk(WCVP native_distribution regions across
  all organism tokens) } U { fixed-scope regional-source tags from all_sources } U
  { existing compounds.region if not global }. The 52-subregion -> 13-macro-region crosswalk
  is an embedded `sub2macro` VALUES table (the planned WCVP crosswalk, applied).
- **Status:** COMPLETE. Multi-valued live: 127,992 compounds with 1 region, 41,548 with 2,
  25,461 with 3, long tail beyond. `build_compound_region_map.py` builds the country-level
  JSON for the browse map (macro-region -> ISO3 expansion).

**Consolidation:** canonical = region_derive_stage2_persist.sql (+ stage1 sampling as its
input study). Archive stage2_backfill variant. Region is not a pending task.

---

## 9. License tiers

- **Canonical: `license_apply_v134.sql`** (most-restrictive-wins per compound across attested
  sources; the erroneous CC0 override discarded). Archives prior labels, staged transaction.
  `reconcile_attestations.sql` reconciles per_source_license_attestation to source_license_ref.
- **Rule (two-level):** per source, most-restrictive applicable license; per compound
  (a structure being a fact reportable by many sources), most-permissive across its sources.
  Tier order (source_license_ref.tier_rank, confirmed): 0 CC0, 1 CC BY 4.0, 2 CC BY-NC 4.0,
  3 CC BY-NC-SA 4.0, 4 CC BY-NC-ND 4.0, 5 Unspecified. CMNPD removed in v32 (share-alike).
- **Live distribution (v1.35):** CC0 1,013,320; CC BY-NC 4.0 84,956; Unspecified 27,051;
  CC BY-NC-ND 4.0 3,758; CC BY 4.0 3,720. Commercial-use (CC0 + CC BY 4.0) = 1,017,040.

**IMPORTANT (deposition):** this distribution differs substantially from the deposited Zenodo
v1.34 doc (which reported CC BY 4.0 994,133 / CC0 12,044). The license re-resolution changed
the headline numbers. The Zenodo v2, HuggingFace, and manuscript updates (done independently,
later) must use the corrected v1.35 distribution. Recorded here so the discrepancy is not lost.

**Consolidation:** canonical = license_apply_v134.sql + reconcile_attestations.sql.

---

## 10. Computed properties and enrichment (specify; heavy, GPU on compute host)

- **ADMET.** `scripts/generate_admet.py` (ADMET-AI, Chemprop D-MPNN ensemble, 104 endpoints;
  the live `admet` table persists a broad set — Zenodo deposit ships 29 ADMET + 12 Tox21 +
  6 RDKit descriptors, 47 columns, 99.999% coverage).
- **Trust score.** `scripts/compute_trust_score.py` (heuristic composite; coefficients
  regression-recovered, R^2 0.93, RMSE 0.06; clamped [0,1]; overwrites corpus-wide).
- **Fingerprints.** `scripts/regen_fingerprints.py` (Morgan 2048-bit radius 2 / ECFP4, MACCS
  167-bit via RDKit; CSR .npz + fingerprint_indices.npy).
- **Embeddings + index.** ChemBERTa (seyonec/ChemBERTa-zinc-base-v1, 768-d mean-pooled);
  FAISS IndexHNSWFlat (M=32, efConstruction=200, efSearch=128), built by
  `scripts/extend_artifacts_v33.py` / `build_nafm_index.py`. NaFM variant scripts present.
- **Scaffolds.** Bemis-Murcko (`scaffolds` table). **SA score.** Ertl-Schuffenhauer.
- **Target predictions.** SEA-style vs ChEMBL v34 actives.

**Consolidation:** these are documented-for-reproducibility (specify, not consolidate). They
run on the GPU compute host; the vast.ai Phase-B bundle (data/tipdb_phase_b) carried the
late-source enrichment. Note (confirmed by audit): scaffold, SA-score, and target-prediction
steps have NO standalone builder script in `scripts/` on the server; they were produced in
the Phase-B bundle and merged at stage 11. ADMET/trust/fingerprints/embeddings do have repo
scripts (named above). Pin the model/version strings above.

---

## 11. Late-source (TipDB) production merge

- **`scripts/assemble_production_merge.py`** combines the local merge output
  (theobroma_v33_merged.csv), the vast.ai Phase-B bundle (ChemBERTa/Morgan/MACCS/scaffolds/
  npc/admet/sa), and the post-cleanup corpus (theobroma_v33_clean.csv) to produce 472 novel
  rows + 7,495 overlap InChIKeys. `scripts/tipdb_production_merge.sql` inserts the novel rows
  (server-assigned THEO_ ids) and enriches overlaps.

**Consolidation:** canonical pair; keep as the late-merge stage.

---

## 12. Serving caches (run last, after data is final)

- `scripts/update_taxonomy_cache.py` -> static/taxonomy_cache.json (+ per-kingdom variants).
- `scripts/build_chem_tree_cache.py` -> static/chem_tree.json (this session; committed).
- `scripts/update_statistics_cache.py` -> static/statistics_cache.json.
- `scripts/build_compound_region_map.py` -> static/compounds_by_country.json.
- static/histograms.json (tracked in git), static/linear_tree_global.json.
- `scripts/update_visitors_map.py` (operational, not corpus).

**Consolidation:** list all cache builders as the final "regenerate caches" stage so a deploy
has them before first request.

---

## 13. Corrections ledger (completeness checklist)

Every data-modifying script found in the audit, with where its result lives. The deposited DB
reflects all of these; this table ensures none is silently dropped in consolidation.

| Correction | Script(s) | Consolidation |
|---|---|---|
| Genus/family joint resolution | build_resolved_taxonomy_v2.py (method) + fix_genus_family_joint_apply.sql (residual) | canonical v2 + archive fix |
| Secondary-kingdom dedup | fix_secondary_kingdoms_apply.sql | archive; DB reflects |
| Lineage backfill | backfill_plant_lineage.py | canonical |
| Lineage coherence stitching | v135_coherence_fix.sql, v135_coherence_fix2.sql | fold or archive |
| APG taxclass relabel | apg_taxclass_fix.sql, rebuild_apg_archive.sql | fold or archive |
| Majority cleanup | cleanup_taxonomy_majority.py | canonical |
| NPClassifier re-run (scramble fix) | npc_load_apply.sql, npc_api_batch.py | canonical |
| Pathway stragglers (27) | fix_pathway_stragglers_apply.sql | append to classification |
| Organism curation + LOTUS | rebuild_curated_lotus_merge_apply.sql | canonical (supersedes 2 earlier) |
| Region multi-valued + crosswalk | region_derive_stage2_persist.sql | canonical |
| License re-resolution | license_apply_v134.sql, reconcile_attestations.sql | canonical |
| Inorganic removal (199) | remove_inorganic.sql | archive; DB reflects |
| Placeholder-name nulling | null_placeholder_names.sql | archive; DB reflects |
| Anne class/organism fixes | anne_fix_class_organism.sql, anne_fix_precise.sql | archive; superseded by npc + curation |
| Filter-options label split | app.py (runtime) | code, committed |
| Superclass search dispatch | app.py (runtime) | code, committed |

---

## 14. Deposition dump scope (for the independent Zenodo/HF update later)

The database contains heavy archival cruft that must be EXCLUDED from any deposition dump.

- **Include (canonical):** compounds, resolved_taxonomy, compound_taxonomy,
  compound_region_map, admet, scaffolds, compound_synonyms, genus_lineage_ref,
  family_lineage_ref, apg_clade_ref, npc_class_parents, npc_super_parents, source_license_ref, per_source_license_attestation; plus the
  precomputed vectors (chemberta_embeddings, faiss_hnsw, morgan/maccs fps, comp_ids).
- **Exclude:** all `*_v134_archive`, `*_v135_preswap_*`, `*_pre_*_2026*`, `*_old_bad`,
  `*_broken_lineage`, `*_buggy`, `resolved_taxonomy_pre_v2`, staging tables (v135_input,
  v135, v135_compid_map, wcvp_distribution_staging, classified_ik, coconut_gate*), and
  `access_log` (operational + privacy-sensitive).
- The existing Zenodo v1.34 deposit shipped six CSV tables + vectors + classification
  artifacts + sources.yaml + dedup_summary.json + MANIFEST + SHA256SUMS; mirror that file set
  with corrected v1.35 content.

---

## 15. Consolidation plan (the deliverable of this phase)

Produce a small set of canonical stage scripts, grouping the scattered one-offs, and an
`archive/` of the raw one-offs with a short index. Proposed canonical set:

1. `01_source_assembly` — convert*.py + dedup_and_merge.py + merge_with_coconut.py (document).
2. `02_load.py` — load_data.py (+ load_sdf.py).
3. `03_compound_taxonomy.py` — populate_compound_taxonomy.py (WCVP download + resolve).
4. `04_resolved_taxonomy.py` — build_resolved_taxonomy_v2.py + genus/family joint fix.
5. `05_lineage_backbone.py` — backfill_plant_lineage.py + coherence + apg taxclass + majority.
6. `06_classification.sql` — npc_load_apply.sql + pathway-straggler fix (+ ontology builder).
7. `07_organism_curation.sql` — rebuild_curated_lotus_merge_apply.sql (canonical).
8. `08_region.sql` — region_derive_stage2_persist.sql.
9. `09_license.sql` — license_apply_v134.sql + reconcile_attestations.sql.
10. `10_enrichment` — generate_admet / trust / fingerprints / embeddings / index (specify).
11. `11_late_merge` — assemble_production_merge.py + tipdb_production_merge.sql.
12. `12_caches` — all static cache builders.

Everything else (backups, dryruns, superseded variants, the fix_help_*/fix_api_* app patchers,
the eval_*/explore_* shell scripts) goes to `archive/` with a one-line index entry each.
The corrections ledger (Section 13) is the checklist: every row must map to a canonical script
or an archived-and-indexed one-off.

---

## 16. Open items (not blocking; capture-and-document goal)

- Optional future nicety: a scoped clean-room rebuild of the data-core stages to confirm the
  consolidated scripts reproduce the checkpoints. Not required, since the deposited DB is the
  ground truth and the scripts + DB together preserve the method.
- (Resolved) The genus_lineage_ref / family_lineage_ref builders are the coherence scripts:
  stage 05b builds genus_lineage_ref, stage 05c builds family_lineage_ref.
- Independent, later: Zenodo v2, HuggingFace refresh, manuscript/preprint update to the
  corrected v1.35 numbers (esp. the license distribution), retire the maintenance banner.

---

*Provenance-and-consolidation manifest. The deposited database is authoritative; scripts
document the method. Section 13 is the completeness checklist; Section 15 is the consolidation
plan.*
