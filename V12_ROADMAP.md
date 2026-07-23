
## Taxonomy resolution issues identified in v1.1

### Gramine-style kingdom misclassification

11 compounds show kingdom mismatch where source_organism indicates plant lineage (family name suffix -aceae) but resolved_taxonomy assigns kingdom=animal. Example: Gramine (THEO_0412143) — source_organism lists Acer rubrum, Acer saccharinum, Arundo donax (all plants) but resolved_taxonomy assigns Stylissa (sponge genus) via Porifera/Demospongiae.

Root cause hypothesis: name collision between species attestations in compound_taxonomy.token field. The same compound name may be reported by both a plant and a sponge organism in different source databases. The three-tier majority-vote resolution selects the sponge attribution despite the source_organism string indicating plants.

Diagnostic query:
```sql
SELECT comp_id FROM compounds c JOIN resolved_taxonomy rt ON rt.comp_id = c.comp_id
WHERE rt.kingdom = 'animal' AND c.source_organism ILIKE '%aceae%';
```

Returns 11 compounds.

Fix path (v1.2): re-run build_resolved_taxonomy_v2.py with name-collision detection: when source_organism contradicts compound_taxonomy.token, prefer source_organism as authoritative. Document the resolution change in v1.2 changelog.

### LOTUS species aggregation

compound_taxonomy schema lacks source_db column, preventing per-source aggregation diagnostics. Anne's reported issue cannot be scoped without first adding source attribution to compound_taxonomy, or by tracing back through ingestion logs to the LOTUS-specific compound list.

Fix path (v1.2): add source_db column to compound_taxonomy (or to a parallel compound_attribution table). Re-ingest with source attribution preserved. Then scope and address LOTUS-specific aggregation.


### Homepage extra-filter row prop-subselect parity

The homepage main-search now has full prop-subselect parity with the Search page (Chemical property type, classification/physicochemical/ADMET optgroups, type-to-jump autocomplete, re-open binding). The addFilterHome function that generates "+ Add filter" rows on the homepage still emits the narrow type list (without "Chemical property") and lacks the prop-dispatch wiring.

Fix path (v1.2): refactor the IIFE's per-row prop-dispatch into a helper function that binds handlers to any (typesel, propsel, input, inputWrap) tuple. Have both the main-search row and addFilterHome-generated rows call the helper. Estimated 100 lines of refactoring.

Workaround for v1.1: users wanting prop-based composition can use the Search page where extra-filter rows already have full parity. The homepage extra-filter rows still work for the existing 12 simple types (name, smiles, kingdom, etc.).
