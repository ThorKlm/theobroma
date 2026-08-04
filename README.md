# THEOBROMA

An open multi-kingdom natural products database with per-compound license tiers for legally-aware drug discovery.

**Live:** [https://theobroma.l3s.uni-hannover.de](https://theobroma.l3s.uni-hannover.de)
**Repository:** [github.com/ThorKlm/theobroma](https://github.com/ThorKlm/theobroma)
**Dataset on Hugging Face:** [huggingface.co/datasets/ThorKl/theobroma](https://huggingface.co/datasets/ThorKl/theobroma)

## Overview

THEOBROMA aggregates 1,132,805 natural product compounds from 29 source databases spanning four biological kingdoms (plant, animal, fungi, bacteria) plus an unresolved category, across 13 resolved geographic regions. Compounds are deduplicated by full 27-character InChIKey, preserving stereoisomers as distinct entries: 486,032 distinct molecular skeletons, of which 314,742 are multi-member stereoisomer families. Each compound carries a per-record license tier so the corpus can be filtered for commercial-use compatibility. The full pipeline integrates ChemBERTa molecular embeddings, Bemis-Murcko scaffold classification, NPClassifier and ClassyFire chemical-class assignments, ADMET-AI predictions (47 endpoints persisted and queryable in the searchable database), synthetic accessibility scores, and SEA-style target predictions against ChEMBL v34 actives.

## Features

Search by compound name, SMILES, InChIKey, source organism, geographic region, kingdom, source database, or chemical class. Multi-filter search combines criteria with AND logic. Browse with kingdom, source, region, and license filters plus a named-only toggle. Each compound has a detail page with structure, properties, classification, provenance across all source databases, and external links. Statistics dashboard shows kingdom, source, region, and license distributions. Bulk export as CSV. JSON API for programmatic access. Similarity search via Morgan or MACCS Tanimoto plus ChemBERTa cosine, the last accelerated by FAISS HNSW indexing. Substructure search via Morgan pre-screen plus RDKit matching. Scaffold browser groups compounds by their Bemis-Murcko core.

## Source databases

COCONUT 2.0, LOTUS v11, FooDB, NPASS 3.0, HERB 2.0, TM-MC 2.0, IMPPAT 2.0, CSIRO Australian NP, ANPDB, NPAtlas, phytochemdb, MicotoXilico, StreptomeDB, MIBiG 4.0, EMNPD, MeFSAT, CyanoMetDB, MycoCentral, NaturAr, LMDB_Lichen, AMDB, TIPdb-3D, Phyto4Health, AfroDb, CMAUPv2, SANCDB, CMDB_Cereals, TMDB_Trichoderma, SMDB_Spice.

The CMNPD database, originally part of the corpus, was removed in v32 due to its CC BY-NC-SA share-alike clause; CMNPD-exclusive compounds were dropped and CMNPD provenance was stripped from multi-source rows. TIPdb-3D was integrated in v33 from the 2015 archive following decommissioning of the original server (structures only; full ethnobotanical metadata pending a successor TIPdb release). LOTUS v11 (April 2026 release) was added in v33 contributing 53,862 novel compounds and enriching 173,436 existing rows with cross-source provenance.

## License tiers

Each compound carries a `license_tier` describing the terms under which its chemical **structure** may be redistributed, resolved by a two-level rule. **First**, each source database is assigned the most restrictive license applicable to its contents: sources that are non-commercial, share-alike, or whose terms cannot be resolved at the per-compound level are treated conservatively at that tier, or recorded as *Unspecified* when no license can be determined. **Second**, because a chemical structure is a fact that may be independently reported by several databases, a compound found in multiple sources takes the **most permissive** license among them, since the structure is genuinely available under that license from at least one source. The assignment is thus conservative where provenance is singular or ambiguous, and permissive only where a concrete more-permissive source exists. Tier order (permissive to restrictive): CC0 < CC BY 4.0 < CC BY-NC 4.0 < CC BY-NC-SA 4.0 < CC BY-NC-ND 4.0 < Unspecified.

Across the 1,132,805 compounds:

- **CC0** (1,013,320; 89.45%): COCONUT 2.0, LOTUS v11, MIBiG 4.0, SANCDB.
- **CC BY 4.0** (3,720; 0.33%): AfroDb, CyanoMetDB, EMNPD, LanaPDB, MeFSAT, TIPdb.
- **CC BY-NC 4.0** (84,956; 7.50%): ANPDB, CMAUP, CSIRO, FooDB, NPASS 3.0, NPAtlas, Phyto4Health, TCMBank, YMDB.
- **CC BY-NC-ND 4.0** (3,758; 0.33%): IMPPAT 2.0.
- **Unspecified** (27,051; 2.39%): HERB 2.0, TM-MC 2.0, StreptomeDB, phytochemdb, MicotoXilico, MycoCentral, AMDB, ConMedNP, NaturAr, LMDB, and specialized/regional collections whose license could not be resolved; excluded from both commercial and non-commercial filters pending source-license confirmation.

**1,017,040 compounds (89.78%) permit commercial use** (CC0 or CC BY 4.0). The `license_tier` reflects **structure** redistribution; reuse of the full integrated annotation record (source organism, references, computed metadata) may additionally require honoring the contributing sources' terms for those fields. The CMNPD database, originally part of the corpus, was removed in v32 due to its CC BY-NC-SA share-alike clause.

## Deployment

```bash
# Clone
git clone https://github.com/ThorKlm/theobroma.git
cd theobroma

# Setup
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Database
sudo -u postgres createuser theobroma
sudo -u postgres createdb -O theobroma theobroma

# Load data (downloads from Hugging Face)
python scripts/download_hf.py --token YOUR_HF_TOKEN
python scripts/load_data.py

# Run
gunicorn -b 0.0.0.0:5000 app:app
```

The L3S Research Center hosts the public instance on a single-host Xeon node (PostgreSQL 17), with a five-year minimum hosting commitment.

## API

```bash
curl "https://theobroma.l3s.uni-hannover.de/api/search?q=curcumin&type=name&limit=10"
```

Full OpenAPI specification at [/api](https://theobroma.l3s.uni-hannover.de/api). Endpoints include compound lookup, search across all field types, autocomplete, similarity, substructure, stereoisomer enumeration, filter-options listing, and streaming bulk export. Source-attribution and license fields ship with every response.

## License

Web application code is MIT. Compound data carries per-record license tiers resolved by the two-level rule described in the License tiers section above (conservative within an ambiguous source, most-permissive across multiple sources); see the per-compound `license_tier` field.

## Citation

```
Klamt, T. et al. (2026). THEOBROMA: an open multi-kingdom natural products
database with per-compound license tiers for legally-aware drug discovery.
Submitted to Nucleic Acids Research.
```

Please also cite the original source databases for any compounds you use; per-source citations and licenses are listed at [/sources](https://theobroma.l3s.uni-hannover.de/sources).

## Hosting

[L3S Research Center](https://www.l3s.de), Leibniz University Hannover.

## Contact

For bug reports and feature requests, open an issue on [GitHub](https://github.com/ThorKlm/theobroma/issues). For other inquiries, contact Thor Klamt at [Thor.Klamt@gmail.com](mailto:Thor.Klamt@gmail.com).
