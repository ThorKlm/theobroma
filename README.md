# THEOBROMA

An open multi-kingdom natural products database with per-compound license tiers for legally-aware drug discovery.

**Live:** [https://theobroma.l3s.uni-hannover.de](https://theobroma.l3s.uni-hannover.de)
**Repository:** [github.com/ThorKlm/theobroma](https://github.com/ThorKlm/theobroma)
**Dataset on Hugging Face:** [huggingface.co/datasets/ThorKl/theobroma](https://huggingface.co/datasets/ThorKl/theobroma)

## Overview

THEOBROMA aggregates 1,133,004 natural product compounds from 29 source databases spanning five biological kingdoms (plant, fungi, bacteria, animal, marine) plus a multi-kingdom category, across 13 resolved geographic regions. Compounds are deduplicated by full 27-character InChIKey (preserving stereoisomers as distinct entries: 486,215 stereoisomer families). Each compound carries a per-record license tier so the corpus can be filtered for commercial-use compatibility. The full pipeline integrates ChemBERTa molecular embeddings, Bemis-Murcko scaffold classification, NPClassifier and ClassyFire chemical-class assignments, ADMET-AI predictions across 104 endpoints (four persisted in the searchable database), synthetic accessibility scores, and SEA-style target predictions against ChEMBL v34 actives.

## Features

Search by compound name, SMILES, InChIKey, source organism, geographic region, kingdom, source database, or chemical class. Multi-filter search combines criteria with AND logic. Browse with kingdom, source, region, and license filters plus a named-only toggle. Each compound has a detail page with structure, properties, classification, provenance across all source databases, and external links. Statistics dashboard shows kingdom, source, region, and license distributions. Bulk export as CSV. JSON API for programmatic access. Similarity search via Morgan or MACCS Tanimoto plus ChemBERTa cosine, the last accelerated by FAISS HNSW indexing. Substructure search via Morgan pre-screen plus RDKit matching. Scaffold browser groups compounds by their Bemis-Murcko core.

## Source databases

COCONUT 2.0, LOTUS v11, FooDB, NPASS 3.0, HERB 2.0, TM-MC 2.0, IMPPAT 2.0, CSIRO Australian NP, ANPDB, NPAtlas, phytochemdb, MicotoXilico, StreptomeDB, MIBiG 4.0, EMNPD, MeFSAT, CyanoMetDB, MycoCentral, NaturAr, LMDB_Lichen, AMDB, TIPdb-3D, Phyto4Health, AfroDb, CMAUPv2, SANCDB, CMDB_Cereals, TMDB_Trichoderma, SMDB_Spice.

The CMNPD database, originally part of the corpus, was removed in v32 due to its CC BY-NC-SA share-alike clause; CMNPD-exclusive compounds were dropped and CMNPD provenance was stripped from multi-source rows. TIPdb-3D was integrated in v33 from the 2015 archive following decommissioning of the original server (structures only; full ethnobotanical metadata pending a successor TIPdb release). LOTUS v11 (April 2026 release) was added in v33 contributing 53,862 novel compounds and enriching 173,436 existing rows with cross-source provenance.

## License tiers

Permissive open (1,031,332 compounds, 91.03%) carry CC BY 4.0 from their source databases. Public domain (12,044, 1.06%) carry CC0. Non-commercial (87,119, 7.69%) carry CC BY-NC 4.0. Unspecified (2,509, 0.22%) await source-license confirmation; the license filter excludes these from both commercial and academic categories pending resolution.

When a compound appears in multiple source databases, the most permissive license applies. Each record carries a `license_tier` field for downstream filtering.

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

Web application code is MIT. Compound data follows per-record license tiers carried forward from source databases; see the License tiers section above and the per-compound `license_tier` field.

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
