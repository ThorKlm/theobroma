# THEOBROMA

An open global multi-kingdom natural products database for virtual screening and drug discovery.

**Live:** [https://theobroma.l3s.uni-hannover.de](https://theobroma.l3s.uni-hannover.de)

## Overview

THEOBROMA aggregates ~960K natural product compounds from 21+ source databases spanning five biological kingdoms (plant, fungi, bacteria, marine, food) and six continents. It provides free, open access to compound structures, physicochemical properties, taxonomic classification, and geographic provenance metadata via a web interface and JSON API.

## Features

- **Search** by compound name, SMILES, InChIKey, source organism, geographic region, kingdom, or source database
- **Browse** with filters by kingdom and source
- **Compound detail** pages with structure, properties, classification, provenance, and external links
- **Statistics** dashboard with kingdom, source, and geographic distribution
- **Bulk download** as CSV and SDF
- **JSON API** for programmatic access

## Source databases

COCONUT 2.0, CMAUP, TCMBank, HERB 2.0, TM-MC 2.0, NPASS 3.0, FooDB, CMNPD, IMPPAT 2.0, Phyto4Health, StreptomeDB, NPAtlas, MeFSAT, MicotoXilico, MIBiG 4.0, phytochemdb, CSIRO Australian NP, VIETHERB, ANPDB, SANCDB, AfroDb, ConMedNP, CamMedNP, LANaPDB, CyanoMetDB, SWMD, YMDB, EMNPD.

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

# Load data
python scripts/download_hf.py --token YOUR_HF_TOKEN
python scripts/load_data.py

# Run
gunicorn -b 0.0.0.0:5000 app:app
```

## API

```bash
curl "https://theobroma.l3s.uni-hannover.de/api/search?q=curcumin&type=name&limit=10"
```

## License

- **Web application code:** MIT License
- **Compound data:** Dual-tier (CC BY 4.0 for open sources, CC BY-NC 4.0 for restricted sources). Each record includes a `license_tier` field.

## Citation

```
Klamt, T. et al. (2026). THEOBROMA: an open global multi-kingdom natural products
database for virtual screening and drug discovery. Nucleic Acids Research (submitted).
```

## Hosted by

[L3S Research Center](https://www.l3s.de), Leibniz University Hannover
