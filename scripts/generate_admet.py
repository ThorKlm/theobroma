"""Generate ADMET predictions for a list of SMILES using ADMET-AI.

ADMET-AI (Swanson et al.) provides 104 endpoint predictions via a Chemprop
D-MPNN ensemble. THEOBROMA persists only four of these endpoints to the
admet table for filter UX (Solubility_AqSolDB, Lipophilicity_AstraZeneca,
DILI, Clearance_Hepatocyte_AZ); the full 104-column CSV is archived under
data/vectors/admet_predictions.csv.

Dependencies (not installed by default in the production venv):
    pip install admet_ai

Run example (CPU; ~5-15 min on L3S for 472 SMILES):
    venv/bin/python scripts/generate_admet.py --input data/tipdb_novel_smiles.csv \\
        --output /tmp/admet_tipdb.csv

Inputs:
    --input  CSV with columns: comp_id, smiles
    --output destination CSV with comp_id plus 104 ADMET-endpoint columns

After running, four columns are inserted into the admet table via:
    INSERT INTO admet (comp_id, "Solubility_AqSolDB", "Lipophilicity_AstraZeneca",
                       "DILI", "Clearance_Hepatocyte_AZ")
    SELECT comp_id, "Solubility_AqSolDB", "Lipophilicity_AstraZeneca",
           "DILI", "Clearance_Hepatocyte_AZ"
    FROM admet_tipdb_staging ON CONFLICT (comp_id) DO UPDATE SET ...
"""
import argparse, csv, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    try:
        from admet_ai import ADMETModel
    except ImportError:
        print("admet_ai not installed; run: pip install admet_ai", file=sys.stderr)
        sys.exit(1)
    model = ADMETModel()
    rows_in = []
    with open(args.input) as f:
        r = csv.DictReader(f)
        for row in r:
            rows_in.append((row["comp_id"], row["smiles"]))
    smiles_list = [s for _, s in rows_in]
    preds = model.predict(smiles=smiles_list)
    # preds is a pandas DataFrame indexed by smiles with 104 columns
    cols = list(preds.columns)
    with open(args.output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["comp_id"] + cols)
        for (comp_id, smi), pred_row in zip(rows_in, preds.itertuples(index=False)):
            w.writerow([comp_id] + list(pred_row))
    print(f"wrote {len(rows_in)} compounds x {len(cols)} endpoints to {args.output}")

if __name__ == "__main__":
    main()
