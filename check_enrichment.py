"""Quick stats after enrichment."""
import pandas as pd

df = pd.read_csv("data/theobroma_final.csv", low_memory=False)
print(f"Total: {len(df):,}")

# Region distribution
regions = df["region"].fillna("unresolved").replace({"": "unresolved", "global": "unresolved"})
print("\nBy region:")
for r, cnt in regions.value_counts().head(15).items():
    print(f"  {r}: {cnt:,}")

# Organism coverage
has_org = (df["source_organism"].fillna("").astype(str).str.strip() != "").sum()
print(f"\nWith source_organism: {has_org:,} ({100*has_org/len(df):.1f}%)")