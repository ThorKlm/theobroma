"""Final dataset summary."""
import pandas as pd
df = pd.read_csv("data/theobroma_final.csv", low_memory=False)
print(f"Total: {len(df):,}")
print(f"\nBy source (top 25):")
for s, cnt in df["source_db"].value_counts().head(25).items():
    print(f"  {s}: {cnt:,}")
print(f"\nMulti-source: {(df['all_sources'].fillna('').str.contains('|', regex=False, na=False)).sum():,}")
print(f"\nWith organism: {(df['source_organism'].fillna('').astype(str).str.strip() != '').sum():,}")
npo = df["source_organism"].fillna("").astype(str).str.contains("NPO", na=False).sum()
print(f"  Real names: {(df['source_organism'].fillna('').astype(str).str.strip() != '').sum() - npo:,}")
print(f"  NPO codes: {npo:,}")