#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

tsv = "/users/avakhari/data/avakhari/archaic_snps/chr22_non_african_archaic_snps.tsv"
meta = "/users/avakhari/data/data/modern_genomes/1000genomes_GRCh38/20130606_g1k_3202_samples_ped_population.txt"
outdir = Path("/users/avakhari/data/avakhari/ld_chr22")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(tsv, sep="\t")
df["is_archaic"] = df["is_archaic"].astype(str).str.lower().isin(
    ["true", "1", "yes", "y"]
)
archaic = df[df["is_archaic"]].copy()
non_archaic = df[~df["is_archaic"]].copy()

archaic[["CHROM", "POS"]].drop_duplicates().to_csv(
    outdir / "archaic_chr22.sites.tsv",
    sep="\t", header=False, index=False
)

non_archaic[["CHROM", "POS"]].drop_duplicates().to_csv(
    outdir / "non_archaic_chr22.sites.tsv",
    sep="\t", header=False, index=False
)

print(f"Archaic SNPs: {len(archaic):,}")
print(f"Non-archaic SNPs: {len(non_archaic):,}")

meta_df = pd.read_csv(meta, sep=r"\s+", engine="python")
sample_col = None
pop_col = None

for c in meta_df.columns:
    lc = c.lower()
    if lc == "sampleid":
        sample_col = c
    if lc == "population":
        pop_col = c

populations = ["CHB", "CEU", "MXL", "PEL", "CLM", "PUR"]
for pop in populations:
    keep = meta_df.loc[meta_df[pop_col] == pop, sample_col].drop_duplicates().sort_values()
    keep_df = pd.DataFrame({"FID": keep.values, "IID": keep.values})
    keep_df.to_csv(outdir / f"{pop}.keep", sep="\t", header=False, index=False)
