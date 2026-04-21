from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("/oscar/data/ehuertas/avakhari/ld_chr22")
INFILE = BASE / "gcta_ldscore_merged" / "chr22.full.PEL.ldscore.merged.tsv"
OUTFILE = BASE / "plots" / "PEL_ldscore_cdf_chr22_maf_gt_0.png"
OUTFILE.parent.mkdir(parents=True, exist_ok=True)

def ecdf(arr):
    x = np.sort(arr)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y

def archaic_mask(series):
    return series.astype(str).str.lower().isin(["true", "1", "yes"])

df = pd.read_csv(INFILE, sep="\t")
df = df[df["ldscore"].notna()].copy()
df = df[df["MAF"] > 0].copy()
mask = archaic_mask(df["is_archaic"])

plt.figure(figsize=(8, 6))

for is_archaic, color, label in [
    (True, "#d62728", "Archaic SNPs"),
    (False, "#1f77b4", "Non-archaic SNPs"),
]:
    sub = df[mask == is_archaic].copy()

    if sub.empty:
        print(f"No rows found for {label}")
        continue

    x, y = ecdf(sub["ldscore"].to_numpy())

    plt.plot(
        x,
        y,
        linewidth=2,
        color=color,
        label=f"{label} (n={len(sub):,})"
    )

plt.xlabel("LD score")
plt.ylabel("Cumulative distribution")
plt.title("chr22 LD-score CDF: PEL (MAF > 0)")
plt.ylim(0, 1)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(OUTFILE, dpi=300)