#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

BASE = Path("/oscar/data/ehuertas/avakhari/ld_chr22")
PLOTDIR = BASE / "plots"
PLOTDIR.mkdir(parents=True, exist_ok=True)

curve_files = {
    "20kb": BASE / "curves" / "ld_decay_curves.tsv",
}

colors = {
    "CHB": "#1f77b4",
    "CEU": "#ff7f0e",
    "MXL": "#2ca02c",
    "PEL": "#d62728",
    "CLM": "#9467bd",
    "PUR": "#8c564b",
}

labels = {
    "CHB": "CHB (non-admixed)",
    "CEU": "CEU (non-admixed)",
    "MXL": "MXL (admixed)",
    "PEL": "PEL (admixed)",
    "CLM": "CLM (admixed)",
    "PUR": "PUR (admixed)",
}

pop_order = ["CHB", "CEU", "MXL", "PEL", "CLM", "PUR"]

def exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

for bin_label, curve_file in curve_files.items():
    df = pd.read_csv(curve_file, sep="\t")

    for snp_set, title in [
        ("archaic", f"LD decay on chr22: archaic SNPs ({bin_label} bins)"),
        ("non_archaic", f"LD decay on chr22: non-archaic SNPs ({bin_label} bins)"),
    ]:
        sub = df[df["snp_set"] == snp_set].copy()

        plt.figure(figsize=(10, 6))

        for pop in pop_order:
            d = sub[sub["population"] == pop].sort_values("bin_mid").copy()
            d = d[d["pair_count"] >= 20]

            if d.empty:
                continue

            x = d["bin_mid"].to_numpy() / 1000.0   # kb
            y = d["mean_r2"].to_numpy()

            plt.scatter(
                x, y,
                color=colors[pop],
                s=28,
                alpha=0.75,
                label=labels[pop]
            )


            a0 = max(y[0] - y[-1], 0.01)
            b0 = 0.01
            c0 = max(min(y[-1], y.min()), 0)

            popt, _ = curve_fit(
                exp_decay,
                x,
                y,
                p0=[a0, b0, c0],
                bounds=(0, [2.0, 5.0, 1.0]),
                maxfev=20000
            )

            xfit = np.linspace(x.min(), x.max(), 400)
            yfit = exp_decay(xfit, *popt)

            plt.plot(
                xfit, yfit,
                color=colors[pop],
                linewidth=2
            )

        plt.xlabel("Distance between SNP pairs (kb)")
        plt.ylabel("Mean $r^2$")
        plt.title(title)
        plt.ylim(bottom=0)
        plt.xlim(left=0)
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(PLOTDIR / f"{snp_set}_ld_decay_points_fit_chr22_{bin_label}.png", dpi=300)
        plt.close()