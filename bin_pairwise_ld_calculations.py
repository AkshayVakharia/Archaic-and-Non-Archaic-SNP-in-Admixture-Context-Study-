#!/usr/bin/env python3
from pathlib import Path
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input PLINK .vcor file")
    parser.add_argument("--population", required=True, help="Population label, e.g. CHB")
    parser.add_argument("--snp-set", required=True, help="archaic or non_archaic")
    parser.add_argument("--output", required=True, help="Output summarized TSV")
    parser.add_argument("--chunk-size", type=int, default=1_000_000, help="Rows per chunk")
    parser.add_argument("--bin-size", type=int, default=20_000, help="Distance bin size in bp")
    parser.add_argument("--max-dist", type=int, default=1_000_000, help="Maximum pair distance in bp")
    args = parser.parse_args()

    infile = Path(args.input)
    outfile = Path(args.output)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    bins = np.arange(0, args.max_dist + args.bin_size, args.bin_size)
    n_bins = len(bins) - 1
    sum_r2 = np.zeros(n_bins, dtype=np.float64)
    count_r2 = np.zeros(n_bins, dtype=np.int64)
    total_rows = 0
    kept_rows = 0
    chunk_idx = 0

    reader = pd.read_csv(
        infile,
        sep=r"\s+",
        engine="c",
        chunksize=args.chunk_size,
        low_memory=True
    )

    for chunk in reader:
        chunk_idx += 1
        total_rows += len(chunk)

        r2_col = None
        for c in chunk.columns:
            if c.upper() in {"UNPHASED_R2", "R2"}:
                r2_col = c
                break
       
        dist = (chunk["POS_B"] - chunk["POS_A"]).abs().to_numpy()
        r2 = pd.to_numeric(chunk[r2_col], errors="coerce").to_numpy()

        valid = np.isfinite(r2) & (dist <= args.max_dist)
        dist = dist[valid]
        r2 = r2[valid]
        kept_rows += len(r2)

        bin_idx = np.searchsorted(bins, dist, side="right") - 1
        valid_bins = (bin_idx >= 0) & (bin_idx < n_bins)
        bin_idx = bin_idx[valid_bins]
        r2 = r2[valid_bins]

        if len(r2) > 0:
            sum_add = np.bincount(bin_idx, weights=r2, minlength=n_bins)
            count_add = np.bincount(bin_idx, minlength=n_bins)
            sum_r2[:len(sum_add)] += sum_add
            count_r2[:len(count_add)] += count_add

    bin_start = bins[:-1]
    bin_end = bins[1:]
    bin_mid = (bin_start + bin_end) / 2
    mean_r2 = np.divide(
        sum_r2,
        count_r2,
        out=np.full_like(sum_r2, np.nan, dtype=np.float64),
        where=count_r2 > 0
    )

    out_df = pd.DataFrame({
        "snp_set": args.snp_set,
        "population": args.population,
        "bin_start": bin_start,
        "bin_end": bin_end,
        "bin_mid": bin_mid,
        "mean_r2": mean_r2,
        "pair_count": count_r2
    })

    out_df.to_csv(outfile, sep="\t", index=False)
