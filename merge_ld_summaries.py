from pathlib import Path
import pandas as pd

base = Path("/oscar/data/ehuertas/avakhari/ld_chr22")
indir = base / "curves" / "by_file_50kb"
outfile = base / "curves" / "ld_decay_curves_50kb.tsv"

files = sorted(indir.glob("*.summary.tsv"))
dfs = []
for f in files:
    df = pd.read_csv(f, sep="\t")
    dfs.append(df)

out = pd.concat(dfs, ignore_index=True)
outfile.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(outfile, sep="\t", index=False)