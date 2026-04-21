from pathlib import Path
import pandas as pd

BASE = Path("/oscar/data/ehuertas/avakhari/ld_chr22")
TSV = Path("/oscar/data/ehuertas/avakhari/archaic_snps/chr22_non_african_archaic_snps.tsv")
LDDIR = BASE / "gcta_ldscore"
OUTDIR = BASE / "gcta_ldscore_merged"
OUTDIR.mkdir(parents=True, exist_ok=True)

f = LDDIR / "chr22.full.PEL.score.ld"
df = pd.read_csv(f, sep=r"\s+", engine="python")
lab = pd.read_csv(TSV, sep="\t")
lab["chrom_clean"] = lab["CHROM"].astype(str).str.replace("^chr", "", regex=True)
lab["snp_id"] = (
    lab["chrom_clean"] + ":" +
    lab["POS"].astype(str) + ":" +
    lab["REF"].astype(str) + ":" +
    lab["ALT"].astype(str)
)
print("First TSV SNP IDs:", lab["snp_id"].head().tolist())
lab = lab[["snp_id", "CHROM", "POS", "REF", "ALT", "is_archaic", "is_non_african"]].copy()
merged = df.merge(lab, left_on="SNP", right_on="snp_id", how="inner")
outfile = OUTDIR / "chr22.full.PEL.ldscore.merged.tsv"
merged.to_csv(outfile, sep="\t", index=False)
