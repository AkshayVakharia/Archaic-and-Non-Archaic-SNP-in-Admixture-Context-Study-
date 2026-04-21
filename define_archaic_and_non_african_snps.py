#!/usr/bin/env python3

import pandas as pd
import numpy as np
from cyvcf2 import VCF

modern_vcf_file = "/users/avakhari/data/avakhari/1000g_qc/chr22.QC.vcf.gz"
metadata_file = "/users/avakhari/data/data/modern_genomes/1000genomes_GRCh38/20130606_g1k_3202_samples_ped_population.txt"
archaic_vcf_file = "/users/avakhari/data/data/archaic_genomes/SkovFiles_GRCh38/individuals_highcov.QC.chr22.vcf.gz"
output_file = "chr22_non_african_archaic_snps.tsv"


def get_ref_alt_freq(genotypes, sample_indices):
    ref_count = 0
    alt_count = 0
    total_alleles = 0

    for idx in sample_indices:
        gt = genotypes[idx]
        a1 = gt[0]
        a2 = gt[1]

        for allele in [a1, a2]:
            if allele == -1:
                continue
            if allele == 0:
                ref_count += 1
                total_alleles += 1
            elif allele == 1:
                alt_count += 1
                total_alleles += 1
            else:
                continue

    if total_alleles == 0:
        return np.nan, np.nan

    ref_freq = ref_count / total_alleles
    alt_freq = alt_count / total_alleles
    return ref_freq, alt_freq

meta = pd.read_csv(metadata_file, sep=r"\s+")

afr_meta = meta[
    (meta["Superpopulation"] == "AFR") &
    (meta["Population"] != "ACB") &
    (meta["Population"] != "ASW")
]
afr_pops = sorted(afr_meta["Population"].unique())

ooa_meta = meta[
    (~meta["Population"].isin(afr_pops)) &
    (meta["Population"] != "ACB") &
    (meta["Population"] != "ASW")
]
ooa_pops = sorted(ooa_meta["Population"].unique())

modern_vcf = VCF(modern_vcf_file)
modern_samples = modern_vcf.samples

sample_to_pop = {}
for i in range(len(meta)):
    sample_id = meta.iloc[i]["SampleID"]
    pop = meta.iloc[i]["Population"]
    sample_to_pop[sample_id] = pop

pop_to_indices = {}
for pop in list(afr_pops) + list(ooa_pops):
    pop_to_indices[pop] = []

for i, sample in enumerate(modern_samples):
    if sample in sample_to_pop:
        pop = sample_to_pop[sample]
        if pop in pop_to_indices:
            pop_to_indices[pop].append(i)

archaic_vcf = VCF(archaic_vcf_file)
archaic_samples = archaic_vcf.samples
archaic_dict = {}
for variant in archaic_vcf:
    if len(variant.ALT) != 1:
        continue
    if len(variant.REF) != 1 or len(variant.ALT[0]) != 1:
        continue
    key = (variant.CHROM, variant.POS, variant.REF, variant.ALT[0])
    info = {}
    info["is_archaic"] = False

    for i, sample in enumerate(archaic_samples):
        gt = variant.genotypes[i]
        a1 = gt[0]
        a2 = gt[1]

        if a1 == 0 and a2 == 0:
            ref_freq = 1.0
            alt_freq = 0.0
        elif a1 == 1 and a2 == 1:
            ref_freq = 0.0
            alt_freq = 1.0
            info["is_archaic"] = True
        elif (a1 == 0 and a2 == 1) or (a1 == 1 and a2 == 0):
            ref_freq = 0.5
            alt_freq = 0.5
        else:
            ref_freq = np.nan
            alt_freq = np.nan

        info[sample + "_REF_FREQ"] = ref_freq
        info[sample + "_ALT_FREQ"] = alt_freq

    archaic_dict[key] = info

results = []
total_snps = 0
num_non_african = 0
num_archaic = 0
for variant in modern_vcf:
    if len(variant.ALT) != 1:
        continue
    if len(variant.REF) != 1 or len(variant.ALT[0]) != 1:
        continue

    total_snps += 1

    chrom = variant.CHROM
    pos = variant.POS
    ref = variant.REF
    alt = variant.ALT[0]
    genotypes = variant.genotypes

    afr_freqs = {}
    for pop in afr_pops:
        pop_ref, pop_alt = get_ref_alt_freq(genotypes, pop_to_indices[pop])
        afr_freqs[pop] = (pop_ref, pop_alt)

    ooa_freqs = {}
    for pop in ooa_pops:
        pop_ref, pop_alt = get_ref_alt_freq(genotypes, pop_to_indices[pop])
        ooa_freqs[pop] = (pop_ref, pop_alt)

    african_ok = True
    for pop in afr_pops:
        pop_alt = afr_freqs[pop][1]
        if not np.isnan(pop_alt) and pop_alt >= 0.01:
            african_ok = False
            break

    non_african_pop_ok = False
    for pop in ooa_pops:
        pop_alt = ooa_freqs[pop][1]
        if not np.isnan(pop_alt) and pop_alt > 0.01:
            non_african_pop_ok = True
            break

    is_non_african = african_ok and non_african_pop_ok

    if not is_non_african:
        continue

    num_non_african += 1

    key = (chrom, pos, ref, alt)

    if key in archaic_dict:
        archaic_info = archaic_dict[key]
        is_archaic = archaic_info["is_archaic"]
    else:
        archaic_info = {
            "AltaiNeandertal_REF_FREQ": np.nan,
            "AltaiNeandertal_ALT_FREQ": np.nan,
            "Vindija33.19_REF_FREQ": np.nan,
            "Vindija33.19_ALT_FREQ": np.nan,
            "Denisova_REF_FREQ": np.nan,
            "Denisova_ALT_FREQ": np.nan,
            "Chagyrskaya-Phalanx_REF_FREQ": np.nan,
            "Chagyrskaya-Phalanx_ALT_FREQ": np.nan,
            "is_archaic": False
        }
        is_archaic = False

    if is_archaic:
        num_archaic += 1

    row = {
        "CHROM": chrom,
        "POS": pos,
        "REF": ref,
        "ALT": alt,
        "AltaiNeandertal_REF_FREQ": archaic_info["AltaiNeandertal_REF_FREQ"],
        "AltaiNeandertal_ALT_FREQ": archaic_info["AltaiNeandertal_ALT_FREQ"],
        "Vindija33.19_REF_FREQ": archaic_info["Vindija33.19_REF_FREQ"],
        "Vindija33.19_ALT_FREQ": archaic_info["Vindija33.19_ALT_FREQ"],
        "Denisova_REF_FREQ": archaic_info["Denisova_REF_FREQ"],
        "Denisova_ALT_FREQ": archaic_info["Denisova_ALT_FREQ"],
        "Chagyrskaya-Phalanx_REF_FREQ": archaic_info["Chagyrskaya-Phalanx_REF_FREQ"],
        "Chagyrskaya-Phalanx_ALT_FREQ": archaic_info["Chagyrskaya-Phalanx_ALT_FREQ"],
        "is_archaic": is_archaic,
        "is_non_african": is_non_african
    }

    for pop in afr_pops:
        row[f"{pop}_REF_FREQ"] = afr_freqs[pop][0]
        row[f"{pop}_ALT_FREQ"] = afr_freqs[pop][1]

    for pop in ooa_pops:
        row[f"{pop}_REF_FREQ"] = ooa_freqs[pop][0]
        row[f"{pop}_ALT_FREQ"] = ooa_freqs[pop][1]

    results.append(row)

results_df = pd.DataFrame(results)
results_df.to_csv(output_file, sep="\t", index=False)

print("Total SNPs analyzed:", total_snps)
print("Number of non-African SNPs:", num_non_african)
print("Number of archaic SNPs:", num_archaic)