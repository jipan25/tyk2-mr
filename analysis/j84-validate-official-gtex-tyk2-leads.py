#!/usr/bin/env python3
"""Direct GTEx-lead follow-up; shared participants and variants remain explicit."""
import csv
import gzip
import math
from pathlib import Path

source = Path(__file__).resolve().parent / "tyk2_focused_rheumatology_immunology/source_results"
root = source / "independent_expression_validation"
official = root / "official_gtex_pipeline_audit"
with (official / "Official_GTEx_TYK2_KEAP1_Lead_Audit.csv").open() as stream:
    leads = [row for row in csv.DictReader(stream) if row["gene_name"] == "TYK2"]
gwas = {}
for endpoint in ("M13_RHEUMA", "RHEUMA_SEROPOS_STRICT"):
    path = source / "clinical_ra_validation/raw_complete_summary_stats" / f"finngen_R13_{endpoint}_TYK2_complete_region.tsv.gz"
    with gzip.open(path, "rt") as stream:
        gwas[endpoint] = {":".join((r["#chrom"],r["pos"],r["ref"],r["alt"])):r
                          for r in csv.DictReader(stream, delimiter="\t")}
comparisons = []
validation = []
for row in leads:
    chrom, pos, ref, alt, build = row["variant_id"].split("_")
    variant = ":".join((chrom.removeprefix("chr"), pos, ref, alt))
    assert build == "b38"
    with (official / f"{row['tissue']}_TYK2_KEAP1_egenes.tsv").open() as stream:
        original = next(r for r in csv.DictReader(stream, delimiter="\t") if r["gene_name"] == "TYK2")
    qualified = (float(row["pval_nominal"]) < 5e-8 and float(row["qval"]) < .05 and float(row["F"]) >= 10
                 and float(original["maf"]) >= .01 and float(original["minor_allele_count"]) >= 10)
    primary = row["tissue"] == "Whole_Blood"
    bx, sx = float(row["slope"]), float(row["slope_se"])
    for endpoint, variants in gwas.items():
        match = variants.get(variant)
        if not match or not qualified:
            continue
        by, sy = float(match["beta"]), float(match["sebeta"])
        b, se = by/bx, sy/abs(bx)
        delta_se = math.sqrt(sy**2/bx**2 + by**2*sx**2/bx**4)
        p = math.erfc(abs(b/se)/math.sqrt(2))
        validation.append({"tissue": row["tissue"], "primary_whole_blood": primary, "variant": variant,
                           "SNP": match["rsids"], "effect_allele": alt, "p_exposure": row["pval_nominal"],
                           "gene_FDR": row["qval"], "F": row["F"], "outcome_id": endpoint,
                           "b": b, "se_NOME": se, "pval_NOME": p, "se_delta": delta_se,
                           "pval_delta": math.erfc(abs(b/delta_se)/math.sqrt(2)),
                           "Bonferroni_two_primary_endpoints": min(1,2*p) if primary else "",
                           "limit": "GTEx/eQTLGen overlap; single-variant follow-up without full official-pipeline regional colocalization; shared variant across tissues"})
    dataset = {"Whole_Blood":"QTD000356", "Spleen":"QTD000326"}.get(row["tissue"])
    if dataset:
        with gzip.open(root / f"{dataset}_TYK2_KEAP1_full_cis.tsv.gz", "rt") as stream:
            rows = [r for r in csv.DictReader(stream, delimiter="\t") if r["variant"] == "chr"+variant.replace(":","_")
                    and r["gene_id"].split(".")[0] == "ENSG00000105397"]
        assert rows and len({(r["beta"],r["se"],r["pvalue"]) for r in rows}) == 1
        q = rows[0]
        comparisons.append({"tissue": row["tissue"], "variant": variant, "official_p": row["pval_nominal"],
                            "official_beta": bx, "official_se": sx, "Catalogue_dataset": dataset,
                            "Catalogue_p": q["pvalue"], "Catalogue_beta": q["beta"], "Catalogue_se": q["se"],
                            "same_ALT_direction": math.copysign(1,bx)==math.copysign(1,float(q["beta"])),
                            "limit": "Same cohort, different quantification/normalization/covariates; neither p-value difference nor slope scale proves superior biology"})
for name, rows in (("Official_GTEx_TYK2_Direct_Lead_RA_Followup.csv",validation),
                   ("Official_GTEx_vs_Catalogue_Lead_Comparison.csv",comparisons)):
    with (official / name).open("w",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
for row in validation:
    if row["primary_whole_blood"]:
        print(row)
print(comparisons)
