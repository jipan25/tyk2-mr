#!/usr/bin/env python3
"""Secondary fixed-marker replication; does not relax discovery selection."""
import csv
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

root = Path(__file__).resolve().parent / "tyk2_focused_rheumatology_immunology/source_results/independent_expression_validation"
ids = ["QTD000021", "QTD000031", "QTD000609", "QTD000612"]
variant = "chr19_10354011_A_G"
plan = {"date_utc": datetime.now(timezone.utc).isoformat(), "role": "Secondary exploratory fixed-marker check added after the primary availability analysis",
        "marker": "rs280497", "variant": variant, "coded_allele": "G",
        "selection_basis": "Official GTEx whole-blood lead selected by expression significance, not by RA outcome",
        "contexts": ids, "replication_rule": "same ALT effect direction, source P<0.05/4, MAF>=0.01, AC>=10, reported imputation r2>=0.8",
        "limits": "Does not change P<5e-8 discovery rule, create a new Catalogue MR estimate, establish cohort independence or RNA mediation"}
(root / "Fixed_Lead_Replication_Plan.json").write_text(json.dumps(plan, indent=2) + "\n")
records = []
for dataset in ids:
    with gzip.open(root / f"{dataset}_TYK2_KEAP1_full_cis.tsv.gz", "rt") as stream:
        matched = [r for r in csv.DictReader(stream, delimiter="\t") if r["variant"] == variant
                   and r["gene_id"].split(".")[0] == "ENSG00000105397"]
    assert matched and len({(r["beta"], r["se"], r["pvalue"]) for r in matched}) == 1
    r = matched[0]
    quality = float(r["maf"]) >= .01 and float(r["ac"]) >= 10 and (r["r2"] == "NA" or float(r["r2"]) >= .8)
    concordant = float(r["beta"]) > 0
    adjusted = min(1, 4*float(r["pvalue"]))
    records.append({"dataset_id": dataset, "SNP": "rs280497", "effect_allele": "G",
                    "beta_expression": r["beta"], "se_expression": r["se"], "pvalue_source": r["pvalue"],
                    "F": (float(r["beta"])/float(r["se"]))**2, "MAF": r["maf"], "AC": r["ac"],
                    "imputation_r2": r["r2"], "quality_pass": quality, "direction_concordant": concordant,
                    "Bonferroni_four_fixed_contexts": adjusted,
                    "strong_instrument_F_gt10": (float(r["beta"])/float(r["se"]))**2 > 10,
                    "fixed_marker_replication_supported": quality and concordant and adjusted < .05,
                    "role": plan["role"], "limit": plan["limits"]})
with (root / "results/Fixed_GTEx_Lead_Expression_Replication.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)
for row in records:
    print(row)
