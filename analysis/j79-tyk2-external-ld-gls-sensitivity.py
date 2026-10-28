#!/usr/bin/env python3
"""Approximate covariance sensitivity; EUR LD is not FinnGen in-sample LD."""
import csv
import math
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parent / "tyk2_focused_rheumatology_immunology/source_results/clinical_ra_validation"
out = root / "instrument_ld_audit"
with (out / "TYK2_1000G_EUR_Signed_LD.csv").open() as stream:
    matrix_rows = list(csv.reader(stream))
labels = matrix_rows[0][1:]
ld = np.array([[float(value) for value in row[1:]] for row in matrix_rows[1:]])
assert np.max(np.abs(ld - ld.T)) < 1e-10 and np.min(np.linalg.eigvalsh(ld)) > 0
with (root / "complete_summary_stats/TYK2_R13_Official_SE_Harmonised.csv").open() as stream:
    harmonized = list(csv.DictReader(stream))
results = []
for endpoint in sorted({row["id.outcome"] for row in harmonized}):
    records = {row["SNP"]: row for row in harmonized if row["id.outcome"] == endpoint and row["mr_keep"] == "TRUE"}
    data, signs = [], []
    for label in labels:
        snp, a1, a2 = label.split("_")
        row = records[snp]
        pair = (row["effect_allele.exposure"], row["other_allele.exposure"])
        assert pair in ((a1, a2), (a2, a1))
        signs.append(1 if pair == (a1, a2) else -1)
        data.append(row)
    aligned_ld = ld * np.outer(signs, signs)
    bx = np.array([float(row["beta.exposure"]) for row in data])
    by = np.array([float(row["beta.outcome"]) for row in data])
    sy = np.array([float(row["se.outcome"]) for row in data])
    covariance = aligned_ld * np.outer(sy, sy)
    precision_bx = np.linalg.solve(covariance, bx)
    denominator = bx @ precision_bx
    beta = float(precision_bx @ by / denominator)
    residual = by - beta * bx
    q = float(residual @ np.linalg.solve(covariance, residual))
    se = float(math.sqrt(max(1, q/(len(bx)-1))/denominator))
    results.append({"outcome_id": endpoint, "method": "EUR-reference generalized IVW sensitivity",
                    "nsnp": len(bx), "b": beta, "se": se, "pval": math.erfc(abs(beta/se)/math.sqrt(2)),
                    "Q": q, "Q_df": len(bx)-1, "min_LD_eigenvalue": float(np.min(np.linalg.eigvalsh(ld))),
                    "limitation": "External 503-person EUR LD, not Finnish in-sample LD; approximate covariance sensitivity"})
with (out / "TYK2_EUR_LD_GLS_Sensitivity.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(results[0]))
    writer.writeheader()
    writer.writerows(results)
for row in results:
    print(row)
