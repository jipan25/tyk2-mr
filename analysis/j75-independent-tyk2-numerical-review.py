#!/usr/bin/env python3
"""Recalculate MR and signal posteriors using only Python's standard library."""
import csv
import gzip
import json
import math
import sys
from pathlib import Path


def read(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t" if ".tsv" in path.name else ","))


def logsumexp(values):
    top = max(values)
    return top + math.log(sum(math.exp(value - top) for value in values))


def key(value):
    return value.removeprefix("chr").replace("_", ":")


def main():
    project = Path(__file__).resolve().parent
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else project / "tyk2_focused_rheumatology_immunology/source_results"
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else project / "tyk2_focused_rheumatology_immunology/submission_v2/scientific_review/numerical_audit"
    out.mkdir(parents=True, exist_ok=True)
    clinical = source / "clinical_ra_validation/complete_summary_stats"
    h = read(clinical / "TYK2_R13_Official_SE_Harmonised.csv")
    published = read(clinical / "TYK2_R13_Official_SE_MR_Methods.csv")
    mr_checks = []
    for endpoint in sorted({r["id.outcome"] for r in h}):
        data = [r for r in h if r["id.outcome"] == endpoint and r["mr_keep"] == "TRUE"]
        assert len(data) == 5
        assert all(r["effect_allele.exposure"] == r["effect_allele.outcome"] and
                   r["other_allele.exposure"] == r["other_allele.outcome"] for r in data)
        bx = [float(r["beta.exposure"]) for r in data]
        by = [float(r["beta.outcome"]) for r in data]
        sy = [float(r["se.outcome"]) for r in data]
        weight_sum = sum(x*x/s**2 for x,s in zip(bx,sy))
        beta = sum(x*y/s**2 for x,y,s in zip(bx,by,sy)) / weight_sum
        q = sum((y-beta*x)**2/s**2 for x,y,s in zip(bx,by,sy))
        se = math.sqrt(max(1, q / 4) / weight_sum)
        p = math.erfc(abs(beta / se) / math.sqrt(2))
        saved = next(r for r in published if r["id.outcome"] == endpoint and r["method"] == "Inverse variance weighted")
        for field, value in (("b", beta), ("se", se), ("pval", p)):
            assert math.isclose(value, float(saved[field]), rel_tol=1e-10, abs_tol=1e-12)
        mr_checks.append({"outcome": endpoint, "beta": beta, "se": se, "p": p,
                          "Q": q, "Q_p_df4": math.exp(-q/2)*(1+q/2),
                          "alleles_aligned": True, "status": "PASS"})

    molecular = source / "molecular_qtl_followup"
    results = read(molecular / "results_signal_coloc/TYK2_Published_SuSiE_Signal_Colocalization_Prior_Sensitivity.csv")
    gwas, qtl = {}, {}
    posterior_checks = []
    for r in results:
        endpoint, dataset = r["outcome_id"], r["dataset_id"]
        if endpoint not in gwas:
            path = molecular / "derived_finemap" / f"finngen_R13_{endpoint}.SUSIE.TYK2_region.tsv.gz"
            gwas[endpoint] = {key(x["v"]): x for x in read(path) if x["region"] == "chr19:8852442-11852442"}
        if dataset not in qtl:
            qtl[dataset] = read(molecular / "derived_finemap" / f"{dataset}.TYK2_KEAP1.lbf_variable.tsv.gz")
        q = {key(x["variant"]): x for x in qtl[dataset] if x["molecular_trait_id"] == r["molecular_trait_id"]}
        g = gwas[endpoint]
        common = sorted(g.keys() & q.keys())
        a = [float(g[v]["lbf_variable" + r["gwas_signal"]]) for v in common]
        b = [float(q[v]["lbf_variable" + r["qtl_signal"]]) for v in common]
        s1, s2, s12 = logsumexp(a), logsumexp(b), logsumexp([x+y for x,y in zip(a,b)])
        separate = s1+s2+math.log1p(-math.exp(s12-s1-s2))
        p1, p2, p12 = (float(r[x]) for x in ("p1", "p2", "p12"))
        logs = [0, math.log(p1)+s1, math.log(p2)+s2,
                math.log(p1*p2)+separate, math.log(p12)+s12]
        total = logsumexp(logs)
        pp = [math.exp(x-total) for x in logs]
        error = max(abs(pp[i]-float(r[f"PP.H{i}.abf"])) for i in range(5))
        assert len(common) == int(r["n_shared_variants"]) and error < 1e-10
        posterior_checks.append({"outcome": endpoint, "dataset": dataset, "molecular_trait": r["molecular_trait_id"],
                                 "gwas_signal": r["gwas_signal"], "qtl_signal": r["qtl_signal"],
                                 "p12": p12, "independent_H4": pp[4], "max_error": error, "status": "PASS"})
    assert len(posterior_checks) == 72
    for name, records in (("Independent_Python_MR.csv", mr_checks), ("Independent_Python_Coloc.csv", posterior_checks)):
        with (out / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    summary = {"MR_outcomes_verified": len(mr_checks), "signal_prior_comparisons_verified": len(posterior_checks),
               "max_posterior_error": max(r["max_error"] for r in posterior_checks),
               "status": "PASS", "limits": "Independent implementation, not an independent human review or new cohort; LD assumptions are not validated by arithmetic agreement"}
    (out / "Independent_Numerical_Audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
