#!/usr/bin/env python3
"""Independently recalculate the unconditioned gene-level ABF results."""
import csv
import gzip
import json
import math
from pathlib import Path

source = Path(__file__).resolve().parent / "tyk2_focused_rheumatology_immunology/source_results"
root = source / "independent_expression_validation"
results = list(csv.DictReader((root / "results/Independent_Expression_Regional_Coloc.csv").open()))


def lse(values):
    top = max(values)
    return top + math.log(sum(math.exp(x-top) for x in values))


def logbf(beta, se, prior):
    variance, w = se*se, prior*prior
    r = w/(variance+w)
    return .5*(math.log1p(-r)+r*(beta/se)**2)


gwas, datasets = {}, {}
qa = []
for row in results:
    endpoint, dataset = row["outcome_id"], row["dataset_id"]
    if endpoint not in gwas:
        path = source / "clinical_ra_validation/raw_complete_summary_stats" / f"finngen_R13_{endpoint}_TYK2_complete_region.tsv.gz"
        with gzip.open(path, "rt") as stream:
            gwas[endpoint] = {":".join((r["#chrom"],r["pos"],r["ref"],r["alt"])):r for r in csv.DictReader(stream,delimiter="\t")}
    if dataset not in datasets:
        path = root / f"{dataset}_TYK2_KEAP1_full_cis.tsv.gz"
        with gzip.open(path, "rt") as stream:
            records = [r for r in csv.DictReader(stream,delimiter="\t") if r["gene_id"].split(".")[0]=="ENSG00000105397"]
        valid = {}
        for r in records:
            quality = (float(r["maf"])>=.01 and float(r["ac"])>=10
                       and (r["r2"]=="NA" or float(r["r2"])>=.8))
            if not quality:
                continue
            key = ":".join((r["chromosome"],r["position"],r["ref"],r["alt"]))
            if key in valid:
                assert (r["beta"],r["se"],r["pvalue"]) == (valid[key]["beta"],valid[key]["se"],valid[key]["pvalue"])
            valid[key] = r
        datasets[dataset] = valid
    q, g = datasets[dataset], gwas[endpoint]
    common = sorted(q.keys() & g.keys())
    a = [logbf(float(q[v]["beta"]),float(q[v]["se"]),.15) for v in common]
    b = [logbf(float(g[v]["beta"]),float(g[v]["sebeta"]),.2) for v in common]
    s1,s2,s12 = lse(a),lse(b),lse([x+y for x,y in zip(a,b)])
    p1,p2,p12 = (float(row[k]) for k in ("p1","p2","p12"))
    weights = [0,math.log(p1)+s1,math.log(p2)+s2,
               math.log(p1*p2)+s1+s2+math.log1p(-math.exp(s12-s1-s2)),math.log(p12)+s12]
    norm = lse(weights)
    pp = [math.exp(value-norm) for value in weights]
    error = max(abs(pp[i]-float(row[f"PP.H{i}.abf"])) for i in range(5))
    assert len(common)==int(row["nsnps"]) and error<1e-10
    qa.append({"dataset":dataset,"outcome":endpoint,"p12":p12,"n_common":len(common),"max_error":error,"status":"PASS"})
with (root / "results/Independent_Expression_Python_Coloc_QA.csv").open("w",newline="") as stream:
    writer=csv.DictWriter(stream,fieldnames=list(qa[0]))
    writer.writeheader()
    writer.writerows(qa)
assert len(qa)==36
print(json.dumps({"comparisons":len(qa),"max_error":max(row["max_error"] for row in qa),"status":"PASS"}))
