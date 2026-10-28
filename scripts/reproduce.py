#!/usr/bin/env python3
"""Recalculate archived core results; outputs never modify the supplied inputs."""
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'


def rows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))


def close(a,b):
    assert math.isclose(float(a),float(b),rel_tol=1e-9,abs_tol=1e-12),(a,b)


def ivw(data):
    assert all(r['effect_allele.exposure']==r['effect_allele.outcome'] and
               r['other_allele.exposure']==r['other_allele.outcome'] for r in data)
    bx=np.array([float(r['beta.exposure']) for r in data])
    by=np.array([float(r['beta.outcome']) for r in data])
    sy=np.array([float(r['se.outcome']) for r in data])
    denominator=float(np.sum(bx*bx/sy**2))
    beta=float(np.sum(bx*by/sy**2)/denominator)
    q=float(np.sum((by-beta*bx)**2/sy**2))
    se=math.sqrt(max(1,q/(len(data)-1))/denominator)
    return {'b':beta,'se':se,'pval':math.erfc(abs(beta/se)/math.sqrt(2)),'Q':q}


def main():
    for record in rows(DATA/'PROVENANCE.csv'):
        assert hashlib.sha256((ROOT/record['path']).read_bytes()).hexdigest()==record['sha256'],record['path']
    harmonized=rows(DATA/'core/harmonized_finngen.csv')
    expected=rows(DATA/'core/finngen_mr_expected.csv')
    loo=rows(DATA/'core/finngen_loo_expected.csv')
    mr=[]
    for endpoint in sorted({r['id.outcome'] for r in harmonized}):
        data=[r for r in harmonized if r['id.outcome']==endpoint and r['mr_keep'].lower()=='true']
        result=ivw(data)
        saved=next(r for r in expected if r['id.outcome']==endpoint and r['method']=='Inverse variance weighted')
        for field in ('b','se','pval'):close(result[field],saved[field])
        mr.append({'outcome':endpoint,'model':'diagonal IVW',**result})
        for omitted in data:
            check=ivw([r for r in data if r['SNP']!=omitted['SNP']])
            saved=next(r for r in loo if r['id.outcome']==endpoint and r['omitted_snp']==omitted['SNP'])
            for field,key in [('b','loo_beta'),('se','loo_se'),('pval','loo_p'),('Q','residual_q')]:close(check[field],saved[key])
        with (DATA/'core/eur_signed_ld.csv').open() as f:matrix=list(csv.reader(f))
        labels=matrix[0][1:];ld=np.array([[float(v) for v in r[1:]] for r in matrix[1:]])
        by_snp={r['SNP']:r for r in data};aligned=[];signs=[]
        for label in labels:
            snp,a1,a2=label.split('_');r=by_snp[snp]
            pair=(r['effect_allele.exposure'],r['other_allele.exposure'])
            assert pair in [(a1,a2),(a2,a1)]
            signs.append(1 if pair==(a1,a2) else -1);aligned.append(r)
        bx=np.array([float(r['beta.exposure']) for r in aligned])
        by=np.array([float(r['beta.outcome']) for r in aligned])
        sy=np.array([float(r['se.outcome']) for r in aligned])
        cov=ld*np.outer(signs,signs)*np.outer(sy,sy)
        precision_x=np.linalg.solve(cov,bx);den=float(bx@precision_x)
        beta=float(precision_x@by/den);residual=by-beta*bx
        q=float(residual@np.linalg.solve(cov,residual));se=math.sqrt(max(1,q/4)/den)
        result={'b':beta,'se':se,'pval':math.erfc(abs(beta/se)/math.sqrt(2)),'Q':q}
        saved=next(r for r in rows(DATA/'core/eur_gls_expected.csv') if r['outcome_id']==endpoint)
        for field in result:close(result[field],saved[field])
        mr.append({'outcome':endpoint,'model':'external EUR covariance sensitivity',**result})
    ukb=rows(DATA/'ukb/UKB_RA_Logged_Route_Estimates.csv')
    for saved in ukb:
        data=rows(DATA/'ukb'/f"UKB_RA_Harmonized_{saved['extraction_route']}_action{saved['harmonisation_action']}.csv")
        data=[r for r in data if r['mr_keep'].lower()=='true']
        result=ivw(data)
        assert len(data)==4 and 'rs144309607' not in {r['SNP'] for r in data}
        for field in ('b','se','pval'):close(result[field],saved[field])
        mr.append({'outcome':'ukb-b-9125','model':saved['extraction_route']+'_action'+saved['harmonisation_action'],**result})
    sufficient=rows(DATA/'core/coloc_sufficient_statistics.csv');expected=rows(DATA/'core/coloc_expected.csv')
    posteriors=[]
    keys=['outcome_id','dataset_id','molecular_trait_id','gwas_signal','qtl_signal','p1','p2','p12']
    for r,saved in zip(sufficient,expected):
        assert all(r[k]==saved[k] for k in keys)
        a,b,j=(float(r[k]) for k in ['log_sum_bf_gwas','log_sum_bf_qtl','log_sum_joint_bf'])
        p1,p2,p12=(float(r[k]) for k in ['p1','p2','p12'])
        weights=[0,math.log(p1)+a,math.log(p2)+b,math.log(p1*p2)+a+b+math.log1p(-math.exp(j-a-b)),math.log(p12)+j]
        m=max(weights);normalizer=m+math.log(sum(math.exp(v-m) for v in weights))
        pp=[math.exp(v-normalizer) for v in weights]
        for i,v in enumerate(pp):close(v,saved[f'PP.H{i}.abf'])
        posteriors.append({**{k:r[k] for k in keys},**{f'PP.H{i}':v for i,v in enumerate(pp)}})
    assert len(posteriors)==72
    models=rows(DATA/'core/event_models.csv')
    selected=next(r for r in models if r['event_id'].endswith('ENST00000525976'))
    comparator=next(r for r in models if r['event_id'].endswith('ENST00000524462'))
    intervals=lambda r:{tuple(map(int,pair.split('-'))) for pair in r['intervals_ascending_genomic'].split(';')}
    assert intervals(comparator)-intervals(selected)=={(10351052,10351162)}
    assert not intervals(selected)-intervals(comparator)
    assert int(selected['n_exons'])==19 and int(selected['exonic_length'])==3074
    features=rows(DATA/'core/transcript_features.csv')
    variant=next(r for r in features if r['transcript_id']=='ENST00000525621' and r['feature']=='exon' and int(r['start'])<=10352442<=int(r['end']))
    assert variant['exon_number']=='23'
    out=ROOT/'outputs';out.mkdir(exist_ok=True)
    for name,records in [('recomputed_mr.csv',mr),('recomputed_posteriors.csv',posteriors)]:
        with (out/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    summary={'status':'PASS','finngen_ivw_endpoints':2,'leave_one_out_checks':10,'external_ld_gls_checks':2,
             'ukb_route_checks':len(ukb),'posterior_checks_from_derived_sufficient_statistics':len(posteriors),
             'event_exons':19,'missing_interval_bases':111,'coding_site_retained_in_exon':23,
             'limits':'Arithmetic checks of archived summaries; no independent cohort, new fine-mapping fit, RNA mediation or raw single-cell reprocessing.'}
    (out/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
