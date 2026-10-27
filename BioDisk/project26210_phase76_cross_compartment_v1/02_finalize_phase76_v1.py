#!/usr/bin/env python3
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np
import pandas as pd

DATASETS = ['GSE207493','GSE264586','GSE294109']
EXTERNAL = ['GSE264586','GSE294109']

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def atomic(frame,path):
    tmp=Path(str(path)+'.part'); frame.to_csv(tmp,sep='\t',index=False,na_rep='NA'); os.replace(tmp,path)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--runs',required=True); p.add_argument('--out',required=True); a=p.parse_args()
    runs=Path(a.runs).resolve(); out=Path(a.out).resolve(); part=Path(str(out)+'.part')
    if out.exists() or part.exists(): raise RuntimeError('Refusing existing output')
    part.mkdir(parents=True)
    prim=[]; sens=[]; units=[]; patients=[]; status=[]; coverage=[]; descriptive=[]; manifest=[]
    for ds in DATASETS:
        d=runs/ds
        files={
          'primary':d/f'{ds}_primary_relations_v1.tsv','sensitivity':d/f'{ds}_sensitivity_relations_v1.tsv',
          'unit':d/f'{ds}_unit_summary_v1.tsv','patient':d/f'{ds}_patient_summary_v1.tsv',
          'status':d/f'{ds}_status_v1.tsv','coverage':d/f'{ds}_program_coverage_v1.tsv',
          'descriptive':d/f'{ds}_unit_relations_descriptive_v1.tsv'}
        for kind,path in files.items(): manifest.append({'dataset':ds,'kind':kind,'path':str(path),'bytes':path.stat().st_size,'sha256':sha256(path)})
        prim.append(pd.read_csv(files['primary'],sep='\t')); sens.append(pd.read_csv(files['sensitivity'],sep='\t'))
        units.append(pd.read_csv(files['unit'],sep='\t')); patients.append(pd.read_csv(files['patient'],sep='\t'))
        status.append(pd.read_csv(files['status'],sep='\t')); coverage.append(pd.read_csv(files['coverage'],sep='\t'))
        descriptive.append(pd.read_csv(files['descriptive'],sep='\t'))
    prim=pd.concat(prim,ignore_index=True); sens=pd.concat(sens,ignore_index=True); units=pd.concat(units,ignore_index=True)
    patients=pd.concat(patients,ignore_index=True); status=pd.concat(status,ignore_index=True); coverage=pd.concat(coverage,ignore_index=True)
    descriptive=pd.concat(descriptive,ignore_index=True)

    decision=[]
    for relation,sub in prim.groupby('relation_id',sort=True):
        states=dict(zip(sub.dataset,sub.status)); supported=[d for d in EXTERNAL if states.get(d)=='COHORT_SUPPORTED']
        formal=[d for d in EXTERNAL if int(sub.loc[sub.dataset==d,'n_patients'].iloc[0])>=8]
        if len(supported)==2: call='TWO_EXTERNAL_SUPPORTED'
        elif len(supported)==1: call='ONE_EXTERNAL_SUPPORTED'
        elif len(formal)<2: call='INSUFFICIENT_TWO_EXTERNAL_FEASIBILITY'
        else: call='NOT_EXTERNALLY_SUPPORTED'
        decision.append({'relation_id':relation,'gse207493_status':states.get('GSE207493','MISSING'),
          'gse264586_status':states.get('GSE264586','MISSING'),'gse294109_status':states.get('GSE294109','MISSING'),
          'formal_external_cohorts':';'.join(formal) if formal else 'None','external_supported_count':len(supported),'decision':call})
    decision=pd.DataFrame(decision)

    meta=[]
    for relation,sub in prim.groupby('relation_id',sort=True):
        s=sub[np.isfinite(sub.spearman_rho)&(sub.n_patients>=4)].copy()
        z=np.arctanh(np.clip(s.spearman_rho.to_numpy(float),-0.999999,0.999999)); w=np.maximum(s.n_patients.to_numpy(float)-3,1)
        mz=float(np.sum(w*z)/np.sum(w)); se=math.sqrt(1/np.sum(w)); zr=mz/se
        pval=float(math.erfc(abs(zr)/math.sqrt(2)))
        meta.append({'relation_id':relation,'n_cohorts':len(s),'total_patients':int(s.n_patients.sum()),
          'weighted_fisher_rho':float(np.tanh(mz)),'z_value':zr,'p_value':pval})
    meta=pd.DataFrame(meta); meta['q_value']=pd.Series(meta.p_value).rank(method='first')
    order=np.argsort(meta.p_value.to_numpy()); ranked=meta.p_value.to_numpy()[order]*len(meta)/np.arange(1,len(meta)+1)
    ranked=np.minimum.accumulate(ranked[::-1])[::-1]; q=np.empty(len(meta)); q[order]=np.minimum(ranked,1); meta['q_value']=q

    chain_required=['P01','P02','P03']
    chain=decision.set_index('relation_id').loc[chain_required]
    if (chain.decision=='TWO_EXTERNAL_SUPPORTED').all(): chain_status='TWO_EXTERNAL_PATIENT_PAIRED_PROXY_CHAIN_SUPPORTED'
    elif (chain.external_supported_count>=1).all(): chain_status='ONE_EXTERNAL_PATIENT_PAIRED_PROXY_CHAIN_SUPPORTED'
    else: chain_status='PATIENT_PAIRED_PROXY_CHAIN_NOT_FULLY_SUPPORTED'

    atomic(prim,part/'phase76_primary_relations_all_v1.tsv'); atomic(sens,part/'phase76_sensitivity_relations_all_v1.tsv')
    atomic(units,part/'phase76_unit_summaries_all_v1.tsv'); atomic(patients,part/'phase76_patient_summaries_all_v1.tsv')
    atomic(status,part/'phase76_dataset_status_all_v1.tsv'); atomic(coverage,part/'phase76_program_coverage_all_v1.tsv')
    atomic(descriptive,part/'phase76_unit_relations_descriptive_all_v1.tsv'); atomic(decision,part/'phase76_cross_cohort_decisions_v1.tsv')
    atomic(meta,part/'phase76_fixed_effect_meta_v1.tsv'); atomic(pd.DataFrame(manifest),part/'INPUT_RESULT_MANIFEST_SHA256_v1.tsv')

    supported=decision.loc[decision.decision.isin(['ONE_EXTERNAL_SUPPORTED','TWO_EXTERNAL_SUPPORTED']),'relation_id'].tolist()
    lines=['# Phase76 patient-paired cross-compartment analysis','',
      f'- Eligible patients: '+', '.join(f"{r.dataset}={int(r.eligible_patients)}" for r in status.itertuples())+'.',
      '- GSE264586 is preregistered as descriptive low-N and cannot provide formal external support.',
      f'- Cross-compartment chain status: `{chain_status}`.','',
      '## Externally supported primary relations','']
    lines += [f'- {x}' for x in supported] if supported else ['- None.']
    lines += ['','## Boundary','',
      'Patient-matched cross-compartment covariation is not measured 5\'-MTA, metabolite transport, or causal mediation. Broad-epithelial and repeated-timepoint analyses are sensitivity/descriptive only.']
    (part/'PHASE76_CROSS_COMPARTMENT_SUMMARY_v1.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (part/'run_metadata_v1.json').write_text(json.dumps({'status':'PASS','chain_status':chain_status,
      'eligible_patients':dict(zip(status.dataset,status.eligible_patients.astype(int))),'inference_unit':'patient'},indent=2)+'\n')
    os.replace(part,out)
    print(json.dumps({'status':'PASS_PHASE76_FINALIZE','chain_status':chain_status}))

if __name__=='__main__': main()

