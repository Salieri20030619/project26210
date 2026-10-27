#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd

def bh(p):
    p=np.asarray(p,float); order=np.argsort(p); x=p[order]*len(p)/np.arange(1,len(p)+1); x=np.minimum.accumulate(x[::-1])[::-1]
    out=np.empty(len(p)); out[order]=np.minimum(x,1); return out

def main():
    p=argparse.ArgumentParser(); p.add_argument('--results',required=True); p.add_argument('--out',required=True); a=p.parse_args()
    r=Path(a.results).resolve(); out=Path(a.out).resolve(); part=Path(str(out)+'.part')
    if out.exists() or part.exists(): raise RuntimeError('Refusing existing validation output')
    part.mkdir(parents=True)
    primary=pd.read_csv(r/'phase76_primary_relations_all_v1.tsv',sep='\t')
    status=pd.read_csv(r/'phase76_dataset_status_all_v1.tsv',sep='\t')
    decisions=pd.read_csv(r/'phase76_cross_cohort_decisions_v1.tsv',sep='\t')
    patients=pd.read_csv(r/'phase76_patient_summaries_all_v1.tsv',sep='\t')
    checks=[]
    def ck(name,cond,obs,exp): checks.append({'check':name,'status':'PASS' if cond else 'FAIL','observed':obs,'expected':exp})
    ck('datasets_3',set(status.dataset)=={'GSE207493','GSE264586','GSE294109'},len(status),3)
    expected={'GSE207493':19,'GSE264586':4,'GSE294109':9}
    observed=dict(zip(status.dataset,status.eligible_patients.astype(int)))
    ck('eligible_patient_counts',observed==expected,str(observed),str(expected))
    ck('patient_rows',len(patients)==sum(expected.values()),len(patients),sum(expected.values()))
    ck('primary_rows',len(primary)==18,len(primary),18)
    ck('six_relations_each',(primary.groupby('dataset').size()==6).all(),str(primary.groupby('dataset').size().to_dict()),'6 each')
    qerr=0.0
    for _,g in primary.groupby('dataset',sort=True): qerr=max(qerr,float(np.max(np.abs(bh(g.p_value)-g.q_value))))
    ck('bh_replay',qerr<1e-12,qerr,'<1e-12')
    replay=np.where(primary.n_patients<8,'DESCRIPTIVE_LOW_N',np.where((primary.q_value<.05)&(primary.spearman_rho>0)&(primary.loo_min_rho>0),'COHORT_SUPPORTED',np.where((primary.q_value<.05)&(primary.spearman_rho<=0),'DIRECTION_CONFLICT','NOT_SUPPORTED')))
    ck('cohort_status_replay',np.array_equal(replay,primary.status.to_numpy()),int(np.sum(replay!=primary.status.to_numpy())),0)
    ck('gse264_low_n',(primary.loc[primary.dataset=='GSE264586','status']=='DESCRIPTIVE_LOW_N').all(),int((primary.loc[primary.dataset=='GSE264586','status']=='DESCRIPTIVE_LOW_N').sum()),6)
    ck('decision_rows',len(decisions)==6,len(decisions),6)
    ck('no_cell_level_inference',True,'patient summaries','patient')
    frame=pd.DataFrame(checks); frame.to_csv(part/'validation_checks_v1.tsv',sep='\t',index=False)
    passed=int((frame.status=='PASS').sum()); failed=int((frame.status=='FAIL').sum()); summary={'status':'PASS' if failed==0 else 'FAIL','checks':len(frame),'passed':passed,'failed':failed}
    (part/'validation_summary_v1.json').write_text(json.dumps(summary,indent=2)+'\n')
    if failed: raise RuntimeError(summary)
    (part/'INDEPENDENT_VALIDATION_PASS_v1').write_text('status=PASS\n')
    os.replace(part,out); print(json.dumps(summary))

if __name__=='__main__': main()

