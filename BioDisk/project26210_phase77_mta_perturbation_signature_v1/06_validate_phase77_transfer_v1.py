#!/usr/bin/env python3
import argparse,json,os
from pathlib import Path
import numpy as np
import pandas as pd
def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);x=p[o]*len(p)/np.arange(1,len(p)+1);x=np.minimum.accumulate(x[::-1])[::-1];r=np.empty(len(p));r[o]=np.minimum(x,1);return r
def main():
 p=argparse.ArgumentParser();p.add_argument('--results',required=True);p.add_argument('--singlecell',required=True);p.add_argument('--spatial',required=True);p.add_argument('--training',required=True);p.add_argument('--out',required=True);a=p.parse_args();r=Path(a.results);sc=Path(a.singlecell);sp=Path(a.spatial);tr=Path(a.training);out=Path(a.out);part=Path(str(out)+'.part')
 if out.exists() or part.exists():raise RuntimeError('refusing existing validation')
 part.mkdir(parents=True);summ=pd.read_csv(r/'singlecell_relation_summary_all_v1.tsv',sep='\t');cov=pd.read_csv(r/'singlecell_signature_coverage_all_v1.tsv',sep='\t');dec=pd.read_csv(r/'singlecell_cross_cohort_decisions_v1.tsv',sep='\t');sprel=pd.read_csv(sp/'relations_by_section_v1.tsv',sep='\t');spmeta=pd.read_csv(sp/'relation_meta_v1.tsv',sep='\t');spqc=pd.read_csv(sp/'section_qc_v1.tsv',sep='\t');spdec=pd.read_csv(sp/'decisions_v1.tsv',sep='\t');qual=pd.read_csv(tr/'phase77_consensus_qualification_v1.tsv',sep='\t')
 checks=[]
 def ck(n,c,o,e):checks.append({'check':n,'status':'PASS' if c else 'FAIL','observed':o,'expected':e})
 ck('consensus_not_qualified',qual.qualification.iloc[0]=='CONSENSUS_SIGNATURE_NOT_QUALIFIED',qual.qualification.iloc[0],'CONSENSUS_SIGNATURE_NOT_QUALIFIED');ck('singlecell_rows',len(summ)==6,len(summ),6);ck('three_datasets',set(summ.dataset)=={'GSE207493','GSE264586','GSE294109'},';'.join(sorted(set(summ.dataset))),'3 datasets');ck('coverage_rows',len(cov)==6,len(cov),6);ck('coverage_eligible',cov.eligible.astype(bool).all(),int(cov.eligible.astype(bool).sum()),6);ck('direction_coverage',((cov.present_up>=8)&(cov.present_down>=8)).all(),f"up_min={cov.present_up.min()},down_min={cov.present_down.min()}",'both >=8')
 qerr=0
 for _,g in summ.groupby('dataset'):qerr=max(qerr,float(np.max(np.abs(bh(g.p_value)-g.q_value))))
 ck('singlecell_bh_replay',qerr<1e-12,qerr,'<1e-12');ck('singlecell_exploratory_only',(summ.evidence_use=='EXPLORATORY_DIRECT_SIGNATURE_ONLY').all(),int((summ.evidence_use=='EXPLORATORY_DIRECT_SIGNATURE_ONLY').sum()),6);ck('decision_rows',len(dec)==2,len(dec),2);ck('decision_no_upgrade',(dec.evidence_use=='EXPLORATORY_ONLY_CONSENSUS_NOT_QUALIFIED').all(),int((dec.evidence_use=='EXPLORATORY_ONLY_CONSENSUS_NOT_QUALIFIED').sum()),2)
 ck('spatial_sections',len(spqc)==24 and (spqc.barcode_match=='PASS').all(),len(spqc),24);ck('spatial_relation_rows',len(sprel)==48,len(sprel),48);ck('spatial_meta_rows',len(spmeta)==4,len(spmeta),4);sqerr=0
 for _,g in spmeta.groupby('preparation'):sqerr=max(sqerr,float(np.max(np.abs(bh(g.p_value)-g.q_value))))
 ck('spatial_bh_replay',sqerr<1e-12,sqerr,'<1e-12');ck('spatial_exploratory_only',(spdec.evidence_use=='EXPLORATORY_ONLY_CONSENSUS_NOT_QUALIFIED').all(),int((spdec.evidence_use=='EXPLORATORY_ONLY_CONSENSUS_NOT_QUALIFIED').sum()),2)
 frame=pd.DataFrame(checks);frame.to_csv(part/'validation_checks_v1.tsv',sep='\t',index=False);failed=int((frame.status=='FAIL').sum());summary={'status':'PASS' if failed==0 else 'FAIL','checks':len(frame),'passed':int((frame.status=='PASS').sum()),'failed':failed};(part/'validation_summary_v1.json').write_text(json.dumps(summary,indent=2)+'\n')
 if failed:raise RuntimeError(summary)
 (part/'INDEPENDENT_VALIDATION_PASS_v1').write_text('status=PASS\n');os.replace(part,out);print(json.dumps(summary))
if __name__=='__main__':main()
