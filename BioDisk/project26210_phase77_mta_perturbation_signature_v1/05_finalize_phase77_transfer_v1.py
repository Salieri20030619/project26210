#!/usr/bin/env python3
import argparse,json,os
from pathlib import Path
import pandas as pd

DATASETS=['GSE207493','GSE264586','GSE294109']
def main():
 p=argparse.ArgumentParser();p.add_argument('--singlecell',required=True);p.add_argument('--spatial',required=True);p.add_argument('--training',required=True);p.add_argument('--out',required=True);a=p.parse_args();sc=Path(a.singlecell);sp=Path(a.spatial);tr=Path(a.training);out=Path(a.out);part=Path(str(out)+'.part')
 if out.exists() or part.exists():raise RuntimeError('refusing existing output')
 part.mkdir(parents=True);summ=[];cov=[];status=[]
 for d in DATASETS:
  x=sc/d;summ.append(pd.read_csv(x/f'{d}_relation_summary_v1.tsv',sep='\t'));cov.append(pd.read_csv(x/f'{d}_signature_coverage_v1.tsv',sep='\t'));status.append(pd.read_csv(x/f'{d}_status_v1.tsv',sep='\t'))
 summ=pd.concat(summ,ignore_index=True);cov=pd.concat(cov,ignore_index=True);status=pd.concat(status,ignore_index=True)
 decisions=[]
 for rid,g in summ.groupby('relation_id',sort=True):
  states=dict(zip(g.dataset,g.statistical_status));n=sum(states.get(x)=='COHORT_SUPPORTED' for x in ['GSE264586','GSE294109']);decisions.append({'relation_id':rid,'gse207493_status':states.get('GSE207493','MISSING'),'gse264586_status':states.get('GSE264586','MISSING'),'gse294109_status':states.get('GSE294109','MISSING'),'external_supported_count':n,'statistical_decision':'TWO_EXTERNAL_SUPPORTED' if n==2 else ('ONE_EXTERNAL_SUPPORTED' if n==1 else 'NOT_SUPPORTED'),'evidence_use':'EXPLORATORY_ONLY_CONSENSUS_NOT_QUALIFIED'})
 decisions=pd.DataFrame(decisions);spdec=pd.read_csv(sp/'decisions_v1.tsv',sep='\t');qual=pd.read_csv(tr/'phase77_consensus_qualification_v1.tsv',sep='\t')
 summ.to_csv(part/'singlecell_relation_summary_all_v1.tsv',sep='\t',index=False);cov.to_csv(part/'singlecell_signature_coverage_all_v1.tsv',sep='\t',index=False);status.to_csv(part/'singlecell_dataset_status_all_v1.tsv',sep='\t',index=False);decisions.to_csv(part/'singlecell_cross_cohort_decisions_v1.tsv',sep='\t',index=False);spdec.to_csv(part/'spatial_decisions_v1.tsv',sep='\t',index=False)
 scsup=decisions.loc[decisions.statistical_decision!='NOT_SUPPORTED','relation_id'].tolist();spsup=spdec.loc[spdec.statistical_decision!='NOT_SUPPORTED','relation_id'].tolist()
 lines=['# Phase77 perturbation-derived MTA signature transfer','',f"- Consensus qualification: `{qual.qualification.iloc[0]}`.",'- The 89-gene direct MTA signature was projected as exploratory sensitivity only.','- No transfer result is eligible to upgrade the mechanism because the consensus signature failed the preregistered bidirectional threshold.','','## Exploratory single-cell statistical support','']
 lines += [f'- {x}' for x in scsup] if scsup else ['- None.'];lines += ['','## Exploratory spatial statistical support',''];lines += [f'- {x}' for x in spsup] if spsup else ['- None.'];lines += ['','## Boundary','','The discovery perturbation is human fibroblast MTA exposure, not macrophage exposure. Statistical transfer cannot establish measured MTA, TYMP causality or macrophage-specific response.']
 (part/'PHASE77_TRANSFER_SUMMARY_v1.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');(part/'run_metadata_v1.json').write_text(json.dumps({'status':'PASS','consensus_qualification':qual.qualification.iloc[0],'evidence_use':'EXPLORATORY_ONLY'},indent=2)+'\n')
 os.replace(part,out);print(json.dumps({'status':'PASS_PHASE77_TRANSFER_FINALIZE'}))
if __name__=='__main__':main()
