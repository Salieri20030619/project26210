#!/usr/bin/env python3
import hashlib,json,os
from pathlib import Path
import pandas as pd

ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_phase71_immune_advanced_v1')
SRC=ROOT/'posthoc_state_audit_v1'; PART=ROOT/'posthoc_state_results_v1.part'; FINAL=ROOT/'posthoc_state_results_v1'
DATASETS=('GSE207493','GSE264586','GSE294109'); EXTERNAL=('GSE264586','GSE294109')
if PART.exists() or FINAL.exists(): raise SystemExit('Refusing existing output')
PART.mkdir()
def read(d,s): return pd.read_csv(SRC/d/f'{d}_{s}',sep='\t')
def write(x,p):
 t=p.with_name(p.name+'.part');x.to_csv(t,sep='\t',index=False,na_rep='NA',lineterminator='\n');os.replace(t,p)
def sha(p):
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
state=pd.concat([read(d,'state_adjusted_summary_v1.tsv') for d in DATASETS],ignore_index=True)
strata=pd.concat([read(d,'state_stratified_summary_v1.tsv') for d in DATASETS],ignore_index=True)
loo=pd.concat([read(d,'response_loo_summary_v1.tsv') for d in DATASETS],ignore_index=True)
status=pd.concat([read(d,'state_audit_status_v1.tsv') for d in DATASETS],ignore_index=True)
state_cross=[]
for rid,g in state.groupby('relation_id'):
 by={r.dataset:r for r in g.itertuples(index=False)}; expected=g.iloc[0].expected
 if expected=='two_sided': n=sum(by[d].decision=='ASSOCIATED' for d in EXTERNAL); dec='TWO_EXTERNAL_ASSOCIATED' if n==2 else ('ONE_EXTERNAL_ASSOCIATED' if n==1 else 'NOT_SUPPORTED')
 else: n=sum(by[d].decision=='SUPPORTED' for d in EXTERNAL); dec='STATE_ADJUSTED_TWO_EXTERNAL' if n==2 else ('STATE_ADJUSTED_ONE_EXTERNAL' if n==1 else 'NOT_SUPPORTED')
 row={'relation_id':rid,'expected':expected,'external_decision':dec,'GSE207493_decision':by['GSE207493'].decision}
 for d in DATASETS: row.update({f'{d}_median_rho':by[d].median_rho,f'{d}_q_value':by[d].q_value,f'{d}_positive_fraction':by[d].positive_fraction,f'{d}_decision':by[d].decision})
 state_cross.append(row)
state_cross=pd.DataFrame(state_cross).sort_values('relation_id')
loo_cross=[]
for (edge,gene),g in loo.groupby(['edge','left_out_gene']):
 by={r.dataset:r for r in g.itertuples(index=False)}; n=sum(by[d].decision=='SUPPORTED' for d in EXTERNAL); dec='TWO_EXTERNAL_SUPPORTED' if n==2 else ('ONE_EXTERNAL_SUPPORTED' if n==1 else 'NOT_SUPPORTED')
 row={'edge':edge,'left_out_gene':gene,'external_decision':dec,'GSE207493_decision':by['GSE207493'].decision}
 for d in DATASETS: row.update({f'{d}_median_rho':by[d].median_rho,f'{d}_q_value':by[d].q_value,f'{d}_decision':by[d].decision})
 loo_cross.append(row)
loo_cross=pd.DataFrame(loo_cross).sort_values(['edge','left_out_gene'])
loo_robust=loo_cross.groupby('edge').agg(variants=('left_out_gene','size'),two_external_supported=('external_decision',lambda x:int((x=='TWO_EXTERNAL_SUPPORTED').sum()))).reset_index()
loo_robust['decision']=loo_robust.apply(lambda r:'LOO_ROBUST' if r.two_external_supported>=3 else 'LOO_NOT_ROBUST',axis=1)
common=[]
for (state_class,rid),g in strata.groupby(['state_class','relation_id']):
 by={r.dataset:r for r in g.itertuples(index=False)}
 available=[d for d in EXTERNAL if d in by]; n=sum(by[d].decision=='SUPPORTED' for d in available)
 dec='TWO_EXTERNAL_SUPPORTED' if len(available)==2 and n==2 else ('ONE_EXTERNAL_SUPPORTED' if n>=1 else 'NOT_SUPPORTED')
 row={'state_class':state_class,'relation_id':rid,'external_decision':dec,'external_datasets_evaluable':len(available)}
 for d in DATASETS:
  if d in by: row.update({f'{d}_n_patients':by[d].n_patients,f'{d}_median_rho':by[d].median_rho,f'{d}_q_value':by[d].q_value,f'{d}_decision':by[d].decision})
 common.append(row)
strata_cross=pd.DataFrame(common).sort_values(['state_class','relation_id'])
write(state,PART/'state_adjusted_summary_all_v1.tsv');write(state_cross,PART/'state_adjusted_cross_cohort_v1.tsv');write(strata,PART/'state_stratified_summary_all_v1.tsv');write(strata_cross,PART/'state_stratified_cross_cohort_v1.tsv');write(loo,PART/'response_loo_summary_all_v1.tsv');write(loo_cross,PART/'response_loo_cross_cohort_v1.tsv');write(loo_robust,PART/'response_loo_robustness_v1.tsv');write(status,PART/'state_audit_status_all_v1.tsv')
lines=['# Phase71 post-hoc state-confounding audit','', '## Fine-state adjusted primary relations','']
for r in state_cross.itertuples(index=False): lines.append(f'- {r.relation_id}: {r.external_decision}; GSE264586 rho={r.GSE264586_median_rho:.4f}, q={r.GSE264586_q_value:.4g}; GSE294109 rho={r.GSE294109_median_rho:.4f}, q={r.GSE294109_q_value:.4g}; discovery={r.GSE207493_decision}.')
lines+=['','## Leave-one-response-gene robustness','']
for r in loo_robust.itertuples(index=False): lines.append(f'- {r.edge}: {r.two_external_supported}/4 variants supported in both external cohorts; {r.decision}.')
lines+=['','## State-stratified diagnostics','']
for r in strata_cross.loc[strata_cross.external_decision.ne('NOT_SUPPORTED')].itertuples(index=False): lines.append(f'- {r.state_class} / {r.relation_id}: {r.external_decision}.')
lines+=['','## Boundary','', 'This post-hoc audit diagnoses frozen cell-state composition and response-score gene dependence. It does not replace the preregistered Phase71 decision or establish measured MTA, receptor activation, lineage conversion, or causality.','']
(PART/'PHASE71_POSTHOC_STATE_AUDIT_SUMMARY_v1.md').write_text('\n'.join(lines),encoding='utf-8')
meta={'status':'PASS_POSTHOC_FINALIZATION_VALIDATION_PENDING','state_rows':len(state),'strata_rows':len(strata),'loo_rows':len(loo),'state_cross_rows':len(state_cross),'loo_cross_rows':len(loo_cross)}
(PART/'run_metadata_v1.json').write_text(json.dumps(meta,indent=2)+'\n')
manifest=[]
for p in sorted(PART.iterdir()):
 if p.is_file():manifest.append({'file':p.name,'size_bytes':p.stat().st_size,'sha256':sha(p)})
write(pd.DataFrame(manifest),PART/'FILE_MANIFEST_SHA256_v1.tsv');(PART/'POSTHOC_FINALIZATION_PASS_v1').write_text('PASS\n');os.replace(PART,FINAL);print(json.dumps(meta))
