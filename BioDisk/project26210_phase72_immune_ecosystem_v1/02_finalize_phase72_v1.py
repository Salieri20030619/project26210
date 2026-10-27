#!/usr/bin/env python3
import hashlib,json,math,os
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_phase72_immune_ecosystem_v1');RUNS=ROOT/'dataset_runs_v1';PART=ROOT/'results_v1.part';FINAL=ROOT/'results_v1'
DIR={'GSE207493':'GSE207493_v2','GSE264586':'GSE264586','GSE294109':'GSE294109'};EXT=('GSE264586','GSE294109')
if PART.exists() or FINAL.exists():raise SystemExit('Refusing existing output')
PART.mkdir()
def read(d,s):return pd.read_csv(RUNS/DIR[d]/f'{d}_{s}',sep='\t')
def write(x,p):t=p.with_name(p.name+'.part');x.to_csv(t,sep='\t',index=False,na_rep='NA',lineterminator='\n');os.replace(t,p)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
summ=pd.concat([read(d,'primary_relations_v1.tsv') for d in DIR],ignore_index=True);sens=pd.concat([read(d,'composition_sensitivity_v1.tsv') for d in DIR],ignore_index=True);status=pd.concat([read(d,'status_v1.tsv') for d in DIR],ignore_index=True);coverage=pd.concat([read(d,'program_coverage_v1.tsv') for d in DIR],ignore_index=True)
if len(summ)!=24 or not status.status.eq('PASS_PHASE72_DATASET_V1').all():raise RuntimeError('Invalid dataset closure')
cross=[]
for rid,g in summ.groupby('relation_id'):
 by={r.dataset:r for r in g.itertuples(index=False)};expected=g.iloc[0].expected;available=[d for d in EXT if d in by];n=sum(by[d].decision=='COHORT_SUPPORTED' for d in available)
 if len(available)==2:decision='TWO_EXTERNAL_SUPPORTED' if n==2 else ('ONE_EXTERNAL_SUPPORTED' if n==1 else 'NOT_SUPPORTED')
 else:decision='GSE294109_SUPPORTED' if 'GSE294109' in by and by['GSE294109'].decision=='COHORT_SUPPORTED' else 'NOT_SUPPORTED'
 row={'relation_id':rid,'predictor':g.iloc[0].predictor,'outcome':g.iloc[0].outcome,'expected':expected,'external_datasets_evaluable':len(available),'external_decision':decision,'GSE207493_decision':by['GSE207493'].decision if 'GSE207493'in by else 'UNEVALUABLE'}
 for d in DIR:
  if d in by:row.update({f'{d}_n_patients':int(by[d].n_patients),f'{d}_rho':by[d].rho,f'{d}_q_value':by[d].q_value,f'{d}_decision':by[d].decision})
 if len(available)==2:
  rs=np.array([by[d].rho for d in available],float);ns=np.array([by[d].n_patients for d in available],float);w=ns-3;z=np.arctanh(np.clip(rs,-.999999,.999999));zm=(w*z).sum()/w.sum();se=1/math.sqrt(w.sum());p=2*stats.norm.sf(abs(zm/se));q=((w*(z-zm)**2).sum());i2=max(0,(q-(len(z)-1))/q)*100 if q>0 else 0;row.update({'meta_rho':math.tanh(zm),'meta_p_value':p,'heterogeneity_Q':q,'heterogeneity_I2_pct':i2})
 cross.append(row)
cross=pd.DataFrame(cross).sort_values('relation_id')
cross['meta_q_value']=np.nan;mask=cross.meta_p_value.notna() if 'meta_p_value'in cross else pd.Series(False,index=cross.index)
if mask.any():cross.loc[mask,'meta_q_value']=stats.false_discovery_control(cross.loc[mask,'meta_p_value'].to_numpy(),method='bh')
long_delta=read('GSE264586','longitudinal_deltas_v1.tsv');long_rel=read('GSE264586','longitudinal_relations_v1.tsv');response=read('GSE294109','response_exploratory_v1.tsv')
phase70=pd.read_csv('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase70_cell2location_spatial_composition_v1/results_v1/primary_decisions_v1.tsv',sep='\t')
def cdec(r):
 x=cross.loc[cross.relation_id.eq(r)];return x.iloc[0].external_decision if len(x) else 'UNEVALUABLE'
tri=pd.DataFrame([
 {'mechanism_segment':'SPP1_state_to_T_dysfunction','phase70_spatial':'ONE_PREPARATION_SUPPORTED','phase72_patient':cdec('E03'),'integrated_call':'MULTIMODAL_SUPPORTED' if cdec('E03')!='NOT_SUPPORTED' else 'SPATIAL_ONLY_PARTIAL'},
 {'mechanism_segment':'SPP1_state_to_cytotoxicity','phase70_spatial':'NOT_DIRECTLY_TESTED','phase72_patient':cdec('E04'),'integrated_call':cdec('E04')},
 {'mechanism_segment':'SPP1_state_to_vascular','phase70_spatial':'TWO_PREPARATION_SUPPORTED','phase72_patient':'NOT_TESTED','integrated_call':'SPATIAL_TWO_PREPARATION_SUPPORTED'},
 {'mechanism_segment':'MTA_response_to_T_dysfunction','phase70_spatial':'NOT_TESTED','phase72_patient':cdec('E01'),'integrated_call':cdec('E01')},
 {'mechanism_segment':'myeloid_to_T_suppression_composite','phase70_spatial':'NOT_TESTED','phase72_patient':cdec('E09'),'integrated_call':cdec('E09')},
])
write(summ.sort_values(['relation_id','dataset']),PART/'phase72_primary_relations_all_v1.tsv');write(cross,PART/'phase72_cross_cohort_decisions_v1.tsv');write(sens.sort_values(['relation_id','dataset']),PART/'phase72_composition_sensitivity_all_v1.tsv');write(status.sort_values('dataset'),PART/'phase72_dataset_status_all_v1.tsv');write(coverage.sort_values(['program','dataset']),PART/'phase72_program_coverage_all_v1.tsv');write(long_delta,PART/'GSE264586_longitudinal_deltas_v1.tsv');write(long_rel,PART/'GSE264586_longitudinal_relations_v1.tsv');write(response,PART/'GSE294109_response_exploratory_v1.tsv');write(tri,PART/'PHASE70_PHASE72_TRIANGULATION_v1.tsv')
lines=['# Phase72 cross-lineage immune ecosystem analysis','',f'- {int(status.immune_cells.sum()):,} frozen immune cells; patient-level inference in 19 discovery, 12 and 13 external patients.','- GSE264586 Treg relations are unevaluable by design because no frozen Treg annotation exists.','', '## Primary relations','']
for r in cross.itertuples(index=False):
 a=f'- {r.relation_id} {r.predictor} -> {r.outcome}: external={r.external_decision}; discovery={r.GSE207493_decision}'
 if hasattr(r,'GSE264586_rho') and pd.notna(r.GSE264586_rho):a+=f'; GSE264586 rho={r.GSE264586_rho:.3f}, q={r.GSE264586_q_value:.3g}'
 if hasattr(r,'GSE294109_rho') and pd.notna(r.GSE294109_rho):a+=f'; GSE294109 rho={r.GSE294109_rho:.3f}, q={r.GSE294109_q_value:.3g}'
 lines.append(a+'.')
lines+=['','## Exploratory layers','',f'- GSE264586 longitudinal: {len(long_delta)} repeated patients; descriptive only.',f'- GSE294109 response comparison: {int(response.n_response.max())} CR/PR versus {int(response.n_no_effect.max())} No_effect; exploratory only.','', '## Boundary','', 'Cross-lineage patient coupling is not direct cell-cell communication, measured MTA, causal immune suppression or a validated clinical predictor. Phase69 remains controlling.','']
(PART/'PHASE72_IMMUNE_ECOSYSTEM_SUMMARY_v1.md').write_text('\n'.join(lines),encoding='utf-8')
meta={'status':'PASS_PHASE72_FINALIZATION_VALIDATION_PENDING','summary_rows':len(summ),'cross_rows':len(cross),'immune_cells':int(status.immune_cells.sum()),'longitudinal_patients':len(long_delta),'response_n':int(response.n_response.max()+response.n_no_effect.max())};(PART/'run_metadata_v1.json').write_text(json.dumps(meta,indent=2)+'\n')
manifest=[]
for p in sorted(PART.iterdir()):
 if p.is_file():manifest.append({'file':p.name,'size_bytes':p.stat().st_size,'sha256':sha(p)})
write(pd.DataFrame(manifest),PART/'FILE_MANIFEST_SHA256_v1.tsv');(PART/'PHASE72_FINALIZATION_PASS_v1').write_text('PASS\n');os.replace(PART,FINAL);print(json.dumps(meta))
