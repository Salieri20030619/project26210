from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.stats import wilcoxon
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'provenance/remaining_analysis_validation_v2_20260915';O.mkdir(exist_ok=False);checks=[]
def ck(name,value):checks.append({'check':name,'pass':bool(value)})
def bh(v):
 a=np.asarray(v,float);q=np.full(len(a),np.nan);ok=np.isfinite(a);x=a[ok];s=np.argsort(x);u=np.minimum.accumulate((x[s]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1];back=np.empty(len(x));back[s]=np.minimum(u,1);q[ok]=back;return q
m=pd.read_csv(R/'07_trajectory/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str},low_memory=False);c4=pd.read_csv(R/'09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz',sep='\t').cell_id
ck('frozen_C4_ids',set(c4)==set(m.index[m.candidate=='4']));ck('unique_input_ids',m.index.is_unique);ck('patient_count',m.curated_patient_id.nunique()==39);ck('library_count',m.sample_id.nunique()==70);ck('all_C4_has_metadata',m.loc[c4].notna().candidate.all())
t=pd.read_csv(R/'07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str,'partition':str},low_memory=False);roots=pd.read_csv(R/'07_trajectory/monocle3_v3/ROOTS.tsv',sep='\t',dtype={'partition':str})
for _,z in roots.iterrows():
 ck(str(z['mode'])+'_'+z.study+'_'+z.partition+'_'+z.root_version+'_root_original_mono',m.loc[z.root_cell,'candidate']=='Monocyte' and m.loc[z.root_cell,'root_eligible'])
for (mode,study,part),z in t.groupby(['mode','study','partition']):
 for version in ['primary','root_q25','root_q50']:
  x=z['pt_'+version].dropna();rr=roots[(roots['mode']==mode)&(roots.study==study)&(roots.partition==part)&(roots.root_version==version)]
  ck(mode+study+part+version+'_no_forced_root',len(x)==0 or len(rr)>0);ck(mode+study+part+version+'_range',len(x)==0 or ((x>=0)&(x<=1)).all())
for p in (R/'07_trajectory/matrices_v1').rglob('HVG.tsv'):ck(str(p.relative_to(R))+'_no_endpoint_genes',not {'SPP1','TYMP'}&set(p.read_text().splitlines()))
for filename,keys,qcol in [('07_trajectory/dynamics_v2/PSEUDOTIME_SUMMARY.tsv',['mode','scope','pt_version'],'FDR'),('07_trajectory/dynamics_v2/DYNAMIC_TREND_SUMMARY.tsv',['mode','scope','pt_version','family','adjustment'],'FDR'),('08_regulons/statistics_v2_20260915/AUC_CONTRAST_SUMMARY.tsv',['scope','endpoint'],'FDR_all_regulons'),('08_regulons/statistics_v2_20260915/AUC_TRAJECTORY_SUMMARY.tsv',['mode','scope','pt_version','adjustment'],'FDR_all_regulons')]:
 df=pd.read_csv(R/filename,sep='\t');ck(filename+'_n_not_cells',df.n_patients.max()<=39)
 for fam,g in df.groupby(keys):ck(filename+str(fam)+'_BH',np.allclose(bh(g.P),g[qcol],atol=1e-12,equal_nan=True))
# Independently aggregate selected AUC comparisons from original all-cell AUC file.
auc=pd.read_csv(R/'08_regulons/network_v2_20260915/AUC_ALL_CELLS.csv',index_col=0);s=pd.read_csv(R/'08_regulons/statistics_v2_20260915/AUC_CONTRAST_SUMMARY.tsv',sep='\t');avail=[]
for tf in ['JUN','JUNB','CEBPD']:
 cols=[c for c in auc.columns if c.split('(')[0]==tf];avail.append({'TF':tf,'available':bool(cols),'regulons':';'.join(cols)})
 for reg in cols:
  for endpoint in ['C4_minus_Monocyte','C4_minus_other_macrophages']:
   values=[]
   for patient,z in m.groupby('curated_patient_id',observed=True):
    A=z.index[z.candidate=='4'];B=z.index[z.candidate=='Monocyte'] if endpoint=='C4_minus_Monocyte' else z.index[~z.candidate.isin(['4','1','Monocyte'])]
    if len(A)>=20 and len(B)>=20:values.append(float(auc.loc[A,reg].mean()-auc.loc[B,reg].mean()))
   row=s[(s.scope=='all_tissues')&(s.endpoint==endpoint)&(s.regulon==reg)].iloc[0];v=np.asarray(values);p=wilcoxon(v,method='approx',correction=True).pvalue if np.any(v) else 1.
   ck(reg+endpoint+'_independent_n',len(v)==row.n_patients);ck(reg+endpoint+'_independent_median',np.isclose(np.median(v),row.median_effect,atol=1e-12));ck(reg+endpoint+'_independent_P',np.isclose(p,row.P,atol=1e-12))
pd.DataFrame(avail).to_csv(O/'PRESPECIFIED_TF_AVAILABILITY.tsv',sep='\t',index=False)
for d in ['06_functional_programs','07_trajectory','08_regulons','09_downstream_crosswalk/remaining_analysis_v2_20260915']:
 bad=[str(p) for p in (R/d).rglob('*') if p.suffix.lower() in ['.png','.svg','.pdf','.jpg','.jpeg','.tiff','.pptx']];ck(d+'_no_figures',not bad)
ck('crosswalk_all_cells',len(pd.read_csv(R/'09_downstream_crosswalk/remaining_analysis_v2_20260915/NEW_TO_LEGACY_CELL_LABELS.tsv.gz',sep='\t'))==len(m))

# Additional independent full-network aggregation and statistical-family checks.
from scipy.stats import false_discovery_control
for fn,keys,qcol in [
 ('AUC_CONTRAST_SUMMARY.tsv',['scope','endpoint'],'FDR_all_regulons'),
 ('AUC_TRAJECTORY_SUMMARY.tsv',['mode','scope','pt_version','adjustment'],'FDR_all_regulons')]:
 df=pd.read_csv(R/'08_regulons/statistics_v2_20260915'/fn,sep='\t')
 for fam,g in df.groupby(keys):
  sel=g[g.TF.isin(['JUN','JUNB','CEBPD'])];valid=np.isfinite(sel.P)
  if valid.any():ck(fn+str(fam)+'_prespecified_BH',np.allclose(false_discovery_control(sel.loc[valid,'P']),sel.loc[valid,'FDR_prespecified_3TF'],atol=1e-12))
 ck(fn+'_small_n_label',df.loc[df.n_patients<10,'reporting_tier'].isin(['descriptive','small_n_exploratory']).all())
p=pd.read_csv(R/'08_regulons/statistics_v2_20260915/PATIENT_AUC_CONTRASTS.tsv',sep='\t')
for (scope,pat,endpoint),g in p.groupby(['scope','patient','endpoint']):
 md=m[m.curated_patient_id==pat]
 if scope=='tumor_only':md=md[md.source_tissue.str.lower()=='tumor']
 a=md.index[md.candidate=='4'];b=md.index[md.candidate=='Monocyte'] if endpoint=='C4_minus_Monocyte' else md.index[~md.candidate.isin(['4','1','Monocyte'])]
 effect=auc.loc[a].mean()-auc.loc[b].mean()
 ck(str((scope,pat,endpoint))+'_all_regulon_raw_means',len(a)>=20 and len(b)>=20 and np.allclose(effect.loc[g.regulon],g.effect,atol=1e-12))
summary=pd.read_csv(R/'08_regulons/statistics_v2_20260915/AUC_CONTRAST_SUMMARY.tsv',sep='\t')
for keys,g in p.groupby(['scope','endpoint','regulon']):
 row=summary[(summary.scope==keys[0])&(summary.endpoint==keys[1])&(summary.regulon==keys[2])].iloc[0]
 v=g.effect.to_numpy();v=v[np.isfinite(v)]
 pv=wilcoxon(v,method='approx',correction=True).pvalue if len(v)>=3 and np.any(v) else (1. if len(v)>=3 else np.nan)
 ck(str(keys)+'_patient_summary',g.patient.is_unique and len(v)==row.n_patients and np.isclose(np.median(v),row.median_effect) and np.isclose(pv,row.P,equal_nan=True))

(O/'CHECKS.json').write_text(json.dumps(checks,indent=2));status={'status':'PASS' if all(x['pass'] for x in checks) else 'FAIL','checks':len(checks),'passed':sum(x['pass'] for x in checks),'failures':[x for x in checks if not x['pass']],'figures_created':False};(O/'STATUS.json').write_text(json.dumps(status,indent=2));print(json.dumps(status),flush=True)
assert status['status']=='PASS'
