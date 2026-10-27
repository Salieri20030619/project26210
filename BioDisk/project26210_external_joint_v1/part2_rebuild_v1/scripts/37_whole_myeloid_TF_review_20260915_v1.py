from pathlib import Path
import json,datetime,hashlib,sys
import pandas as pd,numpy as np
from scipy.stats import spearmanr,rankdata,wilcoxon,false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1'
O=P/'phase150_whole_myeloid_threeTF_20260915_v1';O.mkdir(exist_ok=False)
M=R/'08_regulons/input_v1/CELL_METADATA.tsv.gz';A=R/'08_regulons/network_v2_20260915/AUC_ALL_CELLS.csv';T=R/'07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz'
m=pd.read_csv(M,sep='\t',index_col=0,dtype={'candidate':str},low_memory=False)
header=pd.read_csv(A,nrows=0);tf=['JUN(+)','JUNB(+)','CEBPD(+)']
a=pd.read_csv(A,index_col=0,usecols=[header.columns[0]]+tf);assert a.index.is_unique and set(a.index)==set(m.index)
t=pd.read_csv(T,sep='\t',index_col=0,dtype={'candidate':str,'partition':str},low_memory=False)
plan=dict(question='Whole Monocyte and mature Macrophage three TF elevation; no C4 selection.',state_population='All frozen 12 mature macrophage candidates pooled versus original Monocyte; C1 boundary excluded.',trajectory_population='All eligible partitions regardless of C4; exclude C1 boundary. Additional mixed partitions require >=20 Monocyte and >=20 mature macrophages in the patient partition.',state_threshold='>=20 cells each per patient',trajectory_threshold='>=50 cells and >=10 distinct pseudotime values per patient partition; patient median of partition rho',models='Reuse frozen main/root/DPT/without-stress-cycle models, no refit',statistics='Two-sided approximate Wilcoxon continuity corrected; n<3 no P; n<10 exploratory; three TF BH within endpoint/scope/mode/pseudotime/population/adjustment.',figures=False,created_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'PLAN.json').write_text(json.dumps(plan,indent=2))
meta=m.groupby('candidate').agg(cells=('study','size'),patients=('curated_patient_id','nunique')).reset_index();meta.to_csv(O/'INPUT_COUNTS.tsv',sep='\t',index=False)
pairs=[]
for scope in ['all_tissues','tumor_only']:
 md=m if scope=='all_tissues' else m[m.source_tissue.str.lower()=='tumor']
 for patient,z in md.groupby('curated_patient_id'):
  mono=z.index[z.candidate=='Monocyte'];macro=z.index[~z.candidate.isin(['1','Monocyte'])]
  if min(len(mono),len(macro))<20:continue
  delta=a.loc[macro].mean()-a.loc[mono].mean()
  for reg in tf:pairs.append(dict(scope=scope,patient=patient,study=z.study.iloc[0],regulon=reg,macro_cells=len(macro),mono_cells=len(mono),effect=float(delta[reg])))
p=pd.DataFrame(pairs);p.to_csv(O/'PATIENT_MACRO_MINUS_MONO.tsv',sep='\t',index=False)
def summary(v):
 v=np.asarray(v,float);v=v[np.isfinite(v)];n=len(v)
 pv=np.nan if n<3 else (wilcoxon(v,method='approx',correction=True).pvalue if np.any(v) else 1.)
 return dict(n_patients=n,median_effect=np.median(v) if n else np.nan,positive=int((v>0).sum()),negative=int((v<0).sum()),P=pv)
def annotate(df,keys):
 df['FDR_threeTF']=np.nan
 for _,idx in df.groupby(keys).groups.items():
  valid=df.loc[idx,'P'].dropna()
  if len(valid):df.loc[valid.index,'FDR_threeTF']=false_discovery_control(valid.to_numpy())
 df['reporting_tier']=np.where(df.n_patients<3,'descriptive',np.where(df.n_patients<10,'small_n_exploratory',np.where(df.FDR_threeTF<.05,'FDR_supported',np.where(df.P<.05,'nominal','not_supported'))))
 return df
s=[]
for k,g in p.groupby(['scope','regulon']):s.append(dict(zip(['scope','regulon'],k))|summary(g.effect))
s=annotate(pd.DataFrame(s),['scope']);s.to_csv(O/'MACRO_MINUS_MONO_SUMMARY.tsv',sep='\t',index=False)
tr=[]
for mode,md0 in t.groupby('mode'):
 for scope in ['all_tissues','tumor_only']:
  md=md0[md0.candidate!='1']
  if scope=='tumor_only':md=md[md.source_tissue.str.lower()=='tumor']
  versions=['pt_primary','pt_root_q25','pt_root_q50','dpt_pseudotime'] if mode=='primary' else ['pt_primary']
  for pv in versions:
   for (study,patient,part),z in md.groupby(['study','curated_patient_id','partition'],observed=True):
    z=z[np.isfinite(z[pv])]
    if len(z)<50 or z[pv].nunique()<10:continue
    mono=int((z.candidate=='Monocyte').sum());macro=len(z)-mono
    x=z[pv].to_numpy();depth=np.log1p(z.total_counts_raw.to_numpy());cov=np.column_stack([np.ones(len(z)),rankdata(depth)])
    xr=rankdata(x);xr=xr-cov@np.linalg.lstsq(cov,xr,rcond=None)[0]
    for reg in tf:
     y=a.loc[z.index,reg].to_numpy();rr=spearmanr(x,y).statistic if np.std(y)>0 else np.nan
     yr=rankdata(y);yr=yr-cov@np.linalg.lstsq(cov,yr,rcond=None)[0];adj=np.corrcoef(xr,yr)[0,1] if np.std(xr)>1e-10 and np.std(yr)>1e-10 else np.nan
     for pop in ['all_partitions']+(['mixed_partitions_only'] if mono>=20 and macro>=20 else []):
      tr.append(dict(mode=mode,scope=scope,pt_version=pv,population=pop,study=study,patient=patient,partition=part,regulon=reg,cells=len(z),mono_cells=mono,macro_cells=macro,rho=rr,depth_adjusted_rho=adj))
tr=pd.DataFrame(tr);tr.to_csv(O/'PARTITION_PATIENT_TF_TRENDS.tsv.gz',sep='\t',index=False)
keys=['mode','scope','pt_version','population','study','patient','regulon'];pt=tr.groupby(keys,as_index=False)[['rho','depth_adjusted_rho']].median();pt.to_csv(O/'PATIENT_TF_TRENDS.tsv',sep='\t',index=False)
sums=[]
for k,g in pt.groupby(['mode','scope','pt_version','population','regulon']):
 for adj in ['rho','depth_adjusted_rho']:
  sums.append(dict(zip(['mode','scope','pt_version','population','regulon'],k))|dict(adjustment=adj)|summary(g[adj]))
d=annotate(pd.DataFrame(sums),['mode','scope','pt_version','population','adjustment']);d.to_csv(O/'WHOLE_TRAJECTORY_TF_SUMMARY.tsv',sep='\t',index=False)
studies=[]
for k,g in p.groupby(['scope','study','regulon']):studies.append(dict(zip(['scope','study','regulon'],k))|{'n_patients':g.patient.nunique(),'median_effect':g.effect.median(),'positive':int((g.effect>0).sum()),'negative':int((g.effect<0).sum())})
pd.DataFrame(studies).to_csv(O/'STATE_STUDY_DIRECTIONS.tsv',sep='\t',index=False)
coverage=[]
for mode,z in t[t.candidate!='1'].groupby('mode'):
 for label,mask in [('Monocyte',z.candidate=='Monocyte'),('all_mature_macrophages',z.candidate!='Monocyte')]:
  g=z[mask];coverage.append(dict(mode=mode,cell_group=label,cells=len(g),finite_primary=int(np.isfinite(g.pt_primary).sum()),fraction=float(np.isfinite(g.pt_primary).mean()),patients=g.curated_patient_id.nunique()))
pd.DataFrame(coverage).to_csv(O/'TRAJECTORY_COVERAGE.tsv',sep='\t',index=False)
assert not p.duplicated(['scope','patient','regulon']).any()
assert not pt.duplicated(['mode','scope','pt_version','population','patient','regulon']).any()
manifest=[]
for path in [M,A,T]:
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 manifest.append(dict(path=str(path),bytes=path.stat().st_size,sha256=h.hexdigest()))
pd.DataFrame(manifest).to_csv(O/'SOURCE_MANIFEST.tsv',sep='\t',index=False)
(O/'STATUS.json').write_text(json.dumps(dict(status='COMPUTATION_COMPLETE_AWAITING_REVIEW',time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),figures_created=False,model_refit=False,scope='whole_monocyte_macrophage_no_C4_condition',python=sys.version),indent=2))
print('GROUP_COMPARISON');print(s.to_string(index=False))
print('PRIMARY_OVERALL');print(d[(d['mode']=='primary')&(d.scope=='all_tissues')&(d.pt_version=='pt_primary')].to_string(index=False))
print('SENSITIVITY');print(d[(d.scope=='all_tissues')&(d.adjustment=='rho')&(d.population=='all_partitions')].to_string(index=False))
print('COVERAGE');print(pd.DataFrame(coverage).to_string(index=False))

