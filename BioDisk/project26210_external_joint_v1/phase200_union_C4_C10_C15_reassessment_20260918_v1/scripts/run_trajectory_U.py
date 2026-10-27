from pathlib import Path
import json,warnings
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import wilcoxon,spearmanr,rankdata
warnings.filterwarnings('ignore',category=RuntimeWarning)
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');T=R/'07_trajectory/monocle3_v3';O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase200_union_C4_C10_C15_reassessment_20260918_v1/trajectory/U_v2');O.mkdir(exist_ok=False)
assert json.loads((T/'STATUS.json').read_text())['status']=='COMPLETE_PARTITIONED_TRAJECTORIES'
def bh(v):
 a=np.asarray(v,float);q=np.full(len(a),np.nan);ok=np.isfinite(a);x=a[ok];ii=np.argsort(x);u=np.minimum.accumulate((x[ii]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1];back=np.empty(len(x));back[ii]=np.minimum(u,1);q[ok]=back;return q
def pval(x):
 x=np.asarray(x);x=x[np.isfinite(x)]
 if len(x)<3:return np.nan
 return float(wilcoxon(x,method='approx',correction=True).pvalue) if np.any(x) else 1.
def summarize(v):
 x=np.asarray(v,float);x=x[np.isfinite(x)];return dict(n_patients=len(x),median_effect=float(np.median(x)) if len(x) else np.nan,positive=int((x>0).sum()),negative=int((x<0).sum()),P=pval(x))
def partialrho(pt,y,depth):
 if np.std(y)==0 or np.std(pt)==0:return np.nan
 cov=np.column_stack([np.ones(len(pt)),rankdata(depth)]);a=rankdata(pt);b=rankdata(y);a-=cov@np.linalg.lstsq(cov,a,rcond=None)[0];b-=cov@np.linalg.lstsq(cov,b,rcond=None)[0]
 if np.std(a)<1e-10 or np.std(b)<1e-10:return np.nan
 return float(np.corrcoef(a,b)[0,1])
scores=pd.read_csv(R/'07_trajectory/dynamics_v2/CELL_EXPRESSION_PROGRAM_SCORES.tsv.gz',sep='\t',index_col=0)
m=pd.read_csv(T/'ALL_CELL_PSEUDOTIME.tsv.gz',sep='\t',dtype={'candidate':str,'partition':str},low_memory=False);m=m.set_index('cell_id');assert m.index.isin(scores.index).all();pairs=[];trends=[]
for mode,md0 in m.groupby('mode'):
 for scope in ['all_tissues','tumor_only']:
  md=md0.copy()
  if scope=='tumor_only':md=md[md.source_tissue.astype(str).str.lower().eq('tumor')]
  versions=['pt_primary','pt_root_q25','pt_root_q50','dpt_pseudotime'] if mode=='primary' else ['pt_primary']
  for pv in versions:
   for (study,patient,part),z in md.groupby(['study','curated_patient_id','partition'],observed=True):
    z=z[np.isfinite(z[pv])];c=z.candidate.astype(str);c4=c.isin(['4', '10', '15']);mono=c.eq('Monocyte');other=~c.isin(['4', '10', '15', '1', 'Monocyte'])
    if c4.sum()==0:continue
    for endpoint,mask in [('U_minus_Monocyte',mono),('U_minus_other_macrophages',other)]:
     if c4.sum()>=20 and mask.sum()>=20:
      pairs.append(dict(mode=mode,scope=scope,pt_version=pv,study=study,patient=patient,partition=part,endpoint=endpoint,U_cells=int(c4.sum()),comparator_cells=int(mask.sum()),effect=float(z.loc[c4,pv].median()-z.loc[mask,pv].median())))
    if len(z)<50 or z[pv].nunique()<10:continue
    pt=z[pv].to_numpy();depth=np.log1p(z.total_counts_raw.to_numpy());ss=scores.loc[z.index]
    for feature in scores.columns:
     y=ss[feature].to_numpy();rho=float(spearmanr(pt,y).statistic) if np.std(y)>0 else np.nan
     trends.append(dict(mode=mode,scope=scope,pt_version=pv,study=study,patient=patient,partition=part,feature=feature,family=feature.split(':')[0],cells=len(z),rho=rho,depth_adjusted_rho=partialrho(pt,y,depth)))
 print('DYNAMICS_MODE',mode,flush=True)
p=pd.DataFrame(pairs);tr=pd.DataFrame(trends);p.to_csv(O/'PARTITION_PATIENT_PAIRED_PSEUDOTIME.tsv',sep='\t',index=False);tr.to_csv(O/'PARTITION_PATIENT_DYNAMIC_TRENDS.tsv.gz',sep='\t',index=False)
# Collapse multiple eligible partitions to one equally weighted patient median; never treat partitions as independent n.
pg=p.groupby(['mode','scope','pt_version','study','patient','endpoint'],as_index=False).effect.median();tg=tr.groupby(['mode','scope','pt_version','study','patient','family','feature'],as_index=False)[['rho','depth_adjusted_rho']].median();pg.to_csv(O/'PATIENT_PAIRED_PSEUDOTIME.tsv',sep='\t',index=False);tg.to_csv(O/'PATIENT_DYNAMIC_TRENDS.tsv',sep='\t',index=False)
pr=[];ts=[]
for k,g in pg.groupby(['mode','scope','pt_version','endpoint']):pr.append(dict(zip(['mode','scope','pt_version','endpoint'],k))|summarize(g.effect))
for k,g in tg.groupby(['mode','scope','pt_version','family','feature']):
 for adj in ['rho','depth_adjusted_rho']:ts.append(dict(zip(['mode','scope','pt_version','family','feature'],k))|{'adjustment':adj}|summarize(g[adj]))
pr=pd.DataFrame(pr);ts=pd.DataFrame(ts)
for df,keys in [(pr,['mode','scope','pt_version']),(ts,['mode','scope','pt_version','family','adjustment'])]:
 df['FDR']=np.nan
 for _,ix in df.groupby(keys).groups.items():df.loc[ix,'FDR']=bh(df.loc[ix,'P'])
 df['reporting_tier']=np.where(df.n_patients<3,'descriptive',np.where(df.n_patients<10,'small_n_exploratory',np.where(df.FDR<.05,'FDR_supported',np.where(df.P<.05,'nominal_exploratory','not_supported'))))
pr.to_csv(O/'PSEUDOTIME_SUMMARY.tsv',sep='\t',index=False);ts.to_csv(O/'DYNAMIC_TREND_SUMMARY.tsv',sep='\t',index=False)
studies=[]
for k,g in tg.groupby(['mode','scope','pt_version','study','family','feature']):studies.append(dict(zip(['mode','scope','pt_version','study','family','feature'],k))|summarize(g.rho))
pd.DataFrame(studies).to_csv(O/'DYNAMIC_STUDY_DIRECTIONS.tsv',sep='\t',index=False)
assert not pg.duplicated(['mode','scope','pt_version','patient','endpoint']).any();assert not tg.duplicated(['mode','scope','pt_version','patient','feature']).any()
(O/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_TRAJECTORY_DYNAMIC_STATISTICS','pair_summary_rows':len(pr),'dynamic_summary_rows':len(ts),'figures_created':False,'same_partition_comparisons_only':True,'patient_unique_checks':True},indent=2));print('DYNAMICS_DONE',flush=True)
