from pathlib import Path
import json,warnings
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,spearmanr,rankdata,hypergeom
warnings.filterwarnings('ignore',category=RuntimeWarning)
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');N=R/'08_regulons/network_v2_20260915';O=R/'08_regulons/statistics_v2_20260915';O.mkdir(exist_ok=False)
def bh(x):
 a=np.asarray(x,float);ok=np.isfinite(a);q=np.full(len(a),np.nan);v=a[ok];s=np.argsort(v);qq=np.minimum.accumulate((v[s]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];z=np.empty(len(v));z[s]=np.minimum(qq,1);q[ok]=z;return q

def summ(x):
 a=np.asarray(x);a=a[np.isfinite(a)];n=len(a);p=np.nan if n<3 else (float(wilcoxon(a,method='approx',correction=True).pvalue) if np.any(a) else 1.);return dict(n_patients=n,median_effect=float(np.median(a)) if n else np.nan,positive=int((a>0).sum()),negative=int((a<0).sum()),P=p)

def rho(x,y):return float(spearmanr(x,y).statistic) if np.std(x)>0 and np.std(y)>0 else np.nan

def partial(x,y,d):
 if np.std(y)==0:return np.nan
 cov=np.column_stack([np.ones(len(x)),rankdata(d)]);a=rankdata(x);b=rankdata(y);a-=cov@np.linalg.lstsq(cov,a,rcond=None)[0];b-=cov@np.linalg.lstsq(cov,b,rcond=None)[0];return float(np.corrcoef(a,b)[0,1]) if np.std(a)>1e-10 and np.std(b)>1e-10 else np.nan
auc=pd.read_csv(N/'AUC_ALL_CELLS.csv',index_col=0);meta=pd.read_csv(R/'08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str},low_memory=False);assert set(auc.index)==set(meta.index);auc=auc.loc[meta.index];assert np.isfinite(auc.to_numpy()).all();reg=auc.columns.tolist();pairs=[];tr=[]
for scope in ['all_tissues','tumor_only']:
 md=meta if scope=='all_tissues' else meta[meta.source_tissue.str.lower()=='tumor']
 for patient,z in md.groupby('curated_patient_id',observed=True):
  c=z.candidate;inside=z.index[c=='4']
  for ep,other in [('C4_minus_Monocyte',z.index[c=='Monocyte']),('C4_minus_other_macrophages',z.index[~c.isin(['4','1','Monocyte'])])]:
   if len(inside)<20 or len(other)<20:continue
   v=auc.loc[inside].mean()-auc.loc[other].mean()
   for tf in reg:pairs.append(dict(scope=scope,patient=patient,study=z.study.iloc[0],endpoint=ep,regulon=tf,TF=tf.split('(')[0],C4_cells=len(inside),other_cells=len(other),effect=float(v[tf])))
tm=pd.read_csv(R/'07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str,'partition':str},low_memory=False)
for mode in ['primary','without_stress_cycle']:
 mm=tm[tm['mode']==mode]
 for scope in ['all_tissues','tumor_only']:
  m=mm if scope=='all_tissues' else mm[mm.source_tissue.str.lower()=='tumor']
  for pv in (['pt_primary','pt_root_q25','pt_root_q50','dpt_pseudotime'] if mode=='primary' else ['pt_primary']):
   for (study,pat,part),z in m.groupby(['study','curated_patient_id','partition'],observed=True):
    z=z[np.isfinite(z[pv])]
    if len(z)<50 or z[pv].nunique()<10 or not (z.candidate=='4').any():continue
    x=z[pv].to_numpy();depth=np.log1p(z.total_counts_raw.to_numpy());v=auc.loc[z.index]
    for tf in reg:tr.append(dict(mode=mode,scope=scope,pt_version=pv,study=study,patient=pat,partition=part,regulon=tf,TF=tf.split('(')[0],rho=rho(x,v[tf].to_numpy()),depth_adjusted_rho=partial(x,v[tf].to_numpy(),depth)))
 print('REGULON_TRENDS',mode,flush=True)
p=pd.DataFrame(pairs);t=pd.DataFrame(tr);p.to_csv(O/'PATIENT_AUC_CONTRASTS.tsv',sep='\t',index=False);t.to_csv(O/'PARTITION_PATIENT_AUC_TRENDS.tsv.gz',sep='\t',index=False);tg=t.groupby(['mode','scope','pt_version','study','patient','regulon','TF'],as_index=False)[['rho','depth_adjusted_rho']].median();tg.to_csv(O/'PATIENT_AUC_TRENDS.tsv.gz',sep='\t',index=False)
rr=[];tt=[]
for k,g in p.groupby(['scope','endpoint','regulon','TF']):rr.append(dict(zip(['scope','endpoint','regulon','TF'],k))|summ(g.effect))
for k,g in tg.groupby(['mode','scope','pt_version','regulon','TF']):
 for adj in ['rho','depth_adjusted_rho']:tt.append(dict(zip(['mode','scope','pt_version','regulon','TF'],k))|{'adjustment':adj}|summ(g[adj]))
rr=pd.DataFrame(rr);tt=pd.DataFrame(tt)
for df,fam in [(rr,['scope','endpoint']),(tt,['mode','scope','pt_version','adjustment'])]:
 df['FDR_all_regulons']=np.nan;df['FDR_prespecified_3TF']=np.nan
 for _,ids in df.groupby(fam).groups.items():
  df.loc[ids,'FDR_all_regulons']=bh(df.loc[ids,'P']);sel=df.loc[ids].index[df.loc[ids,'TF'].isin(['JUN','JUNB','CEBPD'])];df.loc[sel,'FDR_prespecified_3TF']=bh(df.loc[sel,'P'])
 df['reporting_tier']=np.where(df.n_patients<3,'descriptive',np.where(df.n_patients<10,'small_n_exploratory',np.where(df.FDR_all_regulons<.05,'FDR_supported',np.where(df.P<.05,'nominal_exploratory','not_supported'))))
rr.to_csv(O/'AUC_CONTRAST_SUMMARY.tsv',sep='\t',index=False);tt.to_csv(O/'AUC_TRAJECTORY_SUMMARY.tsv',sep='\t',index=False)
rr[rr.TF.isin(['JUN','JUNB','CEBPD'])].to_csv(O/'PRESPECIFIED_TF_CONTRASTS.tsv',sep='\t',index=False);tt[tt.TF.isin(['JUN','JUNB','CEBPD'])].to_csv(O/'PRESPECIFIED_TF_TRAJECTORY.tsv',sep='\t',index=False)
# Target overlap is descriptive mechanistic correspondence plus separate ORA family, not TF binding proof.
ed=pd.read_csv(N/'REGULON_TARGET_EDGES.tsv',sep='\t');gm=pd.read_csv(R/'04_markers_annotation/run_v4_macrophage_only/full_marker_tables/cluster_4_all_genes.tsv.gz',sep='\t');eligible=set(gm.loc[gm.test_eligible==True,'gene']) & set(pd.read_csv(R/'08_regulons/input_v1/GENES.tsv',sep='\t').gene);up=set(gm.loc[(gm.patient_wilcoxon_p<.05)&(gm.median_patient_log2CPM_difference>=.25),'gene'])&eligible;down=set(gm.loc[(gm.patient_wilcoxon_p<.05)&(gm.median_patient_log2CPM_difference<=-.25),'gene'])&eligible;ov=[]
for tf,e in ed.groupby('regulon'):
 targets=set(e.target)&eligible
 for direction,sel in [('up',up),('down',down)]:
  hits=sel&targets;pv=float(hypergeom.sf(len(hits)-1,len(eligible),len(targets),len(sel))) if targets else np.nan;ov.append(dict(regulon=tf,direction=direction,targets_measured=len(targets),marker_overlap=len(hits),P=pv,overlap_genes=';'.join(sorted(hits))))
ov=pd.DataFrame(ov);ov['FDR']=bh(ov.P);ov.to_csv(O/'REGULON_TARGET_MARKER_OVERLAP.tsv',sep='\t',index=False)
# Leave-one-patient effects for all regulon state contrasts.
los=[]
for k,g in p.groupby(['scope','endpoint','regulon']):
 if len(g)<3:continue
 vals=g.effect.to_numpy();l=np.asarray([np.median(np.delete(vals,i)) for i in range(len(vals))]);los.append(dict(zip(['scope','endpoint','regulon'],k))|dict(n_patients=len(vals),full_median=float(np.median(vals)),LOPO_min_median=float(l.min()),LOPO_max_median=float(l.max()),direction_stable=bool(np.all(np.sign(l)==np.sign(np.median(vals))))))
pd.DataFrame(los).to_csv(O/'LEAVE_ONE_PATIENT_CONTRAST_STABILITY.tsv',sep='\t',index=False)
assert not p.duplicated(['scope','patient','endpoint','regulon']).any();assert not any(f.suffix in ['.png','.svg','.pdf'] for f in O.rglob('*'))
(O/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_REGULON_PATIENT_ANALYSIS','regulons':len(reg),'all_cells':len(auc),'contrasts':len(rr),'trajectory_rows':len(tt),'figures_created':False},indent=2));print('REGULON_STATS_DONE',flush=True)
