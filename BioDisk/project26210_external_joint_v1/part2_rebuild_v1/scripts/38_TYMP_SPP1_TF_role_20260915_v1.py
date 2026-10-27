from pathlib import Path
import json,datetime,hashlib
import pandas as pd,numpy as np
from scipy.stats import spearmanr,rankdata,wilcoxon,false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase151_TYMP_SPP1_threeTF_role_20260915_v1';O.mkdir(exist_ok=False)
tf=['JUN(+)','JUNB(+)','CEBPD(+)']
M=R/'08_regulons/input_v1/CELL_METADATA.tsv.gz';A=R/'08_regulons/network_v2_20260915/AUC_ALL_CELLS.csv';G=R/'07_trajectory/dynamics_v2/CELL_EXPRESSION_PROGRAM_SCORES.tsv.gz';T=R/'07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz';ED=R/'08_regulons/network_v2_20260915/REGULON_TARGET_EDGES.tsv';B=P/'phase109_targeted_gap_analysis_local_v1/analysis_v3/phase109_patient_summary_v1.tsv'
plan=dict(question='Do JUN/JUNB/CEBPD specifically link tumour TYMP to macrophage SPP1?',hypotheses=['motif-pruned TF target membership for TYMP/SPP1','within-mature-macrophage TF activity versus TYMP or SPP1 expression','paired tumour TYMP versus macrophage TF activity, SPP1 expression and C4 fraction'],scope='All mature macrophage candidates, C1 and Monocyte excluded; tumour-only primary with all-tissue within-mac sensitivity.',unit='Independent patient; cell rho aggregated to patient median; no cell-level significance.',thresholds='Within-patient >=50 macrophages and >=10 target-gene-positive cells; same-partition trajectory-adjusted additionally >=10 distinct pseudotimes. Tumour bridge retains frozen eligible_abundance_primary and >=50 tumour macrophages.',FDR='Within each scope/gene/adjustment threeTF; tumour bridge each frozen exposure by 5 readouts, raw and study-adjusted separate; no mixed families.',causal_claims=False,new_model_fits=False,figures=False,time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'PLAN.json').write_text(json.dumps(plan,indent=2))
m=pd.read_csv(M,sep='\t',index_col=0,dtype={'candidate':str},low_memory=False);h=pd.read_csv(A,nrows=0)
a=pd.read_csv(A,index_col=0,usecols=[h.columns[0]]+tf);g=pd.read_csv(G,sep='\t',index_col=0,usecols=['cell_id','gene:SPP1','gene:TYMP'])
assert set(m.index)==set(a.index)==set(g.index)
ed=pd.read_csv(ED,sep='\t');er=[]
for reg in tf:
 for target in ['TYMP','SPP1']:
  hit=ed[(ed.regulon==reg)&(ed.target==target)]
  er.append(dict(regulon=reg,target=target,present=len(hit)>0,weight=float(hit.weight.iloc[0]) if len(hit) else np.nan))
pd.DataFrame(er).to_csv(O/'THREE_TF_TARGET_MEMBERSHIP.tsv',sep='\t',index=False)
def rho(x,y):
 return float(spearmanr(x,y).statistic) if np.std(x)>0 and np.std(y)>0 else np.nan
def partial(x,y,cols):
 cov=np.column_stack([np.ones(len(x))]+[rankdata(v) for v in cols]);xr=rankdata(x);yr=rankdata(y)
 xr=xr-cov@np.linalg.lstsq(cov,xr,rcond=None)[0];yr=yr-cov@np.linalg.lstsq(cov,yr,rcond=None)[0]
 return float(np.corrcoef(xr,yr)[0,1]) if np.std(xr)>1e-10 and np.std(yr)>1e-10 else np.nan
rows=[]
mac=m[~m.candidate.isin(['1','Monocyte'])]
for scope in ['tumor_only','all_tissues']:
 mm=mac[mac.source_tissue.str.lower()=='tumor'] if scope=='tumor_only' else mac
 for patient,z in mm.groupby('curated_patient_id'):
  if len(z)<50:continue
  for gene in ['SPP1','TYMP']:
   y=g.loc[z.index,'gene:'+gene].to_numpy()
   if (y>0).sum()<10 or np.std(y)==0:continue
   dep=np.log1p(z.total_counts_raw.to_numpy())
   for reg in tf:
    x=a.loc[z.index,reg].to_numpy()
    for adj,r in [('raw',rho(x,y)),('depth_adjusted',partial(x,y,[dep]))]:
     rows.append(dict(scope=scope,patient=patient,study=z.study.iloc[0],gene=gene,regulon=reg,adjustment=adj,cells=len(z),rho=r))
pt=pd.DataFrame(rows);pt.to_csv(O/'PATIENT_TF_GENE_CORRELATIONS.tsv',sep='\t',index=False)
tm=pd.read_csv(T,sep='\t',index_col=0,dtype={'candidate':str,'partition':str},low_memory=False);tm=tm[(tm['mode']=='primary')&(~tm.candidate.isin(['1','Monocyte']))]
pr=[]
for scope in ['tumor_only','all_tissues']:
 mm=tm[tm.source_tissue.str.lower()=='tumor'] if scope=='tumor_only' else tm
 for (study,patient,part),z in mm.groupby(['study','curated_patient_id','partition']):
  z=z[np.isfinite(z.pt_primary)]
  if len(z)<50 or z.pt_primary.nunique()<10:continue
  for gene in ['SPP1','TYMP']:
   y=g.loc[z.index,'gene:'+gene].to_numpy()
   if (y>0).sum()<10 or np.std(y)==0:continue
   for reg in tf:pr.append(dict(scope=scope,patient=patient,study=study,partition=part,gene=gene,regulon=reg,adjustment='depth_pseudotime_adjusted',cells=len(z),rho=partial(a.loc[z.index,reg].to_numpy(),y,[np.log1p(z.total_counts_raw.to_numpy()),z.pt_primary.to_numpy()])))
part=pd.DataFrame(pr);part.to_csv(O/'PARTITION_ADJUSTED_TF_GENE_CORRELATIONS.tsv',sep='\t',index=False)
agg=part.groupby(['scope','patient','study','gene','regulon','adjustment'],as_index=False).rho.median();combined=pd.concat([pt.drop(columns='cells'),agg],ignore_index=True)
combined.to_csv(O/'ALL_PATIENT_TF_GENE_CORRELATIONS.tsv',sep='\t',index=False)
summ=[]
for key,z in combined.groupby(['scope','gene','regulon','adjustment']):
 vals=z.rho.dropna().to_numpy();n=len(vals);pv=np.nan if n<3 else wilcoxon(vals,method='approx',correction=True).pvalue if np.any(vals) else 1.
 summ.append(dict(zip(['scope','gene','regulon','adjustment'],key))|dict(n_patients=n,median_rho=float(np.median(vals)) if n else np.nan,positive=int((vals>0).sum()),negative=int((vals<0).sum()),P=pv))
summ=pd.DataFrame(summ);summ['FDR_threeTF']=np.nan
for _,idx in summ.groupby(['scope','gene','adjustment']).groups.items():
 v=summ.loc[idx,'P'].dropna()
 if len(v):summ.loc[v.index,'FDR_threeTF']=false_discovery_control(v.to_numpy())
summ.to_csv(O/'TF_GENE_SUMMARY.tsv',sep='\t',index=False)
# Frozen malignant-side aggregate versus current tumour macrophages, patient matched.
bridge=pd.read_csv(B,sep='\t');br=[]
for pat,z in mac[mac.source_tissue.str.lower()=='tumor'].groupby('curated_patient_id'):
 if len(z)<50:continue
 row=dict(patient=pat,study=z.study.iloc[0],macrophage_cells=len(z),C4_fraction=float((z.candidate=='4').mean()),SPP1_RNA_mean=float(g.loc[z.index,'gene:SPP1'].mean()))
 for reg in tf:row[reg]=float(a.loc[z.index,reg].mean())
 br.append(row)
br=pd.DataFrame(br).merge(bridge[['patient','study_v1','TYMP_axis_primary','z_TYMP_RNA','malignant_axis_cells','eligible_abundance_primary']],on='patient',validate='one_to_one')
br=br[br.eligible_abundance_primary==True].copy();assert (br.study==br.study_v1).all()
br.to_csv(O/'TUMOUR_TYMP_MACROPHAGE_MATCHED_PATIENTS.tsv',sep='\t',index=False)
outs=[];rng=np.random.default_rng(262151);nperm=9999;responses=tf+['SPP1_RNA_mean','C4_fraction']
for exposure in ['z_TYMP_RNA','TYMP_axis_primary']:
 z=br.dropna(subset=[exposure]+responses).reset_index(drop=True);n=len(z);assert n>=10
 x=rankdata(z[exposure]);Y=np.column_stack([rankdata(z[r]) for r in responses])
 groups=[np.asarray(ids) for ids in z.groupby('study').groups.values()]
 xr=x.copy();yr=Y.copy()
 for ids in groups:xr[ids]-=xr[ids].mean();yr[ids]-=yr[ids].mean(axis=0)
 den=np.sqrt((xr*xr).sum()*((yr*yr).sum(axis=0)));obs=xr@yr/den
 hits=np.zeros(len(responses),int)
 for _ in range(nperm):
  perm=xr.copy()
  for ids in groups:perm[ids]=rng.permutation(perm[ids])
  r=perm@yr/den;hits+=(np.abs(r)>=np.abs(obs)-1e-14)
 adjp=(hits+1)/(nperm+1)
 for j,target in enumerate(responses):
  r,pv=spearmanr(z[exposure],z[target])
  outs.append(dict(exposure=exposure,target=target,adjustment='raw',n_patients=n,effect_rho=r,P=pv))
  outs.append(dict(exposure=exposure,target=target,adjustment='study_adjusted_rank_within_study_permutation',n_patients=n,effect_rho=obs[j],P=adjp[j]))
out=pd.DataFrame(outs);out['FDR_five_readouts']=np.nan
for _,idx in out.groupby(['exposure','adjustment']).groups.items():out.loc[idx,'FDR_five_readouts']=false_discovery_control(out.loc[idx,'P'].to_numpy())
out.to_csv(O/'TUMOUR_TYMP_BRIDGE_SUMMARY.tsv',sep='\t',index=False)
assert not combined.duplicated(['scope','patient','gene','regulon','adjustment']).any()
src=[]
for path in [M,A,G,T,ED,B]:
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 src.append(dict(path=str(path),bytes=path.stat().st_size,sha256=h.hexdigest()))
pd.DataFrame(src).to_csv(O/'SOURCE_MANIFEST.tsv',sep='\t',index=False)
(O/'STATUS.json').write_text(json.dumps(dict(status='COMPUTED_AWAITING_REVIEW',figures_created=False,model_refit=False,time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2))
print('TARGETS');print(pd.DataFrame(er).to_string(index=False))
print('MACRO_TUMOR');print(summ[summ.scope=='tumor_only'].to_string(index=False))
print('BRIDGE');print(out.to_string(index=False))

