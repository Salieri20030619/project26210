from pathlib import Path
import json,sys,hashlib,shutil
import numpy as np,pandas as pd
from scipy.stats import rankdata,wilcoxon
w=Path(__file__).resolve().parents[1];q=w/'qc'
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1');t=b/'09_metabolism_trajectory'
sources=[]
def tab(p,index=None):
 p=Path(p);sources.append(p);return pd.read_csv(p,sep='\t',index_col=index)
auc=['AUC_Pyrimidine','AUC_CysMet','AUC_ArgPro','AUC_Glycerophospholipid']
flux=['Pyrimidine_flux','CysMet_flux','ArgPro_flux','Choline_proxy']
genes=['TYMP','TK1','SLC29A1','SLC29A2','ODC1','SRM','SMS','MTAP','CHKA','PCYT1A','PCYT2','CEPT1']
methods=['Monocle3','Slingshot','scTour'];tcols=[x+'_pseudotime_patient_v1' for x in methods]
m=tab(t/'trajectory_v1/input/cell_metadata_v1.tsv.gz',0)
u=tab(t/'trajectory_v1/input/scvi_umap_v1.tsv.gz',0).loc[m.index]
a=tab(t/'metabolism_v1/source_data/four_pathway_AUCell_cell_scores_v1.tsv.gz',0).loc[m.index]
f=tab(t/'metabolism_v1/source_data/scFEA_four_axis_cell_scores_v1.tsv.gz',0).loc[m.index]
tr=tab(t/'trajectory_v1/source_data/trajectory_all_methods_cell_scores_v1.tsv.gz',0).loc[m.index]
e=tab(w/'data/ENZYME_LOG1P_CP10K.tsv.gz',0).loc[m.index]
eraw=tab(w/'data/ENZYME_RAW_COUNTS.tsv.gz',0).loc[m.index]
assert len(m)==11950 and m.index.is_unique and m.curated_patient_id.nunique()==26
assert m.study.nunique()==5 and np.allclose(m[auc],a[auc]) and np.allclose(m[flux],f[flux])
assert np.array_equal(a.curated_patient_id,m.curated_patient_id)
assert np.array_equal(tr.curated_patient_id,m.curated_patient_id)
assert np.allclose(e.TYMP,m.TYMP_log1p_CP10K_v1,atol=2e-6)
assert np.array_equal(a.TYMP_group=='TYMP_positive',eraw.TYMP>0)
for c in tcols:assert np.isfinite(tr[c]).all() and tr[c].between(0,1).all()
fates=['CellRank_fate_'+x+'_v1' for x in ['Cycle','Hypoxia','TYMP_state']]
assert np.allclose(tr[fates].sum(axis=1),1,atol=1e-8)
assert tr[fates].ge(0).all().all() and tr[fates].le(1).all().all()
patient_ids=[p for p,g in m.groupby('curated_patient_id') if len(g)>=30]
assert len(patient_ids)==20
eligible=m.curated_patient_id.isin(patient_ids)
checks=[]
def rho(x,y):
 ok=np.isfinite(x)&np.isfinite(y);x=np.asarray(x)[ok];y=np.asarray(y)[ok]
 if len(x)<30 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
 return float(np.corrcoef(rankdata(x),rankdata(y))[0,1])
def bh(p):
 p=np.asarray(p,float);order=np.argsort(p);res=np.empty(len(p));res[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return res

def wp(x,approx=True):
 x=np.asarray(x);x=x[np.isfinite(x)]
 if len(x)<2:return np.nan
 if np.all(np.abs(x)<1e-12):return 1.
 return float(wilcoxon(x,method='approx' if approx else 'auto',correction=approx).pvalue)
# All original patient deltas, correlations and BH families are independently checked.
summary={};details={};cors={};cor_summary={}
for key,frame,prefix in [('AUCell',a,'four_pathway_AUCell'),('scFEA',f,'scFEA_four_axis')]:
 folder=t/'metabolism_v1/tables'
 ds=tab(folder/(prefix+'_patient_deltas_v1.tsv'));ss=tab(folder/(prefix+'_paired_statistics_v1.tsv'))
 cs=tab(folder/(prefix+'_correlations_by_patient_v1.tsv'));su=tab(folder/(prefix+'_correlation_summary_v1.tsv'))
 for r in ds.itertuples():
  sub=frame[frame.curated_patient_id==r.curated_patient_id]
  if r.contrast=='TYMP-positive_vs_negative':col,lo,hi='TYMP_group','TYMP_negative','TYMP_positive'
  else:col,lo,hi='clone_group','other_primary_clones','F1_high_primary_clone'
  low=sub.loc[sub[col]==lo,r.pathway];high=sub.loc[sub[col]==hi,r.pathway]
  assert len(low)>=20 and len(high)>=20
  assert np.allclose([low.median(),high.median(),high.median()-low.median()],[r.low_score,r.high_score,r.delta],atol=1e-11)
 for r in cs.itertuples():
  sub=frame[frame.curated_patient_id==r.curated_patient_id]
  rr=rho(sub[r.pathway].values,sub[r.predictor].values)
  assert np.isclose(rr,r.rho,equal_nan=True,atol=1e-10)
 for r in ss.itertuples():
  d=ds[(ds.pathway==r.pathway)&(ds.contrast==r.contrast)].delta.values
  assert len(d)==r.n_patients
  assert np.allclose([np.median(d),np.mean(d>0),wp(d)],[r.median_patient_delta,r.fraction_patients_positive,r.paired_wilcox_PValue],atol=1e-10)
 for r in su.itertuples():
  d=cs[(cs.pathway==r.pathway)&(cs.predictor==r.predictor)].rho.dropna().values
  assert len(d)==r.n_patients
  assert np.allclose([np.median(d),np.mean(d>0),wp(d)],[r.median_patient_rho,r.fraction_positive,r.signed_rank_PValue],atol=1e-10)
 for _,g in ss.groupby('contrast'):assert np.allclose(bh(g.paired_wilcox_PValue),g.FDR,atol=1e-10)
 for _,g in su.groupby('predictor'):assert np.allclose(bh(g.signed_rank_PValue),g.FDR,atol=1e-10)
 summary[key]=ss;details[key]=ds;cors[key]=cs;cor_summary[key]=su
 checks.append({'analysis':key,'patient_deltas_checked':len(ds),'patient_correlations_checked':len(cs),'summary_P_and_BH':'PASS'})
# Rebuild scFEA standardized module aggregates, without refitting the model.
mod=tab(t/'metabolism_v1/source_data/scFEA_M168_cell_flux.tsv.gz',0).loc[m.index]
reg=tab(t/'metabolism_v1/source_data/scFEA_four_axis_module_definitions_v1.tsv')
z=(mod-mod.mean())/mod.std(ddof=1)
constant=mod.max()==mod.min();z.loc[:,constant]=0
flux_errors={}
for feature,sub in reg.groupby('pathway'):
 rebuilt=z[sub.module.tolist()].mean(axis=1);err=float(np.max(np.abs(rebuilt-f[feature])))
 assert err<1e-9;flux_errors[feature]=err
# Recheck all 48 original trajectory summary endpoints with their 12-test/method families.
tc=tab(t/'trajectory_v1/tables/trajectory_all_methods_patient_correlations_v1.tsv')
ts=tab(t/'trajectory_v1/tables/trajectory_all_methods_correlation_summary_v1.tsv')
for r in tc.itertuples():
 idx=m.index[m.curated_patient_id==r.curated_patient_id]
 rr=rho(tr.loc[idx,r.method].values,m.loc[idx,r.target].values)
 assert np.isclose(rr,r.rho,equal_nan=True,atol=1e-9),(r.method,r.target,rr,r.rho)
for r in ts.itertuples():
 d=tc[(tc.method==r.method)&(tc.target==r.target)].rho.dropna().values
 assert len(d)==r.n_patients
 assert np.allclose([np.median(d),np.mean(d>0),wp(d,False)],[r.median_patient_rho,r.fraction_positive,r.signed_rank_PValue],atol=1e-10)
for _,g in ts.groupby('method'):assert np.allclose(bh(g.signed_rank_PValue),g.FDR,atol=1e-10)
# New families fixed before comparison: 85 KEGG pathways, 12 named enzymes vs F1.
kegg=tab(t/'metabolism_v1/source_data/KEGG_AUCell_all_pathways_v1.tsv.gz',0).loc[m.index]
assert kegg.shape==(11950,85)
rows=[];all_delta=[]
for pathway in kegg.columns:
 ds=[]
 for pid in sorted(m.curated_patient_id.unique()):
  ix=m.index[m.curated_patient_id==pid];pos=a.loc[ix,'TYMP_group']=='TYMP_positive'
  if pos.sum()<20 or (~pos).sum()<20:continue
  lo=float(kegg.loc[ix[~pos],pathway].median());hi=float(kegg.loc[ix[pos],pathway].median());ds.append(hi-lo)
  all_delta.append({'patient':pid,'pathway':pathway,'negative_median':lo,'positive_median':hi,'delta':hi-lo})
 rows.append({'pathway':pathway,'n_patients':len(ds),'median_delta':np.median(ds),'positive_patients':int(np.sum(np.array(ds)>0)),'P':wp(ds)})
ks=pd.DataFrame(rows);ks['FDR']=bh(ks.P);ks.to_csv(q/'ALL85_TYMP_PATIENT_STATISTICS.tsv',sep='\t',index=False)
pd.DataFrame(all_delta).to_csv(q/'ALL85_TYMP_PATIENT_DELTAS.tsv',sep='\t',index=False)
# Four-pathway prespecified and 85-pathway exploratory FDRs stay distinct.
kgmap=dict(zip(auc,['Pyrimidine metabolism','Cysteine and methionine metabolism','Arginine and proline metabolism','Glycerophospholipid metabolism']))
for c,name in kgmap.items():assert np.allclose(kegg[name],a[c],atol=1e-12)
ge=[];gsummary=[]
for gene in genes:
 vals=[]
 for pid in patient_ids:
  ix=m.index[m.curated_patient_id==pid];r=rho(m.loc[ix,'F1'].values,e.loc[ix,gene].values)
  ge.append({'patient':pid,'gene':gene,'rho':r,'n_cells':len(ix)})
  if np.isfinite(r):vals.append(r)
 gsummary.append({'gene':gene,'n_patients':len(vals),'median_rho':np.median(vals),'positive_patients':int(np.sum(np.array(vals)>0)),'P':wp(vals)})
gs=pd.DataFrame(gsummary);gs['FDR']=bh(gs.P)
gs.to_csv(q/'ENZYME_F1_STATISTICS.tsv',sep='\t',index=False);pd.DataFrame(ge).to_csv(q/'ENZYME_F1_PATIENT_EFFECTS.tsv',sep='\t',index=False)
# Tied F1 values remain together; require every quartile to contain >=5 cells.
dot=[];inclusion=[]
for pid in patient_ids:
 ix=m.index[m.curated_patient_id==pid];bins=np.clip(np.ceil(rankdata(m.loc[ix,'F1'])/len(ix)*4),1,4).astype(int)
 counts=np.bincount(bins,minlength=5)[1:];keep=bool((counts>=5).all())
 inclusion.append({'patient':pid,'Q1_n':int(counts[0]),'Q2_n':int(counts[1]),'Q3_n':int(counts[2]),'Q4_n':int(counts[3]),'included':keep})
 if not keep:continue
 for k in range(1,5):
  ids=ix[bins==k]
  for gene in genes:dot.append({'patient':pid,'quartile':k,'gene':gene,'mean_log1p':float(e.loc[ids,gene].mean()),'fraction_detected':float((eraw.loc[ids,gene]>0).mean()),'n_cells':len(ids)})
dot=pd.DataFrame(dot);dot.to_csv(q/'ENZYME_DOTPLOT_PATIENT_PROFILES.tsv',sep='\t',index=False)
pd.DataFrame(inclusion).to_csv(q/'ENZYME_QUARTILE_INCLUSION.tsv',sep='\t',index=False)
# Twenty patient-rank bins, fixed Gaussian bandwidth 0.075, fitted within each patient.
features=auc+flux+genes;curve_rows=[];fit_rows=[];state_rows=[];grid=np.linspace(.025,.975,100)
for method,col in zip(methods,tcols):
 for pid in patient_ids:
  ix=m.index[m.curated_patient_id==pid];ranks=rankdata(tr.loc[ix,col])/len(ix)
  ib=np.clip(np.ceil(ranks*20),1,20).astype(int);ub=np.unique(ib);nn=np.array([(ib==k).sum() for k in ub]);x=(ub-.5)/20
  weights=np.exp(-.5*((grid[:,None]-x[None,:])/.075)**2)*nn[None,:];weights/=weights.sum(axis=1,keepdims=True)
  for feature in features:
   v=(e if feature in genes else m).loc[ix,feature].values
   means=np.array([v[ib==k].mean() for k in ub]);fit=weights@means
   for k,n,mu in zip(ub,nn,means):curve_rows.append({'method':method,'patient':pid,'feature':feature,'bin':int(k),'n':int(n),'mean':float(mu)})
   for x0,v0 in zip(grid,fit):fit_rows.append({'method':method,'patient':pid,'feature':feature,'rank':float(x0),'fitted':float(v0)})
  for state in ['Base','Hypoxia','Cycle','TYMP_state']:
   v=(tr.loc[ix,'malignant_state_v1']==state).astype(float).values
   fit=weights@np.array([v[ib==k].mean() for k in ub])
   for x0,v0 in zip(grid,fit):state_rows.append({'method':method,'patient':pid,'state':state,'rank':float(x0),'fraction':float(v0)})
curves=pd.DataFrame(fit_rows);curves.to_csv(w/'data/PATIENT_FITTED_CURVES.tsv.gz',sep='\t',index=False)
pd.DataFrame(curve_rows).to_csv(w/'data/PATIENT_RANK_BINS.tsv.gz',sep='\t',index=False)
states=pd.DataFrame(state_rows);assert np.allclose(states.groupby(['method','patient','rank']).fraction.sum(),1)
states.to_csv(w/'data/PATIENT_STATE_COMPOSITION.tsv.gz',sep='\t',index=False)
for path in set(sources):
 if path.parent.name in ['tables','qc']:shutil.copy2(path,q/('original_'+path.name))
validation={'status':'PASS','cells':len(m),'patients':26,'eligible_correlation_patients':20,'studies':5,'original_metabolism':checks,'scFEA_aggregate_max_errors':flux_errors,'constant_scFEA_modules':list(mod.columns[constant]),'trajectory_patient_rows_checked':len(tc),'trajectory_summary_rows_checked':len(ts),'new_families':{'KEGG_TYMP_paired':85,'enzyme_F1_correlation':12},'dotplot_patients':int(dot.patient.nunique()),'all85_FDR_significant':int((ks.FDR<.05).sum()),'smoothing':'20 patient-rank bins; Gaussian h=0.075; patient median and IQR','fate_probabilities_sum_to_one':True}
(q/'INDEPENDENT_VALIDATION.json').write_text(json.dumps(validation,indent=2))
(q/'SOURCE_PATHS.json').write_text(json.dumps([str(p) for p in sorted(set(sources))],indent=2))
print(json.dumps(validation),flush=True)
