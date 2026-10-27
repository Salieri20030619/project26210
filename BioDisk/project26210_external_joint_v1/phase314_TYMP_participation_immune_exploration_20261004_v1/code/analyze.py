from pathlib import Path
import os,json,gc,subprocess,hashlib,sys
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata,binomtest
P=Path(__file__).resolve().parents[1];R=P.parent;S=R/'phase305_larger_immune_sequential_20261002_v1'
pr=json.loads((P/'resources/ANALYSIS_PROTOCOL.json').read_text())
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
meta=pd.read_csv(S/'tables/OBRADOVIC_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0)
meta.index=meta.index.astype(str)
broad=meta.selection.eq('CD45+')&meta.myeloid_cluster.astype(str).isin(['2','3','5','6','8','9','10','15','17'])
strict=meta.lineage_Macrophage.astype(bool)
assert (~strict|broad).all()
a=ad.read_h5ad(S/'inputs/OBRADOVIC_counts_qc.h5ad')
assert a.obs_names.equals(meta.index)
a=a[broad.to_numpy()].copy();a.obs=meta.loc[broad].copy()
a.obs['original_row']=np.flatnonzero(broad.to_numpy());a.obs['L0_strict']=strict[broad].to_numpy()
a.write_h5ad(P/'inputs/OBRADOVIC_macrophage_raw.h5ad',compression='gzip')
a.obs.to_csv(P/'tables/CELL_MEMBERSHIP.tsv.gz',sep='\t')
x=a.X.tocsr().astype(float);umi=np.asarray(x.sum(axis=1)).ravel()
assert np.allclose(umi,a.obs.total_counts)
norm=x.multiply((1e4/umi)[:,None]).tocsr();norm.data=np.log1p(norm.data)
present={}
for lib in a.obs.library.unique():
 b=ad.read_h5ad(S/f'inputs/OBRADOVIC_{lib}_qc.h5ad',backed='r')
 present[lib]=set(b.var_names);b.file.close()
genes=sorted(set(sum(pr['main_pairs'],[])))
ctx=json.loads((S/'resources/PATIENT_ASSOCIATION_PROTOCOL.json').read_text())['context']
ctx={k:[g for g in v if g not in genes] for k,v in ctx.items()}
(P/'resources/CONTEXT_EXCLUDING_MAIN_GENES.json').write_text(json.dumps(ctx,indent=2))
panel_names=list(pr['immune_panels'])
value_names=genes+panel_names
effects=[];coverage=[];pseudobulk=[];members=[];cache_manifest=[]
def design(n,nums,cats):
 cols=[np.ones(n)]
 for v in nums:
  if np.ptp(v)>0:cols.append(rankdata(v))
 for v in cats:cols.extend(pd.get_dummies(pd.Series(v).astype(str),drop_first=True,dtype=float).to_numpy().T)
 return np.column_stack(cols)
def rho_calc(v,w,D):
 rank=np.linalg.matrix_rank(D)
 z=np.column_stack([rankdata(v),rankdata(w)])
 res=z-D@np.linalg.lstsq(D,z,rcond=None)[0]
 if np.std(res[:,0])<1e-9 or np.std(res[:,1])<1e-9:return np.nan,rank
 return float(np.corrcoef(res.T)[0,1]),int(rank)
for condition,definition in pr['memberships'].items():
 mask=np.ones(a.n_obs,dtype=bool) if condition=='L2_broad' else a.obs.L0_strict.to_numpy(dtype=bool)
 for patient in sorted(a.obs.patient.unique()):
  ii=np.flatnonzero(mask&a.obs.patient.eq(patient).to_numpy());n=len(ii);obs=a.obs.iloc[ii]
  available=set.intersection(*[present[l] for l in obs.library.unique()]) if n else set()
  raw=x[ii];y=norm[ii];det=np.asarray((raw>0).sum(axis=0)).ravel()
  summ=np.asarray(raw.sum(axis=0)).ravel();pb=np.log1p(1e4*summ/max(1,summ.sum()))
  vals={};ok={};why={}
  for gene in genes:
   measured=gene in available
   jj=a.var_names.get_loc(gene) if gene in a.var_names else -1
   vals[gene]=y[:,jj].toarray().ravel() if measured else np.full(n,np.nan)
   detected=int(det[jj]) if measured else 0
   ok[gene]=measured and detected>=10
   why[gene]='' if ok[gene] else ('UNMEASURED' if not measured else 'TARGET_DETECTION_LT10')
   coverage.append(dict(condition=condition,patient=patient,feature=gene,kind='gene',cells=n,genes_expected=1,genes_measured=int(measured),detected_genes=int(detected>0),detected_cells=detected,eligible=bool(ok[gene]),reason=why[gene]))
   pseudobulk.append(dict(condition=condition,patient=patient,feature=gene,cells=n,value=float(pb[jj]) if measured else np.nan,detected_cells=detected,eligible=bool(n>=definition['minimum_cells'] and measured and detected>=5)))
  for name,pan in pr['immune_panels'].items():
   full=pan['genes'];keep=[g for g in full if g in available]
   jj=a.var_names.get_indexer(keep);assert (jj>=0).all()
   vals[name]=np.asarray(y[:,jj].mean(axis=1)).ravel() if len(jj) else np.full(n,np.nan)
   seen=int((det[jj]>=max(3,int(np.ceil(.01*n)))).sum())
   ok[name]=len(keep)>=.8*len(full) and seen>=.2*len(full)
   why[name]='' if ok[name] else ('PROGRAM_MEASUREMENT_LT80PCT' if len(keep)<.8*len(full) else 'PROGRAM_DETECTION_LT20PCT')
   coverage.append(dict(condition=condition,patient=patient,feature=name,kind='program',cells=n,genes_expected=len(full),genes_measured=len(keep),detected_genes=seen,detected_cells=int((vals[name]>0).sum()),eligible=bool(ok[name]),reason=why[name]))
  nums=[obs.total_counts.to_numpy(),obs.n_genes.to_numpy()]
  cats=[obs.library.to_numpy(),obs.myeloid_cluster.to_numpy()]
  D=design(n,nums,cats)
  env=[];env_ok=True
  for name,gs in ctx.items():
   keep=[g for g in gs if g in available];env_ok &= len(keep)>=.8*len(gs)
   jj=a.var_names.get_indexer(keep);env.append(np.asarray(y[:,jj].mean(axis=1)).ravel())
  Dc=design(n,nums+env,cats)
  cache=P/'inputs'/f"{condition}__{patient.replace(':','_')}.npz"
  np.savez_compressed(cache,values=np.column_stack([vals[v] for v in value_names]),value_names=np.array(value_names),D_core=D,D_context=Dc,raw_rows=ii,env=np.column_stack(env),available=np.array(sorted(available)))
  cache_manifest.append(dict(condition=condition,patient=patient,path=str(cache)))
  members.append(dict(condition=condition,patient=patient,cells=n,TYMP_positive=int((vals['TYMP']>0).sum()),SPP1_positive=int((vals['SPP1']>0).sum()),libraries=';'.join(obs.library.unique()),core_rank=int(np.linalg.matrix_rank(D)),context_rank=int(np.linalg.matrix_rank(Dc))))
  jobs=[('participation',s,t,model) for s,t in pr['main_pairs'] for model in ['core','context']]
  jobs += [('immune',s,t,'core') for s in ['TYMP','SPP1'] for t in panel_names]
  for branch,source,target,model in jobs:
   DD=D if model=='core' else Dc;dr=int(np.linalg.matrix_rank(DD));reason=[]
   if n<definition['minimum_cells']:reason.append('CELL_COUNT')
   if not ok[source]:reason.append('SOURCE_'+why[source])
   if not ok[target]:reason.append('TARGET_'+why[target])
   if n-dr<8:reason.append('RESIDUAL_DIMENSION')
   if model=='context' and not env_ok:reason.append('CONTEXT_COVERAGE')
   rho=np.nan
   if not reason:
    rho,dr=rho_calc(vals[source],vals[target],DD)
    if not np.isfinite(rho):reason.append('NO_RESIDUAL_VARIATION')
   effects.append(dict(branch=branch,condition=condition,patient=patient,source=source,target=target,model=model,cells=n,design_rank=dr,rho=rho,reason=';'.join(reason)))
 print('WITHIN_DONE',condition,flush=True)
e=pd.DataFrame(effects);e.to_csv(P/'tables/ALL_PATIENT_EFFECTS.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(P/'tables/FEATURE_COVERAGE.tsv',sep='\t',index=False)
pbdf=pd.DataFrame(pseudobulk);pbdf.to_csv(P/'tables/MACRO_PSEUDOBULK.tsv',sep='\t',index=False)
pd.DataFrame(members).to_csv(P/'tables/PATIENT_MEMBERSHIP_COUNTS.tsv',sep='\t',index=False)
(P/'inputs/CACHE_MANIFEST.json').write_text(json.dumps(cache_manifest,indent=2))
summaries=[]
for keys,g in e.groupby(['branch','condition','source','target','model'],sort=False):
 valid=g[g.rho.notna()];assert valid.patient.nunique()==len(valid)
 n=len(valid);pos=int((valid.rho>0).sum());neg=int((valid.rho<0).sum());nn=pos+neg
 summaries.append(dict(zip(['branch','condition','source','target','model'],keys))|dict(patients=n,positive=pos,negative=neg,zero=n-nn,median_rho=float(valid.rho.median()) if n else np.nan,min_rho=float(valid.rho.min()) if n else np.nan,max_rho=float(valid.rho.max()) if n else np.nan,P=float(binomtest(pos,nn,.5).pvalue) if nn>=10 else np.nan,family=keys[0],family_planned_n=24 if keys[0]=='participation' else 224,status='inferential' if nn>=10 else ('descriptive' if n else 'not_evaluable')))
su=pd.DataFrame(summaries);su['FDR']=np.nan
for branch,gg in su.groupby('branch'):
 ids=gg[gg.P.notna()].sort_values('P').index
 pp=su.loc[ids,'P'].to_numpy();m=len(pp)
 su.loc[ids,'FDR']=np.minimum(1,np.minimum.accumulate((pp*m/np.arange(1,m+1))[::-1])[::-1])
 su.loc[gg.index,'family_tested_n']=m
su['nominal']=su.P.lt(.05);su['FDR_significant']=su.FDR.lt(.05)
su.to_csv(P/'tables/ALL_EFFECT_SUMMARIES.tsv',sep='\t',index=False)
su[su.branch.eq('participation')].to_csv(P/'tables/PRIORITY_PARTICIPATION_RESULTS.tsv',sep='\t',index=False)
su[su.branch.eq('immune')&su.target.isin(pr['immune_display_fixed'])].to_csv(P/'tables/IMMUNE_FEATURES_FIXED_VIEW.tsv',sep='\t',index=False)
old=pd.read_csv(S/'tables/WITHIN_MACRO_PATIENT_EFFECTS.tsv',sep='\t')
old=old[old.dataset.eq('OBRADOVIC')&old.model.eq('depth_complexity_library_subcluster')]
new=e[e.condition.eq('L0_strict')&e.source.eq('TYMP')&e.target.eq('SPP1')&e.model.eq('core')]
check=old[['patient','rho']].merge(new[['patient','rho']],on='patient',suffixes=('_old','_new'),validate='one_to_one')
assert len(check)==11 and np.max(np.abs(check.rho_old-check.rho_new))<1e-10
check.to_csv(P/'qc/PHASE305_EXACT_REPRODUCTION.tsv',sep='\t',index=False)
pm=pd.read_csv(S/'tables/PATIENT_MEASURES.tsv',sep='\t')
tum=pm[pm.dataset.eq('OBRADOVIC')&pm.lineage.eq('Tumour_like')&pm.measure.eq('TYMP')&pm.eligible][['patient','value','technical_batch','pT_category']]
cross=[];cross_inputs=[]
for condition in pr['memberships']:
 for target in ['SPP1','CTSB','TGFB1','RAB37']:
  g=tum.merge(pbdf[pbdf.condition.eq(condition)&pbdf.feature.eq(target)&pbdf.eligible],on='patient',suffixes=('_tumour','_macro'),validate='one_to_one');n=len(g)
  for model in ['unadjusted','batch_pT_adjusted']:
   DD=design(n,[],[] if model=='unadjusted' else [g.technical_batch,g.pT_category]);rho=np.nan;dr=int(np.linalg.matrix_rank(DD)) if n else 0
   if n>=3 and n-dr>=2:rho,dr=rho_calc(g.value_tumour,g.value_macro,DD)
   cross.append(dict(condition=condition,source='Tumour_like:TYMP',target='Macrophage:'+target,model=model,patients=n,rho=rho,design_rank=dr,P=np.nan,FDR=np.nan,status='descriptive' if np.isfinite(rho) else 'not_evaluable'))
   for row in g.to_dict('records'):cross_inputs.append(dict(model=model,target=target,**row))
pd.DataFrame(cross).to_csv(P/'tables/TUMOUR_TYMP_MACRO_PROCESS.tsv',sep='\t',index=False)
pd.DataFrame(cross_inputs).to_csv(P/'tables/TUMOUR_MACRO_PATIENT_PAIRS.tsv',sep='\t',index=False)
qc={'status':'ANALYSIS_COMPLETE_PENDING_INDEPENDENT_VALIDATION','time':datetime.now(ZoneInfo('Asia/Taipei')).isoformat(),'cells_broad':a.n_obs,'cells_strict':int(a.obs.L0_strict.sum()),'patients':a.obs.patient.nunique(),'summary_models':len(su),'patient_effects':len(e),'main_models':24,'immune_models':224,'cross_source_models':len(cross),'old_strict_TYMP_SPP1_matched':len(check),'no_cell_P':True,'no_P_for_cross_source_n_less10':True}
(P/'qc/ANALYSIS_COMPLETE.json').write_text(json.dumps(qc,indent=2))
print(json.dumps(qc),flush=True)
print(su[su.branch.eq('participation')].to_string(index=False),flush=True)
print('IMMUNE_NOMINAL',flush=True)
print(su[su.condition.eq('L2_broad')&su.branch.eq('immune')&su.nominal&su.target.isin(pr['immune_display_fixed'])][['source','target','patients','positive','negative','median_rho','P','FDR']].to_string(index=False),flush=True)
