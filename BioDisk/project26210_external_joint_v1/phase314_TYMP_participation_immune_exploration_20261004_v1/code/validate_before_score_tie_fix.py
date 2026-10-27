from pathlib import Path
import json,hashlib,math,gc,subprocess
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import false_discovery_control
from scipy.linalg import orth
P=Path(__file__).resolve().parents[1];S=P.parent/'phase305_larger_immune_sequential_20261002_v1'
assert (P/'logs/analyze.exit').read_text().strip()=='0'
pr=json.loads((P/'resources/ANALYSIS_PROTOCOL.json').read_text());checks=[]
def check(name,cond,detail=''):
 checks.append(dict(check=name,pass_=bool(cond),detail=str(detail)))
 if not cond:raise AssertionError(name+': '+str(detail))
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
for rec in json.loads((P/'resources/INPUT_MANIFEST.json').read_text()):
 h=hashlib.sha256()
 with open(rec['path'],'rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 check('source_hash_'+Path(rec['path']).name,h.hexdigest()==rec['sha256'])
a=ad.read_h5ad(P/'inputs/OBRADOVIC_macrophage_raw.h5ad')
b=ad.read_h5ad(S/'inputs/OBRADOVIC_counts_qc.h5ad')
meta=pd.read_csv(S/'tables/OBRADOVIC_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0)
mask=meta.selection.eq('CD45+')&meta.myeloid_cluster.astype(str).isin(['2','3','5','6','8','9','10','15','17'])
expected=np.flatnonzero(mask)
check('exact_broad_identity',np.array_equal(expected,a.obs.original_row))
check('exact_broad_barcodes',np.array_equal(meta.index[expected],a.obs_names))
check('exact_strict_flags',np.array_equal(meta.iloc[expected].lineage_Macrophage.astype(bool),a.obs.L0_strict))
check('raw_gene_order',a.var_names.equals(b.var_names))
check('all_selected_raw_counts', (a.X-b.X[expected]).nnz==0)
del b;gc.collect()
x=a.X.tocsr().astype(float);umi=np.asarray(x.sum(axis=1)).ravel()
check('all_raw_depths',np.array_equal(umi,a.obs.total_counts))
check('original_feature_complexity_preserved',np.array_equal(meta.iloc[expected].n_genes.to_numpy(),a.obs.n_genes.to_numpy()))
check('premerge_postmerge_complexity_resolved',json.loads((P/'qc/COMPLEXITY_DEFINITION_VALIDATED.json').read_text())['status']=='PASS')
library_sets={}
for lib in a.obs.library.unique():
 z=ad.read_h5ad(S/f'inputs/OBRADOVIC_{lib}_qc.h5ad',backed='r');library_sets[lib]=set(z.var_names);z.file.close()
ef=pd.read_csv(P/'tables/ALL_PATIENT_EFFECTS.tsv',sep='\t')
cov=pd.read_csv(P/'tables/FEATURE_COVERAGE.tsv',sep='\t')
pb=pd.read_csv(P/'tables/MACRO_PSEUDOBULK.tsv',sep='\t')
ctx=json.loads((P/'resources/CONTEXT_EXCLUDING_MAIN_GENES.json').read_text())
gene_names=sorted(set(sum(pr['main_pairs'],[])))
def ranks(v):return pd.Series(v).rank(method='average').to_numpy()
def independent_design(obs,extra=()):
 cols=[np.ones(len(obs)),ranks(obs.total_counts),ranks(obs.n_genes)]
 for v in extra:
  if np.ptp(v)>0:cols.append(ranks(v))
 for column in ['library','myeloid_cluster']:
  vals=obs[column].astype(str).to_numpy()
  for level in sorted(set(vals))[1:]:cols.append((vals==level).astype(float))
 return np.column_stack(cols)
def icorr(v,w,D):
 U=orth(D)
 vr=ranks(v);wr=ranks(w)
 vr=vr-U@(U.T@vr);wr=wr-U@(U.T@wr)
 if np.std(vr)<1e-9 or np.std(wr)<1e-9:return np.nan
 return float(np.dot(vr,wr)/np.sqrt(np.dot(vr,vr)*np.dot(wr,wr)))
maxerr=0.
for rec in json.loads((P/'inputs/CACHE_MANIFEST.json').read_text()):
 cc=rec['condition'];pt=rec['patient'];cache=np.load(rec['path']);rows=cache['raw_rows'];obs=a.obs.iloc[rows]
 exp=a.obs.patient.eq(pt).to_numpy()&(np.ones(a.n_obs,dtype=bool) if cc=='L2_broad' else a.obs.L0_strict.to_numpy(dtype=bool))
 check('exact_patient_cells_'+cc+pt,np.array_equal(rows,np.flatnonzero(exp)))
 avail=set.intersection(*[library_sets[l] for l in obs.library.unique()])
 check('available_gene_set_'+cc+pt,avail==set(cache['available']))
 raw=x[rows];libsum=np.asarray(raw.sum(axis=1)).ravel()
 y=raw.copy();y.data=np.log1p(y.data*np.repeat(1e4/libsum,np.diff(y.indptr)))
 det=np.asarray((raw>0).sum(axis=0)).ravel();summ=np.asarray(raw.sum(axis=0)).ravel();pv=np.log1p(summ*1e4/summ.sum())
 values={};eligibility={}
 for name in cache['value_names']:
  gs=[name] if name in gene_names else pr['immune_panels'][name]['genes']
  check('target_not_in_program_'+name, name in gene_names or not (set(gs)&{'TYMP','SPP1'}))
  keep=[g for g in gs if g in avail];jj=a.var_names.get_indexer(keep)
  values[name]=np.asarray(y[:,jj].sum(axis=1)).ravel()/len(jj) if len(jj) else np.full(len(rows),np.nan)
  j=list(cache['value_names']).index(name)
  check('score_'+cc+pt+name,np.allclose(values[name],cache['values'][:,j],rtol=1e-10,atol=1e-10,equal_nan=True))
  crow=cov[cov.condition.eq(cc)&cov.patient.eq(pt)&cov.feature.eq(name)].iloc[0]
  seen=int(det[jj[0]]) if name in gene_names and len(jj) else (int((det[jj]>=max(3,math.ceil(.01*len(rows)))).sum()) if name not in gene_names else 0)
  ok=(len(keep)==1 and seen>=10) if name in gene_names else (len(keep)>=.8*len(gs) and seen>=.2*len(gs))
  eligibility[name]=ok
  check('coverage_'+cc+pt+name,bool(crow.eligible)==ok and crow.genes_measured==len(keep) and (crow.detected_cells==seen if name in gene_names else crow.detected_genes==seen))
  if name in gene_names:
   pp=pb[pb.condition.eq(cc)&pb.patient.eq(pt)&pb.feature.eq(name)].iloc[0]
   check('pseudobulk_'+cc+pt+name,np.isclose(pp.value,pv[jj[0]],rtol=1e-12,atol=1e-12) if len(jj) else pd.isna(pp.value))
 env=[]
 for gs in ctx.values():
  jj=a.var_names.get_indexer([g for g in gs if g in avail])
  env.append(np.asarray(y[:,jj].sum(axis=1)).ravel()/len(jj))
 env=np.column_stack(env)
 check('context_'+cc+pt,np.allclose(env,cache['env'],atol=1e-10))
 for model in ['core','context']:
  D=independent_design(obs,extra=env.T if model=='context' else ())
  U=orth(D);stored=cache['D_'+model]
  check('design_space_'+cc+pt+model,np.linalg.norm(stored-U@(U.T@stored))<1e-6 and np.linalg.matrix_rank(stored)==U.shape[1])
  selected=ef[ef.condition.eq(cc)&ef.patient.eq(pt)&ef.model.eq(model)]
  for _,row in selected.iterrows():
   valid=len(rows)>=pr['memberships'][cc]['minimum_cells'] and eligibility[row.source] and eligibility[row.target] and len(rows)-U.shape[1]>=8
   if model=='context':valid &= all(len(set(gs)&avail)>=.8*len(gs) for gs in ctx.values())
   val=icorr(values[row.source],values[row.target],D) if valid else np.nan
   good=(np.isnan(val) and pd.isna(row.rho)) or np.isclose(val,row.rho,rtol=1e-7,atol=1e-9)
   check('effect_'+cc+pt+model+row.source+row.target,good)
   if np.isfinite(val) and pd.notna(row.rho):maxerr=max(maxerr,abs(val-row.rho))
su=pd.read_csv(P/'tables/ALL_EFFECT_SUMMARIES.tsv',sep='\t')
check('planned_model_count',len(su)==248 and len(ef)==248*11)
for _,row in su.iterrows():
 g=ef[ef.branch.eq(row.branch)&ef.condition.eq(row.condition)&ef.source.eq(row.source)&ef.target.eq(row.target)&ef.model.eq(row.model)].dropna(subset=['rho'])
 check('patient_unique_'+str(row.name),len(g)==g.patient.nunique())
 pos=int((g.rho>0).sum());neg=int((g.rho<0).sum());n=pos+neg
 pp=min(1.,2*sum(math.comb(n,k) for k in range(min(pos,neg)+1))/2**n) if n>=10 else np.nan
 check('summary_'+str(row.name),row.patients==len(g) and row.positive==pos and row.negative==neg and (np.isclose(row.median_rho,g.rho.median()) if len(g) else pd.isna(row.median_rho)) and (np.isclose(row.P,pp) if np.isfinite(pp) else pd.isna(row.P)))
for branch,g in su.groupby('branch'):
 v=g[g.P.notna()]
 q=false_discovery_control(v.P.to_numpy(),method='bh')
 check('BH_'+branch,np.allclose(q,v.FDR,rtol=1e-12,atol=1e-12) and (g.family_tested_n==len(v)).all())
cross=pd.read_csv(P/'tables/TUMOUR_TYMP_MACRO_PROCESS.tsv',sep='\t')
ci=pd.read_csv(P/'tables/TUMOUR_MACRO_PATIENT_PAIRS.tsv',sep='\t')
check('cross_source_no_inference',(cross.patients<10).all() and cross.P.isna().all() and cross.FDR.isna().all())
for _,row in cross.iterrows():
 g=ci[ci.condition.eq(row.condition)&ci.model.eq(row.model)&ci.target.eq(row.target.split(':')[1])]
 cols=[np.ones(len(g))]
 if row.model=='batch_pT_adjusted':
  for cn in ['technical_batch','pT_category']:
   vv=g[cn].astype(str).to_numpy()
   for level in sorted(set(vv))[1:]:cols.append((vv==level).astype(float))
 D=np.column_stack(cols);val=icorr(g.value_tumour,g.value_macro,D) if len(g)>=3 and len(g)-np.linalg.matrix_rank(D)>=2 else np.nan
 check('cross_source_'+str(row.name),row.patients==g.patient.nunique() and (np.isclose(val,row.rho,atol=1e-9) if np.isfinite(val) else pd.isna(row.rho)))
pd.DataFrame(checks).to_csv(P/'qc/INDEPENDENT_CHECKS.tsv',sep='\t',index=False)
qc={'status':'PASS','time':datetime.now(ZoneInfo('Asia/Taipei')).isoformat(),'checks':len(checks),'max_effect_difference':maxerr,'raw_selected_matrix_exact':True,'all_scores_rebuilt_from_raw':True,'all_patient_effects_SVD_verified':True,'all_sign_tests_combinatorial_verified':True,'all_BH_scipy_verified':True,'models':248,'cross_source_models':16,'source_SHA_unchanged':True}
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(qc,indent=2))
print(json.dumps(qc),flush=True)
