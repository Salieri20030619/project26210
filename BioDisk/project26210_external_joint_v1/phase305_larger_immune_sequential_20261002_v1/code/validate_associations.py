from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,math,subprocess
import numpy as np,pandas as pd,anndata as ad
from scipy.linalg import orth
from scipy.stats import t as tdist,false_discovery_control
P=Path(__file__).resolve().parents[1];checks=[]
def ck(name,ok):
 checks.append({'check':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
for n in ['OBRADOVIC_preparation','HCRN_preparation','refine_myeloid','qualify_patient_inputs','patient_associations']:ck('exit_'+n,(P/f'logs/{n}.exit').read_text().strip()=='0')
proto=json.loads((P/'resources/PATIENT_ASSOCIATION_PROTOCOL.json').read_text());m=pd.read_csv(P/'tables/PATIENT_MEASURES.tsv',sep='\t');w=pd.read_csv(P/'tables/WITHIN_MACRO_PATIENT_EFFECTS.tsv',sep='\t');ws=pd.read_csv(P/'tables/WITHIN_MACRO_SUMMARY.tsv',sep='\t');e=pd.read_csv(P/'tables/ALL_PATIENT_ASSOCIATIONS.tsv',sep='\t');pairs=pd.read_csv(P/'tables/ASSOCIATION_PATIENT_PAIRS.tsv',sep='\t');cov=pd.read_csv(P/'tables/PROGRAM_GENE_COVERAGE.tsv',sep='\t')
interest=sorted(set(['TYMP','SPP1']+sum(proto['programs'].values(),[])));gi={g:i for i,g in enumerate(interest)};acc={};cellvals={}
for ds in ['OBRADOVIC','HCRN']:
 meta=pd.read_csv(P/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0);parts=[]
 for lib in meta.library.unique():
  subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
  a=ad.read_h5ad(P/f'inputs/{ds}_{lib}_qc.h5ad');mm=meta.loc[a.obs_names];x=a.X.tocsr().astype(float);ck(lib+'_integer_nonnegative',np.isfinite(x.data).all() and np.all(x.data>=0) and np.all(x.data==np.floor(x.data)));jj=a.var_names.get_indexer(interest)
  for lin in m.lineage.unique():
   sel=mm['lineage_'+lin].to_numpy(dtype=bool);n=int(sel.sum());unit=str(mm.analysis_unit.iloc[0]);key=(unit,lin)
   if key not in acc:acc[key]={'n':0,'total':0.,'gene':np.zeros(len(interest)),'det':np.zeros(len(interest)),'available':None}
   q=acc[key];q['n']+=n
   if n:
    q['total']+=float(x[sel].sum());q['gene'][jj>=0]+=np.asarray(x[sel][:,jj[jj>=0]].sum(axis=0)).ravel();q['det'][jj>=0]+=np.asarray((x[sel][:,jj[jj>=0]]>0).sum(axis=0)).ravel();q['available']=set(a.var_names) if q['available'] is None else q['available']&set(a.var_names)
   if lin=='Macrophage' and n:
    v=pd.DataFrame(index=a.obs_names[sel]);total=np.asarray(x[sel].sum(axis=1)).ravel()
    for g in ['TYMP','SPP1']:v[g]=np.log1p(x[sel,a.var_names.get_loc(g)].toarray().ravel()*1e4/total)
    parts.append(v)
 full=pd.concat(parts) if parts else pd.DataFrame(columns=['TYMP','SPP1']);cellvals[ds]=(meta,full)
for i,r in m.iterrows():
 q=acc[(r.analysis_unit,r.lineage)];ck(f'measure{i}_n',q['n']==r.cells);gs=[r.measure] if r.measure in ['TYMP','SPP1'] else proto['programs'][r.measure];keep=[g for g in gs if g in (q['available'] or set())];ck(f'measure{i}_coverage',len(keep)==r.genes_measured);vec=np.log1p(q['gene']*1e4/max(1,q['total']));value=np.mean([vec[gi[g]] for g in keep]) if keep else np.nan;ck(f'measure{i}_original_library_score',np.isclose(value,r.value,atol=1e-10,equal_nan=True))
for i,r in cov.iterrows():
 q=acc[(r.analysis_unit,r.lineage)];ck(f'gene{i}_count',q['det'][gi[r.gene]]==r.detected_cells);v=np.log1p(q['gene'][gi[r.gene]]*1e4/max(1,q['total']));ck(f'gene{i}_value',np.isclose(v,r.pseudobulk_logCP10K,atol=1e-10))
def residual_rho(x,y,D):
 ranks=pd.DataFrame({'x':x,'y':y}).rank().to_numpy();Q=orth(D);res=ranks-Q@(Q.T@ranks);rho=float(res[:,0]@res[:,1]/np.sqrt((res[:,0]@res[:,0])*(res[:,1]@res[:,1])));return rho,Q.shape[1]
caches={r['analysis_unit']:r['path'] for r in json.loads((P/'inputs/WITHIN_CACHE_MANIFEST.json').read_text())}
for i,r in w.iterrows():
 if pd.isna(r.rho):continue
 ds=r.dataset;meta,full=cellvals[ds];ii=meta.analysis_unit.eq(r.analysis_unit)&meta.lineage_Macrophage;v=full.loc[meta.index[ii]];z=np.load(caches[r.analysis_unit]);ck(f'within{i}_TYMP_original',np.allclose(v.TYMP,z['TYMP'],atol=1e-10));ck(f'within{i}_SPP1_original',np.allclose(v.SPP1,z['SPP1'],atol=1e-10));rho,rank=residual_rho(z['TYMP'],z['SPP1'],z[r.model+'_D']);ck(f'within{i}_rho',abs(rho-r.rho)<1e-10);ck(f'within{i}_rank',rank==r.design_rank)
for i,r in ws.iterrows():
 g=w[w.dataset.eq(r.dataset)&w.model.eq(r.model)&w.rho.notna()];g=g[g.donor_first_available] if r.scope=='donor_first_available' else g[g.stratum.eq(r.scope)];ck(f'summary{i}_uniqueN',len(g)==r.patients and g.patient.nunique()==r.patients);ck(f'summary{i}_direction',int((g.rho>0).sum())==r.positive and int((g.rho<0).sum())==r.negative);ck(f'summary{i}_median',np.isclose(g.rho.median(),r.median_rho,equal_nan=True))
 if pd.notna(r.P):
  n=int(r.positive+r.negative);k=min(int(r.positive),int(r.negative));p=min(1.,2*sum(math.comb(n,j) for j in range(k+1))/2**n);ck(f'summary{i}_exactP',abs(p-r.P)<1e-12 and r.patients>=10)
for i,r in e.iterrows():
 g=pairs[pairs.dataset.eq(r.dataset)&pairs.scope.eq(r.scope)&pairs.source.eq(r.source)&pairs.target.eq(r.target)&pairs.model.eq(r.model)];ck(f'edge{i}_N',len(g)==r.patients and g.patient.nunique()==r.patients)
 for side in ['source','target']:
  lin,measure=r[side].split(':');q=m[m.dataset.eq(r.dataset)&m.stratum.eq(r.scope)&m.lineage.eq(lin)&m.measure.eq(measure)&m.eligible].set_index('analysis_unit');ck(f'edge{i}_{side}_qualified',set(g.analysis_unit)<=set(q.index))
  if len(g):ck(f'edge{i}_{side}_values',np.allclose(g['value_'+side],q.loc[g.analysis_unit,'value'],atol=1e-12))
 if pd.notna(r.rho):
  D=np.ones((len(g),1))
  if r.dataset=='OBRADOVIC' and r.model=='prespecified_adjusted':D=np.column_stack([D,pd.get_dummies(g[['technical_batch','pT_category']].astype(str),drop_first=True,dtype=float).to_numpy()])
  rho,rank=residual_rho(g.value_source,g.value_target,D);ck(f'edge{i}_rho',abs(rho-r.rho)<1e-10);ck(f'edge{i}_df',len(g)-rank-1==r.df)
  if pd.notna(r.P):p=float(2*tdist.sf(np.sqrt(r.df)*abs(rho)/np.sqrt(max(1-rho*rho,1e-15)),r.df));ck(f'edge{i}_P',abs(p-r.P)<1e-10 and r.patients>=10 and r.df>=5)
for name,df,groups in [('within',ws,['dataset','scope','family']),('across',e,['dataset','scope','family'])]:
 for key,g in df.groupby(groups):
  a=g[g.P.notna()];ck(name+'_'+str(key)+'_familyN',g.family_tested_n.eq(len(a)).all())
  if len(a):ck(name+'_'+str(key)+'_BH',np.allclose(false_discovery_control(a.P.to_numpy()),a.FDR,atol=1e-12))
ck('no_cell_inference',w.P.isna().all());ck('prespecified240cross_edges',len(e)==240)
(P/'qc/ASSOCIATION_VALIDATION_CHECKS.json').write_text(json.dumps(checks,indent=2));(P/'qc/FINAL_ASSOCIATION_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'original_library_aggregation':True,'independent_orthogonal_projection':True,'sign_P_recomputed':True,'scipy_BH_recomputed':True,'cross_edges':len(e)},indent=2));print('VALIDATED',len(checks),'checks',flush=True)
