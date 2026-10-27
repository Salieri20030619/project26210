from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import json,hashlib,subprocess,gc,math
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata,t
P=Path(__file__).resolve().parents[1];R=P.parent;S=R/'phase305_larger_immune_sequential_20261002_v1';pr=json.loads((P/'resources/ANALYSIS_PROTOCOL.json').read_text());e=pd.read_csv(P/'tables/ALL_ASSOCIATIONS.tsv',sep='\t');pairs=pd.read_csv(P/'tables/ALL_PATIENT_PAIRS.tsv.gz',sep='\t');measures=pd.read_csv(P/'tables/CD8_PATIENT_MEASURES.tsv',sep='\t');src=pd.read_csv(P/'tables/SOURCE_PATIENT_MEASURES.tsv',sep='\t');man=json.loads((P/'inputs/CACHE_MANIFEST.json').read_text());checks=[]
def ck(name,ok,detail=None):
 checks.append({'check':name,'pass':bool(ok),'detail':detail})
 if not ok:raise AssertionError((name,detail))
ck('protocol_hash',hashlib.sha256((P/'resources/ANALYSIS_PROTOCOL.json').read_bytes()).hexdigest()==(P/'resources/ANALYSIS_PROTOCOL.sha256').read_text())
for ds in ['OBRADOVIC','HCRN']:
 q=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(P/'logs'/f'validate_memory_gate_{ds}.log').write_text(q.stdout+q.stderr);ck(ds+'_memory',q.returncode==0)
 a=ad.read_h5ad(S/f'inputs/{ds}_counts_qc.h5ad');m=pd.read_csv(S/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0);ck(ds+'_all_obs_identity',a.obs_names.equals(m.index));x=a.X.tocsr();present={}
 for lib in m.library.unique():
  b=ad.read_h5ad(S/f'inputs/{ds}_{lib}_qc.h5ad',backed='r');present[lib]=set(b.var_names);b.file.close()
 for z in [v for v in man if v['dataset']==ds]:
  unit=z['analysis_unit'];cache=np.load(z['path']);mask=m.analysis_unit.eq(unit).to_numpy()&m.lineage_CD8.to_numpy(bool);n=int(mask.sum());rows=np.flatnonzero(mask);ss=np.asarray(x[rows].sum(axis=0)).ravel();dd=np.asarray((x[rows]>0).sum(axis=0)).ravel();available=set.intersection(*[present[v] for v in m.loc[mask,'library'].unique()]) if n else set()
  ck(unit+'_rows',np.array_equal(rows,cache['raw_rows']));ck(unit+'_all_gene_sums',np.array_equal(ss,cache['raw_sum']));ck(unit+'_all_gene_detection',np.array_equal(dd,cache['detected_cells']));ck(unit+'_gene_order',np.array_equal(a.var_names,cache['genes']));ck(unit+'_measured_genes',available==set(cache['available']))
  for name,v in pr['programs'].items():
   keep=[g for g in v['genes'] if g in available];jj=a.var_names.get_indexer(keep);vals=[math.log1p(float(ss[j])*10000/float(ss.sum())) for j in jj] if ss.sum() else [0]*len(jj);score=math.fsum(vals)/len(vals) if vals else np.nan;mr=measures[(measures.dataset==ds)&(measures.analysis_unit==unit)&(measures.endpoint==name)].iloc[0]
   ck(unit+'_'+name+'_score',np.isclose(score,mr.value,atol=1e-12,equal_nan=True));ok=n>=20 and len(jj)>=.8*len(v['genes']) and int(np.sum(dd[jj]>=3))>=.2*len(v['genes']);strict=ok and n>=100 and int(np.sum(dd[jj]>=max(5,math.ceil(.05*n))))>=.5*len(v['genes']);ck(unit+'_'+name+'_eligibility',bool(mr.eligible_primary)==ok and bool(mr.eligible_strict)==strict)
  ns=int((m.analysis_unit.eq(unit)&m.lineage_Macrophage).sum());mr=measures[(measures.dataset==ds)&(measures.analysis_unit==unit)&(measures.family=='capture_composition')].iloc[0];ck(unit+'_composition',np.isclose(mr.value,n/(n+ns) if n+ns else np.nan,equal_nan=True))
  for source,sr in src[(src.dataset==ds)&(src.analysis_unit==unit)].groupby('source'):
   rr=sr.iloc[0];lin,gene=source.split(':');sel=m.analysis_unit.eq(unit).to_numpy()&m['lineage_'+lin].to_numpy(bool);nn=int(sel.sum());sums=np.asarray(x[sel].sum(axis=0)).ravel();avail=set.intersection(*[present[l] for l in m.loc[sel,'library'].unique()]) if nn else set();j=a.var_names.get_loc(gene);v=math.log1p(sums[j]*10000/sums.sum()) if gene in avail and sums.sum() else np.nan;seen=int((x[sel,j]>0).sum());ck(unit+'_'+source+'_raw_source',nn==rr.cells and np.isclose(v,rr.value,atol=1e-12,equal_nan=True));ck(unit+'_'+source+'_qualification',bool(rr.eligible_primary)==bool(nn>=30 and gene in avail and seen>=3))
 print('VALIDATED_RAW',ds,flush=True);del a,m,x;gc.collect()
def independent(g,model,depth=False):
 n=len(g)
 if n==0:return np.nan,0,-1
 cols=[np.ones(n)]
 if model=='batch_stage_adjusted' and g.dataset.iloc[0]=='OBRADOVIC':
  for name in ['technical_batch','pT_category']:
   vals=g[name].astype(str).to_numpy();levels=sorted(set(vals));cols.extend([(vals==v).astype(float) for v in levels[1:]])
 if depth:
  for name in ['source_median_UMI','receiver_median_UMI']:cols.append(rankdata(np.round(g[name].to_numpy(float),12)))
 D=np.array(cols).T;U,s,V=np.linalg.svd(D,full_matrices=False);rk=int(np.linalg.matrix_rank(D));df=n-rk-1
 if n<3 or df<1:return np.nan,rk,df
 ranks=np.column_stack([rankdata(np.round(g[c].to_numpy(float),12)) for c in ['source_value','receiver_value']]);z=ranks-U[:,:rk]@(U[:,:rk].T@ranks)
 if (np.std(z,axis=0)<1e-9).any():return np.nan,rk,df
 return float(np.dot(z[:,0],z[:,1])/np.linalg.norm(z[:,0])/np.linalg.norm(z[:,1])),rk,df
pg={k:v for k,v in pairs.groupby('edge_id')};sx=pd.read_csv(P/'tables/ALL_SENSITIVITIES.tsv',sep='\t');ll=pd.read_csv(P/'tables/LEAVE_ONE_PATIENT.tsv',sep='\t')
for row in e.itertuples():
 g=pg.get(row.edge_id,pairs.iloc[:0]);ck(row.edge_id+'_unique_patient',g.patient.nunique()==len(g)==row.patients)
 rho,rk,df=independent(g,row.model);ck(row.edge_id+'_rho',np.isclose(rho,row.rho,atol=1e-10,equal_nan=True));ck(row.edge_id+'_rank_df',rk==row.design_rank and df==row.df)
 pv=2*t.sf(abs(rho)*np.sqrt(df/max(1-rho*rho,1e-15)),df) if len(g)>=10 and df>=5 and np.isfinite(rho) else np.nan;ck(row.edge_id+'_P',np.isclose(pv,row.P,atol=1e-10,equal_nan=True))
 for z in sx[sx.edge_id==row.edge_id].itertuples():
  sub=g[g.source_strict&g.receiver_strict] if z.sensitivity=='strict_qualification' else g;v,_,_=independent(sub,row.model,z.sensitivity=='source_receiver_depth_adjusted');ck(row.edge_id+'_'+z.sensitivity,np.isclose(v,z.rho,atol=1e-9,equal_nan=True))
 for z in ll[ll.edge_id==row.edge_id].itertuples():
  v,_,_=independent(g[g.patient!=z.left_out],row.model);ck(row.edge_id+'_leave_'+z.left_out,np.isclose(v,z.rho,atol=1e-9,equal_nan=True))
for key,g in e.groupby(['dataset','stratum','family']):
 valid=g[g.P.notna()].sort_values('P');n=len(valid)
 if n:
  q=np.minimum.accumulate((valid.P.to_numpy()*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1);ck(str(key)+'_BH',np.allclose(q,valid.FDR,atol=1e-12))
 ck(str(key)+'_family_n',all(g.family_tested_n==n))
ck('all_planned_edges',len(e)==918);ck('no_small_sample_P',e[e.patients<10].P.isna().all());ck('all_candidate_rule',(e.nominal_candidate==(e.P<.05)).all())
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'complete_raw_patient_gene_count_validation':True,'independent_SVD_models_sensitivities_leave_patient':True,'all_checks':checks},indent=2));print('FINAL_VALIDATION_PASS',len(checks),flush=True)
