from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,hashlib,subprocess
import pandas as pd,numpy as np,anndata as ad
from scipy.stats import rankdata
P=Path(__file__).resolve().parents[1];checks=[]
def ck(name,ok):
 checks.append({'check':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
for job in ['qc','author','lineages','sequential']:ck('exit_'+job,(P/f'logs/{job}.exit').read_text().strip()=='0')
for step in [1,2,3]:ck('step_marker_'+str(step),json.loads((P/f'qc/STEP{step}_COMPLETE.json').read_text())['status']=='COMPLETE')
e=pd.read_csv(P/'tables/ALL_SEQUENTIAL_EDGES.tsv',sep='\t');pairs=pd.read_csv(P/'tables/ALL_MATCHED_PATIENT_PAIRS.tsv',sep='\t');m=pd.read_csv(P/'tables/PATIENT_LINEAGE_MEASURES.tsv',sep='\t');w=pd.read_csv(P/'tables/STEP1_WITHIN_MACRO_PATIENT.tsv',sep='\t');coverage=pd.read_csv(P/'tables/FUNCTION_GENE_COVERAGE.tsv',sep='\t');proto=json.loads((P/'resources/SEQUENTIAL_ANALYSIS_PROTOCOL.json').read_text())
ck('120_prespecified_edges',len(e)==120);ck('no_small_N_inference',e.P.isna().all() and e.FDR.isna().all() and w.P.isna().all())
for i,r in e.iterrows():
 q=pairs[pairs.stage.eq(r.stage)&pairs.dataset.eq(r.dataset)&pairs.source.eq(r.source)&pairs.target.eq(r.target)]
 ck(f'edge_{i}_unique_patient_N',len(q)==r.n_patients and q.patient.nunique()==r.n_patients)
 if pd.notna(r.rho):ck(f'edge_{i}_rank_correlation',abs(q[['value_source','value_target']].rank().corr().iloc[0,1]-r.rho)<1e-12)
 for role in ['source','target']:
  lin,meas=r[role].split(':');v=m[m.dataset.eq(r.dataset)&m.lineage.eq(lin)&m.measure.eq(meas)&m.main&m.eligible].set_index('sample')
  ck(f'edge_{i}_{role}_qualification',set(q['sample'])<=set(v.index))
  if len(q):ck(f'edge_{i}_{role}_values',np.allclose(q['value_'+role],v.loc[q['sample'],'value']))
for sample,sm in m.groupby('sample'):
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 a=ad.read_h5ad(P/f'inputs/{sample}_qc.h5ad');ds=sm.dataset.iloc[0];meta=pd.read_csv(P/f'tables/{ds}_QUALIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col=0).loc[a.obs_names];x=a.X.tocsr().astype(np.float64)
 ck(sample+'_integer_counts',np.all(x.data==np.floor(x.data)) and np.all(x.data>=0))
 for lin,g in sm.groupby('lineage'):
  sel=meta['lineage_'+lin].to_numpy(dtype=bool);n=int(sel.sum());total=np.asarray(x[sel].sum(axis=0)).ravel();pb=np.log1p(total*1e4/max(1,total.sum()))
  cv=coverage[coverage['sample'].eq(sample)&coverage.lineage.eq(lin)];jj=a.var_names.get_indexer(cv.gene)
  ck(sample+'_'+lin+'_n',g.n_cells.eq(n).all());ck(sample+'_'+lin+'_genes',np.all(jj>=0))
  ck(sample+'_'+lin+'_original_library_pseudobulk',np.allclose(pb[jj],cv.pseudobulk_logCP10K,atol=1e-10))
  for r in g.itertuples():
   panel=[r.measure] if r.measure in ['TYMP','SPP1'] else proto['programs'][r.measure]
   avail=[gene for gene in panel if gene in a.var_names]
   ck(sample+'_'+lin+'_'+r.measure+'_original_coverage',len(avail)==r.genes_measured)
   val=pb[a.var_names.get_indexer(avail)].mean() if avail else np.nan
   ck(sample+'_'+lin+'_'+r.measure+'_score',np.isclose(val,r.value,equal_nan=True))
 for r in w[w['sample'].eq(sample)].itertuples():
  ix=meta.lineage_Macrophage.to_numpy(dtype=bool);xx=x[ix];norm=xx.multiply((1e4/np.asarray(xx.sum(axis=1)).ravel())[:,None]).tocsr();norm.data=np.log1p(norm.data)
  xy=np.column_stack([norm[:,a.var_names.get_loc(g)].toarray().ravel() for g in ['TYMP','SPP1']]);xy=np.apply_along_axis(rankdata,0,xy)
  cs=[meta.loc[ix,'total_counts'].to_numpy(),meta.loc[ix,'n_genes'].to_numpy()]
  if r.model.endswith('environment'):
   for k,gs in proto['context'].items():cs.append(np.asarray(norm[:,a.var_names.isin(gs)].mean(axis=1)).ravel())
  D=np.column_stack([np.ones(sum(ix))]+[rankdata(c) for c in cs if np.ptp(c)>0]);res=xy-D@np.linalg.pinv(D)@xy;rho=np.corrcoef(res.T)[0,1]
  ck(sample+'_'+r.model+'_independent_partial',abs(rho-r.rho)<1e-10)
print('VALIDATED',len(checks),'checks',flush=True)
(P/'qc/VALIDATION_CHECKS.json').write_text(json.dumps(checks,indent=2));(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'independent_library_aggregation':True,'all_edges_recomputed':True,'no_cell_P':True},indent=2))
