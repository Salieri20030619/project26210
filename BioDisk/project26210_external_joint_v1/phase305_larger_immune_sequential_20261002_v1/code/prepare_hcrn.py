from pathlib import Path
import os
os.environ['OMP_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['NUMBA_NUM_THREADS']='4'
import json,subprocess,gc
import numpy as np,pandas as pd,scipy.sparse as sp,scipy.io as sio,anndata as ad,scanpy as sc
P=Path(__file__).resolve().parents[1]
def gate(stage):
 r=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with (P/'logs/memory_gate.log').open('a') as f:f.write('HCRN_'+stage+'\n'+r.stdout+r.stderr)
 if r.returncode:raise RuntimeError('MEMORY_GATE_FAILED')
assert json.loads((P/'qc/HCRN_ARCHIVE_INPUT_VALIDATION.json').read_text())['status']=='PASS'
proto=json.loads((P/'resources/HCRN_QC_PROTOCOL.json').read_text());markers=proto['identity_markers']
plan=pd.read_csv(P/'tables/HCRN_SAMPLE_PLAN.tsv',sep='\t');plan=plan[plan.main_ccRCC];assert len(plan)==14 and plan.patient.nunique()==13
qc=[];objects=[]
for row in plan.itertuples():
 lib=row.library;dest=P/'inputs'/f'HCRN_{lib}_qc.h5ad';qf=P/'qc'/f'HCRN_{lib}.json'
 if dest.exists() and qf.exists():objects.append(ad.read_h5ad(dest));qc.append(json.loads(qf.read_text()));continue
 gate(lib);x=sio.mmread(row.matrix_path).T.tocsr();assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all();x=x.astype(np.int32)
 ft=pd.read_csv(row.genes_path,sep='\t',header=None);genes=ft.iloc[:,1].astype(str).to_numpy();bc=pd.read_csv(row.barcodes_path,sep='\t',header=None).iloc[:,0].astype(str).to_numpy();assert x.shape==(len(bc),len(genes))
 if ft.shape[1]>=3:
  kk=ft.iloc[:,2].eq('Gene Expression').to_numpy();x=x[:,kk].tocsr();genes=genes[kk]
 rawN=len(bc);umi=np.asarray(x.sum(axis=1)).ravel();ng=np.diff(x.indptr);mt=np.asarray(x[:,np.char.startswith(genes.astype(str),'MT-')].sum(axis=1)).ravel()/np.maximum(umi,1)
 keep=(umi>=500)&(ng>=200)&(ng<=8000)&(mt<=.20);x=x[keep].copy();un,inv=np.unique(genes,return_inverse=True)
 if len(un)<len(genes):x=(x@sp.csr_matrix((np.ones(len(inv),dtype=np.int32),(np.arange(len(inv)),inv)),shape=(len(inv),len(un)))).tocsr();genes=un
 obs=pd.DataFrame(index=[lib+':'+b for b in bc[keep]]);obs['dataset']='HCRN';obs['patient']='HCRN:'+row.patient;obs['library']=lib;obs['histology']='ccRCC';obs['timepoint']='Baseline' if row.timepoint=='BL' else 'Progression';obs['treatment']='pre_PD1' if row.timepoint=='BL' else 'post_PD1';obs['clinical_response']=row.response;obs['total_counts']=umi[keep];obs['n_genes']=ng[keep];obs['pct_mito']=mt[keep]
 a=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=genes));assert a.n_obs>=50,'Too few QC cells for fixed preprocessing; review required'
 print('QC',lib,rawN,a.n_obs,flush=True);sc.pp.scrublet(a,expected_doublet_rate=.06,random_state=305,n_prin_comps=30,verbose=False);nd=int(a.obs.predicted_doublet.sum());a=a[~a.obs.predicted_doublet.to_numpy()].copy();a.write_h5ad(dest,compression='gzip');objects.append(a)
 rec=dict(library=lib,patient=row.patient,timepoint=obs.timepoint.iloc[0],raw_cells=rawN,qc_cells=int(keep.sum()),doublets=nd,singlets=a.n_obs,genes=a.n_vars,TYMP_measured='TYMP' in a.var_names,SPP1_measured='SPP1' in a.var_names);qf.write_text(json.dumps(rec,indent=2));qc.append(rec);print('SINGLETS',lib,a.n_obs,flush=True);del x;gc.collect()
pd.DataFrame(qc).to_csv(P/'tables/HCRN_QC_BY_LIBRARY.tsv',sep='\t',index=False);gate('cluster');a=ad.concat(objects,join='outer',merge='same',fill_value=0);a.X=a.X.tocsr();assert a.obs_names.is_unique;del objects;gc.collect();a.write_h5ad(P/'inputs/HCRN_counts_qc.h5ad',compression='gzip')
norm=a.copy();sc.pp.normalize_total(norm,target_sum=1e4);sc.pp.log1p(norm);sc.pp.highly_variable_genes(norm,n_top_genes=2000,batch_key='patient',flavor='seurat');norm.var.loc[norm.var_names.isin(['TYMP','SPP1']),'highly_variable']=False
z=norm[:,norm.var.highly_variable].copy();sc.pp.scale(z,max_value=10);sc.pp.pca(z,n_comps=30,random_state=305);sc.pp.neighbors(z,n_neighbors=15,n_pcs=30,random_state=305)
for res,key in [(.6,'cluster'),(1.,'cluster_high')]:sc.tl.leiden(z,resolution=res,random_state=305,key_added=key,flavor='igraph',n_iterations=2);a.obs[key]=z.obs[key]
a.write_h5ad(P/'inputs/HCRN_counts_qc.h5ad',compression='gzip');a.obs.to_csv(P/'tables/HCRN_CELL_METADATA.tsv.gz',sep='\t');glob=np.asarray(norm.X.mean(axis=0)).ravel();rows=[];detail=[]
for cl in a.obs.cluster.cat.categories:
 ix=a.obs.cluster.eq(cl).to_numpy();mean=np.asarray(norm.X[ix].mean(axis=0)).ravel();pct=np.asarray((a.X[ix]>0).mean(axis=0)).ravel();eff=mean-glob;eff[a.var_names.isin(['TYMP','SPP1'])]=-np.inf
 rr={'cluster':cl,'cells':int(ix.sum()),'patients':a.obs.loc[ix,'patient'].nunique(),'top_genes':','.join(a.var_names[np.argsort(eff)[-25:][::-1]])}
 for k,gs in markers.items():
  jj=np.flatnonzero(a.var_names.isin(gs));rr[k]=float(mean[jj].mean()) if len(jj) else np.nan
  for j in jj:detail.append(dict(cluster=cl,marker_group=k,gene=a.var_names[j],mean_logCP10K=float(mean[j]),detection=float(pct[j])))
 rows.append(rr)
pd.DataFrame(rows).to_csv(P/'tables/HCRN_CLUSTER_OVERVIEW.tsv',sep='\t',index=False);pd.DataFrame(detail).to_csv(P/'tables/HCRN_CLUSTER_MARKERS.tsv',sep='\t',index=False)
(P/'qc/HCRN_PREPARATION_COMPLETE.json').write_text(json.dumps({'status':'READY_FOR_IDENTITY_REVIEW','cells':a.n_obs,'patients':a.obs.patient.nunique(),'libraries':a.obs.library.nunique(),'association_tests_run':False},indent=2));print('READY_FOR_IDENTITY_REVIEW',a.shape,flush=True)
