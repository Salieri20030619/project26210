from pathlib import Path
import os
os.environ['OMP_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['NUMBA_NUM_THREADS']='4'
import json,gzip,hashlib,subprocess,gc,zipfile,shutil,datetime
import numpy as np,pandas as pd,scipy.sparse as sp,scipy.io as sio,anndata as ad,scanpy as sc
P=Path(__file__).resolve().parents[1];R=P.parent
def gate(stage):
 r=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with (P/'logs/memory_gate.log').open('a') as f:f.write(stage+'\n'+r.stdout+r.stderr)
 if r.returncode:raise RuntimeError('MEMORY_GATE_FAILED')
assert (P/'logs/OBRADOVIC_download.exit').read_text().strip()=='0'
assert json.loads((P/'qc/OBRADOVIC_DOWNLOAD_STATE.json').read_text())['status']=='COMPLETE'
rule=json.loads((P/'resources/OBRADOVIC_QC_PROTOCOL.json').read_text());markers=rule['identity_markers']
samples=json.loads((P/'metadata/OBRADOVIC_AUTHOR_SAMPLE_MAP.json').read_text());samples=[s for s in samples if s['tissue']=='Tumor']
assert len(set(s['patient'] for s in samples if s['cd45']=='CD45+'))==11
raw=P/'inputs/OBRADOVIC_raw';raw.mkdir(exist_ok=True)
with zipfile.ZipFile(P/'downloads/OBRADOVIC_singlecell_rawdata.zip') as z:
 names=z.namelist()
 for s in samples:
  s['library']=Path(s['path']).parts[2];s['selected_files']=[]
  for f in s['files']:
   target=f['path'];matches=[n for n in names if n==target or n.endswith('/'+target)]
   if not matches:raise RuntimeError('Missing official archive member '+target)
   assert len(matches)==1;dest=raw/target;assert dest.resolve().is_relative_to(raw.resolve());dest.parent.mkdir(parents=True,exist_ok=True)
   if not dest.exists():
    with z.open(matches[0]) as source,dest.open('wb') as out:shutil.copyfileobj(source,out,8*1024*1024)
   assert dest.stat().st_size==z.getinfo(matches[0]).file_size;s['selected_files'].append(str(dest))
pd.DataFrame([{k:v for k,v in s.items() if k not in ['files','selected_files']} for s in samples]).to_csv(P/'tables/OBRADOVIC_SAMPLE_PLAN.tsv',sep='\t',index=False)
qc=[];objects=[]
for s in samples:
 lib=s['library'];dest=P/'inputs'/f'OBRADOVIC_{lib}_qc.h5ad';qf=P/'qc'/f'OBRADOVIC_{lib}.json'
 if dest.exists() and qf.exists():objects.append(ad.read_h5ad(dest));qc.append(json.loads(qf.read_text()));continue
 gate(lib);d=raw/s['path'];m=next(d.glob('matrix.mtx*'));x=sio.mmread(m).T.tocsr();assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all();x=x.astype(np.int32)
 gf=next(iter(list(d.glob('genes.tsv*'))+list(d.glob('features.tsv*'))));feature=pd.read_csv(gf,sep='\t',header=None);genes=feature.iloc[:,1].astype(str).to_numpy()
 bf=next(d.glob('barcodes.tsv*'));bc=pd.read_csv(bf,sep='\t',header=None).iloc[:,0].astype(str).to_numpy();assert x.shape==(len(bc),len(genes));assert np.isfinite(x.data).all() and (x.data>=0).all()
 if feature.shape[1]>=3:
  keep=feature.iloc[:,2].eq('Gene Expression').to_numpy();x=x[:,keep].tocsr();genes=genes[keep]
 rawN=len(bc);umi=np.asarray(x.sum(axis=1)).ravel();ng=np.diff(x.indptr);mt=np.asarray(x[:,np.char.startswith(genes.astype(str),'MT-')].sum(axis=1)).ravel()/np.maximum(umi,1)
 keep=(umi>=500)&(ng>=200)&(ng<=8000)&(mt<=.20);x=x[keep].copy()
 un,inv=np.unique(genes,return_inverse=True)
 if len(un)<len(genes):x=(x@sp.csr_matrix((np.ones(len(inv),dtype=np.int32),(np.arange(len(inv)),inv)),shape=(len(inv),len(un)))).tocsr();genes=un
 obs=pd.DataFrame(index=[lib+':'+b for b in bc[keep]]);obs['dataset']='OBRADOVIC';obs['patient']='OBRADOVIC:'+s['patient'];obs['library']=lib;obs['tissue']=s['tissue'];obs['selection']=s['cd45'];obs['treatment']='Untreated';obs['total_counts']=umi[keep];obs['n_genes']=ng[keep];obs['pct_mito']=mt[keep]
 a=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=genes));assert a.n_obs>=50,'Too few cells for fixed preprocessing; manual review required'
 print('QC',lib,s['patient'],s['cd45'],rawN,a.n_obs,flush=True)
 sc.pp.scrublet(a,expected_doublet_rate=.06,random_state=305,n_prin_comps=30,verbose=False);nd=int(a.obs.predicted_doublet.sum());a=a[~a.obs.predicted_doublet.to_numpy()].copy()
 a.write_h5ad(dest,compression='gzip');objects.append(a);rec=dict(library=lib,patient=s['patient'],selection=s['cd45'],raw_cells=rawN,qc_cells=int(keep.sum()),doublets=nd,singlets=a.n_obs,genes=a.n_vars,TYMP_measured='TYMP' in a.var_names,SPP1_measured='SPP1' in a.var_names)
 qf.write_text(json.dumps(rec,indent=2));qc.append(rec);print('SINGLETS',lib,a.n_obs,flush=True);del x;gc.collect()
pd.DataFrame(qc).to_csv(P/'tables/OBRADOVIC_QC_BY_LIBRARY.tsv',sep='\t',index=False)
gate('clustering');a=ad.concat(objects,join='outer',merge='same',fill_value=0);a.X=a.X.tocsr();assert a.obs_names.is_unique;del objects;gc.collect();a.write_h5ad(P/'inputs/OBRADOVIC_counts_qc.h5ad',compression='gzip')
norm=a.copy();sc.pp.normalize_total(norm,target_sum=1e4);sc.pp.log1p(norm);sc.pp.highly_variable_genes(norm,n_top_genes=2000,batch_key='patient',flavor='seurat');norm.var.loc[norm.var_names.isin(['TYMP','SPP1']),'highly_variable']=False
z=norm[:,norm.var.highly_variable].copy();sc.pp.scale(z,max_value=10);sc.pp.pca(z,n_comps=30,random_state=305);sc.pp.neighbors(z,n_neighbors=15,n_pcs=30,random_state=305)
for resolution,key in [(.6,'cluster'),(1.,'cluster_high')]:sc.tl.leiden(z,resolution=resolution,random_state=305,key_added=key,flavor='igraph',n_iterations=2);a.obs[key]=z.obs[key]
a.write_h5ad(P/'inputs/OBRADOVIC_counts_qc.h5ad',compression='gzip');a.obs.to_csv(P/'tables/OBRADOVIC_CELL_METADATA.tsv.gz',sep='\t')
glob=np.asarray(norm.X.mean(axis=0)).ravel();rows=[];detail=[]
for cl in a.obs.cluster.cat.categories:
 ix=a.obs.cluster.eq(cl).to_numpy();mean=np.asarray(norm.X[ix].mean(axis=0)).ravel();pct=np.asarray((a.X[ix]>0).mean(axis=0)).ravel();effect=mean-glob;effect[a.var_names.isin(['TYMP','SPP1'])]=-np.inf
 rr={'cluster':cl,'cells':int(ix.sum()),'patients':a.obs.loc[ix,'patient'].nunique(),'CD45_pos_cells':int((a.obs.loc[ix,'selection']=='CD45+').sum()),'top_genes':','.join(a.var_names[np.argsort(effect)[-25:][::-1]])}
 for k,gs in markers.items():
  jj=np.flatnonzero(a.var_names.isin(gs));rr[k]=float(mean[jj].mean()) if len(jj) else np.nan
  for j in jj:detail.append(dict(cluster=cl,marker_group=k,gene=a.var_names[j],mean_logCP10K=float(mean[j]),detection=float(pct[j])))
 rows.append(rr)
pd.DataFrame(rows).to_csv(P/'tables/OBRADOVIC_CLUSTER_OVERVIEW.tsv',sep='\t',index=False);pd.DataFrame(detail).to_csv(P/'tables/OBRADOVIC_CLUSTER_MARKERS.tsv',sep='\t',index=False)
(P/'qc/OBRADOVIC_PREPARATION_COMPLETE.json').write_text(json.dumps({'status':'READY_FOR_IDENTITY_REVIEW','cells':a.n_obs,'patients':a.obs.patient.nunique(),'libraries':a.obs.library.nunique(),'association_tests_run':False},indent=2));print('READY_FOR_IDENTITY_REVIEW',a.shape,flush=True)
