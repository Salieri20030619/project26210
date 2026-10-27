from pathlib import Path
import os
os.environ['OMP_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['NUMBA_NUM_THREADS']='4'
import json,gzip,hashlib,subprocess,re,datetime,gc
import numpy as np,pandas as pd,scipy.sparse as sp,scipy.io as sio
import anndata as ad,scanpy as sc
P=Path(__file__).resolve().parents[1];R=P.parent
def gate(stage):
 z=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with (P/'logs/qc_memory.log').open('a') as f:f.write(stage+'\n'+z.stdout+z.stderr)
 if z.returncode:raise RuntimeError('MEMORY_GATE')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def aggregate_symbols(x,genes):
 genes=np.asarray(genes,dtype=str);unique,inv=np.unique(genes,return_inverse=True)
 if len(unique)==len(genes):return x,genes
 mapping=sp.csr_matrix((np.ones(len(inv),dtype=np.int32),(np.arange(len(inv)),inv)),shape=(len(inv),len(unique)))
 return (x@mapping).tocsr(),unique
idpan=json.loads((R/'phase286_CAF_identity_boundary_audit_20260929_v1/resources/IDENTITY_PANELS_FROZEN.json').read_text())
panels=json.loads((R/'phase292_barrier_execution_and_external_feasibility_20260930_v1/PANELS_FROZEN.json').read_text())
markers={
 'Fibroblast':idpan['fibro_core']+idpan['fibro_support'],
 'Mural':idpan['pericyte']+idpan['smc'],
 'Endothelial':idpan['endothelial'],
 'Macrophage':['C1QA','C1QB','C1QC','CSF1R','CD68','CD163','MRC1','APOE','LST1'],
 'Monocyte':['FCN1','S100A8','S100A9','VCAN','CD14'],
 'Tcell':['CD3D','CD3E','TRAC','CD2','CD247'],
 'NK':['GNLY','NKG7','KLRD1','NCR1','FCGR3A'],
 'Bcell':['MS4A1','CD79A','CD79B','CD37','CD74'],
 'Plasma':['JCHAIN','MZB1','IGHG1','SDC1'],
 'Epithelial':['EPCAM','KRT8','KRT18','KRT19','PAX8'],
 'RCC':['CA9','NDUFA4L2','EGLN3','NNMT','KRT8','KRT18'],
 'NormalTubule':['LRP2','CUBN','SLC34A1','ALDOB','UMOD','SLC12A1'],
 'DC':['CD1C','FCER1A','CLEC10A','CLEC9A','XCR1','LILRA4'],
 'Mast':['TPSAB1','TPSB2','CPA3','KIT'],
 'CD8':['CD8A','CD8B'],'CD4':['CD4','IL7R','LTB'],'Treg':['FOXP3','IL2RA','IKZF2'],
 'CAF_activation':['FAP','POSTN','COL11A1','CTHRC1','CTGF'],
}
rules=dict(created_at=datetime.datetime.now().astimezone().isoformat(),qc=dict(min_UMI=500,min_genes=200,max_genes=8000,max_mito_fraction=.20),doublet='Scanpy native scrublet per library, expected rate .06, automatic threshold, seed303; failure blocks rather than silently omits',cluster=dict(hvg=2000,exclude=['TYMP','SPP1'],pcs=30,neighbors=15,resolution=.6,comparison_resolution=1.0,seed=303),markers=markers,identity_panels=idpan,function_panels=panels,selection='GSE224630 fresh RCC1-5; GSE210038 T1-7 with treated T1/T3 separate; GSE222703 original 3 tumour donors; controls/repeat/thrombus not pooled into patient n',association_tests_run=False)
freeze=P/'resources/QUALIFICATION_PROTOCOL.json'
if freeze.exists():assert json.loads(freeze.read_text())['qc']==rules['qc']
else:freeze.write_text(json.dumps(rules,ensure_ascii=False,indent=2))
samples=[]
for ds in ['GSE210038','GSE224630','GSE222703']:
 for s in json.loads((P/'metadata'/f'{ds}_samples.json').read_text()):
  ch={x['tag'].lower():x['value'].strip() for x in s['characteristics']};title=s['title'];selected=False;main=False
  if ds=='GSE210038':
   selected=title.startswith('Tumour');pt=ch['patient'];tr=ch['prior treatment'];main=selected and tr=='Naive'
  elif ds=='GSE224630':
   selected='frozen' not in title.lower() and 'Thrombus' not in title;pt='Patient '+re.search('patient (\d+)',title).group(1);tr='Untreated (original publication)';main=selected
  else:
   selected='tumor -' in title;pt=ch['patient'];tr='Not reported in GEO sample';main=selected
  samples.append(dict(dataset=ds,sample=s['accession'],patient=ds+':'+pt,treatment=tr,title=title,selected=selected,main=main,raw_files=';'.join(s['files'])))
plan=pd.DataFrame(samples);plan.to_csv(P/'tables/SAMPLE_PLAN.tsv',sep='\t',index=False)
qc=[];hashes=json.loads((P/'inputs/RAW_INPUT_SHA256.json').read_text()) if (P/'inputs/RAW_INPUT_SHA256.json').exists() else {}
for ds in ['GSE210038','GSE224630','GSE222703']:
 dplan=plan[plan.dataset.eq(ds)&plan.selected];objects=[]
 if all((P/f'{folder}/{ds}_{name}').exists() for folder,name in [('inputs','counts_qc.h5ad'),('tables','CLUSTER_IDENTITY_OVERVIEW.tsv'),('tables','CLUSTER_MARKERS.tsv')]):
  qc.extend(json.loads((P/'qc'/f'{sample}.json').read_text()) for sample in dplan['sample']);print('REUSE_VALIDATED_QC_CLUSTER',ds,flush=True);continue
 author=None;authorX=None
 if ds=='GSE222703':
  if not (P/'inputs/GSE222703_AUTHOR_COUNTS.mtx.gz').exists():raise RuntimeError('AUTHOR_EXPORT_MISSING')
  author=pd.read_csv(P/'inputs/GSE222703_AUTHOR_METADATA.tsv.gz',sep='\t')
  authorX=sio.mmread(P/'inputs/GSE222703_AUTHOR_COUNTS.mtx.gz').T.tocsr().astype(np.int32)
  authorGenes=np.array((P/'inputs/GSE222703_AUTHOR_GENES.txt').read_text().splitlines())
  print('AUTHOR_METADATA_FIELDS',author.columns.tolist(),flush=True)
 for row in dplan.itertuples():
  dest=P/'inputs'/f'{row.sample}_qc.h5ad'
  if dest.exists():
   aa=ad.read_h5ad(dest);objects.append(aa);qc.append(json.loads((P/'qc'/f'{row.sample}.json').read_text()));continue
  gate(row.sample);print('READ',row.sample,flush=True)
  if ds=='GSE210038':
   f=next((P/'downloads'/ds).glob(row.sample+'*.gz'));parts=[];genes=[];barcodes=None
   for chunk in pd.read_csv(f,sep='\t',index_col=0,chunksize=1000):
    if barcodes is None:barcodes=chunk.columns.to_numpy()
    genes.extend(chunk.index.astype(str));parts.append(sp.csr_matrix(chunk.to_numpy(dtype=np.int32)))
   x=sp.vstack(parts,format='csr').T.tocsr();del parts
   hashes[str(f)]=sha(f);obs=pd.DataFrame(index=barcodes)
  elif ds=='GSE224630':
   f=next((P/'downloads'/ds).glob(row.sample+'*matrix.mtx.gz'));x=sio.mmread(f).T.tocsr().astype(np.int32)
   ft=next((P/'downloads'/ds).glob(row.sample+'*features.tsv.gz'));genes=pd.read_csv(ft,sep='\t',header=None)[1].astype(str).to_numpy()
   fb=next((P/'downloads'/ds).glob(row.sample+'*barcodes.tsv.gz'));barcodes=pd.read_csv(fb,sep='\t',header=None)[0].astype(str).to_numpy();obs=pd.DataFrame(index=barcodes)
   for path in [f,ft,fb]:hashes[str(path)]=sha(path)
  else:
   # Official object identifiers are resolved from original metadata before running this branch.
   mp=json.loads((P/'resources/AUTHOR_SAMPLE_MAP.json').read_text());field=mp['field'];value=mp['mapping'][row.sample]
   ii=np.flatnonzero(author[field].astype(str).eq(value));assert len(ii)>0
   x=authorX[ii];genes=authorGenes;obs=author.iloc[ii].copy().set_index('cell_id');obs.index=obs.index.astype(str)
   for path in [P/'inputs/GSE222703_AUTHOR_METADATA.tsv.gz',P/'inputs/GSE222703_AUTHOR_COUNTS.mtx.gz']:hashes[str(path)]=sha(path)
  assert x.shape[0]==len(obs) and np.isfinite(x.data).all() and (x.data>=0).all()
  rawN=x.shape[0];umi=np.asarray(x.sum(axis=1)).ravel();ng=np.diff(x.indptr);mito=np.char.startswith(np.asarray(genes,dtype=str),'MT-');mt=np.asarray(x[:,mito].sum(axis=1)).ravel()/np.maximum(umi,1)
  keep=(umi>=500)&(ng>=200)&(ng<=8000)&(mt<=.20)
  x=x[keep].copy();obs=obs.iloc[np.flatnonzero(keep)].copy();obs['total_counts']=umi[keep];obs['n_genes']=ng[keep];obs['pct_mito']=mt[keep]
  x,genes=aggregate_symbols(x,genes);aa=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=genes))
  aa.obs_names=[row.sample+':'+str(b) for b in aa.obs_names]
  for k in ['dataset','sample','patient','treatment','main']:aa.obs[k]=getattr(row,k)
  print('QC',row.sample,rawN,aa.n_obs,flush=True)
  sc.pp.scrublet(aa,expected_doublet_rate=.06,random_state=303,n_prin_comps=30,verbose=False)
  nd=int(aa.obs.predicted_doublet.sum());aa=aa[~aa.obs.predicted_doublet.to_numpy()].copy()
  aa.write_h5ad(dest,compression='gzip');objects.append(aa)
  rec=dict(dataset=ds,sample=row.sample,raw_barcodes=rawN,qc_cells=int(keep.sum()),predicted_doublets=nd,singlets=aa.n_obs,genes=aa.n_vars,total_UMI=int(aa.X.sum()))
  qc.append(rec);(P/'qc'/f'{row.sample}.json').write_text(json.dumps(rec,indent=2));print('SINGLET',rec,flush=True)
  del x,obs,umi,ng,mt;gc.collect()
 gate(ds+'_cluster');a=ad.concat(objects,join='outer',merge='same',fill_value=0);a.X=a.X.tocsr();del objects;gc.collect()
 a.write_h5ad(P/'inputs'/f'{ds}_counts_qc.h5ad',compression='gzip')
 norm=a.copy();sc.pp.normalize_total(norm,target_sum=1e4);sc.pp.log1p(norm)
 sc.pp.highly_variable_genes(norm,n_top_genes=2000,batch_key='sample',flavor='seurat');norm.var.loc[norm.var_names.isin(['TYMP','SPP1']),'highly_variable']=False
 z=norm[:,norm.var.highly_variable].copy();sc.pp.scale(z,max_value=10);sc.pp.pca(z,n_comps=30,random_state=303);sc.pp.neighbors(z,n_neighbors=15,n_pcs=30,random_state=303)
 sc.tl.leiden(z,resolution=.6,random_state=303,key_added='cluster',flavor='igraph',n_iterations=2)
 sc.tl.leiden(z,resolution=1.,random_state=303,key_added='cluster_high',flavor='igraph',n_iterations=2)
 a.obs['cluster']=z.obs.cluster;a.obs['cluster_high']=z.obs.cluster_high
 a.obs.to_csv(P/'tables'/f'{ds}_CELL_CLUSTER_METADATA.tsv.gz',sep='\t')
 a.write_h5ad(P/'inputs'/f'{ds}_counts_qc.h5ad',compression='gzip')
 coverage=[];rows=[];top=[]
 for cl in a.obs.cluster.cat.categories:
  ix=np.flatnonzero(a.obs.cluster.eq(cl));mean=np.asarray(norm.X[ix].mean(axis=0)).ravel();pct=np.asarray((a.X[ix]>0).mean(axis=0)).ravel()
  glob=np.asarray(norm.X.mean(axis=0)).ravel();effects=mean-glob
  topgenes=np.argsort(effects)[-25:][::-1]
  rr=dict(dataset=ds,cluster=cl,cells=len(ix),patients=a.obs.iloc[ix].patient.nunique(),samples=';'.join(a.obs.iloc[ix]['sample'].value_counts().index.astype(str)),top_genes=','.join(a.var_names[topgenes]))
  for label,gs in markers.items():
   jj=np.flatnonzero(a.var_names.isin(gs));rr[label]=float(mean[jj].mean()) if len(jj) else 0
   for j in jj:coverage.append(dict(dataset=ds,cluster=cl,marker_group=label,gene=a.var_names[j],mean_logCP10K=mean[j],detected_fraction=pct[j],cells=len(ix)))
  rows.append(rr)
 pd.DataFrame(rows).to_csv(P/'tables'/f'{ds}_CLUSTER_IDENTITY_OVERVIEW.tsv',sep='\t',index=False)
 pd.DataFrame(coverage).to_csv(P/'tables'/f'{ds}_CLUSTER_MARKERS.tsv',sep='\t',index=False)
 print('CLUSTERS_READY',ds,len(rows),a.n_obs,flush=True)
 del a,norm,z;gc.collect()
 pd.DataFrame(qc).to_csv(P/'tables/QC_BY_LIBRARY.tsv',sep='\t',index=False)
 (P/'inputs/RAW_INPUT_SHA256.json').write_text(json.dumps(hashes,indent=2))
(P/'qc/QC_CLUSTER_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',association_tests_run=False,datasets=3),indent=2))
