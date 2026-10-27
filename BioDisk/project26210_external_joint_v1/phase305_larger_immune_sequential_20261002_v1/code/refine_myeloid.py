from pathlib import Path
import os
os.environ['OMP_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['NUMBA_NUM_THREADS']='4'
import json,subprocess,gc
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
P=Path(__file__).resolve().parents[1]
rule={'OBRADOVIC':[0,2,3,25],'HCRN':[9,10,20],'exclude_hvg':['TYMP','SPP1'],'purpose':'Separate mature macrophages, monocytes and dendritic cells before inspecting target associations','hvg':2000,'pcs':30,'neighbors':15,'resolution':0.7,'seed':305}
f=P/'resources/MYELOID_REFINEMENT_FROZEN.json'
if not f.exists():f.write_text(json.dumps(rule,indent=2))
else:assert json.loads(f.read_text())==rule
markers=json.loads((P/'resources/OBRADOVIC_QC_PROTOCOL.json').read_text())['identity_markers'];markers.update({'cDC':['CD1C','FCER1A','CLEC10A','CLEC9A','XCR1'],'pDC':['LILRA4','IL3RA','CLEC4C','TCF4'],'MatureMacro':['C1QA','C1QB','C1QC','CSF1R','CD163','MRC1','CD68','APOE','TREM2'],'Neutrophil':['CSF3R','FCGR3B','CXCR2','S100A12']})
for ds in ['OBRADOVIC','HCRN']:
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 a=ad.read_h5ad(P/f'inputs/{ds}_counts_qc.h5ad');a=a[a.obs.cluster.astype(str).isin([str(c) for c in rule[ds]])].copy();a.obs['parent_cluster']=a.obs.cluster.astype(str)
 norm=a.copy();sc.pp.normalize_total(norm,target_sum=1e4);sc.pp.log1p(norm);sc.pp.highly_variable_genes(norm,n_top_genes=2000,batch_key='patient',flavor='seurat');norm.var.loc[norm.var_names.isin(rule['exclude_hvg']),'highly_variable']=False
 z=norm[:,norm.var.highly_variable].copy();sc.pp.scale(z,max_value=10);sc.pp.pca(z,n_comps=30,random_state=305);sc.pp.neighbors(z,n_neighbors=15,n_pcs=30,random_state=305);sc.tl.leiden(z,resolution=.7,random_state=305,key_added='myeloid_cluster',flavor='igraph',n_iterations=2)
 a.obs['myeloid_cluster']=z.obs.myeloid_cluster;a.write_h5ad(P/f'inputs/{ds}_myeloid_refined.h5ad',compression='gzip');a.obs.to_csv(P/f'tables/{ds}_MYELOID_METADATA.tsv.gz',sep='\t')
 glob=np.asarray(norm.X.mean(axis=0)).ravel();rows=[];detail=[]
 for cl in a.obs.myeloid_cluster.cat.categories:
  ix=a.obs.myeloid_cluster.eq(cl).to_numpy();mean=np.asarray(norm.X[ix].mean(axis=0)).ravel();pct=np.asarray((a.X[ix]>0).mean(axis=0)).ravel();effect=mean-glob;effect[a.var_names.isin(rule['exclude_hvg'])]=-np.inf
  rr={'cluster':cl,'cells':int(ix.sum()),'patients':a.obs.loc[ix,'patient'].nunique(),'top_genes':','.join(a.var_names[np.argsort(effect)[-20:][::-1]])}
  for k,gs in markers.items():
   jj=np.flatnonzero(a.var_names.isin(gs));rr[k]=float(mean[jj].mean()) if len(jj) else np.nan
   for j in jj:detail.append(dict(cluster=cl,group=k,gene=a.var_names[j],mean=float(mean[j]),detection=float(pct[j])))
  rows.append(rr)
 pd.DataFrame(rows).to_csv(P/f'tables/{ds}_MYELOID_OVERVIEW.tsv',sep='\t',index=False);pd.DataFrame(detail).to_csv(P/f'tables/{ds}_MYELOID_MARKERS.tsv',sep='\t',index=False);print('REFINED',ds,a.shape,flush=True);del a,norm,z;gc.collect()
(P/'qc/MYELOID_REFINEMENT_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','association_tests_run':False},indent=2))
