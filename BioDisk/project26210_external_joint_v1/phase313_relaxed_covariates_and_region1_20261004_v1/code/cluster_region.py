from pathlib import Path
import json,sys,datetime,re
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy import sparse
P=Path(__file__).resolve().parents[1];region=sys.argv[1];assert region in ['region1','region6']
R=P/'resources'/region;R.mkdir(exist_ok=True)
if region=='region1':
 meta=pd.read_csv(R/'CELL_METADATA.tsv',sep='\t',dtype={'cell_id':str});gf=pd.read_csv(R/'FEATURE_AUDIT.tsv',sep='\t');A=sparse.load_npz(R/'assigned_counts.npz');N=sparse.load_npz(R/'nuclear_counts.npz')
else:
 O=P.parent/'phase310_RCC_mapped_spatial_inputs_20261004_v1'
 meta=pd.read_csv(O/'tables/PILOT_CELL_METADATA.tsv',sep='\t',dtype={'cell_id':str});gf=pd.read_csv(O/'tables/PILOT_FEATURE_AUDIT.tsv',sep='\t');A=sparse.load_npz(O/'resources/pilot_assigned_counts_qv20.npz');N=sparse.load_npz(O/'resources/pilot_nuclear_counts_qv20.npz')
keep=(meta.assigned_RNA>=20)&(meta.assigned_genes>=5)&(meta.nuclear_RNA>=5)&np.isfinite(meta[['nucleus_x_proxy_um','nucleus_y_proxy_um']]).all(axis=1)
genes=gf.gene.astype(str).to_numpy();gkeep=np.array([not bool(re.match(r'(?i)^(NegControl|NegativeControl|Blank|Unassigned|antisense)',g)) for g in genes])
a=ad.AnnData(A[keep.to_numpy()][:,gkeep].astype(np.float32),obs=meta.loc[keep].set_index('cell_id'),var=pd.DataFrame(index=genes[gkeep]))
a.layers['assigned_counts']=a.X.copy();a.layers['nuclear_counts']=N[keep.to_numpy()][:,gkeep].astype(np.float32);a.obsm['spatial']=a.obs[['nucleus_x_proxy_um','nucleus_y_proxy_um']].to_numpy();a.obs['region']='tumor1_'+region
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
z=a[:,~a.var_names.isin(['TYMP','SPP1'])].copy();sc.pp.scale(z,max_value=10);sc.tl.pca(z,n_comps=30,svd_solver='arpack',random_state=26210);sc.pp.neighbors(z,n_neighbors=15,n_pcs=30,random_state=26210);sc.tl.leiden(z,resolution=.6,random_state=26210,key_added='identity_cluster',flavor='igraph',directed=False,n_iterations=2)
a.obs['identity_cluster']=z.obs.identity_cluster;a.obsm['X_pca_identity']=z.obsm['X_pca'];a.uns['identity_method']={'exclude_targets':['TYMP','SPP1'],'resolution':.6,'seed':26210,'input':'original assigned QV20 counts; no Proseg gate'}
a.write_h5ad(R/'original_count_identity.h5ad',compression='gzip')
markers=json.loads((P.parent/'phase311_RCC_pathology_and_identity_20261004_v1/resources/SEGMENTATION_IDENTITY_PROTOCOL.json').read_text())['identity_plan']['markers']
L=a.X.toarray();D=a.layers['assigned_counts'].toarray()>0;allmean=L.mean(axis=0);gn=list(a.var_names);summ=[];canonical=[];allgenes=[]
for cl in sorted(a.obs.identity_cluster.unique(),key=int):
 ix=a.obs.identity_cluster.eq(cl).to_numpy();m=L[ix].mean(axis=0);det=D[ix].mean(axis=0);delta=m-allmean
 order=[j for j in np.argsort(-delta) if gn[j] not in ['TYMP','SPP1']][:25]
 row={'region':region,'cluster':cl,'cells':int(ix.sum()),'top25_nontarget_markers':';'.join(gn[j] for j in order)}
 for fam,gs in markers.items():
  js=[gn.index(g) for g in gs if g in gn];row[fam+'_mean_log']=float(m[js].mean()) if js else None
  for g in gs:
   if g in gn:
    j=gn.index(g);canonical.append({'region':region,'cluster':cl,'family':fam,'gene':g,'mean_log':float(m[j]),'fraction_detected':float(det[j])})
 for j,g in enumerate(gn):
  if g not in ['TYMP','SPP1']:allgenes.append({'cluster':cl,'gene':g,'mean_log':float(m[j]),'fraction_detected':float(det[j]),'delta':float(delta[j])})
 summ.append(row)
pd.DataFrame(summ).to_csv(P/'tables'/(region+'_CLUSTER_REVIEW.tsv'),sep='\t',index=False);pd.DataFrame(canonical).to_csv(P/'tables'/(region+'_CANONICAL_MARKERS.tsv'),sep='\t',index=False);pd.DataFrame(allgenes).to_csv(P/'tables'/(region+'_ALL_NONTARGET_MARKERS.tsv.gz'),sep='\t',index=False,compression='gzip')
rec={'status':'READY_FOR_MANUAL_REVIEW','region':region,'QC_cells':len(a),'genes_observed_noncontrol':a.n_vars,'clusters':len(summ),'target_excluded':True,'no_new_effects_yet':True}
(P/'qc'/(region+'_IDENTITY_COMPLETE.json')).write_text(json.dumps(rec,indent=2));print(json.dumps(rec,indent=2));print(pd.DataFrame(summ)[['cluster','cells','top25_nontarget_markers']].to_string(index=False),flush=True)
