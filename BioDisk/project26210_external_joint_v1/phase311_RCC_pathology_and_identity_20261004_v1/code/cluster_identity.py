from pathlib import Path
import json,os,sys,datetime,hashlib
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
P=Path(__file__).resolve().parents[1];q=P/'resources/proseg_region6_v1'
assert (P/'logs/proseg.exit').read_text().strip()=='0'
counts=pd.read_csv(q/'maxpost-counts.csv.gz')
meta=pd.read_csv(q/'cell-metadata.csv.gz',dtype={'original_cell_id':str})
assert len(counts)==len(meta) and np.array_equal(meta.cell.to_numpy(),np.arange(len(meta)))
assert meta.original_cell_id.is_unique
assert np.isfinite(counts.to_numpy()).all() and (counts.to_numpy()>=0).all()
assert np.equal(counts.to_numpy(),np.floor(counts.to_numpy())).all()
meta['counts']=counts.sum(axis=1).to_numpy();meta['genes']=(counts>0).sum(axis=1).to_numpy()
meta['technical_qc']=(meta.counts>=20)&(meta.genes>=5)&np.isfinite(meta[['centroid_x','centroid_y']].to_numpy()).all(axis=1)
meta['identity_cluster']='not_QC'
keep=meta.technical_qc.to_numpy();gene_names=counts.columns.tolist()
a=ad.AnnData(sparse.csr_matrix(counts.loc[keep].to_numpy(dtype=np.float32)),obs=meta.loc[keep].set_index('original_cell_id'),var=pd.DataFrame(index=gene_names))
a.layers['maxpost_counts']=a.X.copy();a.obsm['spatial']=a.obs[['centroid_x','centroid_y']].to_numpy()
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a)
identity_features=[g for g in gene_names if g not in ['TYMP','SPP1']]
z=a[:,identity_features].copy()
sc.pp.scale(z,max_value=10);sc.tl.pca(z,n_comps=30,svd_solver='arpack',random_state=26210)
sc.pp.neighbors(z,n_neighbors=15,n_pcs=30,random_state=26210)
sc.tl.leiden(z,resolution=0.6,random_state=26210,key_added='identity_cluster',flavor='igraph',directed=False,n_iterations=2)
a.obs['identity_cluster']=z.obs.identity_cluster;a.obsm['X_pca_identity']=z.obsm['X_pca']
a.uns['identity_parameters']={'target_excluded':['TYMP','SPP1'],'resolution':0.6,'seed':26210,'n_neighbors':15,'n_pcs':30}
a.write_h5ad(P/'resources/proseg_region6_identity.h5ad',compression='gzip')
meta.loc[keep,'identity_cluster']=z.obs.identity_cluster.astype(str).to_numpy()
meta.to_csv(P/'tables/PROSEG_CELL_IDENTITY_INPUT.tsv',sep='\t',index=False)
L=a.X.toarray();D=a.layers['maxpost_counts'].toarray()>0;all_mean=L.mean(axis=0);families=json.loads((P/'resources/SEGMENTATION_IDENTITY_PROTOCOL.json').read_text())['identity_plan']['markers']
summaries=[];long=[];panel=[]
for cl in sorted(a.obs.identity_cluster.unique(),key=int):
 ix=a.obs.identity_cluster.eq(cl).to_numpy();m=L[ix].mean(axis=0);d=D[ix].mean(axis=0);delta=m-all_mean
 order=[j for j in np.argsort(-delta) if gene_names[j] not in ['TYMP','SPP1']][:25]
 row={'cluster':cl,'cells':int(ix.sum()),'median_counts':float(np.median(a.obs.loc[ix,'counts'])),'top25_mean_difference_markers':';'.join(gene_names[j] for j in order)}
 for fam,gs in families.items():
  js=[gene_names.index(g) for g in gs if g in gene_names]
  row[fam+'_mean_log']=float(m[js].mean()) if js else None
  row[fam+'_detect2_fraction']=float((D[ix][:,js].sum(axis=1)>=2).mean()) if js else None
  for g in gs:
   if g in gene_names:
    j=gene_names.index(g);panel.append({'cluster':cl,'family':fam,'gene':g,'mean_log':float(m[j]),'fraction_detected':float(d[j])})
 for j,g in enumerate(gene_names):
  if g not in ['TYMP','SPP1']:long.append({'cluster':cl,'gene':g,'mean_log':float(m[j]),'difference_from_all_mean':float(delta[j]),'fraction_detected':float(d[j])})
 summaries.append(row)
pd.DataFrame(summaries).to_csv(P/'tables/IDENTITY_CLUSTER_REVIEW.tsv',sep='\t',index=False)
pd.DataFrame(panel).to_csv(P/'tables/IDENTITY_CANONICAL_MARKERS.tsv',sep='\t',index=False)
pd.DataFrame(long).to_csv(P/'tables/IDENTITY_ALL_NONTARGET_GENES.tsv.gz',sep='\t',index=False,compression='gzip')
rec={'status':'READY_FOR_MANUAL_IDENTITY_REVIEW','all_cells':len(meta),'QC_cells':int(keep.sum()),'observed_genes':len(gene_names),'clusters':len(summaries),'target_genes_excluded_from_identity':['TYMP','SPP1'],'no_biological_tests':True,'output':str(P/'resources/proseg_region6_identity.h5ad'),'software':{'python':sys.version,'scanpy':sc.__version__,'anndata':ad.__version__}}
(P/'qc/IDENTITY_CLUSTERING_COMPLETE.json').write_text(json.dumps(rec,indent=2))
print(json.dumps(rec,indent=2));print(pd.DataFrame(summaries)[['cluster','cells','top25_mean_difference_markers']].to_string(index=False),flush=True)
