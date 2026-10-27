from common import *
import h5py,scanpy as sc,anndata as ad,harmonypy,importlib.metadata as md
from scipy import sparse
def main():
 gate('stromal_start')
 audit=R/'phase286_CAF_identity_boundary_audit_20260929_v1'
 m=read(audit/'inputs/CELL_IDENTITY_AUDIT.tsv.gz')
 wanted=m.cell_type_broad_v1.isin(['Fibroblast','Pericyte/SMC'])|m.fibro_candidate
 take=wanted&m.qc_pass;rows=np.flatnonzero(take.to_numpy());obs=m.loc[take].copy().reset_index(drop=True)
 m.loc[wanted&~m.qc_pass].to_csv(P/'tables/EXCLUDED_QC_CELLS.tsv.gz',sep='\t',index=False)
 obs['original_row']=rows;obs.to_csv(P/'inputs/STROMAL_INPUT_METADATA.tsv.gz',sep='\t',index=False)
 print('INPUT',len(obs),obs.curated_patient_id.nunique(),flush=True)
 chunks=[]
 with h5py.File(H,'r') as h:
  genes=h['var/gene'].asstr()[:];ids=h['obs/cell_id'].asstr()[:];assert ids[rows].tolist()==obs.cell_id.tolist()
  old=h['obsm/X_scVI'][:][rows];x=h['layers/counts'];ptr=x['indptr'][:]
  for start in range(0,len(m),4096):
   end=min(start+4096,len(m));local=np.flatnonzero(take.iloc[start:end].to_numpy())
   if len(local)==0:continue
   lo,hi=int(ptr[start]),int(ptr[end])
   a=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(genes)))
   chunks.append(a[local])
 raw=sparse.vstack(chunks,format='csr');assert np.array_equal(np.asarray(raw.sum(1)).ravel(),obs.total_counts_raw)
 sparse.save_npz(P/'inputs/STROMAL_COUNTS.npz',raw)
 np.save(P/'inputs/GENES.npy',genes.astype(str))
 np.save(P/'inputs/ORIGINAL_SCVI.npy',old)
 norm=raw.astype(np.float32).multiply((1e4/obs.total_counts_raw.to_numpy())[:,None]).tocsr()
 log=norm.copy();log.data=np.log1p(log.data)
 o=obs[['cell_id','study','curated_patient_id','sample_id','source_tissue','cell_type_broad_v1']].copy();o.index=o.cell_id
 a=ad.AnnData(log,obs=o,var=pd.DataFrame(index=genes))
 sc.pp.highly_variable_genes(a,n_top_genes=3000,flavor='seurat',batch_key='study')
 allowed=np.array([g not in ['TYMP','SPP1'] and not g.startswith(('MT-','RPL','RPS')) and g not in ['HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ'] for g in genes])
 hv=a.var.highly_variable.to_numpy()&allowed
 a.var.to_csv(P/'tables/HVG_SELECTION.tsv',sep='\t')
 v=a[:,hv].copy();sc.pp.scale(v,max_value=10);sc.tl.pca(v,n_comps=30,svd_solver='arpack',random_state=287)
 print('PCA',v.shape,flush=True)
 hm=harmonypy.run_harmony(v.obsm['X_pca'],v.obs,'study',max_iter_harmony=20,random_state=287,verbose=True)
 zz=np.asarray(hm.Z_corr)
 if zz.shape!=v.obsm['X_pca'].shape:zz=zz.T
 assert zz.shape==v.obsm['X_pca'].shape
 v.obsm['X_harmony']=zz
 sc.pp.neighbors(v,n_neighbors=15,use_rep='X_harmony',random_state=287)
 for key,res in [('cluster_primary',.4),('cluster_sensitive',.8)]:
  sc.tl.leiden(v,resolution=res,key_added=key,random_state=287,flavor='igraph',n_iterations=2)
  obs[key]=v.obs[key].astype(str).to_numpy()
 olda=ad.AnnData(sparse.csr_matrix((len(obs),1)),obs=o.copy());olda.obsm['X_old']=old
 sc.pp.neighbors(olda,n_neighbors=15,use_rep='X_old',random_state=287)
 sc.tl.leiden(olda,resolution=.4,key_added='old_scvi_cluster',random_state=287,flavor='igraph',n_iterations=2)
 obs['cluster_original_scvi']=olda.obs.old_scvi_cluster.astype(str).to_numpy()
 obs.to_csv(P/'tables/STROMAL_CLUSTER_MEMBERSHIP.tsv.gz',sep='\t',index=False)
 np.savez_compressed(P/'inputs/REDUCED_COORDINATES.npz',cell_id=obs.cell_id.to_numpy(str),pca=v.obsm['X_pca'],harmony=zz)
 panels=json.loads((audit/'resources/IDENTITY_PANELS_FROZEN.json').read_text())
 union=sorted(set(sum(list(panels.values()),[])));idx={g:i for i,g in enumerate(genes)}
 agg=[];markers=[];summ=[]
 for mode in ['cluster_primary','cluster_sensitive','cluster_original_scvi']:
  for cluster,gg in obs.groupby(mode):
   ii=gg.index.to_numpy();other=np.setdiff1d(np.arange(len(obs)),ii)
   mean=np.asarray(log[ii].mean(0)).ravel();pct=np.asarray((raw[ii]>0).mean(0)).ravel()
   bg=np.asarray(log[other].mean(0)).ravel();pctbg=np.asarray((raw[other]>0).mean(0)).ravel()
   diff=mean-bg;order=np.argsort(-diff)
   for j in [j for j in order if allowed[j] and pct[j]>=.1][:60]:
    markers.append(dict(mode=mode,cluster=cluster,gene=genes[j],mean_logCP10K=float(mean[j]),mean_other=float(bg[j]),difference=float(diff[j]),detection=float(pct[j]),detection_other=float(pctbg[j])))
   for g in union:
    if g in idx:agg.append(dict(mode=mode,cluster=cluster,gene=g,cells=len(ii),mean_logCP10K=float(mean[idx[g]]),detection=float(pct[idx[g]])))
   pc=gg.curated_patient_id.value_counts();ss=gg.study.value_counts()
   summ.append(dict(mode=mode,cluster=cluster,cells=len(ii),patients=len(pc),studies=len(ss),largest_patient_fraction=pc.iloc[0]/len(ii),largest_study_fraction=ss.iloc[0]/len(ii),tumor_cells=int(gg.tumor.sum()),original_labels=json.dumps(gg.cell_type_broad_v1.value_counts().to_dict()),study_counts=json.dumps(ss.to_dict()),high_priority_review=int(gg.high_priority_review.sum())))
 tab(agg,'CLUSTER_LINEAGE_MARKERS.tsv');tab(markers,'CLUSTER_TOP_MARKERS_DESCRIPTIVE.tsv');tab(summ,'CLUSTER_COMPOSITION.tsv')
 pat=[]
 for (cluster,study,patient,tissue),g in obs.groupby(['cluster_primary','study','curated_patient_id','source_tissue']):
  for gene in union:
   if gene in idx:pat.append(dict(cluster=cluster,study=study,patient=patient,tissue=tissue,cells=len(g),gene=gene,mean_logCP10K=float(np.asarray(log[g.index,idx[gene]].mean()).item()),detection=float(np.asarray((raw[g.index,idx[gene]]>0).mean()).item())))
 tab(pat,'PATIENT_CLUSTER_LINEAGE_MARKERS.tsv')
 for other in ['cluster_sensitive','cluster_original_scvi']:
  pd.crosstab(obs.cluster_primary,obs[other]).to_csv(P/'tables'/f'PRIMARY_VS_{other}.tsv',sep='\t')
 dump(P/'qc/STROMAL_CLUSTER_COMPLETE.json',dict(status='COMPUTED',cells=len(obs),patients=obs.curated_patient_id.nunique(),hvg=int(hv.sum()),source_path=str(H),source_bytes=H.stat().st_size,source_mtime_ns=H.stat().st_mtime_ns,counts_match=True,excluded_TYMP_SPP1_from_features=True,original_labels_unchanged=True,seed=287,versions={x:md.version(x) for x in ['scanpy','anndata','harmonypy','numpy','scipy','pandas','igraph']}))
 print('STROMAL_CLUSTER_COMPLETE',obs.cluster_primary.value_counts().to_dict(),flush=True)
if __name__=='__main__':main()
