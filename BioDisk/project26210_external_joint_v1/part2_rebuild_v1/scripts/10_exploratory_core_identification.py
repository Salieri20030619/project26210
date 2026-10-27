from pathlib import Path
import os,json,datetime,concurrent.futures,multiprocessing

def partition_task(args):
 path,res,seed=args
 import numpy as np
 from scipy import sparse
 import igraph as ig,leidenalg
 A=sparse.load_npz(path).tocoo();g=ig.Graph(n=A.shape[0],edges=list(zip(A.row.tolist(),A.col.tolist())),directed=False);g.es['weight']=A.data.astype(float).tolist()
 part=leidenalg.find_partition(g,leidenalg.RBConfigurationVertexPartition,weights=g.es['weight'],resolution_parameter=res,seed=seed,n_iterations=-1)
 return np.asarray(part.membership,dtype=str)

def main():
 import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy,shutil,gc
 from scipy import sparse
 from sklearn.metrics import adjusted_rand_score
 R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');PRE=R/'03_macrophage_clustering/run_v3';O=R/'03_macrophage_clustering/run_v4_exploratory';O.mkdir();seed=262147;sc.settings.n_jobs=8
 def log(t):print(datetime.datetime.now().isoformat(),t,flush=True)
 meta=pd.read_csv(R/'02_myeloid_qc/run_v2/PURIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col='cell_id',low_memory=False);meta=meta.drop(columns=[c for c in meta if c=='cluster' or c.startswith('r0.') or c.startswith('r1.')]);met=pd.read_csv(PRE/'RESOLUTION_STABILITY.tsv',sep='\t');selected=float(met[met.clusters>=2].sort_values(['minimum_between_seed_ARI','mean_between_seed_ARI','resolution'],ascending=[False,False,True]).iloc[0].resolution);assert selected==.8
 for name in ['FEATURE_FILTER.tsv','HVG_ALL_CANDIDATE_STATISTICS.tsv','PCA_LOADINGS.tsv.gz','PCA_VARIANCE_RATIO.npy','RESOLUTION_STABILITY.tsv']:shutil.copy2(PRE/name,O/name)
 (O/'STATUS.json').write_text(json.dumps({'status':'EXPLORATORY_BUILD_RUNNING','global_gate_pass':False,'selected_resolution':selected}))
 rng=np.random.default_rng(seed+1);idx=np.sort(np.concatenate([rng.choice(np.flatnonzero(meta.study.astype(str).to_numpy()==s),size=int(.8*sum(meta.study.astype(str).to_numpy()==s)),replace=False) for s in sorted(meta.study.astype(str).unique())]))
 log('reproduce three seed partitions and independent subsample partition on frozen adjacency')
 for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
 with concurrent.futures.ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn')) as pool:
  jobs=[pool.submit(partition_task,(str(PRE/'primary_adjacency.npz'),selected,s)) for s in [seed,42,31415]];sj=pool.submit(partition_task,(str(PRE/'subsample_adjacency.npz'),selected,seed));labs=[f.result() for f in jobs];sl=sj.result()
 aris=[float(adjusted_rand_score(labs[i],labs[j])) for i in range(3) for j in range(i+1,3)];assert abs(min(aris)-float(met[met.resolution==selected].minimum_between_seed_ARI.iloc[0]))<1e-8
 core=[];matches=[]
 for c in sorted(set(labs[0]),key=int):
  jj=[]
  for tag,aa,bb in [('seed42',labs[0],labs[1]),('seed31415',labs[0],labs[2]),('subsample80',labs[0][idx],sl)]:
   tab=pd.crosstab(pd.Series(aa),pd.Series(bb));row=tab.loc[c];jac=row/(row.sum()+tab.sum(axis=0)-row);best=jac.idxmax();score=float(jac[best]);matches.append({'cluster':c,'comparison':tag,'best_matching_cluster':str(best),'Jaccard':score,'intersection':int(row[best]),'reference_cells_compared':int(row.sum())});jj.append(score)
  core.append({'cluster':c,'minimum_seed_Jaccard':min(jj[:2]),'subsample_Jaccard':jj[2],'cluster_core_reporting_gate':min(jj[:2])>=.75 and jj[2]>=.70})
 pd.DataFrame(matches).to_csv(O/'CLUSTER_MATCHES.tsv',sep='\t',index=False);pd.DataFrame(core).to_csv(O/'CLUSTER_CORE_STABILITY.tsv',sep='\t',index=False)
 # Recover PCA display basis using saved loadings and the same raw-count normalization/scaling.
 log('recover display basis; clustering graph remains frozen')
 a=ad.read_h5ad(R/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad');a=a[meta.index].copy();load=pd.read_csv(PRE/'PCA_LOADINGS.tsv.gz',sep='\t',index_col=0);hv=load.index.astype(str).tolist();total=np.asarray(a.X.sum(axis=1)).ravel();xx=a[:,hv].X.astype(np.float32).multiply(10000/total[:,None]).tocsr();xx.data=np.log1p(xx.data);h=ad.AnnData(X=xx);sc.pp.scale(h,max_value=10);centered=h.X-h.X.mean(axis=0);pc=centered@load.to_numpy(np.float32);vr=np.var(pc,axis=0,ddof=1)/np.var(centered,axis=0,ddof=1).sum();old=np.load(PRE/'PCA_VARIANCE_RATIO.npy');err=float(np.max(abs(vr-old)));assert err<1e-4;del a,h,xx,centered;gc.collect();hm=harmonypy.run_harmony(pc,meta,['study'],theta=2,nclust=50,max_iter_harmony=20,random_state=seed,ncores=8);hp=np.asarray(hm.Z_corr);hp=hp.T if hp.shape==pc.T.shape else hp
 z=ad.AnnData(X=sparse.csr_matrix((len(meta),0),dtype=np.float32),obs=meta.copy());z.obsm['X_analysis']=hp[:,:30];z.obsm['X_pca_uncorrected']=pc;z.obsm['X_harmony50']=hp;sc.pp.neighbors(z,n_neighbors=50,use_rep='X_analysis',random_state=seed);fresh=z.obsp['connectivities'];frozen=sparse.load_npz(PRE/'primary_adjacency.npz');delta=fresh-frozen;relative=float(sparse.linalg.norm(delta)/sparse.linalg.norm(frozen));z.obsp['connectivities']=frozen;z.obs['cluster']=pd.Categorical(labs[0])
 for s,l in zip([seed,42,31415],labs):z.obs['r'+str(selected)+'_s'+str(s)]=pd.Categorical(l)
 z.uns['selected_resolution']=selected;z.uns['selection_status']='EXPLORATORY_NOT_GLOBALLY_STABLE';z.uns['neighbors_k']=50
 log('UMAP from frozen purified graph');sc.tl.umap(z,min_dist=.35,random_state=seed)
 z.obs.to_csv(O/'CELL_CLUSTER_ASSIGNMENTS.tsv.gz',sep='\t',compression='gzip',index_label='cell_id');pd.DataFrame(z.obsm['X_umap'],index=z.obs_names,columns=['UMAP1','UMAP2']).to_csv(O/'UMAP_COORDINATES.tsv.gz',sep='\t',compression='gzip',index_label='cell_id');pd.DataFrame({'gene':hv}).to_csv(O/'HVG_SELECTED.tsv',sep='\t',index=False);z.write_h5ad(O/'CLUSTERING_GRAPH.h5ad',compression='gzip')
 for field,name in [('study','CLUSTER_BY_STUDY.tsv'),('curated_patient_id','CLUSTER_BY_PATIENT.tsv'),('lineage_review_flag','CLUSTER_BY_REVIEW_FLAG.tsv')]:pd.crosstab(z.obs.cluster,z.obs[field]).to_csv(O/name,sep='\t')
 status={'status':'READY_FOR_EXPLORATORY_IDENTIFICATION_REVIEW','selection_status':'EXPLORATORY_NOT_GLOBALLY_STABLE','global_gate_pass':False,'selected_resolution':selected,'n_cells':len(z),'clusters':len(set(labs[0])),'minimum_seed_ARI':min(aris),'subsample_regraph_ARI':float(adjusted_rand_score(labs[0][idx],sl)),'clusters_with_core_support':int(pd.DataFrame(core).cluster_core_reporting_gate.sum()),'PCA_variance_ratio_reconstruction_max_error':err,'display_reconstructed_graph_relative_difference':relative,'actual_clustering_uses_frozen_graph':True,'core_stability':core};(O/'STATUS.json').write_text(json.dumps(status,indent=2));print(json.dumps(status),flush=True)
if __name__=='__main__':main()
