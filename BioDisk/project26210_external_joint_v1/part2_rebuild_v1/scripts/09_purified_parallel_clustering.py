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
 import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy,re,gc
 from scipy import sparse
 from sklearn.metrics import adjusted_rand_score
 R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'03_macrophage_clustering/run_v3';O.mkdir();seed=262147;sc.settings.n_jobs=8
 def log(t):print(datetime.datetime.now().isoformat(),t,flush=True)
 def graph(emb,obs,k=50):
  a=ad.AnnData(X=sparse.csr_matrix((len(obs),0),dtype=np.float32),obs=obs.copy());a.obsm['X_analysis']=np.asarray(emb,np.float32);sc.pp.neighbors(a,n_neighbors=k,use_rep='X_analysis',random_state=seed);return a
 def mix(a):
  t=pd.crosstab(a.obs.cluster,a.obs.study);f=t.div(t.sum(axis=1),axis=0);e=-(f*np.log(f.replace(0,np.nan))).sum(axis=1)/np.log(5);return float(np.average(e,weights=t.sum(axis=1))),float(t.max(axis=1).sum()/len(a))
 spec=json.loads((R/'02_myeloid_qc/run_v2/STATUS.json').read_text());assert spec['status']=='PASS_PURIFIED_INPUT_REVIEW'
 plan={'started':datetime.datetime.now().isoformat(),'reason':'Repeat unsupervised fine clustering after explicit complexity and lineage purification; no SPP1-based selection','input':spec,'HVG':3000,'PCA':50,'Harmony_study_theta':2,'k':50,'PCs':30,'resolutions_descending':[.8,.6,.4,.2,.15,.1,.05],'seeds':[seed,42,31415],'gates':{'minimum_seed_ARI':.85,'subsample_regraph_ARI':.8,'small_cluster_fraction':.01,'min_nontrivial_clusters':2},'parallel_workers':3,'selection':'highest tested resolution passing unchanged stability gates; no marker inspection until selection complete'}
 (R/'00_plan/PURIFIED_CLUSTERING_PLAN_v3.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2))
 log('load purified macrophage raw counts')
 meta=pd.read_csv(R/'02_myeloid_qc/run_v2/PURIFIED_CELL_METADATA.tsv.gz',sep='\t',index_col='cell_id',low_memory=False);a=ad.read_h5ad(R/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad');a=a[meta.index].copy();a.obs=meta.copy();a.obs=a.obs.drop(columns=[c for c in a.obs if c=='cluster' or c.startswith('r0.') or c.startswith('r1.')]);assert a.n_obs==spec['retained_macrophages']
 countkeep=np.asarray((a.X>0).sum(axis=0)).ravel()>=20;legacy_mt={'COX1','COX2','COX3','ND1','ND2','ND3','ND4','ND4L','ND5','ND6','ATP6','ATP8','CYTB'};hb={'HBA1','HBA2','HBB','HBD','HBG1','HBG2','HBE1','HBM','HBQ1'}
 technical=np.array([g.startswith('MT-') or g in legacy_mt or g in hb or bool(re.fullmatch(r'RP[SL]\d+[A-Z]?(?:P\d+)?',g)) or g in {'RPLP0','RPLP1','RPLP2','RPS4Y1','RPS4Y2'} for g in a.var_names])
 pd.DataFrame({'gene':a.var_names,'at_least20_cells':countkeep,'excluded_technical_feature':technical}).to_csv(O/'FEATURE_FILTER.tsv',sep='\t',index=False)
 h=a[:,countkeep&~technical].copy();log('new study-aware HVGs');sc.pp.highly_variable_genes(h,flavor='seurat_v3',n_top_genes=3000,batch_key='study');h.var.to_csv(O/'HVG_ALL_CANDIDATE_STATISTICS.tsv',sep='\t');hv=h.var_names[h.var.highly_variable].tolist();assert len(hv)==3000
 total=np.asarray(a.X.sum(axis=1)).ravel();x=a[:,hv].X.astype(np.float32).multiply(10000/total[:,None]).tocsr();x.data=np.log1p(x.data);h=ad.AnnData(X=x,obs=a.obs.copy(),var=pd.DataFrame(index=hv));del a,x;gc.collect();log('new scaled PCA and Harmony');sc.pp.scale(h,max_value=10);sc.tl.pca(h,n_comps=50,svd_solver='arpack',random_state=seed);pc=h.obsm['X_pca'].copy();obs=h.obs.copy();pd.DataFrame(h.varm['PCs'],index=h.var_names).to_csv(O/'PCA_LOADINGS.tsv.gz',sep='\t',compression='gzip');np.save(O/'PCA_VARIANCE_RATIO.npy',h.uns['pca']['variance_ratio']);del h;gc.collect()
 hm=harmonypy.run_harmony(pc,obs,['study'],theta=2,nclust=50,max_iter_harmony=20,random_state=seed,ncores=8);hp=np.asarray(hm.Z_corr);hp=hp.T if hp.shape==pc.T.shape else hp;assert hp.shape==pc.shape
 log('primary and 80% study-stratified regraphs');z=graph(hp[:,:30],obs);rng=np.random.default_rng(seed+1);idx=np.sort(np.concatenate([rng.choice(np.flatnonzero(obs.study.astype(str).to_numpy()==s),size=int(.8*sum(obs.study.astype(str).to_numpy()==s)),replace=False) for s in sorted(obs.study.astype(str).unique())]));sub=graph(hp[idx,:30],obs.iloc[idx]);sparse.save_npz(O/'primary_adjacency.npz',z.obsp['connectivities']);sparse.save_npz(O/'subsample_adjacency.npz',sub.obsp['connectivities']);metrics=[];selected=None
 # Spawned seed workers are single-threaded; they do not use labels or expression markers.
 for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
 import numba
 numba.set_num_threads(4);sc.settings.n_jobs=4
 with concurrent.futures.ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn')) as pool:
  for r in plan['resolutions_descending']:
   log('parallel seeds resolution '+str(r));labs=list(pool.map(partition_task,[(str(O/'primary_adjacency.npz'),r,s) for s in plan['seeds']]));aris=[adjusted_rand_score(labs[i],labs[j]) for i in range(3) for j in range(i+1,3)];cnt=pd.Series(labs[0]).value_counts();small=float(cnt[cnt<100].sum()/len(z));sa=np.nan
   for sd,lab in zip(plan['seeds'],labs):z.obs['r'+str(r)+'_s'+str(sd)]=pd.Categorical(lab)
   if min(aris)>=.85 and small<=.01 and len(cnt)>=2:sl=pool.submit(partition_task,(str(O/'subsample_adjacency.npz'),r,seed)).result();sa=adjusted_rand_score(labs[0][idx],sl)
   passed=min(aris)>=.85 and sa>=.8 and small<=.01 and len(cnt)>=2;z.obs['cluster']=pd.Categorical(labs[0]);e,d=mix(z);metrics.append({'resolution':r,'clusters':len(cnt),'min_cluster_cells':int(cnt.min()),'small_cluster_cell_fraction':small,'minimum_between_seed_ARI':float(min(aris)),'mean_between_seed_ARI':float(np.mean(aris)),'subsample_regraph_ARI':sa,'weighted_study_entropy':e,'weighted_dominant_study_fraction':d,'structural_gate_pass':passed});pd.DataFrame(metrics).to_csv(O/'RESOLUTION_STABILITY.tsv',sep='\t',index=False)
   if passed:selected=r;break
  if selected is None:(O/'STATUS.json').write_text(json.dumps({'status':'STABILITY_UNRESOLVED','metrics':metrics}));raise RuntimeError('No purified resolution passed unchanged nontrivial stability gates')
  log('selected '+str(selected)+' clusters '+str(z.obs.cluster.nunique()));sens_jobs=[];sense_meta=[]
  for name,emb,ids,k in [('original_k30',hp[:,:30],np.arange(len(z)),30),('uncorrected_pca30',pc[:,:30],np.arange(len(z)),50),('harmony50',hp,np.arange(len(z)),50),('drop_lineage_review_flags',hp[~obs.lineage_review_flag.to_numpy(bool),:30],np.flatnonzero(~obs.lineage_review_flag.to_numpy(bool)),50)]:
   log('sensitivity graph '+name);v=graph(emb,obs.iloc[ids],k);path=O/(name+'_adjacency.npz');sparse.save_npz(path,v.obsp['connectivities']);sens_jobs.append(pool.submit(partition_task,(str(path),selected,seed)));sense_meta.append((name,ids,v))
  sens=[]
  for future,(name,ids,v) in zip(sens_jobs,sense_meta):
   lab=future.result();v.obs['cluster']=pd.Categorical(lab);e,d=mix(v);sens.append({'sensitivity':name,'n_cells':len(ids),'clusters':len(set(lab)),'ARI_to_primary':float(adjusted_rand_score(z.obs.cluster.astype(str).to_numpy()[ids],lab)),'weighted_study_entropy':e,'weighted_dominant_study_fraction':d});pd.DataFrame({'cell_id':obs.index[ids],'cluster':lab}).to_csv(O/(name+'_labels.tsv.gz'),sep='\t',index=False,compression='gzip')
 pd.DataFrame(sens).to_csv(O/'GRAPH_SENSITIVITY.tsv',sep='\t',index=False);log('UMAP and final checkpoint');sc.tl.umap(z,min_dist=.35,random_state=seed);z.obsm['X_pca_uncorrected']=pc;z.obsm['X_harmony50']=hp;z.uns['selected_resolution']=selected;z.uns['selection_status']='STRUCTURAL_GATE_PASS';z.uns['neighbors_k']=50;z.obs.to_csv(O/'CELL_CLUSTER_ASSIGNMENTS.tsv.gz',sep='\t',compression='gzip',index_label='cell_id');pd.DataFrame(z.obsm['X_umap'],index=z.obs_names,columns=['UMAP1','UMAP2']).to_csv(O/'UMAP_COORDINATES.tsv.gz',sep='\t',compression='gzip',index_label='cell_id');pd.DataFrame({'gene':hv}).to_csv(O/'HVG_SELECTED.tsv',sep='\t',index=False)
 for field,name in [('study','CLUSTER_BY_STUDY.tsv'),('curated_patient_id','CLUSTER_BY_PATIENT.tsv'),('lineage_review_flag','CLUSTER_BY_REVIEW_FLAG.tsv')]:pd.crosstab(z.obs.cluster,z.obs[field]).to_csv(O/name,sep='\t')
 z.write_h5ad(O/'CLUSTERING_GRAPH.h5ad',compression='gzip');s={'status':'CLUSTERING_COMPLETE_PENDING_REVIEW','selection_status':'STRUCTURAL_GATE_PASS','selected_resolution':selected,'n_cells':len(z),'clusters':z.obs.cluster.nunique(),'neighbors_k':50,'resolution_metrics':metrics,'sensitivities':sens,'input_purification':str(R/'02_myeloid_qc/run_v2/STATUS.json'),'SPP1_threshold_used':False};(O/'STATUS.json').write_text(json.dumps(s,indent=2));print(json.dumps(s),flush=True)
if __name__=='__main__':main()
