from pathlib import Path
import json,datetime,gc,shutil
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');PRE=R/'03_macrophage_clustering/run_v1';O=R/'03_macrophage_clustering/run_v2';O.mkdir();sc.settings.n_jobs=8;seed=262147
p=ad.read_h5ad(PRE/'CLUSTERING_GRAPH.h5ad');obs=p.obs.copy();obs=obs.drop(columns=[c for c in obs if c.startswith('r0.') or c.startswith('r1.') or c=='cluster']);pc=p.obsm['X_pca_uncorrected'];hp=p.obsm['X_harmony50'];z=None
for name in ['FEATURE_FILTER.tsv','HVG_ALL_CANDIDATE_STATISTICS.tsv','HVG_SELECTED.tsv','PCA_LOADINGS.tsv.gz','PCA_VARIANCE_RATIO.npy','HARMONY_OBJECTIVE.json']:
 if (PRE/name).exists():shutil.copy2(PRE/name,O/name)
def log(t):print(datetime.datetime.now().isoformat(),t,flush=True)
def graph(emb,meta,k=50):
 a=ad.AnnData(X=sparse.csr_matrix((len(meta),0),dtype=np.float32),obs=meta.copy());a.obsm['X_analysis']=np.asarray(emb,np.float32);sc.pp.neighbors(a,n_neighbors=k,use_rep='X_analysis',random_state=seed,method='umap');return a
def cl(a,r,s,k):sc.tl.leiden(a,resolution=r,random_state=s,key_added=k,flavor='leidenalg',n_iterations=-1,directed=False);return a.obs[k].astype(str).to_numpy()
def mix(a):
 t=pd.crosstab(a.obs.cluster,a.obs.study);f=t.div(t.sum(axis=1),axis=0);e=-(f*np.log(f.replace(0,np.nan))).sum(axis=1)/np.log(5);return float(np.average(e,weights=t.sum(axis=1))),float(t.max(axis=1).sum()/len(a))
log('k50 graph and study-stratified 80% regraph')
z=graph(hp[:,:30],obs);rng=np.random.default_rng(seed+1);idx=np.sort(np.concatenate([rng.choice(np.flatnonzero(obs.study.astype(str).to_numpy()==s),size=int(.8*sum(obs.study.astype(str).to_numpy()==s)),replace=False) for s in sorted(obs.study.astype(str).unique())]));sub=graph(hp[idx,:30],obs.iloc[idx]);metrics=[];selected=None
for r in [.2,.15,.1,.05]:
 log('recovery candidate resolution '+str(r));labs=[];aris=[];failed_early=False
 for sd in [seed,42,31415]:
  lab=cl(z,r,sd,'r'+str(r)+'_s'+str(sd));labs.append(lab)
  if len(labs)>=2:
   aris=[adjusted_rand_score(labs[i],labs[j]) for i in range(len(labs)) for j in range(i+1,len(labs))]
   if min(aris)<.85:failed_early=True;break
 counts=pd.Series(labs[0]).value_counts();small=float(counts[counts<100].sum()/len(z));sa=np.nan
 if not failed_early and small<=.01:
  sl=cl(sub,r,seed,'r'+str(r));sa=adjusted_rand_score(labs[0][idx],sl)
 passed=(not failed_early) and min(aris)>=.85 and sa>=.8 and small<=.01
 z.obs['cluster']=pd.Categorical(labs[0]);e,d=mix(z)
 metrics.append({'resolution':r,'clusters':len(counts),'min_cluster_cells':int(counts.min()),'small_cluster_cell_fraction':small,'minimum_between_seed_ARI':float(min(aris)),'mean_between_seed_ARI':float(np.mean(aris)),'subsample_regraph_ARI':sa,'weighted_study_entropy':e,'weighted_dominant_study_fraction':d,'seed_runs_completed':len(labs),'structural_gate_pass':passed,'early_failure_skip':failed_early})
 pd.DataFrame(metrics).to_csv(O/'RESOLUTION_STABILITY.tsv',sep='\t',index=False)
 if passed:selected=r;break
if selected is None:
 (O/'STATUS.json').write_text(json.dumps({'status':'STABILITY_UNRESOLVED','metrics':metrics},indent=2));raise RuntimeError('Recovery grid did not pass unchanged gates; do not promote stable fine clusters.')
log('selected stable resolution '+str(selected));sens=[]
for label,emb,ids,k in [('original_k30',hp[:,:30],np.arange(len(obs)),30),('uncorrected_pca30',pc[:,:30],np.arange(len(obs)),50),('harmony50',hp,np.arange(len(obs)),50),('drop_lineage_review_flags',hp[~obs.lineage_review_flag.to_numpy(),:30],np.flatnonzero(~obs.lineage_review_flag.to_numpy()),50)]:
 log('sensitivity '+label);a=graph(emb,obs.iloc[ids],k);ll=cl(a,selected,seed,'cluster');e,d=mix(a);sens.append({'sensitivity':label,'n_cells':len(ids),'clusters':len(set(ll)),'ARI_to_primary':float(adjusted_rand_score(z.obs.cluster.astype(str).to_numpy()[ids],ll)),'weighted_study_entropy':e,'weighted_dominant_study_fraction':d});pd.DataFrame({'cell_id':obs.index[ids],'cluster':ll}).to_csv(O/(label+'_labels.tsv.gz'),sep='\t',index=False,compression='gzip');del a;gc.collect()
pd.DataFrame(sens).to_csv(O/'GRAPH_SENSITIVITY.tsv',sep='\t',index=False);log('UMAP and graph checkpoint');sc.tl.umap(z,min_dist=.35,random_state=seed);z.obsm['X_pca_uncorrected']=pc;z.obsm['X_harmony50']=hp;z.uns['selected_resolution']=selected;z.uns['selection_status']='STRUCTURAL_GATE_PASS';z.uns['neighbors_k']=50;z.uns['expression_source']=str(R/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad')
z.obs.to_csv(O/'CELL_CLUSTER_ASSIGNMENTS.tsv.gz',sep='\t',compression='gzip',index_label='cell_id');pd.DataFrame(z.obsm['X_umap'],index=z.obs_names,columns=['UMAP1','UMAP2']).to_csv(O/'UMAP_COORDINATES.tsv.gz',sep='\t',compression='gzip',index_label='cell_id')
for field,name in [('study','CLUSTER_BY_STUDY.tsv'),('curated_patient_id','CLUSTER_BY_PATIENT.tsv'),('lineage_review_flag','CLUSTER_BY_REVIEW_FLAG.tsv')]:pd.crosstab(z.obs.cluster,z.obs[field]).to_csv(O/name,sep='\t')
f=O/'CLUSTERING_GRAPH.h5ad';z.write_h5ad(str(f)+'.part',compression='gzip');Path(str(f)+'.part').replace(f)
status={'status':'CLUSTERING_COMPLETE_PENDING_REVIEW','selection_status':'STRUCTURAL_GATE_PASS','selected_resolution':selected,'n_cells':len(z),'clusters':z.obs.cluster.nunique(),'neighbors_k':50,'legacy_or_spp1_labels_used':False,'resolution_metrics':metrics,'sensitivities':sens,'initial_failed_grid_preserved':str(PRE)};(O/'STATUS.json').write_text(json.dumps(status,indent=2));print(json.dumps(status),flush=True)
