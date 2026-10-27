from pathlib import Path
import os,json,datetime,hashlib,re,gc
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'03_macrophage_clustering/run_v1';O.mkdir();plan=json.loads((R/'00_plan/CLUSTERING_AND_IDENTIFICATION_PLAN_v1.json').read_text());seed=plan['seed'];sc.settings.n_jobs=8
assert json.loads((R/'01_input_audit/run_v1/STATUS.json').read_text())['status']=='PASS_INPUT_SCOPE_AND_LINEAGE_REVIEW'
def log(t):print(datetime.datetime.now().isoformat(),t,flush=True)
def graph(emb,obs):
 z=ad.AnnData(X=sparse.csr_matrix((len(obs),0),dtype=np.float32),obs=obs.copy());z.obsm['X_analysis']=np.asarray(emb,dtype=np.float32);sc.pp.neighbors(z,n_neighbors=30,use_rep='X_analysis',random_state=seed,method='umap');return z
def cl(z,r,sd,key):sc.tl.leiden(z,resolution=r,random_state=sd,key_added=key,flavor='leidenalg',n_iterations=-1,directed=False);return z.obs[key].astype(str).to_numpy()
def mix(z,key):
 labels=z.obs[key].astype(str);study=z.obs.study.astype(str);tab=pd.crosstab(labels,study);fra=tab.div(tab.sum(axis=1),axis=0);ent=-(fra*np.log(fra.replace(0,np.nan))).sum(axis=1)/np.log(study.nunique());return float(np.average(ent,weights=tab.sum(axis=1))),float(tab.max(axis=1).sum()/len(labels))
log('load frozen myeloid raw counts and retain macrophages only')
a=ad.read_h5ad(R/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad');a=a[a.obs.cell_type_broad_v1.astype(str)=='Macrophage'].copy();assert a.n_obs==42546
qc=pd.read_csv(R/'02_myeloid_qc/run_v1/cell_lineage_and_qc_audit.tsv.gz',sep='\t').set_index('cell_id');a.obs['lineage_review_flag']=qc.reindex(a.obs_names).mixed_lineage_review_flag.to_numpy(bool)
keep=np.asarray((a.X>0).sum(axis=0)).ravel()>=20
genes=a.var_names.astype(str);legacy_mt={'COX1','COX2','COX3','ND1','ND2','ND3','ND4','ND4L','ND5','ND6','ATP6','ATP8','CYTB'};hb={'HBA1','HBA2','HBB','HBD','HBG1','HBG2','HBE1','HBM','HBQ1'}
technical=np.array([g.startswith('MT-') or g in legacy_mt or g in hb or bool(re.fullmatch(r'RP[SL]\d+[A-Z]?(?:P\d+)?',g)) or g in {'RPLP0','RPLP1','RPLP2','RPS4Y1','RPS4Y2'} for g in genes])
featurefilter=pd.DataFrame({'gene':genes,'at_least20_cells':keep,'excluded_technical_feature':technical});featurefilter.to_csv(O/'FEATURE_FILTER.tsv',sep='\t',index=False)
log('batch-aware HVG selection on integer raw counts')
h=a[:,keep&~technical].copy();sc.pp.highly_variable_genes(h,flavor='seurat_v3',n_top_genes=3000,batch_key='study',span=0.3)
h.var.to_csv(O/'HVG_ALL_CANDIDATE_STATISTICS.tsv',sep='\t');hvg=h.var_names[h.var.highly_variable].tolist();assert len(hvg)==3000
# Normalize by total library RNA, not by only HVGs.
tot=np.asarray(a.X.sum(axis=1)).ravel();x=a[:,hvg].X.astype(np.float32).multiply(10000/tot[:,None]).tocsr();x.data=np.log1p(x.data)
h=ad.AnnData(X=x,obs=a.obs.copy(),var=pd.DataFrame(index=pd.Index(hvg)));del a,x;gc.collect()
log('scaled HVG PCA')
sc.pp.scale(h,max_value=10);sc.tl.pca(h,n_comps=50,svd_solver='arpack',random_state=seed)
pc=h.obsm['X_pca'].copy();pd.DataFrame(h.varm['PCs'],index=h.var_names).to_csv(O/'PCA_LOADINGS.tsv.gz',sep='\t',compression='gzip');np.save(O/'PCA_VARIANCE_RATIO.npy',h.uns['pca']['variance_ratio']);obs=h.obs.copy();del h;gc.collect()
log('Harmony correction by study; no subtype labels')
ho=harmonypy.run_harmony(pc,obs,['study'],theta=2,nclust=50,max_iter_harmony=20,random_state=seed,ncores=8,verbose=True)
hpc=np.asarray(ho.Z_corr)
if hpc.shape==pc.T.shape:hpc=hpc.T
assert hpc.shape==pc.shape and np.isfinite(hpc).all()
if hasattr(ho,'objective_harmony'):(O/'HARMONY_OBJECTIVE.json').write_text(json.dumps([float(v) for v in ho.objective_harmony]))
log('primary 30-PC Harmony neighborhood graph')
z=graph(hpc[:,:30],obs);z.obsm['X_pca_uncorrected']=pc;z.obsm['X_harmony50']=hpc.astype(np.float32)
rng=np.random.default_rng(seed+1);subidx=np.sort(np.concatenate([rng.choice(np.flatnonzero(obs.study.astype(str).to_numpy()==s),size=int(.8*sum(obs.study.astype(str).to_numpy()==s)),replace=False) for s in sorted(obs.study.astype(str).unique())]));log('80% study-stratified subsample regraph')
sub=graph(hpc[subidx,:30],obs.iloc[subidx]);metrics=[]
for r in plan['resolutions']:
 log('Leiden candidate resolution '+str(r));labs=[]
 for sd in plan['leiden_seeds']:labs.append(cl(z,r,sd,'r'+str(r)+'_s'+str(sd)))
 sublab=cl(sub,r,seed,'r'+str(r));counts=pd.Series(labs[0]).value_counts();aris=[adjusted_rand_score(labs[i],labs[j]) for i in range(3) for j in range(i+1,3)];sa=adjusted_rand_score(labs[0][subidx],sublab);small=float(counts[counts<100].sum()/len(obs));entropy,dom=mix(z,'r'+str(r)+'_s'+str(seed))
 metrics.append({'resolution':r,'clusters':len(counts),'min_cluster_cells':int(counts.min()),'small_cluster_cell_fraction':small,'minimum_between_seed_ARI':min(aris),'mean_between_seed_ARI':float(np.mean(aris)),'subsample_regraph_ARI':sa,'weighted_study_entropy':entropy,'weighted_dominant_study_fraction':dom,'structural_gate_pass':min(aris)>=.85 and sa>=.8 and small<=.01})
m=pd.DataFrame(metrics);m.to_csv(O/'RESOLUTION_STABILITY.tsv',sep='\t',index=False)
valid=m[m.structural_gate_pass]
if len(valid):selected=float(valid.resolution.max());selection_status='STRUCTURAL_GATE_PASS'
else:selected=float(m.sort_values(['minimum_between_seed_ARI','subsample_regraph_ARI','resolution'],ascending=[False,False,True]).iloc[0].resolution);selection_status='PROVISIONAL_NO_RESOLUTION_PASSED_ALL_GATES'
z.obs['cluster']=z.obs['r'+str(selected)+'_s'+str(seed)].astype('category');log('selected resolution '+str(selected)+' '+selection_status)
sens=[]
for label,emb,ids in [('uncorrected_pca30',pc[:,:30],np.arange(len(obs))),('harmony50',hpc,np.arange(len(obs))),('drop_lineage_review_flags',hpc[~obs.lineage_review_flag.to_numpy(),:30],np.flatnonzero(~obs.lineage_review_flag.to_numpy()))]:
 log('sensitivity '+label);v=graph(emb,obs.iloc[ids]);ll=cl(v,selected,seed,'cluster');ari=adjusted_rand_score(z.obs.cluster.astype(str).to_numpy()[ids],ll);e,d=mix(v,'cluster');sens.append({'sensitivity':label,'n_cells':len(ids),'clusters':len(set(ll)),'ARI_to_primary':ari,'weighted_study_entropy':e,'weighted_dominant_study_fraction':d});pd.DataFrame({'cell_id':obs.index[ids],'cluster':ll}).to_csv(O/(label+'_labels.tsv.gz'),sep='\t',index=False,compression='gzip');del v;gc.collect()
pd.DataFrame(sens).to_csv(O/'GRAPH_SENSITIVITY.tsv',sep='\t',index=False)
log('UMAP display only after structural resolution selection')
sc.tl.umap(z,min_dist=.35,spread=1,random_state=seed)
z.uns['expression_source']=str(R/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad');z.uns['selected_resolution']=selected;z.uns['selection_status']=selection_status
z.obs.to_csv(O/'CELL_CLUSTER_ASSIGNMENTS.tsv.gz',sep='\t',compression='gzip',index_label='cell_id');pd.DataFrame(z.obsm['X_umap'],index=z.obs_names,columns=['UMAP1','UMAP2']).to_csv(O/'UMAP_COORDINATES.tsv.gz',sep='\t',compression='gzip',index_label='cell_id');pd.DataFrame({'gene':hvg}).to_csv(O/'HVG_SELECTED.tsv',sep='\t',index=False)
pd.crosstab(z.obs.cluster,z.obs.study).to_csv(O/'CLUSTER_BY_STUDY.tsv',sep='\t');pd.crosstab(z.obs.cluster,z.obs.curated_patient_id).to_csv(O/'CLUSTER_BY_PATIENT.tsv',sep='\t');pd.crosstab(z.obs.cluster,z.obs.lineage_review_flag).to_csv(O/'CLUSTER_BY_REVIEW_FLAG.tsv',sep='\t')
f=O/'CLUSTERING_GRAPH.h5ad';z.write_h5ad(str(f)+'.part',compression='gzip');Path(str(f)+'.part').replace(f)
status={'status':'CLUSTERING_COMPLETE_PENDING_REVIEW','selection_status':selection_status,'selected_resolution':selected,'n_cells':z.n_obs,'clusters':z.obs.cluster.nunique(),'legacy_or_spp1_labels_used':False,'seed':seed,'resolution_metrics':metrics,'sensitivities':sens};(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2));print(json.dumps(status,ensure_ascii=False),flush=True)
