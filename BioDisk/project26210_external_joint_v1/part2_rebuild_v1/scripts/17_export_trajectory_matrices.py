from pathlib import Path
import json,gzip
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy.io import mmwrite
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');I=R/'07_trajectory/input_v1';O=R/'07_trajectory/matrices_v1';O.mkdir(exist_ok=False)
a=ad.read_h5ad(I/'TRAJECTORY_RAW_COUNTS.h5ad');plan={'status':'FROZEN_BEFORE_TRAJECTORY','study_separate':True,'fit_excludes':['SPP1','TYMP'],'features':'study-local raw-count seurat_v3 HVG2500; exclude mitochondrial/ribosomal/Hb; log1p CP10K full library and scaled PCA30','primary':'All eligible mono and macrophage candidates; study-local graph; use_partition TRUE; close_loop FALSE; >=20 credible original monocytes per rooted partition','roots':'Monocyte closest to partition mono median PCA; sensitivity distance quantiles0.25 and0.50','sensitivity':'DPT on study PCA k30; independent refit excluding C7 C11 C15 stress/cycle; tumor-restricted patient summaries','statistics':'Patient paired C4 vs Monocyte and C4 vs other mature macrophages >=20 cells each. Normalize pseudotime0-1 within rooted partition; no comparison between disconnected partitions. Rootless C4 fractions reported.','figures':False};(O/'PLAN.json').write_text(json.dumps(plan,indent=2));registry=[]
for mode in ['primary','without_stress_cycle']:
 for study in sorted(a.obs.study.unique()):
  mask=a.obs.study.eq(study)
  if mode!='primary':mask&=~a.obs.candidate.astype(str).isin(['7','11','15'])
  z=a[mask].copy();p=O/mode/study;p.mkdir(parents=True)
  ids=z.obs_names.tolist();keep=np.asarray((z.X>0).sum(axis=0)).ravel()>=20;names=z.var_names;keep&=~names.str.match(r'^(MT-|RPL|RPS|HB[ABDEGMQZ])');keep&=~names.isin(['SPP1','TYMP']);h=z[:,keep].copy();sc.pp.highly_variable_genes(h,flavor='seurat_v3',n_top_genes=min(2500,h.n_vars));sel=h.var_names[h.var.highly_variable];z.var['fit_hvg']=z.var_names.isin(sel)
  raw=z.X.copy();tot=np.asarray(raw.sum(axis=1)).ravel();h=z[:,sel].copy();h.X=h.X.astype(float).multiply(10000/tot[:,None]).tocsr();h.X.data=np.log1p(h.X.data);sc.pp.scale(h,max_value=10);sc.tl.pca(h,n_comps=min(30,h.n_obs-2),random_state=262148,svd_solver='arpack');sc.pp.neighbors(h,n_neighbors=30,n_pcs=30,random_state=262148)
  sc.tl.diffmap(h,n_comps=15,random_state=262148);h.obs=z.obs.copy();mono=np.flatnonzero(h.obs.root_eligible.to_numpy());dpt=np.full(h.n_obs,np.nan)
  if len(mono)>=20:
   pos=h.obsm['X_pca'];center=np.median(pos[mono],axis=0);iroot=mono[np.argmin(((pos[mono]-center)**2).sum(axis=1))];h.uns['iroot']=int(iroot);sc.tl.dpt(h,n_dcs=10);dpt=h.obs.dpt_pseudotime.to_numpy()
  md=z.obs.copy();md['dpt_pseudotime']=dpt;md.to_csv(p/'metadata.tsv.gz',sep='\t',index_label='cell_id');pd.DataFrame(h.obsm['X_pca'],index=ids).to_csv(p/'pca.tsv.gz',sep='\t',index_label='cell_id');pd.Series(z.var_names).to_csv(p/'genes.tsv',index=False,header=False);pd.Series(sel).to_csv(p/'HVG.tsv',index=False,header=False)
  with gzip.open(p/'counts.mtx.gz','wb') as f:mmwrite(f,raw.T)
  registry.append({'mode':mode,'study':study,'cells':z.n_obs,'root_candidates':len(mono),'C4':int((z.obs.candidate.astype(str)=='4').sum()),'DPT_finite':int(np.isfinite(dpt).sum())});print(registry[-1],flush=True)
(O/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_MATRIX_EXPORT','studies_modes':registry,'figures_created':False},indent=2))
