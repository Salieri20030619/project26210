from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['NUMBA_NUM_THREADS']='4';os.environ['OMP_NUM_THREADS']='4'
import json,datetime,subprocess,gc
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
P=Path(__file__).resolve().parents[1];plan=pd.read_csv(P/'tables/SAMPLE_PLAN.tsv',sep='\t');assert json.loads((P/'qc/INPUT_ARCHIVES_VALIDATED.json').read_text())['status']=='PASS';rows=[];objects=[]
panels={
'T':['Cd3d','Cd3e','Trac','Cd247'], 'CD8':['Cd8a','Cd8b1'],'CD4_Treg':['Cd4','Foxp3','Il2ra'],
'NK':['Nkg7','Klrd1','Ncr1','Klrb1c'],'B':['Cd79a','Cd79b','Ms4a1','Cd19'],'Plasma':['Jchain','Mzb1','Sdc1'],
'Myeloid':['Lyz2','Tyrobp','Lst1','Fcer1g'],'Macrophage':['C1qa','C1qb','C1qc','Adgre1','Csf1r','Cd68','Mrc1'],
'Monocyte':['Ly6c2','Ccr2','Plac8','S100a8','S100a9'],'Neutrophil':['Csf3r','Retnlg','S100a8','S100a9','Ly6g'],
'DC':['Flt3','Zbtb46','Xcr1','Clec9a','Cd209a'],'pDC':['Siglech','Gzmb','Tcf4','Bst2'],'Mast':['Kit','Cpa3','Mcpt4','Tpsb2'],
'Fibroblast':['Col1a1','Col1a2','Col3a1','Dcn','Lum','Pdgfra'],'Mural':['Rgs5','Cspg4','Mcam','Myh11','Des'],
'Endothelial':['Pecam1','Cdh5','Vwf','Kdr','Esam'],'Epithelial':['Epcam','Krt8','Krt18','Krt19','Krt7'],
'Mesothelial':['Msln','Upk3b','Wt1','Krt19'],'MHCII':['Cd74','H2-Aa','H2-Ab1','H2-Eb1'],'Cycling':['Mki67','Top2a','Stmn1'],
'Erythroid':['Hba-a1','Hba-a2','Hbb-bs','Alas2'],'Acinar':['Prss1','Prss2','Cpa1','Cela3b'],'Lymphatic':['Prox1','Lyve1','Pdpn','Flt4']}
(P/'resources/IDENTITY_MARKER_PANELS.json').write_text(json.dumps(panels,indent=2))
for r in plan.itertuples():
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 a=sc.read_10x_mtx(P/'inputs'/r.library,var_names='gene_symbols',make_unique=True,cache=False);a.var_names_make_unique();a.obs_names=[r.library+':'+b for b in a.obs_names];x=a.X.tocsr();assert np.isfinite(x.data).all() and (x.data>=0).all() and (x.data==np.floor(x.data)).all();total=np.asarray(x.sum(1)).ravel();ng=np.asarray((x>0).sum(1)).ravel();mito=np.asarray(x[:,a.var_names.str.lower().str.startswith('mt-')].sum(1)).ravel()/np.maximum(total,1)*100;upper=float(np.quantile(total,.995));mask=(ng>=200)&(ng<=6000)&(total>=500)&(total<=upper)&(mito<=20);a.obs['total_counts']=total;a.obs['n_genes']=ng;a.obs['pct_mito']=mito;a.obs['library']=r.library;a.obs['group']=r.group;a.obs['pool_tumors']=2;raw_n=a.n_obs;a=a[mask].copy()
 sc.pp.scrublet(a,expected_doublet_rate=.06,random_state=26210,verbose=False);before=a.n_obs;nd=int(a.obs.predicted_doublet.sum());a=a[~a.obs.predicted_doublet].copy();a.X=a.X.tocsr().astype(np.int32);a.write_h5ad(P/f'inputs/{r.library}_counts_qc.h5ad',compression='gzip');objects.append(a);rows.append({'library':r.library,'group':r.group,'raw_cells':raw_n,'QC_before_doublets':before,'doublets_removed':nd,'final_cells':a.n_obs,'genes_measured':a.n_vars,'upper_UMI_995':upper,'median_UMI':float(a.obs.total_counts.median()),'median_genes':float(a.obs.n_genes.median()),'Spp1_measured':'Spp1' in a.var_names,'Tymp_measured':'Tymp' in a.var_names});print('QC_LIBRARY',r.library,raw_n,before,nd,a.n_obs,flush=True)
a=ad.concat(objects,join='outer',merge='first');a.obs['library']=a.obs.library.astype('category');del objects;gc.collect();a.layers['counts']=a.X.copy();sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);sc.pp.highly_variable_genes(a,n_top_genes=3000,batch_key='library',flavor='seurat');a.var.loc[a.var_names.isin(['Spp1','Tymp']),'highly_variable']=False;b=a[:,a.var.highly_variable].copy();sc.pp.scale(b,max_value=10);sc.tl.pca(b,n_comps=40,svd_solver='arpack',random_state=26210);a.obsm['X_pca']=b.obsm['X_pca'];del b;gc.collect();sc.pp.neighbors(a,n_neighbors=20,n_pcs=40,random_state=26210);sc.tl.leiden(a,resolution=.7,key_added='cluster',flavor='igraph',n_iterations=2,directed=False,random_state=26210)
x=a.X.tocsr();counts=a.layers['counts'];overall=np.asarray(x.mean(0)).ravel();markers=[];profiles=[];coverage=[]
for cl in sorted(a.obs.cluster.unique(),key=int):
 mask=a.obs.cluster.eq(cl).to_numpy();n=int(mask.sum());avg=np.asarray(x[mask].mean(0)).ravel();det=np.asarray((counts[mask]>0).mean(0)).ravel();enrich=avg-overall;allowed=(det>=.20)&(~a.var_names.isin(['Spp1','Tymp']));ids=np.where(allowed)[0];top=ids[np.argsort(enrich[ids])[::-1][:35]]
 for j in top:markers.append({'cluster':cl,'gene':a.var_names[j],'mean_logCP10K':avg[j],'detected_fraction':det[j],'enrichment_vs_all':enrich[j]})
 rec={'cluster':cl,'cells':n,'median_UMI':float(a.obs.loc[mask,'total_counts'].median()),'median_genes':float(a.obs.loc[mask,'n_genes'].median())}
 for lin,gs in panels.items():
  jj=a.var_names.get_indexer([g for g in gs if g in a.var_names]);rec[lin+'_score']=float(avg[jj].mean()) if len(jj) else np.nan
  for gene,j in zip([g for g in gs if g in a.var_names],jj):coverage.append({'cluster':cl,'panel':lin,'gene':gene,'mean_logCP10K':avg[j],'detected_fraction':det[j]})
 profiles.append(rec)
pd.DataFrame(rows).to_csv(P/'tables/QC_BY_LIBRARY.tsv',sep='\t',index=False);pd.DataFrame(markers).to_csv(P/'tables/CLUSTER_MARKERS.tsv',sep='\t',index=False);pd.DataFrame(profiles).to_csv(P/'tables/CLUSTER_IDENTITY_PROFILE.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(P/'tables/CLUSTER_MARKER_COVERAGE.tsv',sep='\t',index=False);a.obs.to_csv(P/'tables/CELL_METADATA_PRE_IDENTITY.tsv.gz',sep='\t');a.layers['lognorm']=a.X.copy();a.X=a.layers.pop('counts');a.write_h5ad(P/'inputs/counts_qc_clustered.h5ad',compression='gzip');(P/'qc/QC_CLUSTER_COMPLETE.json').write_text(json.dumps({'status':'READY_FOR_IDENTITY_REVIEW','cells':a.n_obs,'genes':a.n_vars,'clusters':a.obs.cluster.nunique(),'time':datetime.datetime.now().astimezone().isoformat()},indent=2));print('QC_CLUSTER_COMPLETE',a.shape,a.obs.cluster.nunique(),flush=True)
