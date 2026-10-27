from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import wilcoxon
from sklearn.metrics import adjusted_rand_score
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');C=R/(R/'provenance/CURRENT_CLUSTERING_RUN.txt').read_text().strip();M=R/(R/'provenance/CURRENT_MARKER_RUN.txt').read_text().strip();V=R/(R/'provenance/CURRENT_VALIDATION_RUN.txt').read_text().strip();O=R/'provenance/three_step_validation_v1';O.mkdir()
expected_n=json.loads((R/'02_myeloid_qc/run_v2/STATUS.json').read_text())['retained_macrophages']
checks=[]
def check(k,x,y):
 ok=bool(np.isclose(x,y,atol=1e-8,rtol=1e-7,equal_nan=True)) if isinstance(x,(float,int,np.number)) else x==y
 checks.append({'check':k,'actual':str(x),'expected':str(y),'status':'PASS' if ok else 'FAIL'})
z=ad.read_h5ad(C/'CLUSTERING_GRAPH.h5ad');raw=ad.read_h5ad(R/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad');a=raw[z.obs_names].copy();del raw;x=a.X.tocsr();clusters=sorted(z.obs.cluster.astype(str).unique(),key=int)
check('unique_cells',z.obs_names.nunique(),expected_n);check('all_parent_broad_labels_macrophage',bool((a.obs.cell_type_broad_v1.astype(str)=='Macrophage').all()),True);check('all_current_QC_pass',bool(a.obs.analysis_keep.astype(bool).all()),True);check('fine_analysis_complexity',bool(((a.obs.n_genes_raw>=500)&(a.obs.total_counts_raw>=1000)).all()),True);check('patients',a.obs.curated_patient_id.nunique(),40);check('libraries',a.obs.sample_id.nunique(),75);check('genes',a.n_vars,37795);check('finite_nonnegative_integer_counts',bool(np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.rint(x.data)).all()),True)
check('graph_shape',z.obsp['connectivities'].shape,(expected_n,expected_n));dif=z.obsp['connectivities']-z.obsp['connectivities'].T;check('graph_symmetric',float(abs(dif).max()) if dif.nnz else 0,0);check('finite_UMAP',bool(np.isfinite(z.obsm['X_umap']).all()),True);check('finite_Harmony',bool(np.isfinite(z.obsm['X_harmony50']).all()),True)
met=pd.read_csv(C/'RESOLUTION_STABILITY.tsv',sep='\t');passed=met[met.structural_gate_pass];expected_resolution=float(met[met.clusters>=2].sort_values(['minimum_between_seed_ARI','mean_between_seed_ARI','resolution'],ascending=[False,False,True]).iloc[0].resolution) if z.uns['selection_status']=='EXPLORATORY_NOT_GLOBALLY_STABLE' else float(passed.resolution.max());check('predefined_resolution_selection',float(z.uns['selected_resolution']),expected_resolution);check('HVG_count',len(pd.read_csv(C/'HVG_SELECTED.tsv',sep='\t')),3000)
ann=pd.read_csv(M/'CLUSTER_ANNOTATIONS.tsv',sep='\t',dtype={'cluster':str});check('one_annotation_per_cluster',sorted(ann.cluster,key=int),clusters);check('annotation_unique',ann.cluster.nunique(),len(clusters))
comp=pd.read_csv(V/'CLUSTER_COMPOSITION_QC.tsv',sep='\t',dtype={'cluster':str});check('composition_conservation',int(comp.n_cells.sum()),expected_n)
for c in clusters:
 tab=pd.read_csv(M/'full_marker_tables'/('cluster_'+c+'_all_genes.tsv.gz'),sep='\t',float_precision='round_trip');qfile=tab.BH_FDR_cluster_gene_family.to_numpy();p=tab.patient_wilcoxon_p.to_numpy();ids=np.flatnonzero(np.isfinite(p));order=ids[np.argsort(p[ids])];q=np.minimum.accumulate((p[order]*len(ids)/np.arange(1,len(ids)+1))[::-1])[::-1];check('BH_cluster_'+c,bool(np.allclose(qfile[order],np.minimum(q,1),rtol=1e-7,atol=1e-12)),True);check('all_genes_cluster_'+c,len(tab),37795)
 d=np.load(M/'paired_patient_effects'/('cluster_'+c+'.npz'));df=d['effects'];pt=d['patient'];check('unique_paired_patients_'+c,len(set(pt)),len(pt));check('gene_effect_dimensions_'+c,df.shape,(int(tab.eligible_patients.iloc[0]),37795));check('median_effects_'+c,bool(np.allclose(np.median(df,axis=0),tab.median_patient_log2CPM_difference.to_numpy(),equal_nan=True)),True)
 selected=['SPP1']+tab.sort_values('annotation_priority',ascending=False).gene.head(2).tolist()
 for g in dict.fromkeys(selected):
  gi=a.var_names.get_loc(g);r=tab[tab.gene==g].iloc[0];m=z.obs.cluster.astype(str).to_numpy()==c;vals=x[m,gi].toarray().ravel();tot=np.asarray(x[m].sum(axis=1)).ravel();check('pct_'+c+'_'+g,float((vals>0).mean()),r.pct_in);check('mean_'+c+'_'+g,float(np.log1p(vals/tot*10000).mean()),r.mean_log1pCP10K_in)
  if r.test_eligible:
   pv=wilcoxon(df[:,gi],method='approx',correction=True,zero_method='wilcox').pvalue;check('P_'+c+'_'+g,float(pv),r.patient_wilcoxon_p)
 # Direct raw-count re-aggregation for two eligible patients and three audit genes.
 for patient in pt[:2]:
  pm=z.obs.curated_patient_id.astype(str).to_numpy()==str(patient);inside=pm&(z.obs.cluster.astype(str).to_numpy()==c);outside=pm&~(z.obs.cluster.astype(str).to_numpy()==c);check('eligible_pair_min_cells_'+c+'_'+str(patient),bool(inside.sum()>=20 and outside.sum()>=20),True)
  for g in dict.fromkeys(selected):
   gi=a.var_names.get_loc(g);u=float(x[inside,gi].sum());v=float(x[outside,gi].sum());effect=np.log2(u/float(x[inside].sum(dtype=np.float64))*1e6+1)-np.log2(v/float(x[outside].sum(dtype=np.float64))*1e6+1);check('raw_paired_effect_'+c+'_'+str(patient)+'_'+g,float(effect),float(df[list(pt).index(patient),gi]))
pd.DataFrame(checks).to_csv(O/'CHECKS.tsv',sep='\t',index=False);status={'status':'PASS' if all(c['status']=='PASS' for c in checks) else 'FAIL','checks':len(checks),'passed':sum(c['status']=='PASS' for c in checks),'failures':[c for c in checks if c['status']!='PASS'],'global_partition_stability_pass':z.uns['selection_status']!='EXPLORATORY_NOT_GLOBALLY_STABLE','interpretation':'Technical/numerical verification of executed discovery workflow; not external biological validation.'};(O/'STATUS.json').write_text(json.dumps(status,indent=2));print(json.dumps(status))
