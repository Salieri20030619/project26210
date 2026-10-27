from pathlib import Path
import json,datetime
import numpy as np,pandas as pd,anndata as ad
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'02_myeloid_qc/run_v2';O.mkdir()
plan={'date':datetime.datetime.now().isoformat(),'reason':'First-pass post-cluster identification revealed low-complexity cohort-dominated cluster C1 and vascular/pericyte mixed cluster C2. These are not accepted macrophage biological subtypes.','frozen_filter':{'minimum_genes':500,'minimum_UMIs':1000,'exclude_first_pass_vascular_cluster':'2','exclude_nonmyeloid_dominant_cells':'At least3 detected genes in a nonmyeloid panel, whose mean log1p CP10K exceeds myeloid backbone score. Panels and all flags retained.'},'selection_not_based_on_SPP1_or_C1QC_expression':True,'new_HVG_PCA_Harmony_required':True,'retain_all_failed_and_excluded_data_as_audit':True}
(O/'PURIFICATION_PLAN.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2))
z=ad.read_h5ad(R/'03_macrophage_clustering/run_v2/CLUSTERING_GRAPH.h5ad');a=ad.read_h5ad(R/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad');a=a[z.obs_names].copy();x=a.X.tocsr();tot=np.asarray(x.sum(axis=1)).ravel();norm=x.astype(float).multiply(10000/tot[:,None]).tocsr();norm.data=np.log1p(norm.data)
panels={'myeloid_backbone':['LST1','TYROBP','FCER1G','AIF1','CTSS','LYZ','CSF1R','CD68','CTSB'],'epithelial':['EPCAM','KRT8','KRT18','KRT19','CD24'],'T_cell':['CD3D','CD3E','TRAC'],'NK_cell':['NKG7','GNLY','KLRD1','PRF1'],'B_cell':['CD79A','MS4A1','CD37'],'endothelial':['PECAM1','VWF','KDR','RAMP2','EMCN'],'pericyte':['RGS5','PDGFRB','CSPG4','MCAM','NOTCH3','ACTA2','TAGLN','MYL9'],'fibroblast':['DCN','COL1A1','COL1A2','COL3A1','LUM'],'mast':['TPSAB1','TPSB2','KIT']}
m=z.obs.copy();m['first_pass_cluster']=m.cluster.astype(str)
for k,g in panels.items():
 idx=a.var_names.get_indexer([s for s in g if s in a.var_names]);m[k+'_score']=np.asarray(norm[:,idx].mean(axis=1)).ravel();m[k+'_detected']=np.asarray((x[:,idx]>0).sum(axis=1)).ravel()
m['low_complexity_for_fine_analysis']=(m.n_genes_raw<500)|(m.total_counts_raw<1000);m['vascular_mixed_cluster']=m.first_pass_cluster.eq('2');m['nonmyeloid_dominant']=False
for k in panels:
 if k!='myeloid_backbone':m['nonmyeloid_dominant']|=(m[k+'_detected']>=3)&(m[k+'_score']>m.myeloid_backbone_score)
m['keep_purified']=~(m.low_complexity_for_fine_analysis|m.vascular_mixed_cluster|m.nonmyeloid_dominant)
m.to_csv(O/'CELL_PURIFICATION_AUDIT.tsv.gz',sep='\t',compression='gzip',index_label='cell_id')
pd.crosstab(m.study,m.keep_purified).to_csv(O/'RETAINED_BY_STUDY.tsv',sep='\t');pd.crosstab(m.curated_patient_id,m.keep_purified).to_csv(O/'RETAINED_BY_PATIENT.tsv',sep='\t');pd.crosstab(m.sample_id,m.keep_purified).to_csv(O/'RETAINED_BY_LIBRARY.tsv',sep='\t');pd.crosstab(m.first_pass_cluster,m.keep_purified).to_csv(O/'RETAINED_BY_INITIAL_CLUSTER.tsv',sep='\t')
kept=m[m.keep_purified];kept.to_csv(O/'PURIFIED_CELL_METADATA.tsv.gz',sep='\t',compression='gzip',index_label='cell_id')
status={'status':'PASS_PURIFIED_INPUT_REVIEW','source_macrophages':len(m),'retained_macrophages':len(kept),'patients':kept.curated_patient_id.nunique(),'libraries':kept.sample_id.nunique(),'studies':kept.study.nunique(),'low_complexity_flag':int(m.low_complexity_for_fine_analysis.sum()),'vascular_cluster_flag':int(m.vascular_mixed_cluster.sum()),'nonmyeloid_dominant_flag':int(m.nonmyeloid_dominant.sum()),'excluded_union':int((~m.keep_purified).sum()),'no_SPP1_selection':True,'panels':panels};(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2));print(json.dumps(status,ensure_ascii=False),flush=True)
