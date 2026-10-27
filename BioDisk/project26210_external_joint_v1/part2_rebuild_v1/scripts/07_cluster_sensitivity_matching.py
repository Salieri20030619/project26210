from pathlib import Path
import pandas as pd,numpy as np
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');C=R/(R/'provenance/CURRENT_CLUSTERING_RUN.txt').read_text().strip()
a=pd.read_csv(C/'CELL_CLUSTER_ASSIGNMENTS.tsv.gz',sep='\t',usecols=['cell_id','cluster'],dtype={'cluster':str}).set_index('cell_id');out=[]
for name in ['uncorrected_pca30','harmony50','drop_lineage_review_flags']:
 b=pd.read_csv(C/(name+'_labels.tsv.gz'),sep='\t',dtype={'cluster':str}).set_index('cell_id');m=a.join(b,how='inner',lsuffix='_primary',rsuffix='_sensitivity');t=pd.crosstab(m.cluster_primary,m.cluster_sensitivity)
 for c in t.index:
  jacc=t.loc[c]/(t.loc[c].sum()+t.sum(axis=0)-t.loc[c]);best=jacc.idxmax();out.append({'sensitivity':name,'primary_cluster':c,'best_matching_cluster':best,'n_primary_compared':int(t.loc[c].sum()),'n_shared_cells':int(t.loc[c,best]),'best_match_Jaccard':float(jacc[best]),'primary_retained_fraction_in_match':float(t.loc[c,best]/t.loc[c].sum())})
pd.DataFrame(out).to_csv(C/'CLUSTER_LEVEL_SENSITIVITY_MATCHES.tsv',sep='\t',index=False);print(pd.DataFrame(out).to_string(index=False))
