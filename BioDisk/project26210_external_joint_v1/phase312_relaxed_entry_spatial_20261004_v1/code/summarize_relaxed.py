from pathlib import Path
import pandas as pd
p=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase312_relaxed_entry_spatial_20261004_v1')
d=pd.read_csv(p/'tables/COMPLETE_DIRECTIONAL_MATCH.tsv',sep='\t')
e=['SPP1_vs_fibroblast_fraction','SPP1_vs_nearest_CD8_distance','fibroblast_fraction_vs_nearest_CD8_distance']
print('L2 PRIMARY50 COUNTS');print(d[(d.condition=='L2_cluster_broad')&(d.center_set=='condition_members')&(d.geometry=='Proseg_centroid')&(d.radius_um==50)][['assay','model']+e+['all_three_positive']].to_string(index=False))
print('FULL DIRECTION COUNTS');print(d.groupby(['center_set','condition','assay','model']).all_three_positive.agg(['sum','count']).to_string())
print('L2 MATCHED50');print(d[(d.condition=='L2_cluster_broad')&(d.center_set=='matched_strict_centers')&(d.geometry=='Proseg_centroid')&(d.assay=='Proseg')&(d.radius_um==50)][['model']+e].to_string(index=False))
