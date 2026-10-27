from pathlib import Path
import pandas as pd,json
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase197_top20_network_dynamics_20260917_v4'
a=pd.read_csv(r/'phase131_sections2_5_2_6_reference_correction_remote_v1/data/PATHWAY_PATIENT_GROUP_DOT_VALUES.tsv.gz',sep='\t')
v=a.pivot(index=['patient','factor','pathway'],columns='group',values='mean_AUC')
delta=(v.TYMP_positive-v.TYMP_negative).rename('delta').reset_index()
patient=delta.groupby(['patient','pathway']).delta.mean().unstack()
tab=pd.DataFrame({'median_patient_delta':patient.median(),'patients':patient.notna().sum(),'positive':(patient>0).sum(),'negative':(patient<0).sum()})
tab['absolute_median_delta']=abs(tab.median_patient_delta)
tab=tab.sort_values(['absolute_median_delta','median_patient_delta'],ascending=[False,False],kind='stable')
tab.insert(0,'rank',range(1,len(tab)+1));tab['selected']=tab['rank']<=20
tab.to_csv(d/'tables/TOP20_RANKING_ALL85.tsv',sep='\t')
patient.to_csv(d/'tables/PATIENT_BALANCED_PATHWAY_DELTA.tsv',sep='\t')
