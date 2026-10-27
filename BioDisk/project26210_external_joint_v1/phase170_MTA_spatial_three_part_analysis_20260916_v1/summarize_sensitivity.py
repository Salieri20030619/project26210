from pathlib import Path
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase170_MTA_spatial_three_part_analysis_20260916_v1')
s=pd.read_csv(O/'PATIENT_SUMMARY.tsv',sep='\t');p=pd.read_csv(O/'PATIENT_EFFECTS.tsv',sep='\t');d=pd.read_csv(O/'SECTION_EFFECTS.tsv',sep='\t');l=pd.read_csv(O/'LEAVE_ONE_PATIENT.tsv',sep='\t')
print('C2 policies')
print(s[(s.stratum=='all')&(s.policy!='equal_fisher_mean')][['role','candidate','policy','n_patients','positive','median_rho','P','FDR']].to_string(index=False))
print('LOO median ranges')
print(l.groupby(['role','candidate']).median_rho.agg(['min','max']).to_string())
print('NEW_EXCLUSIONS')
print(d[(d.model=='coordinate_adjusted')&(d.spots>=50)&(~d.eligible)][['section','role','candidate','spots','condition_number','reason']].to_string(index=False))
a=p[(p.model=='extended_primary')&(p.stratum=='all')&(p.policy=='equal_fisher_mean')]
b=p[(p.model=='coordinate_adjusted')&(p.stratum=='all')&(p.policy=='equal_fisher_mean')]
j=a.merge(b,on=['patient_id','role','candidate'],suffixes=('_primary','_coordinate'))
j.to_csv(O/'COORDINATE_ADJUSTMENT_PAIRED.tsv',sep='\t',index=False)
print('MATCHED_COORD')
print(j.groupby(['role','candidate']).agg(n=('patient_id','size'),primary=('rho_primary','median'),coordinate=('rho_coordinate','median')).to_string())
pair=pd.read_csv(O/'MTAP_ADJUSTMENT_PAIRED.tsv',sep='\t')
print('MTAP CHANGE')
print(pair.groupby(['role','candidate']).agg(n=('patient_id','size'),change=('rho_change','median')).to_string())
print('EXCLUDE_COHORT')
print(s[s.stratum.str.startswith('exclude_cohort_')][['role','candidate','stratum','n_patients','positive','median_rho']].to_string(index=False))
