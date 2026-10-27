from pathlib import Path
import pandas as pd,json
D=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase157_ordered_bridge_reassessment_20260915_v1/04_other_X_screen')
d=pd.read_csv(D/'CANDIDATE_EVIDENCE_RANKING.tsv',sep='\t')
print('BOTH_LEGS',d[d.both_legs_nominal][['candidate','a_rho','a_P','b_partial_rho','b_P','study_a_rho','study_b_rho','same_path_directions','LOPO_ab_positive']].to_dict('records'))
sel=d[d.candidate.isin(['ligand:CXCL16','ligand:IGFBP3','ligand:GRN','program:Glycolysis','program:Cholesterol Homeostasis'])]
cols=['candidate','XY_rho','XY_P','a_P_FDR_family','b_P_FDR_family','indirect_ab','ab_low','ab_high','study_ab','study_ab_low','study_ab_high','receptor_label','receptor_genes','n_C4_patients','patients_all_subunits_detected_10pct','median_min_subunit_detection','median_limiting_subunit_C4_minus_other','LOPO_ab_min','LOPO_ab_max','LOSO_ab_positive']
print('DETAIL',sel[cols].fillna('NA').to_dict('records'))
print('EXCLUDED_METABOLIC',pd.read_csv(D/'CANDIDATE_EXCLUSIONS.tsv',sep='\t').query("label.str.startswith('metabolism:')",engine='python').to_dict('records'))
