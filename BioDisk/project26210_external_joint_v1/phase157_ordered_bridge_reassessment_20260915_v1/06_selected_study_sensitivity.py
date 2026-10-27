from pathlib import Path
import pandas as pd
from bridge_stats import bridge
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase157_ordered_bridge_reassessment_20260915_v1';D=O/'04_other_X_screen'
b=pd.read_csv(P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t');X=pd.read_csv(D/'PATIENT_CANDIDATE_VALUES.tsv',sep='\t',index_col=0).loc[b.patient]
selected=['program:Cholesterol Homeostasis','program:Glycolysis','ligand:IGFBP3','ligand:CXCL16','ligand:GRN'];rows=[]
for name in selected:
 for study in b.study.unique():
  keep=b.study!=study
  r=bridge(b.TYMP_pseudobulk[keep],X[name].to_numpy()[keep],b.C4_fraction[keep],b.study[keep],'raw')
  rows.append(dict(candidate=name,omitted_study=study,n_patients=int(keep.sum()),a_rho=r['a_rho'],b_partial_rho=r['b_partial_rho'],indirect_ab=r['indirect_ab']))
pd.DataFrame(rows).to_csv(D/'SELECTED_CANDIDATE_LEAVE_STUDY_DETAILS.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_dict('records'))
