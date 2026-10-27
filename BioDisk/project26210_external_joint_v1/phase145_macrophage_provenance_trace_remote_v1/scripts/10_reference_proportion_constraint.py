from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase145_macrophage_provenance_trace_remote_v1');f=Path('/mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1/reference_v2/reference_v2_cell_assignments.tsv.gz')
a=pd.read_csv(f,sep='\t');states=['SPP1-state Macrophage','C1QC-state Macrophage','Other Macrophage'];a=a[a.state_v2.isin(states)];c=pd.crosstab([a.study_v1,a.patient],a.state_v2).reindex(columns=states,fill_value=0);c['total_macro']=c.sum(axis=1);c['spp1_fraction']=c[states[0]]/c.total_macro;c.to_csv(O/'identity/CURRENT_RULE_STATE_PATIENT_PROPORTIONS.tsv',sep='\t')
s={'n_patients':len(c),'n_macrophages':int(c.total_macro.sum()),'n_spp1_state':int(c[states[0]].sum()),'spp1_fraction_min':float(c.spp1_fraction.min()),'spp1_fraction_median':float(c.spp1_fraction.median()),'spp1_fraction_max':float(c.spp1_fraction.max()),'interpretation':'Proportions are influenced by within-patient quantile state definition; not unconditioned discovered-cluster prevalence.'}
(O/'identity/CURRENT_RULE_STATE_PROPORTION_STATUS.json').write_text(json.dumps(s,indent=2));print(json.dumps(s))
