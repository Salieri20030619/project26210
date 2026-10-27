from pathlib import Path
import json,datetime
import pandas as pd,numpy as np
from scipy.stats import spearmanr
from stats_shared import effect,residual_ranks
R=Path.home()/'project26210_external_joint_v1';O=R/'phase213_immune_restart_execution_20260919_v1';assert json.loads((O/'I2_main_statistics/STATUS.json').read_text())['status']=='COMPLETE_VALIDATED';D=O/'I3_robustness';D.mkdir(exist_ok=False)
d=pd.read_csv(O/'I1_patient_inputs/PATIENT_HYPOTHESIS_INPUTS.tsv',sep='\t');reg=pd.read_csv(R/'phase212_immune_restart_plan_20260919_v1/HYPOTHESIS_REGISTER.tsv',sep='\t');loo=[];los=[];sensitivity=[];activation=[];common=[]
for h in reg.itertuples():
 a=d[(d.hypothesis_id==h.hypothesis_id)&(d.variant=='primary')];base=effect(a.x,a.y,a.study)
 for patient in a.patient:
  z=a[a.patient!=patient];r=effect(z.x,z.y,z.study);loo.append(dict(hypothesis_id=h.hypothesis_id,removed_patient=patient,n_patients=len(z),rho=r,primary_rho=base,direction_preserved=bool(np.isfinite(r) and np.sign(r)==np.sign(base))))
 for study in a.study.unique():
  z=a[a.study!=study];r=effect(z.x,z.y,z.study);los.append(dict(hypothesis_id=h.hypothesis_id,removed_study=study,n_patients=len(z),rho=r,primary_rho=base,direction_preserved=bool(np.isfinite(r) and np.sign(r)==np.sign(base)),reporting='descriptive_only'))
 for variant in ['exclude_concentrated','high_confidence_only']:
  z=d[(d.hypothesis_id==h.hypothesis_id)&(d.variant==variant)];r=effect(z.x,z.y,z.study) if len(z) else np.nan;sensitivity.append(dict(hypothesis_id=h.hypothesis_id,check=variant,n_patients=len(z),rho=r,primary_rho=base,estimable=np.isfinite(r),residual_df=np.nan,reporting='derived receiver population; descriptive effect only; not replacement'))
 if h.endpoint!='Treg_CD4':
  _,rankZ=residual_ranks(a.x,a.study,a.receiver_mean_n_genes);df=len(a)-rankZ-1;r=effect(a.x,a.y,a.study,a.receiver_mean_n_genes) if df>=5 else np.nan;sensitivity.append(dict(hypothesis_id=h.hypothesis_id,check='receiver_gene_complexity_adjusted',n_patients=len(a),rho=r,primary_rho=base,estimable=np.isfinite(r),residual_df=df,reporting='one prespecified technical covariate; effect only'))
 if h.endpoint=='CD8_dysfunction':
  _,rankZ=residual_ranks(a.x,a.study,a.activation);df=len(a)-rankZ-1;r=effect(a.x,a.y,a.study,a.activation) if df>=5 else np.nan;sensitivity.append(dict(hypothesis_id=h.hypothesis_id,check='activation_adjusted_dysfunction',n_patients=len(a),rho=r,primary_rho=base,estimable=np.isfinite(r),residual_df=df,reporting='secondary explanation; not main endpoint replacement'))
  activation.append(dict(hypothesis_id=h.hypothesis_id,n_patients=len(a),predictor_vs_activation_rho=effect(a.x,a.activation,a.study),dysfunction_vs_activation_rho=effect(a.y,a.activation,a.study),dysfunction_adjusted_for_activation_rho=r))
 z=d[(d.hypothesis_id==h.hypothesis_id)&(d.variant=='common_C4_libraries')];common.append(dict(hypothesis_id=h.hypothesis_id,predictor=h.predictor,endpoint=h.endpoint,n_patients=len(z),rho=effect(z.x,z.y,z.study),scope='identical C4-intrinsic-qualified libraries for all three predictors',patients=';'.join(sorted(z.patient))))
for name,rows in [('LEAVE_ONE_PATIENT.tsv',loo),('LEAVE_ONE_STUDY.tsv',los),('PRESET_SENSITIVITY_EFFECTS.tsv',sensitivity),('ACTIVATION_DYSFUNCTION_COMPARISON.tsv',activation),('COMMON_C4_LIBRARY_EFFECTS.tsv',common)]:pd.DataFrame(rows).to_csv(D/name,sep='\t',index=False)
# summarize influence without multiplicity of opportunistic P values.
a=pd.DataFrame(loo);b=pd.DataFrame(los);summary=[]
for h in reg.hypothesis_id:
 x=a[a.hypothesis_id==h];y=b[b.hypothesis_id==h];summary.append(dict(hypothesis_id=h,leave_patient_n=len(x),leave_patient_estimable=int(x.rho.notna().sum()),leave_patient_direction_preserved=int(x.direction_preserved.sum()),leave_patient_rho_min=x.rho.min(),leave_patient_rho_max=x.rho.max(),leave_study_n=len(y),leave_study_estimable=int(y.rho.notna().sum()),leave_study_direction_preserved=int(y.direction_preserved.sum()),leave_study_rho_min=y.rho.min(),leave_study_rho_max=y.rho.max()))
pd.DataFrame(summary).to_csv(D/'INFLUENCE_SUMMARY.tsv',sep='\t',index=False)
cc=pd.DataFrame(common)
for _,g in cc.groupby('endpoint'):assert g.n_patients.nunique()==1 and g.patients.nunique()==1
assert len(sensitivity)==36;assert not any('p_value' in row for row in sensitivity)
(D/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_VALIDATED','time':datetime.datetime.now().astimezone().isoformat(),'leave_patient_rows':len(loo),'leave_study_rows':len(los),'preset_sensitivity_rows':len(sensitivity),'common_library_predictors_exact_patient_match':True,'no_new_p_values':True,'primary_inputs_unchanged':True},indent=2));print((D/'STATUS.json').read_text(),flush=True)
