from pathlib import Path
import json,datetime
import pandas as pd,numpy as np
from scipy.stats import spearmanr,rankdata,false_discovery_control
from stats_shared import *
R=Path.home()/'project26210_external_joint_v1';O=R/'phase213_immune_restart_execution_20260919_v1';assert json.loads((O/'I1_patient_inputs/STATUS.json').read_text())['status']=='COMPLETE_VALIDATED';D=O/'I2_main_statistics';D.mkdir(exist_ok=False);(D/'distributions').mkdir()
inputs=pd.read_csv(O/'I1_patient_inputs/PATIENT_HYPOTHESIS_INPUTS.tsv',sep='\t');reg=pd.read_csv(R/'phase212_immune_restart_plan_20260919_v1/HYPOTHESIS_REGISTER.tsv',sep='\t');rows=[];studyrows=[]
for h in reg.itertuples():
 d=inputs[(inputs.hypothesis_id==h.hypothesis_id)&(inputs.variant=='primary')].sort_values(['study','patient']);assert len(d)==h.frozen_patient_n;x=d.x.to_numpy();y=d.y.to_numpy();s=d.study.to_numpy();res,dist=main_test(x,y,s);lo,hi,nb,boot=bootstrap_effect(x,y,s);raw=spearmanr(x,y)
 gs=groups(s);info=[g for g in gs if len(g)>1 and np.ptp(x[g])>0 and np.ptp(y[g])>0]
 z=dict(hypothesis_id=h.hypothesis_id,predictor=h.predictor,endpoint=h.endpoint,tier=h.analysis_tier,expected_direction=h.expected_direction,family=h.family,family_size=h.family_size,adjustment=h.adjustment,n_patients=len(d),n_studies=d.study.nunique(),study_counts=json.dumps(d.study.value_counts().to_dict()),informative_studies=len(info),informative_patients=sum(len(g) for g in info),rho=res['rho'],p_value=res['p'],permutation_method=res['method'],permutations=res['permutations'],possible_permutations=res['possible_permutations'],residual_df=res['residual_df'],bootstrap_CI_low=lo,bootstrap_CI_high=hi,bootstrap_valid=nb,bootstrap_total=5000,unadjusted_rho=float(raw.statistic),unadjusted_p_asymptotic=float(raw.pvalue),inference='patient-level; library/technical replication excluded')
 rows.append(z);np.savez_compressed(D/f'distributions/{h.hypothesis_id}.npz',permutation_rhos=dist,bootstrap_rhos=boot)
 for st,g in d.groupby('study'):
  rho=float(spearmanr(g.x,g.y).statistic) if len(g)>=2 and g.x.nunique()>1 and g.y.nunique()>1 else np.nan
  studyrows.append(dict(hypothesis_id=h.hypothesis_id,study=st,n_patients=len(g),rho_descriptive=rho,reporting='descriptive_only; n<10; no study-level P or FDR',median_predictor=float(g.x.median()),median_outcome=float(g.y.median())))
 print(h.hypothesis_id,'n',len(d),'rho',res['rho'],'P',res['p'],flush=True)
a=pd.DataFrame(rows);a['FDR']=np.nan
for fam,g in a.groupby('family'):
 if len(g)>1:
  p=g.p_value.fillna(1).to_numpy();q=bh(p);a.loc[g.index,'FDR']=np.where(g.p_value.notna(),q,np.nan)
a['expected_direction_observed']=np.where(a.expected_direction=='positive',a.rho>0,a.rho<0)
def grade(row):
 if not np.isfinite(row.rho) or not np.isfinite(row.p_value):return 'not_estimable'
 value=row.FDR if row.adjustment=='BH' else row.p_value
 if value<.05:return ('expected_direction_supported' if row.expected_direction_observed else 'opposite_direction_supported')+('_FDR' if row.adjustment=='BH' else '_prespecified_single_endpoint_P')
 if row.p_value<.05:return 'nominal_only_'+('expected_direction' if row.expected_direction_observed else 'opposite_direction')
 return 'not_supported_'+('expected_direction' if row.expected_direction_observed else 'opposite_direction')
a['evidence_grade']=a.apply(grade,axis=1);a.to_csv(D/'PRIMARY_12_HYPOTHESES.tsv',sep='\t',index=False);pd.DataFrame(studyrows).to_csv(D/'WITHIN_STUDY_DESCRIPTIVE_EFFECTS.tsv',sep='\t',index=False)
# Independent design-matrix residualization and stored-distribution P/BH checks.
checks=[]
for row in a.itertuples():
 d=inputs[(inputs.hypothesis_id==row.hypothesis_id)&(inputs.variant=='primary')].sort_values(['study','patient']);Z=np.column_stack([np.ones(len(d)),pd.get_dummies(d.study,drop_first=True,dtype=float).to_numpy()]);Q=np.linalg.qr(Z,mode='reduced')[0];rx=rankdata(d.x);ry=rankdata(d.y);xx=rx-Q@(Q.T@rx);yy=ry-Q@(Q.T@ry);rho=np.dot(xx,yy)/(np.linalg.norm(xx)*np.linalg.norm(yy));assert abs(rho-row.rho)<1e-10
 dist=np.load(D/f'distributions/{row.hypothesis_id}.npz');v=dist['permutation_rhos'];hits=int((np.abs(v)>=abs(row.rho)-1e-12).sum());p=hits/len(v) if row.permutation_method.startswith('exact') else (hits+1)/(len(v)+1);assert abs(p-row.p_value)<1e-15
 b=dist['bootstrap_rhos'];valid=np.isfinite(b);assert valid.sum()==row.bootstrap_valid
 if valid.mean()>=.8:assert np.allclose(np.quantile(b[valid],[.025,.975]),[row.bootstrap_CI_low,row.bootstrap_CI_high],atol=1e-12)
 checks.append(dict(hypothesis_id=row.hypothesis_id,QR_effect_match=True,stored_permutation_p_match=True,bootstrap_match=True))
for _,g in a.groupby('family'):
 if len(g)>1:assert np.allclose(false_discovery_control(g.p_value.fillna(1).to_numpy()),g.FDR,equal_nan=True)
pd.DataFrame(checks).to_csv(D/'INDEPENDENT_VALIDATION.tsv',sep='\t',index=False);(D/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_VALIDATED','time':datetime.datetime.now().astimezone().isoformat(),'hypotheses':len(a),'independent_checks':len(checks)*3+sum(len(g)>1 for _,g in a.groupby('family')),'seed_permutation':262212,'seed_bootstrap':263212,'bootstrap_n':5000,'primary_p_not_selected_from_raw_comparison':True},indent=2));print(a[['hypothesis_id','n_patients','rho','p_value','FDR','evidence_grade']].to_string(index=False),flush=True)
