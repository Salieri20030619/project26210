from pathlib import Path
import json,datetime
import pandas as pd,numpy as np
from scipy.stats import rankdata
R=Path.home()/'project26210_external_joint_v1';O=R/'phase213_immune_restart_execution_20260919_v1';D=O/'I5_external_description';F=R/'phase211_immune_receiver_refreeze_20260919_v1/freeze_v2';assert json.loads((O/'I4_directed_communication/STATUS.json').read_text())['status']=='COMPLETE_VALIDATED';assert (D/'R_EXTRACTION_STATUS.txt').read_text().strip()=='PASS_SAMPLE_METRICS_NO_INFERENCE'
a=pd.read_csv(D/'EXTERNAL_SAMPLE_METRICS.tsv',sep='\t');f=pd.read_csv(F/'EXTERNAL_PATIENT_ENDPOINT_FREEZE.tsv',sep='\t');f=f[(f.dataset=='GSE294109')&f.eligible]
keys=['dataset','patient','timepoint','treatment','endpoint'];p=a.groupby(keys).agg(x=('x_mature_SPP1','mean'),y=('y_outcome','mean'),n_samples=('sample_id','size'),samples=('sample_id',lambda x:';'.join(sorted(x)))).reset_index();assert not p.duplicated(keys).any();assert len(p)==len(f)
for z in p.itertuples():
 row=f[(f.patient==z.patient)&(f.timepoint==z.timepoint)&(f.treatment==z.treatment)&(f.endpoint==z.endpoint)].iloc[0];assert set(z.samples.split(';'))==set(row.eligible_samples.split(';')) and z.n_samples==row.n_samples
 assert z.C4_used==False if hasattr(z,'C4_used') else True
p.to_csv(D/'EXTERNAL_PATIENT_METRICS.tsv',sep='\t',index=False)
rows=[]
for k,g in p.groupby(['timepoint','treatment','endpoint']):
 n=len(g);rho=np.nan
 if n>=3 and g.x.nunique()>1 and g.y.nunique()>1:rho=float(np.corrcoef(rankdata(g.x),rankdata(g.y))[0,1])
 direction='not_assessable_n1' if n==1 else ('descriptive_only_no_test')
 delta=np.nan
 if n>=2:
  high=g[g.x==g.x.max()].y.mean();low=g[g.x==g.x.min()].y.mean();delta=float(high-low) if g.x.nunique()>1 else np.nan
 rows.append(dict(timepoint=k[0],treatment=k[1],endpoint=k[2],n_patients=n,n_samples=int(g.n_samples.sum()),rho_descriptive_only=rho,outcome_difference_highest_vs_lowest_SPP1=delta,predictor_min=float(g.x.min()),predictor_max=float(g.x.max()),outcome_min=float(g.y.min()),outcome_max=float(g.y.max()),interpretation=direction,inferential_P_not_computed=True,FDR_not_computed=True,formal_external_validation=False))
s=pd.DataFrame(rows);s.to_csv(D/'EXTERNAL_STRATIFIED_DESCRIPTIVE_SUMMARY.tsv',sep='\t',index=False);assert len(s)==12 and set(s.n_patients)=={1,2,4};assert p.patient.nunique()==9
# Direct independent scalar verification of pseudobulk and library aggregation.
assert np.allclose(a.x_mature_SPP1,np.log2(1+1e6*a.mature_SPP1_counts/a.mature_total_counts),rtol=0,atol=1e-12)
for z in p.itertuples():
 x=a[(a.patient==z.patient)&(a.endpoint==z.endpoint)&(a.timepoint==z.timepoint)&(a.treatment==z.treatment)];assert abs(sum(x.x_mature_SPP1)/len(x)-z.x)<1e-12 and abs(sum(x.y_outcome)/len(x)-z.y)<1e-12
checks={'exact_frozen_samples_and_patient_sets':True,'pseudobulk_formula_independent_match':True,'equal_sample_patient_aggregation':True,'strata_remain_separate':len(s)==12,'no_n1_correlation':s.loc[s.n_patients==1,'rho_descriptive_only'].isna().all(),'no_inferential_tests':True,'no_C4_proxy_substitution':not a.C4_used.any()};checks={k:bool(v) for k,v in checks.items()};assert all(checks.values());(D/'VALIDATION.json').write_text(json.dumps(checks,indent=2));(D/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_VALIDATED','time':datetime.datetime.now().astimezone().isoformat(),'n_unique_patients':9,'n_patient_endpoint_rows':len(p),'strata':len(s),'no_inference':True,'GSE264586':'not_added; frozen identity/tumour-status gap preserved','GSE207493':'not independent external validation; overlap preserved','checks':checks},indent=2));print(s[['treatment','endpoint','n_patients','rho_descriptive_only','outcome_difference_highest_vs_lowest_SPP1']].to_string(index=False),flush=True)
