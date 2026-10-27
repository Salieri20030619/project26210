from pathlib import Path
import pandas as pd,numpy as np,json,warnings
from statsmodels.duration.hazard_regression import PHReg
p=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase270_SPP1_literature_clinical_revalidation_20260927_v1')
r=pd.read_csv(p/'tables/ALL_SURVIVAL_MODELS.tsv',sep='\t'); mem=pd.read_csv(p/'tables/MODEL_MEMBERSHIP.tsv',sep='\t')
checks=[]
for _,row in r[r.model=='continuous_adjusted'].iterrows():
 d=pd.read_csv(p/f'tables/{row.dataset.upper()}_ANALYSIS_INPUT.tsv',sep='\t')
 z=mem[mem.model_id==row.model_id].merge(d,on='sample_id')
 cols=['x','age','sex']+(['trial','risk'] if row.dataset=='CheckMate' else [])
 ex=pd.get_dummies(z[cols],drop_first=True,dtype=float)
 ex=ex.loc[:,ex.nunique()>1]
 with warnings.catch_warnings(record=True) as ws:
  warnings.simplefilter('always');f=PHReg(z.time,ex,status=z.event,ties='efron').fit(disp=0,maxiter=200)
 # confirm parameter scores are near zero even where optimizer convergence flags depend on tolerance
 grad=f.model.score(f.params); checks.append(dict(model_id=int(row.model_id),max_abs_score=float(np.max(abs(grad))),warnings=';'.join(str(w.message) for w in ws),parameters_finite=bool(np.isfinite(f.params).all())))
pd.DataFrame(checks).to_csv(p/'qc/ADJUSTED_OPTIMIZER_CHECK.tsv',sep='\t',index=False)
out=[]
it=pd.read_csv(p/'tables/TREATMENT_INTERACTIONS.tsv',sep='\t')
for _,rr in it.iterrows():
 d=pd.read_csv(p/f'tables/{rr.dataset.upper()}_ANALYSIS_INPUT.tsv',sep='\t')
 if rr.dataset=='CheckMate':d=d[d.Cohort=='CM-025'].copy()
 d['x']=(d[rr.score]-d[rr.score].mean())/d[rr.score].std()
 d['tx']=(d.arm==('NIVOLUMAB' if rr.dataset=='CheckMate' else 'Avelumab+Axitinib')).astype(float)
 d['x_tx']=d.x*d.tx
 cols=['x','tx','x_tx','age','sex']+(['risk'] if rr.dataset=='CheckMate' else [])
 d=d.dropna(subset=cols+[rr.endpoint+'_time',rr.endpoint+'_event'])
 ex=pd.get_dummies(d[cols],drop_first=True,dtype=float)
 f=PHReg(d[rr.endpoint+'_time'],ex,status=d[rr.endpoint+'_event'],ties='efron').fit(disp=0,maxiter=200)
 i=ex.columns.get_loc('x_tx');hr=np.exp(f.params[i]);pv=f.pvalues[i]
 assert abs(hr-rr.HR_ratio)<1e-7 and abs(pv-rr.P)<1e-7
 out.append(dict(dataset=rr.dataset,endpoint=rr.endpoint,score=rr.score,n=len(d),HR_ratio=hr,P=pv,max_abs_score=float(np.max(abs(f.model.score(f.params))))))
pd.DataFrame(out).to_csv(p/'qc/INDEPENDENT_INTERACTION_CHECK.tsv',sep='\t',index=False)
print('ADJUSTED_CHECK',pd.DataFrame(checks).to_string(index=False));print('INTERACTIONS_VERIFIED',len(out))
