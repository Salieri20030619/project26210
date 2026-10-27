from pathlib import Path
import pandas as pd,numpy as np,json,warnings
from scipy.stats import chi2
from statsmodels.duration.hazard_regression import PHReg
p=Path(__file__).resolve().parents[1]
m=pd.read_csv(p/'tables/C4_MODEL_MEMBERSHIP.tsv',sep='\t')
r=pd.read_csv(p/'tables/C4_PFS_TIME_VARYING_MODELS.tsv',sep='\t')
im=pd.read_csv(p/'tables/C4_INTERACTION_MEMBERSHIP.tsv',sep='\t')
ir=pd.read_csv(p/'tables/C4_TIME_VARYING_TREATMENT_INTERACTIONS.tsv',sep='\t')
def fit_risksets(d,cols,ttcols):
 base=pd.get_dummies(d[cols],drop_first=True,dtype=float);base=base.loc[:,base.nunique()>1]
 xx=[];yy=[];ss=[]
 for j,t in enumerate(np.unique(d.loc[d.event==1,'time'])):
  risk=d.time>=t;at=base.loc[risk].copy()
  for v in ttcols:at[v+'_logtime']=d.loc[risk,v]*np.log1p(t)
  xx.append(at);yy.extend(((d.loc[risk,'time']==t)&(d.loc[risk,'event']==1)).astype(int));ss.extend([j]*sum(risk))
 X=pd.concat(xx,ignore_index=True);ev=np.asarray(yy);strat=np.asarray(ss)
 with warnings.catch_warnings(record=True) as ws:
  warnings.simplefilter('always');fit=PHReg(np.ones(len(X)),X,status=ev,strata=strat,ties='efron').fit(disp=0,maxiter=200)
 return fit,X.columns,';'.join(str(w.message) for w in ws)
out=[]
for _,q in r.iterrows():
 d=m[m.model_id==q.model_id].copy().reset_index(drop=True);cv=['age','sex','trial','risk']+(['lit','mac'] if q.model=='plus_Jiang5_macro' else [])
 f,cols,w=fit_risksets(d,['x']+cv,['x']);f0,_,_=fit_risksets(d,cv,[]);fc,_,_=fit_risksets(d,['x']+cv,[])
 b0=f.params[cols.get_loc('x')];b1=f.params[cols.get_loc('x_logtime')]
 jp=chi2.sf(2*(f.llf-f0.llf),2);tp=chi2.sf(2*(f.llf-fc.llf),1)
 assert max(abs(b0-q.beta_0),abs(b1-q.beta_log1p_month),abs(jp-q.joint_P),abs(tp-q.time_dependence_P))<1e-6
 out.append({'type':'within_arm','id':int(q.model_id),'beta0_diff':abs(b0-q.beta_0),'beta1_diff':abs(b1-q.beta_log1p_month),'joint_P_diff':abs(jp-q.joint_P),'warnings':w})
 print('VERIFIED_TIME_MODEL',int(q.model_id),flush=True)
for _,q in ir.iterrows():
 d=im[im.interaction_id==q.interaction_id].copy().reset_index(drop=True);d['xt']=d.x*d.tx;d['lt']=d.lit*d.tx
 cv=['age','sex','risk']+(['lit','lt','mac'] if q.model!='clinical' else [])
 f,cols,w=fit_risksets(d,['x','tx','xt']+cv,['x','xt']);f0,_,_=fit_risksets(d,['x','tx']+cv,['x'])
 b0=f.params[cols.get_loc('xt')];b1=f.params[cols.get_loc('xt_logtime')];jp=chi2.sf(2*(f.llf-f0.llf),2)
 assert max(abs(b0-q.beta_interaction_0),abs(b1-q.beta_interaction_log1p_month),abs(jp-q.joint_interaction_P))<1e-6
 out.append({'type':'interaction','id':int(q.interaction_id),'beta0_diff':abs(b0-q.beta_interaction_0),'beta1_diff':abs(b1-q.beta_interaction_log1p_month),'joint_P_diff':abs(jp-q.joint_interaction_P),'warnings':w})
 print('VERIFIED_TIME_INTERACTION',int(q.interaction_id),flush=True)
pd.DataFrame(out).to_csv(p/'qc/INDEPENDENT_TIME_VARIATION_CHECK.tsv',sep='\t',index=False)
(p/'qc/INDEPENDENT_TIME_VARIATION_COMPLETE.json').write_text(json.dumps({'independent_method':'explicit event-time risk-set expansion with Efron ties; statsmodels PHReg vs R coxph tt','models_verified':len(out),'tolerance':1e-6},indent=2))
print('ALL_TIME_MODELS_VERIFIED',len(out))
