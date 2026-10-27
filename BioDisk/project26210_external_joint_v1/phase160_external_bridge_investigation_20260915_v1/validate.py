from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import rankdata,t,norm
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase160_external_bridge_investigation_20260915_v1'
plan=json.loads((O/'PLAN_BEFORE_ASSOCIATIONS.json').read_text());R=pd.read_csv(O/'COHORT_BRIDGE_RESULTS.tsv',sep='\t');M=pd.read_csv(O/'META_BRIDGE_RESULTS.tsv',sep='\t')
checks=[]
for cohort,g in R.groupby('cohort'):
 f=pd.read_csv(O/('PATIENT_VALUES_'+cohort+'.tsv'),sep='\t',index_col=0)
 for ix,r in g.iterrows():
  cov=plan['adjustments'][r.adjustment]
  d=f[['TYMP','SPP1',r.candidate]+cov].dropna()
  for e,cols in [('a',['TYMP',r.candidate]+cov),('XY',[r.candidate,'SPP1']+cov),('b',[r.candidate,'SPP1','TYMP']+cov)]:
   mat=np.array([rankdata(d[c]) for c in cols]);precision=np.linalg.inv(np.corrcoef(mat));rho=-precision[0,1]/np.sqrt(precision[0,0]*precision[1,1]);k=len(cols)-2;p=2*t.sf(abs(rho)*np.sqrt((len(d)-k-2)/(1-rho*rho)),len(d)-k-2)
   checks.append(dict(check='precision_'+str(ix)+'_'+e,passed=abs(rho-r[e+'_rho'])<1e-9 and abs(p-r[e+'_P'])<1e-9))
for (family,adj,e),ix in M.groupby(['family','adjustment','edge']).groups.items():
 for model in ['fixed','random']:
  q=multipletests(M.loc[ix,model+'_P'],method='fdr_bh')[1]
  checks.append(dict(check='metaBH_'+family+adj+e+model,passed=np.allclose(q,M.loc[ix,model+'_FDR'])))
# Reproduce the historical marginal measured-MTA result, same source universe.
r=M[(M.candidate=='MET:5-methylthioadenosine (MTA)')&(M.adjustment=='raw')&(M.edge=='a')].iloc[0]
checks.append(dict(check='historical_MTA_fixed_rho',passed=abs(r.fixed_rho-0.039467)<1e-5))
for _,row in M.iterrows():
 g=R[(R.candidate==row.candidate)&(R.adjustment==row.adjustment)]
 e=row.edge;z=np.arctanh(g[e+'_rho'].to_numpy());v=1/(g.n.to_numpy()-g[e+'_k'].to_numpy()-3)
 w=1/v;fixed=np.average(z,weights=w);q=np.sum(w*(z-fixed)**2);tau=max(0,(q-(len(g)-1))/(w.sum()-np.sum(w*w)/w.sum())) if len(g)>1 else 0.
 wr=1/(v+tau);est=np.average(z,weights=wr);p=2*norm.sf(abs(est)*np.sqrt(wr.sum()))
 checks.append(dict(check='independent_meta_'+row.candidate+row.adjustment+e,passed=abs(np.tanh(est)-row.random_rho)<1e-10 and abs(p-row.random_P)<1e-10))
V=pd.DataFrame(checks);assert V.passed.all();V.to_csv(O/'INDEPENDENT_VALIDATION.tsv',sep='\t',index=False)
status=json.loads((O/'STATUS.json').read_text());status.update(independent_checks=len(V),independent_passed=int(V.passed.sum()),total_checks=status['checks']+len(V));(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print('SLC16A3 COHORTS\n'+R[R.candidate=='SLC16A3'][['cohort','adjustment','n','a_rho','a_P','a_FDR','XY_rho','XY_P','b_rho','b_P','b_FDR','indirect_ab']].to_string(index=False))
L=pd.read_csv(O/'LEAVE_ONE_COHORT_META.tsv',sep='\t')
print('SLC16A3 LOO\n'+L[(L.candidate=='SLC16A3')&(L.adjustment=='myeloid_hypoxia')][['omitted','edge','random_rho','random_P','positive_cohorts']].to_string(index=False))
print('MEASURED LAC/CHOL\n'+R[R.candidate.isin(['MET:lactate','MET:cholesterol'])][['cohort','candidate','adjustment','n','a_rho','a_P','XY_rho','XY_P','b_rho','b_P']].to_string(index=False))
print('VALIDATION',json.dumps(status))

