from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import rankdata,t,norm
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase163_polyamine_MTA_serial_receiver_20260915_v1';S=P/'phase162_MTA_polyamine_part1_reassessment_20260915_v1'
plan=json.loads((O/'PLAN_BEFORE_FIT.json').read_text());R=pd.read_csv(O/'COHORT_SERIAL_DECOMPOSITION.tsv',sep='\t');E=pd.read_csv(O/'COHORT_CONDITIONAL_EDGES.tsv',sep='\t');M=pd.read_csv(O/'META_CONDITIONAL_EDGES.tsv',sep='\t');draws=np.load(O/'PATIENT_BOOTSTRAP_SERIAL_DRAWS.npz')
checks=[]
for ci,cohort in enumerate(['RC18','RC20','CPTAC','CPTAC_val']):
 f=pd.read_csv(S/('PUBLIC_PATIENTS_'+cohort+'.tsv'),sep='\t',index_col=0)
 for xi,x in enumerate(plan['candidates']):
  for ai,(adj,cn) in enumerate(plan['adjustments'].items()):
   d=f[['TYMP',x,'MET:5-methylthioadenosine (MTA)','SPP1']+cn].dropna()
   rng=np.random.default_rng(262163+ci*100+xi*10+ai);ii=rng.integers(0,len(d),(2000,len(d)))
   expected=draws[cohort+'__'+x+'__'+adj]
   for rep in [0,1,17,1999]:
    dd=d.iloc[ii[rep]]
    def z(v):
     a=rankdata(v);return (a-a.mean())/a.std(ddof=1)
    a=z(dd.TYMP);xx=z(dd[x]);zz=z(dd['MET:5-methylthioadenosine (MTA)']);y=z(dd.SPP1)
    C=np.column_stack([np.ones(len(dd))]+[z(dd[c]) for c in cn])
    aa=np.linalg.lstsq(np.column_stack([C,a]),xx,rcond=None)[0][-1]
    dt=np.linalg.lstsq(np.column_stack([C,a,xx]),zz,rcond=None)[0][-1]
    b=np.linalg.lstsq(np.column_stack([C,a,xx,zz]),y,rcond=None)[0][-1]
    checks.append(dict(check='bootstrap_'+cohort+x+adj+str(rep),passed=abs(aa*dt*b-expected[rep])<1e-9))
for i,r in E.iterrows():
 df=r.n-r.k-2;p=2*t.sf(abs(r.rho)*np.sqrt(df/(1-r.rho**2)),df)
 checks.append(dict(check='partial_P_'+str(i),passed=abs(p-r.P)<1e-10))
for i,r in M.iterrows():
 d=E[(E.candidate==r.candidate)&(E.adjustment==r.adjustment)&(E.edge==r.edge)]
 z=np.arctanh(d.rho.to_numpy());v=1/(d.n.to_numpy()-d.k.to_numpy()-3);w=1/v;m=np.average(z,weights=w);q=np.sum(w*(z-m)**2);tau=max(0,(q-len(d)+1)/(w.sum()-(w*w).sum()/w.sum()));wr=1/(v+tau);mu=np.average(z,weights=wr)
 checks.append(dict(check='meta_'+str(i),passed=abs(np.tanh(mu)-r.random_rho)<1e-10))
for _,ix in M.groupby(['adjustment','edge']).groups.items():checks.append(dict(check='metaBH_'+str(min(ix)),passed=np.allclose(multipletests(M.loc[ix,'random_P'],method='fdr_bh')[1],M.loc[ix,'random_FDR'])))
V=pd.DataFrame(checks);assert V.passed.all();V.to_csv(O/'INDEPENDENT_VALIDATION.tsv',sep='\t',index=False)
status=json.loads((O/'STATUS.json').read_text());status.update(independent_checks=len(V),independent_passed=int(V.passed.sum()),total_checks=status['checks']+len(V));(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print('SERIAL CI SUMMARY\n'+R.assign(excludes_zero=(R.serial_low>0)|(R.serial_high<0),all_three_positive=(R.a>0)&(R.d>0)&(R.b>0)).groupby('adjustment')[['excludes_zero','all_three_positive']].sum().to_string())
print('META SENSITIVITY\n'+M[(M.edge.isin(['b','d','MTA_without_X','X_direct']))&(M.adjustment.isin(['raw','MTAP']))][['candidate','adjustment','edge','random_rho','random_P','random_FDR','positive_cohorts']].to_string(index=False))
print('STATUS',json.dumps(status))

