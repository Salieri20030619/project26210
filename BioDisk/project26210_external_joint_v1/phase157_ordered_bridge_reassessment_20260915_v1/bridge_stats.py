import numpy as np,pandas as pd
from scipy.stats import rankdata,t
def rankz(v):
 a=rankdata(np.asarray(v,float));s=a.std(ddof=1)
 return (a-a.mean())/s if s>0 else a*0
def residual(v,C):return v-C@np.linalg.lstsq(C,v,rcond=None)[0]
def corrp(x,y,C):
 a=residual(x,C);b=residual(y,C)
 if a.std()<1e-12 or b.std()<1e-12:return np.nan,np.nan
 r=float(np.corrcoef(a,b)[0,1]);df=len(x)-np.linalg.matrix_rank(C)-1
 p=float(2*t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-15)),df))
 return r,p
def bridge(A,M,Y,studies,adjustment):
 A=rankz(A);M=rankz(M);Y=rankz(Y);n=len(A)
 C=np.ones((n,1)) if adjustment=='raw' else pd.get_dummies(pd.Series(studies),dtype=float).to_numpy()
 a,pa=corrp(A,M,C);my,pmy=corrp(M,Y,C);b,pb=corrp(M,Y,np.column_stack([C,A]))
 ca=np.linalg.lstsq(np.column_stack([C,A]),M,rcond=None)[0][-1]
 cb=np.linalg.lstsq(np.column_stack([C,A,M]),Y,rcond=None)[0]
 total=np.linalg.lstsq(np.column_stack([C,A]),Y,rcond=None)[0][-1]
 return dict(n_patients=n,a_rho=a,a_P=pa,XY_rho=my,XY_P=pmy,b_partial_rho=b,b_P=pb,a_beta=float(ca),b_beta=float(cb[-1]),indirect_ab=float(ca*cb[-1]),direct_cprime=float(cb[-2]),total_c=float(total),identity_error=float(total-cb[-2]-ca*cb[-1]))
def bootstrap_ab(A,M,Y,studies,adjustment,seed=262157,B=2000):
 A=rankz(A);M=rankz(M);Y=rankz(Y);n=len(A);rng=np.random.default_rng(seed)
 groups=[np.flatnonzero(np.asarray(studies)==s) for s in sorted(set(studies))]
 ix=np.tile(np.arange(n),(B,1))
 for ids in groups:ix[:,ids]=rng.choice(ids,size=(B,len(ids)),replace=True)
 aa=A[ix].copy();mm=M[ix].copy();yy=Y[ix].copy()
 center=[np.arange(n)] if adjustment=='raw' else groups
 for ids in center:
  aa[:,ids]-=aa[:,ids].mean(axis=1)[:,None];mm[:,ids]-=mm[:,ids].mean(axis=1)[:,None];yy[:,ids]-=yy[:,ids].mean(axis=1)[:,None]
 av=(aa*aa).sum(axis=1);mv=(mm*mm).sum(axis=1);am=(aa*mm).sum(axis=1);ay=(aa*yy).sum(axis=1);my=(mm*yy).sum(axis=1);det=av*mv-am*am
 ok=(av>1e-12)&(det>1e-10);ab=(am[ok]/av[ok])*(my[ok]*av[ok]-am[ok]*ay[ok])/det[ok]
 if not len(ab):return dict(bootstrap_valid=0,ab_low=np.nan,ab_high=np.nan,bootstrap_positive_fraction=np.nan)
 return dict(bootstrap_valid=len(ab),ab_low=float(np.quantile(ab,.025)),ab_high=float(np.quantile(ab,.975)),bootstrap_positive_fraction=float((ab>0).mean()))
