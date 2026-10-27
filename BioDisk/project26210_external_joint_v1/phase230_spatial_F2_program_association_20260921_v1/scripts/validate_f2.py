from pathlib import Path
import json,math
import pandas as pd,numpy as np
from scipy import stats
P=Path(__file__).resolve().parents[1]
s=pd.read_csv(P/'tables/COHORT_SUMMARY.tsv',sep='\t');p=pd.read_csv(P/'tables/PATIENT_EFFECTS.tsv',sep='\t');a=pd.read_csv(P/'tables/SECTION_EFFECTS.tsv',sep='\t')
checks=[]
for row in s.itertuples():
 v=p[(p.dataset==row.dataset)&(p.target==row.target)&(p.model==row.model)].rho.dropna().to_numpy()
 assert len(v)==row.n and (v>0).sum()==row.positive and (v<0).sum()==row.negative
 assert abs(np.median(v)-row.median_rho)<1e-13
 if len(v)>=10:
  assert len(np.unique(abs(v)))==len(v) and (v!=0).all()
  ranks=stats.rankdata(abs(v)).astype(int);observed=int(ranks[v>0].sum());total=int(ranks.sum())
  dp=[0]*(total+1);dp[0]=1
  for w in ranks:
   for j in range(total,int(w)-1,-1):dp[j]+=dp[j-int(w)]
  lo=min(observed,total-observed);exact=min(1.,2*sum(dp[:lo+1])/2**len(v))
  sign=min(1.,2*sum(math.comb(len(v),k) for k in range(min(row.positive,row.negative)+1))/2**len(v))
  assert abs(exact-row.P)<1e-12 and abs(sign-row.P_sign)<1e-12
  checks.append(dict(dataset=row.dataset,target=row.target,model=row.model,P_exact_DP=exact,sign_exact=sign))
 else:assert pd.isna(row.P)
for (dataset,model),g in s.groupby(['dataset','model']):
 if g.P.notna().all():
  order=np.argsort(g.P.to_numpy());val=g.P.to_numpy()[order];q=np.minimum.accumulate((val*len(val)/np.arange(1,len(val)+1))[::-1])[::-1]
  np.testing.assert_allclose(g.FDR.to_numpy()[order],np.minimum(q,1),rtol=1e-12)
for row in p.itertuples():
 v=a[(a.patient==row.patient)&(a.dataset==row.dataset)&(a.target==row.target)&(a.model==row.model)].rho.dropna().to_numpy()
 assert row.n_sections==len(v)
 if len(v):assert abs(np.tanh(np.arctanh(v).mean())-row.rho)<1e-13
pd.DataFrame(checks).to_csv(P/'audit/PATIENT_STATISTICS_EXACT_CHECK.tsv',sep='\t',index=False)
(P/'FINAL_VALIDATION.json').write_text(json.dumps(dict(status='PASS',section_models=152,patient_model_rows=len(p),summary_rows=len(s),independent_exact_patient_tests=len(checks),FDR_verified=True,repeat_patient_aggregation_verified=True),indent=2))
print('PASS',len(checks))
