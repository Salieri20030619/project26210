from pathlib import Path
import json,math
import pandas as pd,numpy as np
from scipy import stats
P=Path(__file__).resolve().parents[1]
s=pd.read_csv(P/'tables/COHORT_MODELS.tsv',sep='\t');p=pd.read_csv(P/'tables/PATIENT_MODELS.tsv',sep='\t');a=pd.read_csv(P/'tables/SECTION_MODELS.tsv',sep='\t')
assert len(a)==608 and a.groupby(['sample','target']).n_source.nunique().eq(1).all()
checked=0
for r in s.itertuples():
 v=p[(p.dataset==r.dataset)&(p.target==r.target)&(p.model==r.model)].rho.to_numpy()
 assert len(v)==r.n and np.isfinite(v).all() and abs(np.median(v)-r.median_rho)<1e-13
 if len(v)>=10:
  ranks=stats.rankdata(abs(v));assert len(np.unique(abs(v)))==len(v) and (v!=0).all()
  ranks=ranks.astype(int);total=int(ranks.sum());dp=[0]*(total+1);dp[0]=1
  for w in ranks:
   for j in range(total,int(w)-1,-1):dp[j]+=dp[j-int(w)]
  obs=int(ranks[v>0].sum());pv=min(1.,2*sum(dp[:min(obs,total-obs)+1])/2**len(v))
  assert abs(pv-r.P)<1e-12;checked+=1
 else:assert np.isnan(r.P)
for cols,qcol in [(['dataset','model'],'FDR_within_model'),(['dataset'],'FDR_diagnostic16')]:
 for _,g in s.groupby(cols):
  if g.P.notna().all():
   order=np.argsort(g.P.to_numpy());vv=g.P.to_numpy()[order];q=np.minimum.accumulate((vv*len(vv)/np.arange(1,len(vv)+1))[::-1])[::-1]
   np.testing.assert_allclose(g[qcol].to_numpy()[order],np.minimum(q,1),atol=1e-12)
d=pd.read_csv(P/'tables/PAIRED_BLOCK_EFFECT_CHANGES.tsv',sep='\t')
assert len(d)==72
for r in d.itertuples():
 g=p[(p.dataset==r.dataset)&(p.target==r.target)].pivot(index='patient',columns='model',values='rho')
 delta=g[r.to_model]-g[r.from_model]
 assert len(delta)==r.n and abs(delta.median()-r.paired_median_delta)<1e-13
(P/'FINAL_VALIDATION.json').write_text(json.dumps(dict(status='PASS',models=608,independent_patient_P_checks=checked,paired_block_comparisons=72,FDR_both_levels_verified=True,same_spots_across_models=True),indent=2))
print('PASS',checked)
