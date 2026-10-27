from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import chi2,spearmanr
from statsmodels.duration.hazard_regression import PHReg
p=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase270_SPP1_literature_clinical_revalidation_20260927_v1')
h=p.parent.parent
r=pd.read_csv(p/'tables/ALL_SURVIVAL_MODELS.tsv',sep='\t')
m=pd.read_csv(p/'tables/MODEL_MEMBERSHIP.tsv',sep='\t')
inputs={s:pd.read_csv(p/f'tables/{s.upper()}_ANALYSIS_INPUT.tsv',sep='\t') for s in ['CheckMate','JAVELIN']}
checks=[]
for ds,d in inputs.items():
 e=pd.read_csv(p/f'inputs/{ds.lower()}_genes.tsv',sep='\t',index_col=0)
 if ds=='JAVELIN':e=np.log2(e+1)
 s=e.loc[['SPP1','FN1','CLEC5A','OLR1','IL1RN']].mean()
 np.testing.assert_allclose(d.Jiang5,s[d.sample_id],rtol=1e-12)
 np.testing.assert_allclose(d.SPP1,e.loc['SPP1',d.sample_id],rtol=1e-12)
 for _,rr in r[r.dataset==ds].iterrows():
  z=m[m.model_id==rr.model_id].merge(d,on='sample_id',validate='one_to_one')
  assert len(z)==rr.n and z.event.sum()==rr.events
  dd=d[d.arm==rr.arm]
  expected=((dd[rr.score]>=rr.cutoff).astype(float) if rr.model.startswith('median') else (dd[rr.score]-dd[rr.score].mean())/dd[rr.score].std())
  np.testing.assert_allclose(z.x,pd.Series(expected.values,index=dd.sample_id)[z.sample_id],atol=1e-12)
  cols=['x']
  if rr.model=='continuous_adjusted':
   cols+=['age','sex']+(['trial','risk'] if ds=='CheckMate' else [])
  ex=pd.get_dummies(z[cols],drop_first=True,dtype=float);ex=ex.loc[:,ex.nunique()>1]
  fit=PHReg(z.time,ex,status=z.event,ties='efron').fit(disp=0)
  j=ex.columns.get_loc('x');diff=abs(np.exp(fit.params[j])-rr.HR);pdiff=abs(fit.pvalues[j]-rr.P)
  assert diff<1e-7 and pdiff<1e-7,(rr.model_id,diff,pdiff)
  ld=None
  if rr.model.startswith('median'):
   oe=0.;v=0.
   for t in np.unique(z.loc[z.event==1,'time']):
    risk=z[z.time>=t];ev=z[(z.time==t)&(z.event==1)]
    n=len(risk);nh=risk.x.sum();de=len(ev);dh=ev.x.sum()
    oe+=dh-de*nh/n
    if n>1:v+=de*(nh/n)*(1-nh/n)*(n-de)/(n-1)
   lp=chi2.sf(oe*oe/v,1);ld=abs(lp-rr.logrank_P);assert ld<1e-8
  checks.append(dict(model_id=int(rr.model_id),hr_difference=diff,p_difference=pdiff,logrank_difference=ld))
pd.DataFrame(checks).to_csv(p/'qc/INDEPENDENT_MODEL_CHECK.tsv',sep='\t',index=False)
corr=[]
old1=pd.read_csv(h/'project26210_phase81_bulk_ici_transcriptome_rebuild_v1/results_v1/phase81_sample_scores_v1.tsv',sep='\t')
old3=pd.read_csv(h/'project26210_phase83_bulk_state_projection_v1/results_v1/clinical_v1/phase83_sample_scores_v1.tsv',sep='\t')
for ds,d in inputs.items():
 z=d[['sample_id','SPP1','Jiang5','TYMP']].merge(old1,on='sample_id',validate='one_to_one').merge(old3.drop(columns=['dataset','tymp_z']),on='sample_id',validate='one_to_one')
 assert len(z)==len(d)
 for a in ['SPP1','Jiang5']:
  for b in ['TYMP','spp1_macrophage','full_chain_proxy','spp1_state_weight','spp1_macrophage_fraction','total_macrophage_weight','official_macro']:
   zz=z[[a,b]].dropna();rho=spearmanr(zz[a],zz[b]).statistic
   corr.append(dict(dataset=ds,score=a,other=b,n=len(zz),spearman_rho=rho))
pd.DataFrame(corr).to_csv(p/'tables/CONSTRUCT_CORRELATIONS.tsv',sep='\t',index=False)
(p/'qc/INDEPENDENT_CHECK_COMPLETE.json').write_text(json.dumps({'models_verified':len(checks),'max_hr_difference':max(x['hr_difference'] for x in checks),'max_p_difference':max(x['p_difference'] for x in checks),'score_formulas_verified':True},indent=2))
print('VERIFIED',len(checks));print(pd.DataFrame(corr).to_string(index=False))
