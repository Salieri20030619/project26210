from pathlib import Path
import pandas as pd,numpy as np,json,warnings
from scipy.stats import chi2,spearmanr
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
p=Path(__file__).resolve().parents[1];r=p.parent
s=pd.read_csv(p/'tables/C4_SURVIVAL_MODELS.tsv',sep='\t');m=pd.read_csv(p/'tables/C4_MODEL_MEMBERSHIP.tsv',sep='\t')
checks=[]
for _,q in s.iterrows():
 d=m[m.model_id==q.model_id].copy()
 cv=['age','sex']+(['trial','risk'] if q.dataset=='CheckMate' else [])
 ext={'clinical':[],'plus_SPP1':['spp'],'plus_Jiang5':['lit'],'plus_Jiang5_macro':['lit','mac']}[q.model]
 cols=['x']+cv+ext;ex=pd.get_dummies(d[cols],drop_first=True,dtype=float);ex=ex.loc[:,ex.nunique()>1]
 with warnings.catch_warnings(record=True) as ws:
  warnings.simplefilter('always');f=PHReg(d.time,ex,status=d.event,ties='efron').fit(disp=0,maxiter=200)
  f0=PHReg(d.time,ex.drop(columns='x'),status=d.event,ties='efron').fit(disp=0,maxiter=200)
 i=ex.columns.get_loc('x');hr=np.exp(f.params[i]);pv=f.pvalues[i];lp=chi2.sf(2*(f.llf-f0.llf),1)
 assert len(d)==q.n and d.event.sum()==q.events
 assert max(abs(hr-q.HR),abs(pv-q.P),abs(lp-q.increment_P))<1e-7
 checks.append({'model_id':q.model_id,'HR_difference':abs(hr-q.HR),'P_difference':abs(pv-q.P),'LRT_P_difference':abs(lp-q.increment_P),'max_gradient':float(np.max(abs(f.model.score(f.params)))),'warnings':';'.join(str(w.message) for w in ws)})
pd.DataFrame(checks).to_csv(p/'qc/INDEPENDENT_COX_CHECK.tsv',sep='\t',index=False)
it=pd.read_csv(p/'tables/C4_TREATMENT_INTERACTIONS.tsv',sep='\t');im=pd.read_csv(p/'tables/C4_INTERACTION_MEMBERSHIP.tsv',sep='\t');oc=[]
for _,q in it.iterrows():
 d=im[im.interaction_id==q.interaction_id].copy();d['xt']=d.x*d.tx;d['lt']=d.lit*d.tx
 cols=['x','tx','xt','age','sex']+(['risk'] if q.dataset=='CheckMate' else [])+(['lit','lt','mac'] if q.model!='clinical' else [])
 ex=pd.get_dummies(d[cols],drop_first=True,dtype=float)
 f=PHReg(d.time,ex,status=d.event,ties='efron').fit(disp=0,maxiter=200);f0=PHReg(d.time,ex.drop(columns='xt'),status=d.event,ties='efron').fit(disp=0,maxiter=200)
 i=ex.columns.get_loc('xt');hr=np.exp(f.params[i]);pv=f.pvalues[i];lp=chi2.sf(2*(f.llf-f0.llf),1)
 assert max(abs(hr-q.HR_ratio),abs(pv-q.P),abs(lp-q.increment_P))<1e-7
 oc.append({'interaction_id':q.interaction_id,'HR_difference':abs(hr-q.HR_ratio),'P_difference':abs(pv-q.P),'LRT_P_difference':abs(lp-q.increment_P)})
pd.DataFrame(oc).to_csv(p/'qc/INDEPENDENT_INTERACTION_CHECK.tsv',sep='\t',index=False)
frozen=json.loads((p/'qc/CLINICAL_SCORE_FREEZE.json').read_text());genes=frozen['genes_used'];cor=[]
for ds in ['CheckMate','JAVELIN']:
 d=pd.read_csv(p/f'tables/{ds.upper()}_ANALYSIS_INPUT.tsv',sep='\t')
 x=pd.read_csv(p/f'inputs/{ds}_C4_genes.tsv',sep='\t',index_col=0).loc[genes]
 if ds=='JAVELIN':x=np.log2(x+1)
 np.testing.assert_allclose(x.mean()[d.sample_id],d.C4_state,rtol=1e-12)
 for other in ['SPP1','Jiang5','TYMP','official_macro']:
  z=d[['C4_state',other]].dropna()
  cor.append({'dataset':ds,'comparison':other,'n':len(z),'rho':spearmanr(z.C4_state,z[other]).statistic})
pd.DataFrame(cor).to_csv(p/'tables/C4_CONSTRUCT_CORRELATIONS.tsv',sep='\t',index=False)
for group,df in s.groupby('family'):
 np.testing.assert_allclose(multipletests(df.P,method='fdr_bh')[1],df.FDR,atol=1e-12)
z=np.load(r/'phase200_union_C4_C10_C15_reassessment_20260918_v1/core/contrasts/tumor_C4.npz',allow_pickle=True)
pan=pd.read_csv(p/'tables/LOSO_TRAINING_GENES.tsv',sep='\t');pr=pd.read_csv(p/'tables/LOSO_PAIRED_STATE_VALIDATION.tsv',sep='\t');ex={'SPP1','TYMP','FN1','CLEC5A','OLR1','IL1RN'}
assert not ex.intersection(genes)
for _,q in pr.iterrows():
 i=np.where(z['patients']==q.patient)[0][0];gs=pan.loc[pan.heldout_study==q.study,'gene'];ii=[np.where(z['genes']==g)[0][0] for g in gs]
 assert abs(float(z['effects'][i,ii].mean())-q.delta)<1e-10
 # training gene tables contain training-only n; heldout study patients excluded
 t=pd.read_csv(p/f'tables/MARKERS_{q.study}.tsv.gz',sep='\t');assert t.patients.unique().tolist()==[int((z['study']!=q.study).sum())]
 sel=t[t.selected];assert (sel.P<.05).all() and not ex.intersection(sel.gene)
(p/'qc/INDEPENDENT_CHECK.json').write_text(json.dumps({'cox_models':len(checks),'interactions':len(oc),'score_formulas':'PASS','LOSO_patient_effects':'PASS','excluded_genes':'PASS','BH_families':'PASS'},indent=2))
print('INDEPENDENT_PASS',len(checks),len(oc));print(pd.DataFrame(cor).to_string(index=False))
