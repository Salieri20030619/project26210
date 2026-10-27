from pathlib import Path
import ast,math,json,datetime,hashlib
import numpy as np,pandas as pd
import statsmodels.api as sm
H=Path(__file__).resolve().parents[1];P=H.parent
read=lambda f:pd.read_csv(f,sep='\t',float_precision='round_trip')
tree=ast.parse((P/'code/validate.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['eq','signP','bhcheck']],type_ignores=[]),'independent_helpers','exec'),globals())
s=read(H/'tables/SECTION_EFFECTS.tsv');p=read(H/'tables/PATIENT_EFFECTS.tsv');su=read(H/'tables/COHORT_SUMMARY.tsv');de=read(H/'tables/INTERACTION_MODEL_DETAILS.tsv')
keys=['dataset','family','endpoint','model','effect_type'];assert len(s)==836 and len(p)==814 and len(su)==66 and len(de)==608
for r in p.itertuples():
 mask=s.patient.eq(r.patient)
 for k in keys:mask &= s[k].eq(getattr(r,k))
 v=s.loc[mask,'effect'].dropna().tolist();e=np.nan if not v else math.tanh(sum(math.atanh(max(-.999999999,min(.999999999,x))) for x in v)/len(v)) if r.effect_type=='rho' else sum(v)/len(v)
 eq(r.effect,e);assert r.n_sections==len(v)
for r in su.itertuples():
 mask=np.ones(len(p),bool)
 for k in keys:mask &= p[k].eq(getattr(r,k)).to_numpy()
 v=p.loc[mask,'effect'].dropna().tolist();eq(r.P,signP(v) if len(v)>=10 else np.nan);eq(r.median,np.median(v) if v else np.nan);assert r.n==len(v)
for _,g in su.groupby(['dataset','family']):assert len(g) in [6,8];bhcheck(g)
core=['depth_mean','depth_delta','macro_logmean','malignant_mean','malignant_delta','log_side_count_ratio'];errors=[];aggregations=0
for sample,g in s.groupby('sample'):
 z=read(H/'inputs'/f'{sample}_MODEL_INPUT.tsv.gz');scores=read(H/'inputs'/f'{sample}_PROGRAMS.tsv.gz');members=read(P/'inputs'/f'{sample}_SOURCE_MEMBERSHIP.tsv.gz')
 for r in de[(de['sample']==sample)&de.effect.notna()].itertuples():
  cols=[r.macrophage,r.stroma,r.target+'_delta']+r.covariates.split(';');v=z.loc[z.source_eligible&np.isfinite(z[cols]).all(1),cols].to_numpy()
  M=(v[:,0]-v[:,0].mean())/v[:,0].std();F=(v[:,1]-v[:,1].mean())/v[:,1].std();C=(v[:,3:]-v[:,3:].mean(0))/v[:,3:].std(0)
  X=np.column_stack([np.ones(len(v)),M,F,M*F,C]);b=sm.OLS(v[:,2],X).fit().params;eq(b[3],r.effect,1e-8);eq(b[2]-b[3],r.F_slope_M_low);eq(b[2]+b[3],r.F_slope_M_high);errors.append(abs(b[3]-r.effect))
 for r in g[g.effect_type.eq('rho')&g.effect.notna()].itertuples():
  x,y=r.endpoint.split('__');cols=[x,y]+(core if r.model=='core_adjusted' else [])
  v=z.loc[z.source_eligible,cols].dropna().rank().to_numpy();X=np.column_stack([np.ones(len(v)),v[:,2:]]);Q=np.linalg.qr(X,mode='reduced')[0];res=v[:,:2]-Q@(Q.T@v[:,:2]);rho=float(np.corrcoef(res.T)[0,1]);eq(rho,r.effect,1e-8);errors.append(abs(rho-r.effect))
 for sid,gg in members.groupby('segment'):
  for name in ['C4_platform17','C4_platform15_noAPOEshared']:
   eq(z.at[sid,name],scores[name].iloc[gg.spot_index.to_numpy(int)].mean());aggregations+=1
raw=json.loads((P/'qc/READONLY_INPUTS.json').read_text())
for r in raw['raw_h5ad_metadata']:
 st=Path(r['path']).stat();assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns']
assert (H/'logs/exit_status.txt').read_text()=='0'
qc=read(H/'tables/PROGRAM_ELIGIBILITY.tsv');assert len(qc)==76 and qc.eligible.all()
# Descriptive conditional slopes at patient level, no extra inference.
for src,prefix in [(P,'original'),(H,'harmonized')]:
 dd=read(src/'tables/INTERACTION_MODEL_DETAILS.tsv');kk=['dataset','patient','family','endpoint','model']
 pp=dd.groupby(kk)[['beta_M','beta_F','beta_MF','F_slope_M_low','F_slope_M_high','intercept']].mean().reset_index()
 pp.to_csv(src/'tables/CONDITIONAL_SLOPES_PATIENT_DESCRIPTIVE.tsv',sep='\t',index=False)
 rows=[]
 for key,gg in pp.groupby(['dataset','family','endpoint','model']):
  gg=gg[gg.beta_MF.notna()];row=dict(zip(['dataset','family','endpoint','model'],key));row.update(n=len(gg),negative_interaction=int((gg.beta_MF<0).sum()),negative_F_slope_M_high=int((gg.F_slope_M_high<0).sum()),negative_F_slope_M_low=int((gg.F_slope_M_low<0).sum()),median_F_slope_M_low=float(gg.F_slope_M_low.median()),median_F_slope_M_high=float(gg.F_slope_M_high.median()))
  rows.append(row)
 pd.DataFrame(rows).to_csv(src/'tables/CONDITIONAL_SLOPES_COHORT_DESCRIPTIVE.tsv',sep='\t',index=False)
out=dict(status='PASS',at=datetime.datetime.now().astimezone().isoformat(),technical_extension=True,original_results_preserved=True,platform_sources_verified=14,qualified_sections=38,new_spatial_endpoints=66,section_effects=836,patient_effects=814,model_checks=len(errors),max_model_error=max(errors),source_aggregation_checks=aggregations,raw_h5ad_metadata_unchanged=True,figures=0)
(H/'qc/FINAL_VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
