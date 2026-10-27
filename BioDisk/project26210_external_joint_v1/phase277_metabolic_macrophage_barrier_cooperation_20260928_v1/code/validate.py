from pathlib import Path
import numpy as np,pandas as pd,math,json,hashlib,datetime
import statsmodels.api as sm
from scipy import sparse
P=Path(__file__).resolve().parents[1];R=P.parent;OLD=R/'phase234_part4_replan_execution_20260921_v1';GEO=R/'phase276_interface_topology_and_localization_20260928_v1'
read=lambda f:pd.read_csv(f,sep='\t',float_precision='round_trip')
s=read(P/'tables/SECTION_EFFECTS.tsv');pat=read(P/'tables/PATIENT_EFFECTS.tsv');su=read(P/'tables/COHORT_SUMMARY.tsv');details=read(P/'tables/INTERACTION_MODEL_DETAILS.tsv')
keys=['dataset','family','endpoint','model','effect_type']
def eq(a,b,tol=1e-9):
 if pd.isna(a) and pd.isna(b):return
 assert np.isfinite(a) and np.isfinite(b) and abs(a-b)<tol,(a,b)
assert len(s)==38*34 and len(pat)==37*34 and len(su)==102 and len(details)==38*24
for r in pat.itertuples():
 mask=s.patient.eq(r.patient)
 for k in keys:mask &= s[k].eq(getattr(r,k))
 v=s.loc[mask,'effect'].dropna().tolist()
 e=np.nan if not v else math.tanh(sum(math.atanh(max(-.999999999,min(.999999999,x))) for x in v)/len(v)) if r.effect_type=='rho' else sum(v)/len(v)
 eq(r.effect,e);assert len(v)==r.n_sections
def signP(v):
 vv=[x for x in v if x!=0];n=len(vv)
 if not n:return 1.
 k=min(sum(x>0 for x in vv),sum(x<0 for x in vv))
 return min(1.,2*sum(math.comb(n,j) for j in range(k+1))/2**n)
for r in su.itertuples():
 mask=np.ones(len(pat),bool)
 for k in keys:mask &= pat[k].eq(getattr(r,k)).to_numpy()
 v=pat.loc[mask,'effect'].dropna().tolist();assert r.n==len(v)
 eq(r.P,signP(v) if len(v)>=10 else np.nan);eq(r.median,np.median(v) if v else np.nan)
def bhcheck(frame):
 for pc,qc in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:
  order=sorted(frame.index,key=lambda j:frame.at[j,pc] if pd.notna(frame.at[j,pc]) else 1.);last=1.
  for rank,j in reversed(list(enumerate(order,1))):
   v=frame.at[j,pc] if pd.notna(frame.at[j,pc]) else 1.;last=min(last,len(frame)*v/rank,1.)
   eq(frame.at[j,qc],last if pd.notna(frame.at[j,pc]) else np.nan)
for _,g in su.groupby(['dataset','family']):
 assert len(g) in [4,6,8];bhcheck(g)
ip=read(P/'tables/C4_IDENTITY_PATIENT_EFFECTS.tsv');isu=read(P/'tables/C4_IDENTITY_SUMMARY.tsv')
for r in isu.itertuples():
 v=ip.loc[ip.program==r.program,'effect'].to_numpy();assert len(v)==r.n==24
 eq(r.P,signP(v));eq(r.median,np.median(v))
bhcheck(isu)
core=['depth_mean','depth_delta','macro_logmean','malignant_mean','malignant_delta','log_side_count_ratio']
errors=[];aggregation_errors=[];geometry_count=0
for sample,group in s.groupby('sample'):
 z=read(P/'inputs'/f'{sample}_MODEL_INPUT.tsv.gz');a=read(P/'inputs'/f'{sample}_SOURCE_MEMBERSHIP.tsv.gz');scores=read(P/'inputs'/f'{sample}_PROGRAMS.tsv.gz');d=read(OLD/'inputs'/f'{sample}_INPUTS.tsv.gz')
 for row in details[(details['sample']==sample)&details.effect.notna()].itertuples():
  cols=[row.macrophage,row.stroma,row.target+'_delta']+row.covariates.split(';');mask=z.source_eligible&np.isfinite(z[cols]).all(1)
  v=z.loc[mask,cols].to_numpy();assert len(v)==row.n>=50
  xx=(v[:,0]-v[:,0].mean())/v[:,0].std();ff=(v[:,1]-v[:,1].mean())/v[:,1].std();cc=(v[:,3:]-v[:,3:].mean(0))/v[:,3:].std(0)
  X=np.column_stack([np.ones(len(v)),xx,ff,xx*ff,cc]);b=sm.OLS(v[:,2],X).fit().params
  eq(b[3],row.effect,1e-8);eq(b[1],row.beta_M);eq(b[2],row.beta_F);eq(b[2]-b[3],row.F_slope_M_low);eq(b[2]+b[3],row.F_slope_M_high);errors.append(abs(b[3]-row.effect))
 for row in group[group.effect_type.eq('rho')&group.effect.notna()].itertuples():
  x,y=row.endpoint.split('__');cols=[x,y]+(core if row.model=='core_adjusted' else [])
  v=z.loc[z.source_eligible,cols].dropna().rank().to_numpy();X=np.column_stack([np.ones(len(v)),v[:,2:]]);Q=np.linalg.qr(X,mode='reduced')[0];res=v[:,:2]-Q@(Q.T@v[:,:2]);rho=float(np.corrcoef(res.T)[0,1])
  eq(rho,row.effect,1e-8);errors.append(abs(rho-row.effect))
 for sid,g in a.groupby('segment'):
  ids=g.spot_index.to_numpy(int);assert len(ids)==len(set(ids))
  for name in ['C4_19','C4_17','APOEMph_noC4','cholesterol_disjoint','polyamine_MTA_context']:
   v=scores[name].to_numpy()[ids];e=float(v.mean());eq(z.at[sid,name],e,1e-9)
   if np.isfinite(e):aggregation_errors.append(abs(z.at[sid,name]-e))
  eq(z.at[sid,'SPP1_share'],float(d.V1_SPP1_q05.iloc[ids].sum()/d.total_macro_q05.iloc[ids].sum()),1e-9)
 nodes=read(GEO/'inputs'/f'{sample}_nodes.tsv.gz');mem=read(GEO/'inputs'/f'{sample}_membership.tsv.gz');D=sparse.load_npz(OLD/'graphs'/f'{sample}_distance800.npz').tocsr()
 for sid in sorted(set([0,int(z.segment.max())])):
  ends=nodes.iloc[mem.loc[mem.segment==sid,'node'].to_numpy(int)][['spot_a','spot_b']].to_numpy(int).ravel()
  expected=set(ends.tolist())
  for j in ends:
   ii=D.indices[D.indptr[j]:D.indptr[j+1]];vv=D.data[D.indptr[j]:D.indptr[j+1]]
   expected.update(ii[vv<=150+1e-6].tolist())
  observed=set(a.loc[a.segment==sid,'spot_index'].astype(int));assert expected==observed;geometry_count+=1
assert (P/'logs/exit_status.txt').read_text()=='0'
freeze=json.loads((P/'FREEZE.json').read_text());assert hashlib.sha256((P/'PROTOCOL_BEFORE_TESTING.md').read_bytes()).hexdigest()==freeze['protocol_sha256']
assert hashlib.sha256((P/'PROGRAMS_FROZEN.json').read_bytes()).hexdigest()==freeze['programs_sha256']
qc=read(P/'tables/SECTION_QC.tsv');assert qc.marker_reconstruction_error.max()<1e-9
assert not list(P.rglob('*.png')) and not list(P.rglob('*.pdf'))
out=dict(status='PASS',at=datetime.datetime.now().astimezone().isoformat(),sections=38,patients=37,spatial_endpoints=len(su),identity_endpoints=4,section_effects=len(s),patient_effects=len(pat),model_and_correlation_checks=len(errors),max_model_error=max(errors) if errors else None,source_aggregation_checks=len(aggregation_errors),max_source_error=max(aggregation_errors),independent_source_geometry_checks=geometry_count,checks=['statsmodels OLS target interaction and conditional slopes','QR partial rank correlation','patient aggregation','exact combinatorial sign test','manual BH','source200um graph membership','new score aggregation','cached score reconstruction and read-only inputs'],figures=0,manuscript_changes=0)
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
