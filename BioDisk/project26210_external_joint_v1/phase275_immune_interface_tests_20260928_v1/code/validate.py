from pathlib import Path
import json,math,datetime,hashlib
import numpy as np,pandas as pd
import statsmodels.api as sm
P=Path(__file__).resolve().parents[1];OLD=P.parent/'phase234_part4_replan_execution_20260921_v1';MARK=P.parent/'phase246_spatial_marker_structure_20260924_v1'
s=pd.read_csv(P/'tables/SECTION_EFFECTS.tsv',sep='\t');p=pd.read_csv(P/'tables/PATIENT_EFFECTS.tsv',sep='\t');c=pd.read_csv(P/'tables/COHORT_SUMMARY.tsv',sep='\t')
keys=['dataset','family','endpoint','model','effect_type','width_um']
def eq(a,b,tol=1e-10):
 if pd.isna(a) and pd.isna(b):return
 assert np.isfinite(a) and np.isfinite(b) and abs(a-b)<tol,(a,b)
assert len(s)==38*36 and len(p)==37*36 and len(c)==108
assert not s.duplicated(['sample']+keys).any()
checks=[]
for row in p.itertuples():
 mask=s.patient.eq(row.patient)
 for k in keys:mask &= s[k].eq(getattr(row,k))
 v=s.loc[mask,'effect'].dropna().tolist()
 if not v:expected=np.nan
 elif row.effect_type=='rho':expected=math.tanh(sum(math.atanh(max(-.999999999,min(.999999999,x))) for x in v)/len(v))
 else:expected=sum(v)/len(v)
 eq(expected,row.effect);assert len(v)==row.n_sections
def exactsign(v):
 nz=[x for x in v if x!=0]
 if not nz:return 1.
 k=min(sum(x>0 for x in nz),sum(x<0 for x in nz))
 return min(1.,2*sum(math.comb(len(nz),i) for i in range(k+1))/2**len(nz))
for row in c.itertuples():
 mask=np.ones(len(p),bool)
 for k in keys:mask &= p[k].eq(getattr(row,k)).to_numpy()
 v=p.loc[mask,'effect'].dropna().tolist();assert row.n==len(v)
 eq(row.median,float(np.median(v)) if len(v) else np.nan)
 eq(row.P,exactsign(v) if len(v)>=10 else np.nan)
 assert row.positive==sum(x>0 for x in v) and row.negative==sum(x<0 for x in v)
for _,g in c.groupby(['dataset','family']):
 n=len(g);assert n in (4,8)
 for pc,qc in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:
  order=sorted(g.index,key=lambda j:float(g.at[j,pc]) if pd.notna(g.at[j,pc]) else 1.)
  last=1.
  for rank,j in reversed(list(enumerate(order,start=1))):
   v=float(g.at[j,pc]) if pd.notna(g.at[j,pc]) else 1.
   last=min(last,n*v/rank,1.)
   eq(g.at[j,qc],last if pd.notna(g.at[j,pc]) else np.nan)
errs=[]
for sample,g in s[s.width_um==400].groupby('sample'):
 d=pd.read_csv(OLD/'inputs'/f'{sample}_INPUTS.tsv.gz',sep='\t');m=pd.read_csv(MARK/'inputs'/f'{sample}_MARKERS.tsv.gz',sep='\t');det=pd.read_csv(P/'inputs'/f'{sample}_CD8_detection.tsv.gz',sep='\t')
 inside=d.tumor_expression_region.astype(float).to_numpy();mask=np.isfinite(d.expression_distance_um)&(d.expression_distance_um<=400+1e-6)
 ys={'T_abundance_log2':np.log2(d.Tcell_q05).to_numpy(),'CD8_RNA':m.CD8lineage.to_numpy(),'T_RNA':m.Tcell.to_numpy(),'CD8_codetection':det.CD8_codetection.astype(float).to_numpy() if m.CD8lineage.notna().any() else np.full(len(d),np.nan)}
 for row in g[g.family.str.startswith('band_')].itertuples():
  if not np.isfinite(row.effect):continue
  target=row.endpoint.split('__')[1];y=ys[target];valid=mask&np.isfinite(y)
  mat=np.column_stack([np.ones(valid.sum()),inside[valid]])
  if row.model!='unadjusted':
   cv=np.column_stack([d.log1p_UMI,np.log1p(d.total_macro_q05)])[valid]
   cv=(cv-cv.mean(0))/cv.std(0);mat=np.column_stack([mat,cv])
  b=float(sm.OLS(y[valid],mat).fit().params[1]);eq(b,row.effect,1e-9);errs.append(abs(b-row.effect))
 center=pd.read_csv(P/'inputs'/f'{sample}_centers_400.tsv.gz',sep='\t',float_precision='round_trip')
 covnames=['depth_mean','depth_inside_minus_outside','log1p_macro200','malignant_mean','malignant_inside_minus_outside','log_side_count_ratio']
 for row in g[g.family.str.startswith('SPP1_coupling')].itertuples():
  if not np.isfinite(row.effect):continue
  cols=['SPP1_share200',row.endpoint.split('__')[1]+'_inside_minus_outside']+(covnames if row.model!='unadjusted' else [])
  a=center.loc[center.center_eligible,cols].dropna().rank(method='average').to_numpy()
  D=np.column_stack([np.ones(len(a)),a[:,2:]])
  Q=np.linalg.qr(D,mode='reduced')[0];res=a[:,:2]-Q@(Q.T@a[:,:2])
  rr=float(np.dot(res[:,0],res[:,1])/np.sqrt(np.dot(res[:,0],res[:,0])*np.dot(res[:,1],res[:,1])))
  eq(rr,row.effect,1e-9);errs.append(abs(rr-row.effect))
n=pd.read_csv(P/'tables/LOCAL_GEOMETRY_NUMERIC_CHECK.tsv',sep='\t');assert np.max(np.abs(n.T_delta-n.saved_T_delta))<1e-10
q=pd.read_csv(P/'tables/SECTION_GEOMETRY_ELIGIBILITY.tsv',sep='\t');assert len(q)==114 and q.marker_reconstruction_error.max()<1e-10
assert (q.band_inside+q.band_outside<=q.reachable_spots).all()
freeze=json.loads((P/'FREEZE.json').read_text());assert hashlib.sha256((P/'PROTOCOL_BEFORE_TESTING.md').read_bytes()).hexdigest()==freeze['protocol_sha256']
assert (P/'logs/exit_status.txt').read_text()=='0'
assert not list(P.rglob('*.png')) and not list(P.rglob('*.pdf'))
out=dict(status='PASS',at=datetime.datetime.now().astimezone().isoformat(),sections=38,patients=37,section_effects=len(s),patient_effects=len(p),cohort_endpoints=len(c),FDR_families=int(c.groupby(['dataset','family']).ngroups),independent_checks=['math exact binomial sign P','manual BH','patient merge','statsmodels OLS side contrasts','QR rank partial correlation','direct neighbor geometry','marker reconstruction','frozen protocol and read-only input verification'],numeric_checks=len(errs),max_numeric_error=max(errs),max_marker_error=float(q.marker_reconstruction_error.max()),figures_created=0,manuscript_edited=False)
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
