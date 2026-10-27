from pathlib import Path
import json,hashlib,subprocess
import pandas as pd,numpy as np
from scipy import stats
P=Path(__file__).resolve().parents[1];R=P.parent;A=R/'phase291_spatial_null_TYMP_matched_regions_20260930_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gate(s):
 with (P/'logs/memory_gate.log').open('a') as f:
  f.write(s+'\n');f.flush();subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,check=True)
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=lambda v:v.item() if hasattr(v,'item') else str(v)))
dump(P/'qc/FREEZE.json',dict(protocol_sha256=sha(P/'PROTOCOL_FROZEN.md'),source_table_sha256=sha(A/'tables/EXACT_NULL_SECTION.tsv')))
old=pd.read_csv(A/'tables/EXACT_NULL_SECTION.tsv',sep='\t')
old=old[old.endpoint.isin(['joint_to_T_distribution','joint_to_CD8_distribution'])]
rows=[];checks=[];manifest={}
for (sample,radius),grp in old.groupby(['sample','radius']):
 gate(sample+'_'+str(radius))
 path=A/'inputs'/f'{sample}_{radius}_REGIONS.tsv.gz';e=pd.read_csv(path,sep='\t').set_index('edge');manifest[str(path)]=sha(path)
 for r in grp.itertuples():
  target='T' if r.endpoint=='joint_to_T_distribution' else 'CD8'
  for program in ['ECM_organization','ECM_remodeling']:
   out=dict(sample=sample,patient=r.patient,dataset=r.dataset,radius=radius,model=r.model,program=program,endpoint=target,effect=np.nan,spatial_P=np.nan,reason=r.reason,regions=r.regions)
   if r.reason!='OK':rows.append(out);continue
   f=A/'inputs'/f'{sample}_{radius}_{r.model}_{r.endpoint}_SPECTRAL.npz';g=A/'inputs'/f'{sample}_{radius}_{r.model}_{r.endpoint}_EXACT_NULL.npz'
   z=np.load(f);v=np.load(g);D=z['D'];W=z['W'];yr=z['yr'];V=v['V']
   zz=e.loc[z['edge']];raw=zz['out_p_'+program].to_numpy()
   if not np.isfinite(raw).all():out['reason']='NONFINITE_EXPOSURE';rows.append(out);continue
   ranked=stats.rankdata(raw);xr=ranked-D@np.linalg.lstsq(D,ranked,rcond=None)[0]
   if np.linalg.norm(xr)<1e-8:out['reason']='CONSTANT_EXPOSURE_RESIDUAL';rows.append(out);continue
   coef=V.T@xr;assert np.max(abs(V@coef-xr))<1e-7
   seed=int(hashlib.sha256((sample+str(radius)+r.model+target+program).encode()).hexdigest()[:8],16)
   signs=np.random.default_rng(seed).choice([-1.,1.],size=(len(coef),999))
   denom=np.linalg.norm(xr)*np.linalg.norm(yr);null=(V.T@yr)@(coef[:,None]*signs)/denom
   obs=float(xr@yr/denom);p=float((1+(abs(null)>=abs(obs)-1e-12).sum())/1000)
   # Recompute observed endpoint residual from actual region readout.
   ycol='T_delta' if target=='T' else 'CD8_delta'
   yy=stats.rankdata(zz[ycol]);yr2=yy-D@np.linalg.lstsq(D,yy,rcond=None)[0]
   assert np.allclose(yr,yr2,atol=1e-8)
   assert abs(np.corrcoef(xr,yr2)[0,1]-obs)<1e-10
   mor=(xr@W@xr)/(xr@xr);mc=nc=cc=0.
   for k in [0,499,998]:
    sur=V@(coef*signs[:,k]);mc=max(mc,abs((sur@W@sur)/(sur@sur)-mor));nc=max(nc,abs(sur@sur/(xr@xr)-1));cc=max(cc,abs(D.T@sur).max())
    assert abs(yr@sur/denom-null[k])<1e-10
   assert mc<1e-7 and nc<1e-8 and cc<1e-6
   out.update(effect=obs,spatial_P=p,reason='OK',outer_T_positive_fraction=float(zz.T_outside.gt(0).mean()),inner_T_mean=float(zz.T_inside.mean()),outer_T_mean=float(zz.T_outside.mean()),outside_ECM_median=float(np.median(raw)),null_q025=float(np.quantile(null,.025)),null_q975=float(np.quantile(null,.975)))
   checks.append(dict(sample=sample,radius=radius,model=r.model,program=program,endpoint=target,max_moran_error=mc,max_norm_error=nc,max_cov_error=cc))
   rows.append(out)
 print('SECTION_RADIUS_COMPLETE',sample,radius,flush=True)
sec=pd.DataFrame(rows);sec['spatial_FDR']=np.nan
for _,ix in sec.groupby(['sample','radius','model']).groups.items():
 vals=sec.loc[ix,'spatial_P'];sec.loc[ix,'spatial_FDR']=np.where(vals.notna(),stats.false_discovery_control(vals.fillna(1)),np.nan)
sec.to_csv(P/'tables/SECTION_RESULTS.tsv',sep='\t',index=False)
keys=['dataset','radius','model','program','endpoint'];ps=[]
for key,z in sec.groupby(keys+['patient']):
 vals=z.effect.dropna().to_numpy();val=float(np.tanh(np.arctanh(np.clip(vals,-.999999999,.999999999)).mean())) if len(vals) else np.nan
 ps.append(dict(zip(keys+['patient'],key),effect=val,sections=len(vals)))
pa=pd.DataFrame(ps);pa.to_csv(P/'tables/PATIENT_RESULTS.tsv',sep='\t',index=False)
rs=[]
for key,z in pa.groupby(keys):
 a=z.effect.dropna().to_numpy();n=len(a);p=float(stats.wilcoxon(a,alternative='two-sided',method='auto').pvalue) if n>=10 and np.any(a!=0) else (1. if n>=10 else np.nan)
 rs.append(dict(zip(keys,key),n=n,negative=int((a<0).sum()),positive=int((a>0).sum()),median_rho=float(np.median(a)) if n else np.nan,P=p))
out=pd.DataFrame(rs);out['FDR']=np.nan
for _,ix in out.groupby(['dataset','radius','model']).groups.items():
 p=out.loc[ix,'P'];assert len(ix)==4
 out.loc[ix,'FDR']=np.where(p.notna(),stats.false_discovery_control(p.fillna(1)),np.nan)
out.to_csv(P/'tables/COHORT_RESULTS.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(P/'qc/MODEL_VALIDATION.tsv',sep='\t',index=False)
assert all(sha(Path(p))==h for p,h in manifest.items())
dump(P/'qc/COMPLETION.json',dict(status='COMPUTED_VALIDATED',section_models=len(sec),valid_models=len(checks),source_regions_hashes=manifest,patient_tests='n>=10 only;4 hypotheses per cohort radius model',source_unchanged=True))
print(out.to_string(index=False));print('PHASE298_COMPUTED_VALIDATED')
