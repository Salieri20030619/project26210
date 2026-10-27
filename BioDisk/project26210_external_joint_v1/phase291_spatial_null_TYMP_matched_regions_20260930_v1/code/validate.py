from analyze import *
def main():
 assert (P/'logs/analyze.exit').read_text().strip()=='0'
 se=read(P/'tables/SPECTRAL_SECTION.tsv');maxrho=0.;maxnull=0.;maxmi=0.;checked=0
 for r in se[se.reason.eq('OK')].itertuples():
  f=P/'inputs'/f'{r.sample}_{r.radius}_{r.model}_{r.endpoint}_SPECTRAL.npz'
  z=np.load(f);D=z['D'];U,sv,_=np.linalg.svd(D,full_matrices=False);U=U[:,sv>np.finfo(float).eps*max(D.shape)*sv[0]]
  xr=z['xr'];yr=z['yr'];rho=np.dot(xr,yr)/(np.linalg.norm(xr)*np.linalg.norm(yr));maxrho=max(maxrho,abs(rho-r.effect))
  coeff=z['a'][:,None]*z['signs'];V=z['V'];W=z['W']
  assert np.max(abs(W-W.T))<1e-12 and np.all(np.diag(W)==0)
  numer=(yr@V)@coeff
  proj=(U.T@V)@coeff
  den=np.sqrt(np.sum(coeff**2,axis=0)-np.sum(proj**2,axis=0))*np.linalg.norm(yr)
  null=numer/den;maxnull=max(maxnull,float(abs(null-z['null']).max()))
  assert abs((1+np.sum(abs(null)>=abs(rho)))/1000-r.spatial_p)<1e-12
  for j in [0,499,998]:
   surr=V@coeff[:,j]
   original_mi=len(xr)/W.sum()*xr@W@xr/np.dot(xr,xr)
   surrogate_mi=len(xr)/W.sum()*surr@W@surr/np.dot(surr,surr)
   maxmi=max(maxmi,abs(original_mi-surrogate_mi))
  assert np.max(abs(D.T@xr))<1e-6 and np.max(abs(D.T@yr))<1e-6
  checked+=1
 assert maxrho<1e-10 and maxnull<1e-8 and maxmi<1e-8
 # Read-only identity and direct reaggregation of every exported T readout.
 for rec in json.loads((P/'qc/INPUT_MANIFEST.json').read_text()):assert sha(rec['path'])==rec['sha256']
 terr=0.;regions=0
 for f in (P/'inputs').glob('*_REGIONS.tsv.gz'):
  sample,radius=f.name.replace('_REGIONS.tsv.gz','').rsplit('_',1);radius=int(radius)
  z=read(f);d=read(O/'inputs'/f'{sample}_INPUTS.tsv.gz')
  mm=json.loads((A/'inputs'/f'{sample}_morphology_supported_MEMBERSHIPS.json').read_text())
  mm={x['edge']:x for x in mm if x['radius']==radius}
  v=d.Tcell_q05.to_numpy()
  for r in z.itertuples():
   vals=[]
   for side in ['inside','outside']:
    ids=mm[r.edge][side];vv=v[ids];vv=vv[np.isfinite(vv)];value=np.mean(vv) if len(vv)>=3 else np.nan
    got=getattr(r,'T_'+side);assert np.isfinite(value)==np.isfinite(got)
    if np.isfinite(value):terr=max(terr,abs(value-got))
    vals.append(value)
   delta=np.log1p(vals[0])-np.log1p(vals[1]);assert np.isclose(delta,r.T_delta,equal_nan=True)
   regions+=1
 assert terr<1e-10
 # Match identities, no reuse, complete spatial disjointness, calipers and balance.
 pairs=read(P/'tables/MATCHED_PAIRS.tsv');qual=read(P/'tables/MATCH_QUALIFICATION.tsv');ma=read(P/'tables/MATCH_SECTION.tsv')
 for (sample,radius,model),g in pairs.groupby(['sample','radius','model']):
  z=read(P/'inputs'/f'{sample}_{radius}_REGIONS.tsv.gz').set_index('edge')
  ids=g.high_edge.tolist()+g.low_edge.tolist();assert len(ids)==len(set(ids))
  mm=json.loads((A/'inputs'/f'{sample}_morphology_supported_MEMBERSHIPS.json').read_text());mm={m['edge']:m for m in mm if m['radius']==radius}
  used=set()
  for ed in ids:
   spots=set(mm[ed]['inside']+mm[ed]['outside']);assert not used&spots;used|=spots
  cov=CORE+(VASC if model=='vascular' else []);b=z[np.isfinite(z[['in_g_TYMP']+cov]).all(1)]
  sd=b[cov].std(ddof=0);ac=sd>1e-12
  delta=(z.loc[g.high_edge,cov].to_numpy()-z.loc[g.low_edge,cov].to_numpy())[:,ac]/sd[ac].to_numpy()
  assert abs(delta).max()<=1+1e-10 and np.sqrt((delta**2).mean(1)).max()<=.75+1e-10
  smd=abs(delta.mean(0)).max()
  qr=qual[(qual['sample']==sample)&(qual.radius==radius)&(qual.model==model)].iloc[0]
  assert abs(qr.max_abs_SMD-smd)<1e-10 and bool(qr.eligible)==(len(g)>=3 and smd<=.25)
  assert np.all(z.loc[g.high_edge,'component'].to_numpy()==z.loc[g.low_edge,'component'].to_numpy())
  for fam,ep,col in MATCH:
   dif=z.loc[g.high_edge,col].to_numpy()-z.loc[g.low_edge,col].to_numpy()
   assert np.allclose(dif,g[ep+'_difference'],equal_nan=True)
   rr=ma[(ma['sample']==sample)&(ma.radius==radius)&(ma.model==model)&(ma.endpoint==ep)].iloc[0]
   finite=dif[np.isfinite(dif)]
   expected=float(finite.mean()) if qr.eligible and len(finite)>=3 else np.nan
   assert np.isclose(rr.effect,expected,equal_nan=True)
   # Endpoint-complete subset balance must match qualification; whole-slice CD8 missing is allowed.
   if 0<len(finite)<len(dif):
    dd=delta[np.isfinite(dif)];assert abs(dd.mean(0)).max()<=.25 or not np.isfinite(rr.effect)
 # Aggregate independence and separate FDR families.
 for name in ['SPECTRAL','MATCH']:
  s=read(P/'tables'/('SPECTRAL_SECTION.tsv' if name=='SPECTRAL' else 'MATCH_SECTION.tsv'));pa=read(P/'tables'/f'{name}_PATIENT.tsv');su=read(P/'tables'/f'{name}_COHORT.tsv')
  keys=['method','dataset','radius','model','family','endpoint']
  for r in pa.itertuples():
   mask=np.ones(len(s),bool)
   for k in keys+['patient']:mask &= s[k].eq(getattr(r,k))
   v=s.loc[mask,'effect'].dropna().to_numpy()
   val=(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean()) if name=='SPECTRAL' else v.mean()) if len(v) else np.nan
   assert np.isclose(val,r.effect,equal_nan=True)
  for r in su.itertuples():
   mask=np.ones(len(pa),bool)
   for k in keys:mask &= pa[k].eq(getattr(r,k))
   v=pa.loc[mask,'effect'].dropna().to_numpy();assert len(v)==r.patients and (v>0).sum()==r.positive and (v<0).sum()==r.negative
   pv=stats.wilcoxon(v,method='auto').pvalue if len(v)>=10 and np.any(v!=0) else (1. if len(v)>=10 else np.nan)
   assert np.isclose(pv,r.p_patient,equal_nan=True)
  for _,g in su.groupby(keys[:-1]):
   g=g[g.p_patient.notna()]
   if len(g):assert np.allclose(stats.false_discovery_control(g.p_patient.to_numpy()),g.FDR_patient)
 for _,g in se.groupby(['sample','radius','model','family']):
  g=g[g.spatial_p.notna()]
  if len(g):assert np.allclose(stats.false_discovery_control(g.spatial_p.to_numpy()),g.spatial_FDR)
 dump(P/'qc/INDEPENDENT_VALIDATION.json',dict(status='PASS',spectral_models_checked=checked,null_replicates_each=999,max_rho_error=maxrho,max_null_error=maxnull,max_moran_error=maxmi,regions_reaggregated=regions,max_T_reaggregation_error=terr,matched_pairs_checked=len(pairs),input_hashes_unchanged=True,patient_units_and_FDR_checked=True))
 print('VALIDATION_PASS',checked,regions,len(pairs),flush=True)
if __name__=='__main__':main()
