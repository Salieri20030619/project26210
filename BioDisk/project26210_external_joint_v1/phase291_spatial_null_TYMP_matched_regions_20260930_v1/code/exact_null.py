from analyze import *
def main():
 sec=read(P/'tables/SPECTRAL_SECTION.tsv');rows=[];qc=[]
 for r in sec.itertuples():
  out=r._asdict();out.pop('Index',None);out['null_method']='covariate_constrained_spectral'
  out['spatial_p']=np.nan;out['spatial_FDR']=np.nan
  if r.reason!='OK':rows.append(out);continue
  z=np.load(P/'inputs'/f'{r.sample}_{r.radius}_{r.model}_{r.endpoint}_SPECTRAL.npz')
  D=z['D'];xr=z['xr'];yr=z['yr'];W=z['W'];n=len(xr)
  Q=linalg.null_space(D.T);ev,U=linalg.eigh(Q.T@W@Q);V=Q@U;a=V.T@xr
  assert np.max(abs(V@a-xr))<1e-7
  seed=int(hashlib.sha256((r.sample+str(r.radius)+r.model+r.endpoint+'exact_cov').encode()).hexdigest()[:8],16)
  signs=np.random.default_rng(seed).choice([-1.,1.],size=(len(a),999));sur=V@(a[:,None]*signs)
  null=yr@sur/(np.linalg.norm(xr)*np.linalg.norm(yr));obs=float(xr@yr/(np.linalg.norm(xr)*np.linalg.norm(yr)))
  mi=n/W.sum()*np.sum(sur*(W@sur),axis=0)/np.sum(sur*sur,axis=0)
  mor=n/W.sum()*(xr@W@xr)/(xr@xr)
  err=float(abs(mi-mor).max());cov=float(abs(D.T@sur).max());norm=float(abs(np.sum(sur*sur,axis=0)/(xr@xr)-1).max())
  assert err<1e-8 and cov<1e-6 and norm<1e-8 and abs(obs-r.effect)<1e-10
  out.update(spatial_p=float((1+(abs(null)>=abs(obs)-1e-12).sum())/1000),null_q025=float(np.quantile(null,.025)),null_q975=float(np.quantile(null,.975)))
  rows.append(out);qc.append(dict(sample=r.sample,radius=r.radius,model=r.model,endpoint=r.endpoint,regions=n,null_dimensions=len(a),max_moran_error=err,max_cov_error=cov,max_norm_error=norm))
  np.savez_compressed(P/'inputs'/f'{r.sample}_{r.radius}_{r.model}_{r.endpoint}_EXACT_NULL.npz',V=V,a=a,signs=signs.astype(np.int8),null=null)
 e=pd.DataFrame(rows)
 for _,ix in e.groupby(['sample','radius','model','family']).groups.items():
  jj=[i for i in ix if pd.notna(e.at[i,'spatial_p'])]
  if jj:e.loc[jj,'spatial_FDR']=bh(e.loc[jj,'spatial_p'])
 tab(e,'EXACT_NULL_SECTION.tsv');tab(qc,'EXACT_NULL_QC.tsv')
 # Independent check: reconstruct every surrogate correlation, check covariance and Moran on 3 reps.
 maxrho=0.;maxmi=0.
 for r in e[e.reason.eq('OK')].itertuples():
  orig=np.load(P/'inputs'/f'{r.sample}_{r.radius}_{r.model}_{r.endpoint}_SPECTRAL.npz')
  z=np.load(P/'inputs'/f'{r.sample}_{r.radius}_{r.model}_{r.endpoint}_EXACT_NULL.npz')
  xr=orig['xr'];yr=orig['yr'];D=orig['D'];W=orig['W']
  nco=(z['V'].T@yr)@ (z['a'][:,None]*z['signs'])/np.sqrt((xr@xr)*(yr@yr))
  maxrho=max(maxrho,float(abs(nco-z['null']).max()))
  assert abs((1+(abs(nco)>=abs(r.effect)-1e-12).sum())/1000-r.spatial_p)<1e-12
  for j in [0,499,998]:
   v=z['V']@(z['a']*z['signs'][:,j]);assert abs(D.T@v).max()<1e-6
   er=abs((v@W@v)/(v@v)-(xr@W@xr)/(xr@xr))*len(xr)/W.sum();maxmi=max(maxmi,float(er))
 for _,g in e.groupby(['sample','radius','model','family']):
  g=g[g.spatial_p.notna()]
  if len(g):assert np.allclose(stats.false_discovery_control(g.spatial_p.to_numpy()),g.spatial_FDR)
 assert maxrho<1e-10 and maxmi<1e-8
 dump(P/'qc/EXACT_NULL_VALIDATION.json',dict(status='PASS',models=len(qc),null_replicates_each=999,max_null_rho_error=maxrho,max_moran_error=maxmi,observed_effects_and_patient_statistics_unchanged=True))
 print('EXACT_NULL_COMPLETE_VALIDATED',len(qc),flush=True)
if __name__=='__main__':main()
