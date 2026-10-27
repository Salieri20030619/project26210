from source import *
def independent_rho(x,y,cov=None,mask=None,min_n=50):
 x=np.asarray(x,float);y=np.asarray(y,float)
 if cov is None:cov=np.empty((len(x),0))
 mask=np.ones(len(x),bool) if mask is None else np.asarray(mask,bool)
 good=mask&np.isfinite(x)&np.isfinite(y)&np.isfinite(cov).all(1);A=np.column_stack([x,y,cov])[good];n=len(A)
 if n<min_n or min(len(np.unique(A[:,0])),len(np.unique(A[:,1])))<3:return np.nan
 ranks=pd.DataFrame(A).rank(method='average').to_numpy();c=ranks[:,2:];sd=c.std(0)
 if (sd==0).any():return np.nan
 D=np.column_stack([np.ones(n),(c-c.mean(0))/sd]) if c.shape[1] else np.ones((n,1))
 if np.linalg.matrix_rank(D)<D.shape[1] or np.linalg.cond(D)>=30:return np.nan
 u,s,vt=np.linalg.svd(D,full_matrices=False);res=ranks[:,:2]-u@(u.T@ranks[:,:2])
 if (res.std(0)<1e-10).any():return np.nan
 return float(stats.pearsonr(res[:,0],res[:,1]).statistic)
def main():
 for name in ['source','spatial']:assert (P/'logs'/f'{name}_exit.txt').read_text().strip()=='0'
 checks=[];maxerr=0.;count=0
 cache=read(P/'inputs/CELL_CACHE.tsv.gz');ef=read(P/'tables/WITHIN_MACRO_PATIENT_EFFECTS.tsv')
 for r in ef[ef.effect.notna()].itertuples():
  z=cache[cache.patient==r.patient];mask=z.candidate.astype(str).eq('4').to_numpy() if r.group=='C4' else np.ones(len(z),bool)
  cov=np.column_stack([z.log_depth,z.n_genes])
  if r.group=='WholeMacro':cov=np.column_stack([cov,pd.get_dummies(z.candidate.astype(str),dtype=float,drop_first=True).to_numpy()])
  rho=independent_rho(z.TYMP,z[r.endpoint],None if r.model=='unadjusted' else cov,mask,50 if r.group=='C4' else 100)
  assert np.isfinite(rho) and abs(rho-r.effect)<1e-10,(r,rho)
  maxerr=max(maxerr,abs(rho-r.effect));count+=1
 checks.append(dict(check='within_patient_all_effects_SVD_recomputed',n=count,max_error=maxerr))
 from spatial import CORE,OUT,HP,H
 sfe=read(P/'tables/SPATIAL_SECTION_EFFECTS.tsv');maxerr=0.;count=0;aggerr=0.;aggcount=0
 reg=read(H/'tables/INPUT_SECTION_REGISTER.tsv')
 for r in reg.itertuples():
  z=read(P/'inputs'/f'{r.sample}_SPLIT_INTERFACE.tsv.gz');d=read(OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz');sc=read(HP/'inputs'/f'{r.sample}_PROGRAMS.tsv.gz');mem=read(H/'inputs'/f'{r.sample}_SOURCE_MEMBERSHIP.tsv.gz');inside=d.malignant_fraction.ge(.2).to_numpy()
  for sid,g in mem.groupby('segment'):
   ids=g.spot_index.to_numpy(int);a=ids[inside[ids]];b=ids[~inside[ids]]
   assert not set(a)&set(b)
   assert z.at[sid,'inside_n']==len(a) and z.at[sid,'outside_n']==len(b)
   if len(a) and len(b):
    for col,val in [('inside_TYMP',d.TYMP_RNA.iloc[a].mean()),('outside_TYMP',d.TYMP_RNA.iloc[b].mean()),('outside_C4',sc.C4_platform17.iloc[b].mean())]:
     assert abs(z.at[sid,col]-val)<1e-10;aggerr=max(aggerr,abs(z.at[sid,col]-val));aggcount+=1
  for rr in sfe[(sfe['sample']==r.sample)&(sfe.family=='spatial_cross_compartment')&sfe.effect.notna()].itertuples():
   cov={'unadjusted':[],'core_adjusted':CORE,'outside_TYMP_adjusted':CORE+['outside_TYMP'],'vascular_sensitivity':CORE+['outside_TYMP','inside_endothelial','outside_endothelial']}[rr.model]
   rho=independent_rho(z.inside_TYMP,z[rr.endpoint],z[cov].to_numpy() if cov else None,z.split_eligible.to_numpy(bool))
   assert np.isfinite(rho) and abs(rho-rr.effect)<1e-10,(rr,rho)
   maxerr=max(maxerr,abs(rho-rr.effect));count+=1
 checks.append(dict(check='spatial_all_correlations_SVD',n=count,max_error=maxerr));checks.append(dict(check='disjoint_source_aggregates',n=aggcount,max_error=aggerr))
 # Verify all sign tests, BH and participant aggregation.
 for name,keys,eff_file,effkeys in [('WITHIN_MACRO_SUMMARY.tsv',['group'],'WITHIN_MACRO_PATIENT_EFFECTS.tsv',['group','endpoint','model']),('SPATIAL_COHORT_SUMMARY.tsv',['dataset','family'],'SPATIAL_PATIENT_EFFECTS.tsv',['dataset','family','endpoint','model','effect_type'])]:
  su=read(P/'tables'/name);e=read(P/'tables'/eff_file)
  for r in su.to_dict('records'):
   sel=np.ones(len(e),bool)
   for k in effkeys:sel &= e[k].eq(r[k]).to_numpy()
   v=e.loc[sel,'effect'].dropna().to_numpy();assert len(v)==r['n'] and int((v>0).sum())==r['positive']
   assert np.isclose(np.median(v),r['median']) if len(v) else np.isnan(r['median'])
   if len(v)>=10:
    nz=np.count_nonzero(v);pp=stats.binomtest(np.count_nonzero(v>0),nz,.5).pvalue if nz else 1.
    assert abs(pp-r['P'])<1e-12
   else:assert pd.isna(r['P'])
  for _,g in su.groupby(keys):
   ref=stats.false_discovery_control(g.P.fillna(1).to_numpy());ix=g.P.notna().to_numpy();assert np.allclose(ref[ix],g.FDR.to_numpy()[ix],atol=1e-12)
  checks.append(dict(check='patient_sign_and_BH_'+name,n=len(su),pass_all=True))
 for name,keys in [('PAIRED_SOURCE_ASSOCIATIONS.tsv',['source']),('TYMP_CD8_FUNCTION_SUMMARY.tsv',['source'])]:
  su=read(P/'tables'/name)
  for _,g in su.groupby(keys):
   ref=stats.false_discovery_control(g.P.fillna(1).to_numpy());ix=g.P.notna().to_numpy();assert np.allclose(ref[ix],g.FDR.to_numpy()[ix],atol=1e-12)
  assert su.loc[su.n<10,'P'].isna().all()
 # Independently calculate all cross-patient target correlations with study dummies.
 b=read(P/'tables/PAIRED_SOURCE_PATIENT_INPUTS.tsv');su=read(P/'tables/PAIRED_SOURCE_ASSOCIATIONS.tsv');maxerr=0
 for r in su.itertuples():
  other='macro_TYMP' if r.source=='tumor_TYMP' else 'tumor_TYMP';g=b[b[r.endpoint].notna()].copy()
  A=g[[r.source,r.endpoint]].rank().to_numpy();cols=[np.ones(len(g))]+[g.study.eq(s).to_numpy(float) for s in sorted(g.study.unique())[1:]]
  if r.model=='study_other_source_adjusted':cols.append(g[other].rank().to_numpy())
  D=np.column_stack(cols);u=np.linalg.svd(D,full_matrices=False)[0];res=A-u@(u.T@A);rho=stats.pearsonr(res[:,0],res[:,1]).statistic
  assert abs(rho-r.rho)<1e-10;maxerr=max(maxerr,abs(rho-r.rho))
 checks.append(dict(check='all_cross_source_rhos_independent_SVD',n=len(su),max_error=maxerr))
 for name in ['SOURCE_INPUTS.json','SPATIAL_INPUTS.json']:
  for rec in json.loads((P/'qc'/name).read_text()):assert sha(rec['path'])==rec['sha256']
 assert not list(P.rglob('*.png')) and not list(P.rglob('*.pdf'))
 dump(P/'qc/INTERNAL_VALIDATION.json',dict(status='PASS',at=now(),checks=checks,original_inputs_unchanged=True,source_exit=0,spatial_exit=0))
 print('INTERNAL_VALIDATION_PASS',checks,flush=True)
if __name__=='__main__':main()
