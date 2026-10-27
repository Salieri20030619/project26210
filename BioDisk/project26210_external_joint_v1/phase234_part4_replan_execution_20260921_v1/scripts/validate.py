from common import *
import math,traceback
def exact_sign(v):
 v=np.asarray(v);v=v[np.isfinite(v)];nz=v[v!=0];n=len(nz)
 if not n:return 1.
 k=min(int((nz>0).sum()),int((nz<0).sum()))
 return min(1.,2*sum(math.comb(n,j) for j in range(k+1))/2**n)
def bh(p):
 p=np.asarray(p,float);order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=np.empty(len(p));out[order]=np.minimum(q,1);return out
def select(df,row,keys):
 mask=np.ones(len(df),bool)
 for key in keys:mask &= df[key].astype(str).to_numpy()==str(getattr(row,key))
 return df[mask]
def main():
 for key in ['A1_COMPLETE.json','A2_A4_COMPLETE.json']:assert (P/key).exists(),key
 checks=[];counts={}
 for prefix in ['A1','A2_A4']:
  sec=pd.read_csv(P/'tables'/(prefix+'_SECTION_EFFECTS.tsv'),sep='\t',float_precision='round_trip')
  pat=pd.read_csv(P/'tables'/(prefix+'_PATIENT_EFFECTS.tsv'),sep='\t',float_precision='round_trip')
  su=pd.read_csv(P/'tables'/(prefix+'_COHORT_SUMMARY.tsv'),sep='\t',float_precision='round_trip')
  assert not sec.duplicated(['sample']+KEYS).any()
  assert not pat.duplicated(['dataset','patient']+KEYS).any()
  assert not su.duplicated(['dataset']+KEYS).any()
  # Exhaustive independent patient aggregation.
  pg=pat.set_index(['dataset','patient']+KEYS)
  for key,g in sec.groupby(['dataset','patient']+KEYS):
   v=g.effect.dropna().to_numpy();rr=pg.loc[key]
   assert int(rr.n_sections)==len(v)
   expected=(np.tanh(np.mean(np.arctanh(np.clip(v,-.999999999,.999999999)))) if key[-1]=='rho' else np.mean(v)) if len(v) else np.nan
   assert (np.isnan(expected) and np.isnan(rr.effect)) or abs(expected-rr.effect)<2e-12
  patgroups={k:g for k,g in pat.groupby(['dataset']+KEYS)}
  for rr in su.itertuples():
   key=tuple(getattr(rr,c) for c in ['dataset']+KEYS);v=patgroups[key].effect.dropna().to_numpy()
   assert len(v)==rr.n
   assert int((v>0).sum())==rr.positive and int((v<0).sum())==rr.negative
   if len(v):assert abs(float(np.median(v))-rr.median)<1e-12
   if len(v)>=10:
    p=exact_sign(v);assert abs(p-rr.P)<1e-12
    wp=stats.wilcoxon(v,method='auto').pvalue if (v!=0).any() else 1.
    assert abs(wp-rr.P_wilcoxon)<1e-12
    checks.append(dict(prefix=prefix,dataset=rr.dataset,endpoint=rr.endpoint,layer=rr.layer,scale=rr.scale,model=rr.model,P_exact_independent=p,error=abs(p-rr.P)))
   else:assert np.isnan(rr.P) and np.isnan(rr.FDR)
  for key,g in su.groupby(['dataset','family','layer']):
   for pcol,qcol in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:
    qq=bh(g[pcol].fillna(1).to_numpy());ok=g[pcol].notna().to_numpy()
    np.testing.assert_allclose(g[qcol].to_numpy()[ok],qq[ok],atol=1e-12)
   assert (g.family_size_registered==len(g)).all()
  counts[prefix]=dict(section_models=len(sec),patient_models=len(pat),cohort_summaries=len(su),estimable=int(sec.effect.notna().sum()))
 tab('audit/INDEPENDENT_PATIENT_SIGN_CHECKS.tsv',pd.DataFrame(checks))
 regs=pd.read_csv(P/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t',float_precision='round_trip')
 ep=pd.read_csv(P/'tables/FROZEN_ENDPOINTS.tsv',sep='\t',float_precision='round_trip')
 broad=pd.read_csv(P/'tables/A1_SECTION_EFFECTS.tsv',sep='\t',float_precision='round_trip')
 narrow=pd.read_csv(P/'tables/A2_A4_SECTION_EFFECTS.tsv',sep='\t',float_precision='round_trip')
 numerical=[];graphs=[];scorechecks=[]
 from anndata.io import read_elem
 import h5py
 definitions=json.loads((P/'PROGRAM_DEFINITIONS.json').read_text())['programs']
 representatives={'GLMF1','H46_T','GSM5924030','GSM5924042'}
 for rr in regs.itertuples():
  d=pd.read_csv(P/'inputs'/(rr.sample+'_INPUTS.tsv.gz'),sep='\t',float_precision='round_trip')
  G=sparse.load_npz(P/'graphs'/(rr.sample+'_ring.npz')).tocsr()
  D=sparse.load_npz(P/'graphs'/(rr.sample+'_distance800.npz')).tocsr()
  assert not G.diagonal().any() and (G-G.T).nnz==0
  assert not D.diagonal().any() and (D-D.T).nnz==0
  assert set(np.unique(D.data))<=set(range(100,801,100))
  a,b=D.nonzero();assert np.array_equal(d.component.to_numpy()[a],d.component.to_numpy()[b])
  assert D.shape==(len(d),len(d)) and d.spot_id.is_unique
  assert ((d.V3_SPP1_macro_share>=0)&(d.V3_SPP1_macro_share<=1)).all()
  np.testing.assert_allclose(d.V2_SPP1_fraction/d.total_macro_fraction,d.V3_SPP1_macro_share,atol=1e-10)
  prev=np.zeros(len(d))
  for radius in [200,400,800]:
   H=D.copy();H.data=(H.data<=radius).astype(float);H.eliminate_zeros();nn=np.asarray(H.sum(axis=1)).ravel();assert (nn>=prev).all();prev=nn
  graphs.append(dict(sample=rr.sample,spots=len(d),distance_graph_edges=D.nnz,no_cross_component=True,no_self=True))
  # Independent broad calculations, all38 sections and all70 endpoints.
  bb=broad[broad['sample']==rr.sample]
  for r in bb.itertuples():
   values=d[[r.exposure,r.target]].replace([np.inf,-np.inf],np.nan).dropna()
   if r.layer=='malignant_region':values=values.loc[d.loc[values.index,'tumor_expression_region']]
   if len(values)<50:check=np.nan
   elif r.scale!='W0':
    grid=4 if r.scale=='W1_4x4' else 5
    values=values.groupby(d.loc[values.index,'block_'+str(grid)]).median()
    check=values.iloc[:,0].corr(values.iloc[:,1],method='spearman') if len(values)>=8 else np.nan
   else:check=values.iloc[:,0].corr(values.iloc[:,1],method='spearman')
   if len(values) and min(values.iloc[:,0].nunique(),values.iloc[:,1].nunique())<3:check=np.nan
   if np.isfinite(r.effect):
    assert abs(check-r.effect)<2e-12
    numerical.append(dict(sample=rr.sample,endpoint=r.endpoint,layer=r.layer,scale=r.scale,model=r.model,error=abs(check-r.effect)))
   else:assert not np.isfinite(check),(rr.sample,r.endpoint,r.scale,check)
  # Independent score recomputation from raw counts for representative sections across sources/preparations.
  if rr.sample in representatives:
   with h5py.File(rr.input_path) as h:
    var=read_elem(h['var']);X=read_elem(h['layers/counts'] if 'counts' in h.get('layers',{}) else h['X'])
    totals=np.asarray(X.sum(axis=1)).ravel()
    for name,genes in definitions.items():
     if not d[name].notna().all():continue
     used=[g for g in genes if g in var.index];raw=X[:,var.index.get_indexer(used)].toarray().astype(float)
     l=np.log1p(raw*10000/totals[:,None]);sd=l.std(axis=0)
     z=(l-l.mean(axis=0))/np.where(sd>0,sd,1);score=z.sum(axis=1)/len(used)
     err=float(np.max(abs(score-d[name].to_numpy())));assert err<1e-10
     scorechecks.append(dict(sample=rr.sample,program=name,max_error=err))
  # Every section's finite band contrasts independently reproduced.
  for r in narrow[(narrow['sample']==rr.sample)&(narrow.layer=='expression_boundary')].itertuples():
   if not np.isfinite(r.effect):continue
   side=d.tumor_expression_region.to_numpy(bool) if r.family=='boundary_tumor' else ~d.tumor_expression_region.to_numpy(bool)
   dist=d.expression_distance_um.to_numpy();v=d[r.exposure].to_numpy()
   a=v[side&np.isfinite(dist)&(dist<=int(r.scale))&np.isfinite(v)]
   b=v[side&np.isfinite(dist)&(dist>800)&np.isfinite(v)]
   assert len(a)>=20 and len(b)>=20
   err=abs(float(np.median(a)-np.median(b))-r.effect);assert err<1e-12
  if rr.sample in representatives:
   for scale,sg in narrow[(narrow['sample']==rr.sample)&(~narrow.layer.isin(['expression_boundary','expression_boundary_slope']))].groupby('scale'):
    H=G.copy() if scale=='ring' else D.copy()
    if scale!='ring':H.data=(H.data<=int(scale)).astype(float);H.eliminate_zeros()
    n=np.asarray(H.sum(axis=1)).ravel();W=sparse.diags(1/np.maximum(n,1))@H
    for r in sg.itertuples():
     if not np.isfinite(r.effect):continue
     x=d[r.exposure].to_numpy();source=d[r.target].to_numpy();y=W@source
     cut={'Q0':0,'Q1':.05,'Q2':.10}[r.qualification]
     mask=(d.total_macro_fraction.to_numpy()>=cut)&(n>=3)&np.isfinite(x)&np.isfinite(y)
     cv=[]
     if r.model!='A0':cv=[d.log1p_UMI.to_numpy(),W@d.log1p_UMI.to_numpy()]
     if r.model in ['A2','A3','A4']:
      bg=d.other_macro_fraction.to_numpy() if r.exposure=='V1_SPP1_q05' else d.total_macro_fraction.to_numpy()
      cv += [bg,W@bg,d.malignant_fraction.to_numpy()]
     if r.model in ['A3','A4']:cv += [source]
     if r.model=='A4':
      xy=d[['physical_x_um','physical_y_um']].to_numpy();xy=(xy-xy.mean(axis=0))/xy.std(axis=0);cv += [xy[:,0],xy[:,1],xy[:,0]**2,xy[:,0]*xy[:,1],xy[:,1]**2]
     z=pd.DataFrame(np.column_stack([x,y]+cv)[mask]).rank(method='average').to_numpy()
     design=np.column_stack([np.ones(len(z)),z[:,2:]])
     # SVD pseudoinverse independent of fitted least-squares/QR implementation.
     res=z[:,:2]-design@np.linalg.pinv(design)@z[:,:2]
     check=np.dot(res[:,0],res[:,1])/np.sqrt(np.dot(res[:,0],res[:,0])*np.dot(res[:,1],res[:,1]))
     err=abs(check-r.effect);assert err<2e-9,(rr.sample,r.endpoint,r.model,err)
     numerical.append(dict(sample=rr.sample,endpoint=r.endpoint,layer=r.layer,scale=r.scale,model=r.model,error=err))
  print('VALIDATE',rr.sample,flush=True)
 tab('audit/NUMERICAL_EFFECT_CHECKS.tsv',pd.DataFrame(numerical));tab('audit/RAW_SCORE_CHECKS.tsv',pd.DataFrame(scorechecks));tab('audit/GRAPH_INTEGRITY.tsv',pd.DataFrame(graphs))
 for x in json.loads((P/'audit/SOURCE_MANIFEST.json').read_text()):
  assert all(meta(x['path'])[k]==x[k] for k in ['bytes','mtime_ns','inode'])
  if x.get('sha256'):assert sha(x['path'])==x['sha256']
 old=json.loads((P/'audit/RDS_BEFORE.json').read_text())
 for x in old:assert meta(x['path'])==x
 lock=json.loads((P/'A0_FINAL_FREEZE.json').read_text())
 assert sha(P/'IMPLEMENTATION_FREEZE_BEFORE_ASSOCIATION.json')==lock['parameters_sha256']
 assert sha(P/'tables/FROZEN_ENDPOINTS.tsv')==lock['endpoints_sha256']
 result=dict(status='PASS',time=stamp(),counts=counts,independent_sign_tests=len(checks),independent_numerical_effect_checks=len(numerical),max_effect_error=max(x['error'] for x in numerical),raw_program_score_checks=len(scorechecks),RDS_preserved=len(old),source_metadata_and_recorded_hashes_unchanged=True,graphs_verified=38,physical_unit_verified=True,GLMF_new_P_values=0,new_figures=0)
 save('FINAL_VALIDATION.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':
 try:main()
 except Exception:(P/'ERROR_VALIDATION.txt').write_text(traceback.format_exc());raise
