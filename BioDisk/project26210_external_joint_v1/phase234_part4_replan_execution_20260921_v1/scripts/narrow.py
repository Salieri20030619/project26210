from common import *
import traceback
from scipy.sparse.csgraph import dijkstra
def matrix(d,sid,scale):
 if scale=='ring':
  H=sparse.load_npz(P/'graphs'/(sid+'_ring.npz')).tocsr()
 else:
  f=P/'graphs'/(sid+'_distance800.npz')
  if not f.exists():return None,None
  H=sparse.load_npz(f).tocsr();H.data=(H.data<=int(scale)).astype(float);H.eliminate_zeros()
 nn=np.asarray(H.sum(axis=1)).ravel()
 return sparse.diags(1/np.maximum(nn,1))@H,nn
def main():
 assert (P/'A0_FINAL_FREEZE.json').exists()
 regs=pd.read_csv(P/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t')
 ep=pd.read_csv(P/'tables/FROZEN_ENDPOINTS.tsv',sep='\t')
 core=ep[((ep.exposure.isin(['V1_SPP1_q05','V3_SPP1_macro_share']))&ep.family.isin(['SPP1_vascular','SPP1_stroma_program','SPP1_function']))|((ep.family=='MTA_context_SPP1')&ep.target.isin(['V1_SPP1_q05','V3_SPP1_macro_share']))|ep.family.isin(['KEGG_four_axes','custom_four_axes'])]
 lr=pd.read_csv(P/'tables/FROZEN_LR_CANDIDATES.tsv',sep='\t')
 rows=[];regions=[];lrcover=[]
 for r in regs.itertuples():
  d=pd.read_csv(P/'inputs'/(r.sample+'_INPUTS.tsv.gz'),sep='\t');N=len(d);macro=d.total_macro_fraction.to_numpy()
  def add(e,layer,scale,q,model,x,y,mask,cov=None,effect_type='rho'):
   value,n,reason=fit(x,y,cov=cov,mask=mask)
   rows.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family=e.family,endpoint=e.endpoint,exposure=e.exposure,target=e.target,layer=layer,scale=str(scale),qualification=q,model=model,effect_type=effect_type,effect=value,n_obs=n,reason=reason))
  for scale in ['800','400','200','ring']:
   W,nn=matrix(d,r.sample,scale)
   if W is None:
    for e in core.itertuples():
     rows.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family=e.family,endpoint=e.endpoint,exposure=e.exposure,target=e.target,layer='scale',scale=scale,qualification='Q1',model='A0',effect_type='rho',effect=np.nan,n_obs=0,reason='PHYSICAL_SCALE_NOT_ELIGIBLE'))
    continue
   neighbor={t:W@d[t].fillna(0).to_numpy() if d[t].notna().all() else np.full(N,np.nan) for t in set(core.target)}
   for e in core.itertuples():
    add(e,'scale',scale,'Q1','A0',d[e.exposure].to_numpy(),neighbor[e.target],(macro>=.05)&(nn>=3))
   # Dedicated LR list: independent family; full receptor complex must have passed expression gate.
   for rr in lr.itertuples():
    x=d[rr.ligand+'_RNA'].to_numpy();rec=d['LR_'+rr.receptor].to_numpy()
    y=W@np.nan_to_num(rec) if np.isfinite(rec).all() else np.full(N,np.nan)
    eligible=(nn>=3)
    class E:pass
    e=E();e.family=rr.family;e.endpoint=rr.interaction;e.exposure=rr.ligand+'_RNA';e.target='LR_'+rr.receptor
    add(e,'LR_compatibility',scale,'Q0','A0',x,y,eligible)
    lrcover.append(dict(sample=r.sample,interaction=rr.interaction,scale=scale,eligible_sources=int(eligible.sum()),ligand_positive_sources=int((eligible&(x>0)).sum()),sources_with_receptor_neighbor=int((eligible&(y>0)).sum()),ligand_positive_with_receptor_neighbor=int((eligible&(x>0)&(y>0)).sum()),complex_eligible=bool(np.isfinite(rec).all()),not_celltype_specific=True))
   if scale!='400':continue
   targets=core[core.exposure.isin(['V1_SPP1_q05','V3_SPP1_macro_share'])]
   for e in targets.itertuples():
    x=d[e.exposure].to_numpy();sc=d[e.target].to_numpy();y=neighbor[e.target];mask=(macro>=.05)&(nn>=3)
    depth=np.column_stack([d.log1p_UMI,W@d.log1p_UMI.to_numpy()])
    bg=d.other_macro_fraction.to_numpy() if e.exposure=='V1_SPP1_q05' else macro
    composition=np.column_stack([bg,W@bg,d.malignant_fraction])
    blocks={'A0':None,'A1':depth,'A2':np.column_stack([depth,composition]),'A3':np.column_stack([depth,composition,sc])}
    if e.target in ['POSTNCAF','ECM4']:
     xy=d[['physical_x_um','physical_y_um']].to_numpy();xy=(xy-xy.mean(axis=0))/xy.std(axis=0)
     blocks['A4']=np.column_stack([depth,composition,sc,xy,xy[:,0]**2,xy[:,0]*xy[:,1],xy[:,1]**2])
    for mod,cov in blocks.items():add(e,'adjustment',scale,'Q1',mod,x,y,mask,cov)
    if e.exposure=='V3_SPP1_macro_share':
     for q,cut in [('Q0',0),('Q1',.05),('Q2',.10)]:
      add(e,'qualification',scale,q,'A0',x,y,(macro>=cut)&(nn>=3))
  # Expression-boundary branch only; no AI histology upgrade. All targets use independent frozen malignant partition.
  distance=d.expression_distance_um.to_numpy();tumor=d.tumor_expression_region.to_numpy(bool)
  featurelist=['V1_SPP1_q05','V2_SPP1_fraction','V3_SPP1_macro_share','V4_SPP1_program','total_macro_fraction','POSTNCAF','ECM4']
  for side,side_mask in [('tumor',tumor),('exterior',~tumor)]:
   finite=side_mask&np.isfinite(distance)
   coremask=finite&(distance>800)
   for feature in featurelist:
    value=d[feature].to_numpy()
    for width in [800,400,200]:
     band=finite&(distance<=width);a=value[band&np.isfinite(value)];b=value[coremask&np.isfinite(value)]
     effect=float(np.median(a)-np.median(b)) if len(a)>=20 and len(b)>=20 else np.nan
     reason='' if np.isfinite(effect) else 'BAND_OR_FIXED_GT800_REFERENCE_LT20'
     rows.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family='boundary_'+side,endpoint=feature+'__band_vs_core',exposure=feature,target=side+'_band_minus_gt800',layer='expression_boundary',scale=str(width),qualification='Q0_region_ge20',model='A0',effect_type='median_difference',effect=effect,n_obs=len(a)+len(b),n_band=len(a),n_reference=len(b),reason=reason))
    # Continuous distance bins are a separate descriptive slope; direction not pooled with band differences.
    dd=pd.DataFrame(dict(distance=distance[finite],value=value[finite]))
    dd['bin']=np.floor(dd.distance/200).astype(int)
    bb=dd.groupby('bin').agg(value=('value','median'),n=('value','count'),distance=('distance','median')).reset_index()
    bb=bb[bb.n>=20].copy()
    for br in bb.itertuples():regions.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,side=side,feature=feature,bin_start_um=br.bin*200,median_distance_um=br.distance,median_value=br.value,spots=br.n))
    slope=float(np.polyfit(bb.distance.to_numpy()/200,bb.value.to_numpy(),1)[0]) if len(bb)>=3 and np.isfinite(bb.value).all() else np.nan
    rows.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family='boundary_slope_'+side,endpoint=feature+'__distance_slope',exposure=feature,target=side+'_distance',layer='expression_boundary_slope',scale='200um_bins',qualification='atleast3_bins_ge20',model='A0',effect_type='slope_per200um',effect=slope,n_obs=len(bb),reason='' if np.isfinite(slope) else 'FEWER_THAN3_VALID_BINS'))
  print('A2_A4',r.sample,flush=True)
 sec=pd.DataFrame(rows);tab('tables/A2_A4_SECTION_EFFECTS.tsv',sec)
 tab('tables/EXPRESSION_BOUNDARY_BIN_PROFILES.tsv',pd.DataFrame(regions));tab('tables/LR_SPATIAL_COVERAGE.tsv',pd.DataFrame(lrcover))
 pat,su=summarize(sec,'A2_A4')
 # Same-patient effect changes: never infer a change from P crossing0.05.
 pairs=[]
 for (dataset,epid,exposure,target,layer),g in pat[pat.layer.isin(['scale','qualification','adjustment'])].groupby(['dataset','endpoint','exposure','target','layer']):
  col={'scale':'scale','qualification':'qualification','adjustment':'model'}[layer]
  piv=g.pivot(index='patient',columns=col,values='effect')
  contrasts={'scale':[('800','400'),('400','200'),('200','ring')],'qualification':[('Q0','Q1'),('Q1','Q2')],'adjustment':[('A0','A1'),('A1','A2'),('A2','A3'),('A3','A4')]}[layer]
  for a,b in contrasts:
   if a not in piv or b not in piv:continue
   gg=piv[[a,b]].dropna();v=(gg[b]-gg[a]).to_numpy()
   z=infer(v);z={('change_'+k):val for k,val in z.items()}
   pairs.append(dict(dataset=dataset,endpoint=epid,exposure=exposure,target=target,layer=layer,from_level=a,to_level=b,n_common=len(gg),median_from=gg[a].median(),median_to=gg[b].median(),**z))
 paired=pd.DataFrame(pairs);paired['family']=paired.endpoint.map(ep.set_index('endpoint').family);paired['change_FDR']=np.nan
 for k,g in paired.groupby(['dataset','family','layer']):
  pp=g.change_P
  paired.loc[g.index,'change_FDR']=np.where(pp.notna(),stats.false_discovery_control(pp.fillna(1).to_numpy()),np.nan)
 tab('tables/PAIRED_LAYER_CHANGES.tsv',paired)
 save('A2_A4_COMPLETE.json',dict(time=stamp(),section_models=len(sec),estimable=int(sec.effect.notna().sum()),patient_models=len(pat),summaries=len(su),paired_changes=len(paired),new_figures=0))
if __name__=='__main__':
 try:main()
 except Exception:(P/'ERROR_NARROW.txt').write_text(traceback.format_exc());raise
