from common import *
import traceback
def main():
 assert (P/'A0_FINAL_FREEZE.json').exists()
 regs=pd.read_csv(P/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t')
 ep=pd.read_csv(P/'tables/FROZEN_ENDPOINTS.tsv',sep='\t');rows=[]
 for r in regs.itertuples():
  d=pd.read_csv(P/'inputs'/(r.sample+'_INPUTS.tsv.gz'),sep='\t')
  for e in ep.itertuples():
   x=d[e.exposure].to_numpy();y=d[e.target].to_numpy()
   for scale in ['W0','W1_4x4','W1_5x5']:
    mask=np.isfinite(x)&np.isfinite(y);nspot=int(mask.sum())
    if scale=='W0':effect,n,reason=fit(x,y,mask=mask)
    elif nspot<50:effect,n,reason=np.nan,0,'SOURCE_LT50'
    else:
     k=4 if scale=='W1_4x4' else 5
     b=pd.DataFrame(dict(x=x[mask],y=y[mask],block=d['block_'+str(k)].to_numpy()[mask])).groupby('block')[['x','y']].median()
     effect,n,reason=fit(b.x,b.y,min_n=8)
    rows.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family=e.family,endpoint=e.endpoint,exposure=e.exposure,target=e.target,layer='broad',scale=scale,qualification='Q0',model='A0',effect_type='rho',effect=effect,n_obs=n,n_spots=nspot,reason=reason))
   if e.family in ['KEGG_four_axes','custom_four_axes','MTA_context_SPP1']:
    mask=np.isfinite(x)&np.isfinite(y)&d.tumor_expression_region.to_numpy(bool)
    effect,n,reason=fit(x,y,mask=mask)
    rows.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family=e.family,endpoint=e.endpoint,exposure=e.exposure,target=e.target,layer='malignant_region',scale='W0',qualification='malignant_fraction_ge0.2',model='A0',effect_type='rho',effect=effect,n_obs=n,n_spots=n,reason=reason))
  print('A1',r.sample,flush=True)
 sec=pd.DataFrame(rows);tab('tables/A1_SECTION_EFFECTS.tsv',sec);pat,su=summarize(sec,'A1')
 save('A1_COMPLETE.json',dict(time=stamp(),section_models=len(sec),estimable=int(sec.effect.notna().sum()),patient_models=len(pat),summaries=len(su)))
 print(su[(su.P<.05)&(su.scale=='W0')][['dataset','family','exposure','target','n','median','positive','P','FDR']].to_string(index=False),flush=True)
if __name__=='__main__':
 try:main()
 except Exception:(P/'ERROR_A1.txt').write_text(traceback.format_exc());raise
