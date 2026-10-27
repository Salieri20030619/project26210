from common import *
CORE=['inside_depth','outside_depth','inside_macro','outside_macro','inside_malignant','outside_malignant','log_side_count_ratio_new']
ORG=['outside_C4','outside_POSTN','outside_ECM','SPP1_POSTN_coverage']
IMM=['CD8_RNA_delta','T_abundance_log2_delta']
def main():
 gate('spatial_start')
 reg=read(R/'phase277_metabolic_macrophage_barrier_cooperation_20260928_v1/tables/INPUT_SECTION_REGISTER.tsv')
 configs=[('unadjusted',[]),('core_adjusted',CORE),('inside_TYMP_adjusted',CORE+['inside_TYMP']),('vascular_stromal_sensitivity',CORE+['inside_TYMP','inside_endothelial','outside_endothelial','fibroblast_mean'])]
 effects=[];quals=[];manifest=[]
 for r in reg.itertuples():
  f=PRE/'inputs'/f'{r.sample}_SPLIT_INTERFACE.tsv.gz';manifest.append(dict(path=str(f),sha256=sha(f)))
  z=read(f);mask=z.split_eligible.to_numpy(bool)
  for family,ends in [('organization',ORG),('immune_distribution',IMM)]:
   for y in ends:
    for model,covs in configs:
     e,n,reason=fit(z.outside_TYMP.to_numpy(),z[y].to_numpy(),z[covs].to_numpy() if covs else None,mask,min_n=50)
     effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family=family,endpoint=y,model=model,effect=e,segments=n,reason=reason))
  for y in IMM:
   covs=CORE+['inside_TYMP','SPP1_POSTN_coverage']
   e,n,reason=fit(z.outside_TYMP.to_numpy(),z[y].to_numpy(),z[covs].to_numpy(),mask,min_n=50)
   effects.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,family='niche_conditioned_sensitivity',endpoint=y,model='inside_TYMP_and_coverage',effect=e,segments=n,reason=reason))
  quals.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,segments=len(z),split_eligible=int(mask.sum()),minimum_segments=50,source_interpretation='outer_compartment_mixed_RNA_not_cell_specific'))
  print('OUTER_INTERFACE',r.sample,flush=True)
 e=pd.DataFrame(effects);tab(e,'OUTER_TYMP_SECTION_EFFECTS.tsv');tab(quals,'SPATIAL_QUALIFICATION.tsv')
 keys=['dataset','family','endpoint','model'];pa=[]
 for key,z in e.groupby(keys+['patient']):
  v=z.effect.dropna().to_numpy();value=float(np.tanh(np.arctanh(np.clip(v,-.99999999,.99999999)).mean())) if len(v) else np.nan
  pa.append(dict(zip(keys+['patient'],key),effect=value,sections=';'.join(z.loc[z.effect.notna(),'sample']),n_sections=len(v)))
 pa=pd.DataFrame(pa);tab(pa,'OUTER_TYMP_PATIENT_EFFECTS.tsv');summary(pa,keys,['dataset','family'],'OUTER_TYMP_COHORT_SUMMARY.tsv')
 loo=[]
 for key,z in pa.groupby(keys):
  z=z[z.effect.notna()]
  if len(z)<10:continue
  for pat in z.patient:
   v=z.loc[z.patient!=pat,'effect'];loo.append(dict(zip(keys,key),omitted_patient=pat,n=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median=v.median()))
 tab(loo,'OUTER_TYMP_LEAVE_PATIENT.tsv')
 old=read(PRE/'tables/SPATIAL_COHORT_SUMMARY.tsv');old['inheritance']='Phase280 original inner-TYMP analysis, original P/FDR unchanged';tab(old,'INHERITED_INNER_TYMP_ALL_RESULTS.tsv')
 for rec in manifest:assert sha(rec['path'])==rec['sha256']
 dump(P/'qc/SPATIAL_INPUTS.json',manifest);dump(P/'qc/SPATIAL_COMPLETE.json',dict(status='COMPUTED',section_effects=len(e),patient_effects=len(pa),sections=len(reg)))
 print('SPATIAL_COMPLETE',flush=True)
if __name__=='__main__':main()
