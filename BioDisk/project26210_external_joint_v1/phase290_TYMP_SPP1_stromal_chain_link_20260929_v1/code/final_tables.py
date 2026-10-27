from common import *
def main():
 pt=read(P/'tables/CHAIN_PATIENT_EFFECTS.tsv');rows=[]
 keys=['dataset','readout','model','family','endpoint']
 for key,g in pt.groupby(keys):
  w=g.pivot(index='patient',columns='radius',values='effect').dropna()
  for radius in w.columns:
   v=w[radius].to_numpy();rows.append(dict(zip(keys,key),radius=int(radius),common_patients=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median_effect=float(np.median(v)) if len(v) else np.nan))
 tab(rows,'COMMON_PATIENT_RADIUS_COMPARISON.tsv')
 s=read(P/'tables/CHAIN_COHORT_SUMMARY.tsv');s['family_size_registered']=2;s['family_size_estimable']=s.groupby(['dataset','radius','readout','model','family']).P.transform(lambda z:int(z.notna().sum()))
 tab(s,'CHAIN_COHORT_SUMMARY.tsv')
 d=read(P/'tables/SAME_PATIENT_CHAIN_DETAILS.tsv')
 d[(d.dataset=='TJRCC')&(d.radius==200)&(d.readout=='SPP1_mapped_share')&(d.model=='core')&(d.source=='tumor')&(d.immune=='CD8')].to_csv(P/'tables/TJRCC_ALL_EVALUABLE_CD8_CHAIN_CASES.tsv',sep='\t',index=False)
 definitions=[R/'phase277_metabolic_macrophage_barrier_cooperation_20260928_v1/harmonized_platform_v1/ADDENDUM_BEFORE_HARMONIZED_MODELS.json',R/'phase289_state_source_boundary_process_20260929_v1/PANELS_FROZEN.json',R/'phase289_state_source_boundary_process_20260929_v1/tables/BOUNDARY_QUALIFICATION.tsv']
 dump(P/'qc/DEFINITION_DEPENDENCIES.json',[dict(path=str(f),sha256=sha(f)) for f in definitions])
 print('FINAL_TABLES_COMPLETE',flush=True)
if __name__=='__main__':main()
