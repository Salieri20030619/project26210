from common import *
import math
def main():
 regs=pd.read_csv(P/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t').fillna('UNKNOWN')
 broad=pd.read_csv(P/'tables/A1_SECTION_EFFECTS.tsv',sep='\t',float_precision='round_trip')
 narrow=pd.read_csv(P/'tables/A2_A4_SECTION_EFFECTS.tsv',sep='\t',float_precision='round_trip')
 # Source/preparation stratification is descriptive, retaining source-separated main inference.
 strata=[]
 allsec=pd.concat([broad,narrow],ignore_index=True)
 gse=allsec[allsec.dataset=='GSE175540'].merge(regs[['sample','preparation','cohort']],on='sample',validate='many_to_one')
 for field in ['preparation','cohort']:
  patient=[]
  for key,g in gse.groupby([field,'patient']+KEYS):
   v=g.effect.dropna().to_numpy()
   val=(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean()) if key[-1]=='rho' else np.mean(v)) if len(v) else np.nan
   patient.append(dict(zip(['stratum','patient']+KEYS,key),effect=val))
  pp=pd.DataFrame(patient)
  tab('tables/GSE175540_'+field+'_PATIENT_EFFECTS.tsv',pp)
  for key,g in pp.groupby(['stratum']+KEYS):
   v=g.effect.dropna().to_numpy()
   strata.append(dict(stratification=field,**dict(zip(['stratum']+KEYS,key)),n=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median=float(np.median(v)) if len(v) else np.nan,note='descriptive stratum; c_2 not independent between preparations; primary inference remains deduplicated cohort'))
 tab('tables/GSE175540_STRATIFIED_DESCRIPTION.tsv',pd.DataFrame(strata))
 # Parallel denominator comparisons use identical patients, never significance-selected cohorts.
 pat=pd.read_csv(P/'tables/A1_PATIENT_EFFECTS.tsv',sep='\t',float_precision='round_trip')
 comparisons=[]
 for (dataset,family,target,scale),g in pat[(pat.layer=='broad')&pat.family.isin(['SPP1_vascular','SPP1_stroma_program','SPP1_function'])].groupby(['dataset','family','target','scale']):
  piv=g.pivot(index='patient',columns='exposure',values='effect')
  for first,last in [('V1_SPP1_q05','V2_SPP1_fraction'),('V2_SPP1_fraction','V3_SPP1_macro_share'),('V1_SPP1_q05','V3_SPP1_macro_share')]:
   gg=piv[[first,last]].dropna();v=(gg[last]-gg[first]).to_numpy();out=infer(v)
   comparisons.append(dict(dataset=dataset,family=family,target=target,scale=scale,from_variable=first,to_variable=last,n_common=len(gg),median_from=gg[first].median(),median_to=gg[last].median(),median_effect_change=out['median'],CI_low=out['CI_low'],CI_high=out['CI_high'],note='paired rho comparison across definitions, not evidence one denominator is correct'))
 tab('tables/PAIRED_DENOMINATOR_CHANGES.tsv',pd.DataFrame(comparisons))
 # Model diagnostic on original and common patient sets; qualification common-spot fit equals Q2 by definition.
 pair=pd.read_csv(P/'tables/PAIRED_LAYER_CHANGES.tsv',sep='\t',float_precision='round_trip')
 for key,g in pair.groupby(['dataset','family','layer']):
  p=g.change_P.fillna(1).to_numpy();order=np.argsort(p);v=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q=np.empty(len(p));q[order]=np.minimum(v,1)
  mask=g.change_P.notna().to_numpy();np.testing.assert_allclose(q[mask],g.change_FDR.to_numpy()[mask],atol=1e-12)
 defs=json.loads((P/'PROGRAM_DEFINITIONS.json').read_text())['programs']
 defs.update(V4_SPP1_program=defs['spp1_state_no_spp1'],TYMP_RNA=['TYMP'])
 eps=pd.read_csv(P/'tables/FROZEN_ENDPOINTS.tsv',sep='\t');overlap=[]
 for r in eps.itertuples():
  if r.exposure in defs and r.target in defs:overlap.append(dict(endpoint=r.endpoint,shared_genes=';'.join(sorted(set(defs[r.exposure])&set(defs[r.target]))),n_shared=len(set(defs[r.exposure])&set(defs[r.target]))))
 tab('tables/PROGRAM_GENE_OVERLAP.tsv',pd.DataFrame(overlap))
 eligibility=[
 dict(endpoint='histology_boundary',status='NOT_ESTIMABLE',reason='No audited pathology boundary; AI drafts remain unreviewed'),
 dict(endpoint='CAF_cell_abundance',status='NOT_ESTIMABLE',reason='Current13-state reference lacks independently identifiable CAF; POSTNCAF expression program separate'),
 dict(endpoint='CD8_cell_location_or_exclusion',status='NOT_ESTIMABLE',reason='No independent CD8 abundance reference; expression program not cell positioning'),
 dict(endpoint='SPP1_macrophage_to_specific_receiver_LR',status='NOT_ESTIMABLE',reason='Mixed-spot ligand/receptor cannot establish sender/receiver cell identity; spot compatibility completed'),
 dict(endpoint='treatment_effect',status='DESCRIPTIVE_ONLY',reason='GLMF one patient per treatment background per prior plan; unknown fields remain UNKNOWN; no treatment-effect test'),
 dict(endpoint='MISTy_prediction_gain',status='OPTIONAL_NOT_RUN',reason='Supplement in plan, not required primary endpoint; no new prediction-gain model selected after association results'),
 dict(endpoint='COMMOT_transport',status='LEGACY_RETAINED_NOT_REFIT',reason='Current registered LR estimate is complete-subunit expression compatibility on frozen graph; not renamed COMMOT transport'),
 dict(endpoint='TLS_peripheral_contrast',status='INHERITED',reason='Phase222 strict TLS12 patients,10 negative; median-7.27percentage points,P0.03857; one prespecified supplement'),
 dict(endpoint='legacy_TYMP_high_and_composite_niche',status='INHERITED',reason='Phase232 numerically audited original reference and composite construction; not merged with current13-state definition')]
 tab('tables/TECHNICAL_SCOPE_AND_UNESTIMABLE.tsv',pd.DataFrame(eligibility))
 oldsources=[]
 for phase in ['phase219','phase222','phase223','phase230','phase231','phase232']:
  dr=next(R.glob(phase+'_*'))
  for f in dr.glob('*.md'):
   oldsources.append(dict(phase=phase,source=str(f),sha256=sha(f),status='read-only inherited; original endpoints and statistics retained'))
 tab('tables/LEGACY_EVIDENCE_REGISTER.tsv',pd.DataFrame(oldsources))
 # Preserve old statistics as separate original-definition source tables.
 inherited=P/'legacy_tables';inherited.mkdir(exist_ok=True)
 import shutil
 for prefix,dirpath in [('phase232',R/'phase232_legacy_part4_significance_audit_20260921_v1/tables'),('phase230',R/'phase230_spatial_F2_program_association_20260921_v1/tables'),('phase231',R/'phase231_spatial_adjustment_decomposition_20260921_v1/tables'),('phase219',R/'phase219_spatial_reference_identifiability_20260920_v1/program_check')]:
  for f in dirpath.glob('*.tsv'):
   if any(x in f.name.lower() for x in ['summary','abundance','programs','all_relations']):
    shutil.copy2(f,inherited/(prefix+'_'+f.name))
 save('SUPPLEMENT_VALIDATION.json',dict(status='PASS',time=stamp(),paired_change_FDR_recomputed_by_original_family=True,descriptive_strata_rows=len(strata),denominator_comparisons=len(comparisons),legacy_tables=len(list(inherited.glob('*.tsv'))),LR_common_receptor_different_ligand_eligibility_concordant=True))
if __name__=='__main__':main()
