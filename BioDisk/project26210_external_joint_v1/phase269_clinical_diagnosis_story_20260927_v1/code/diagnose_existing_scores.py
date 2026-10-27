from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
P=Path(__file__).resolve().parents[1];H=P.parent.parent
paths={
's81':H/'project26210_phase81_bulk_ici_transcriptome_rebuild_v1/results_v1/phase81_sample_scores_v1.tsv',
's83':H/'project26210_phase83_bulk_state_projection_v1/results_v1/clinical_v1/phase83_sample_scores_v1.tsv',
'cm':H/'project26210_external_joint_v1_phase64_large_immune_validation_v1/inputs_v1/CheckMate_Braun2020/processed_v1/checkmate_clinical_rna_subset_v1.tsv'}
a=pd.read_csv(paths['s81'],sep='\t');b=pd.read_csv(paths['s83'],sep='\t');quality=[];correlations=[]
for dataset,d in b.groupby('dataset'):
 row=dict(dataset=dataset,n=len(d),other_macrophage_zero_percent=float((d.other_macrophage_weight==0).mean()*100),c1qc_zero_percent=float((d.c1qc_state_weight==0).mean()*100),SPP1_fraction_ceiling_percent=float(np.isclose(d.spp1_macrophage_fraction,1).mean()*100),SPP1_C1QC_fraction_sum_max_error=float(np.max(abs(d.spp1_macrophage_fraction+d.c1qc_macrophage_fraction-1))),SPP1_fraction_sd=float(d.spp1_macrophage_fraction.std()),SPP1_fraction_IQR=float(d.spp1_macrophage_fraction.quantile(.75)-d.spp1_macrophage_fraction.quantile(.25)))
 quality.append(row)
 for outcome in ['spp1_state_weight','total_macrophage_weight','spp1_macrophage_fraction','official_macro']:
  v=d[['tymp_z',outcome]].dropna();correlations.append(dict(dataset=dataset,predictor='tymp_z',outcome=outcome,n=len(v),descriptive_spearman=float(v.corr(method='spearman').iloc[0,1]),new_P_test=False))
for dataset,d in a.groupby('dataset'):
 for x,y in [('tymp_z','mta_response'),('tymp_z','spp1_macrophage'),('mta_response','spp1_macrophage')]:
  v=d[[x,y]].dropna();correlations.append(dict(dataset=dataset,predictor=x,outcome=y,n=len(v),descriptive_spearman=float(v.corr(method='spearman').iloc[0,1]),new_P_test=False))
pd.DataFrame(quality).to_csv(P/'tables/PROJECTION_DIAGNOSTICS.tsv',sep='\t',index=False)
pd.DataFrame(correlations).to_csv(P/'tables/EXISTING_SCORE_DESCRIPTIVE_CORRELATIONS.tsv',sep='\t',index=False)
modelpaths=[]
for pattern in ['project26210_phase81_*/results_v1/*models*.tsv','project26210_phase81_*/results_v1/*sensitivity*.tsv','project26210_phase82_*/results_v1/*models*.tsv','project26210_phase83_*/results_v1/clinical_v1/*models*.tsv','project26210_phase83_*/results_v1/clinical_v1/*sensitivity*.tsv','project26210_phase84_*/results_v1/coupling_v1/*models*.tsv','project26210_phase84_*/results_v1/coupling_v1/*outcomes*.tsv','project26210_phase85_*/results_v1/decomposition_v1/*models*.tsv']:
 modelpaths.extend(H.glob(pattern))
ph=[];models=[]
for f in sorted(set(modelpaths)):
 d=pd.read_csv(f,sep='\t');d['source_path']=str(f);models.append(d)
 if 'ph_p_value' in d:
  ph.extend(d[d.ph_p_value<.05][[c for c in ['cohort','model_id','feature','n','n_event','p_value','q_value','ph_p_value','source_path'] if c in d]].to_dict('records'))
allmodels=pd.concat(models,ignore_index=True);allmodels.to_csv(P/'tables/CORRECTED_STAGE81_85_EXISTING_MODELS.tsv',sep='\t',index=False)
pd.DataFrame(ph).to_csv(P/'tables/EXISTING_PH_DIAGNOSTIC_FLAGS.tsv',sep='\t',index=False)
inputs=list(paths.values())+list(set(modelpaths))+[H/'project26210_phase83_bulk_state_projection_v1/scripts_v1/run_phase83_clinical_models_v1.R',H/'project26210_phase81_bulk_ici_transcriptome_rebuild_v1/PHASE81_SIGNATURE_MEMBERS_v1.tsv']
(P/'qc/SOURCE_MANIFEST.json').write_text(json.dumps([dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in inputs],indent=2))
print('PROJECTION',pd.DataFrame(quality).to_string(index=False));print('CORRELATIONS',pd.DataFrame(correlations).to_string(index=False));print('PH_FLAGS',len(ph),'EXISTING_MODELS',len(allmodels))
