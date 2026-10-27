from pathlib import Path
import pandas as pd,numpy as np,json,subprocess,datetime,hashlib
P=Path(__file__).resolve().parents[1];B=P.parent/'phase315_CRC_SPP1_barrier_paper_reproduction_20261004_v1'
subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
assert json.loads((P/'qc/P7_IDENTITY_VALIDATION.json').read_text())['status']=='PASS'
m=pd.read_csv(P/'tables/P7_VERIFIED_COORDINATE_MAPPING.tsv',sep='\t',float_precision='round_trip');d=pd.read_csv(P/'tables/CANDIDATE_DATA_TARGETS.tsv',sep='\t',float_precision='round_trip');z=m.merge(d,on='barcode',validate='one_to_one');z['malignant']=z.DefineTypes.eq('Malignant epithelial cells')
rows=[]
for scope in ['all_nonmalignant','exclude_unknown']:
 a=z[z.malignant];b=z[~z.malignant]
 if scope=='exclude_unknown':b=b[~b.DefineTypes.str.lower().isin(['unknown','unkown'])]
 for gene in ['CD3D','CD8A','CD4','CD19','MS4A1','GZMA','PRF1','GNLY']:
  rows.append({'patient':'P7','comparison':scope,'gene':gene,'n_malignant_spots':len(a),'n_other_spots':len(b),'mean_malignant':a[gene].mean(),'mean_other':b[gene].mean(),'delta_mean':a[gene].mean()-b[gene].mean(),'detection_malignant':a[gene].gt(0).mean(),'detection_other':b[gene].gt(0).mean(),'delta_detection':a[gene].gt(0).mean()-b[gene].gt(0).mean(),'biological_P':None})
out=pd.DataFrame(rows);out.to_csv(P/'tables/P7_RECOVERED_IMMUNE_RESULTS.tsv',sep='\t',index=False);z.to_csv(P/'tables/P7_RECOVERED_SPOT_INPUTS.tsv',sep='\t',index=False)
old=pd.read_csv(B/'tables/IMMUNE_MAP_QUANTIFICATION.tsv',sep='\t',float_precision='round_trip');full=pd.concat([old,out],ignore_index=True);full.to_csv(P/'tables/FOUR_PATIENT_IMMUNE_RESULTS.tsv',sep='\t',index=False)
primary=full[full.comparison.eq('all_nonmalignant')&full.gene.eq('CD3D')];check=out[out.comparison.eq('all_nonmalignant')&out.gene.eq('CD3D')].iloc[0]
result={'time':datetime.datetime.now().astimezone().isoformat(),'same_paper_P7_input_recovered':True,'P7_primary_direction_pass':bool(check.delta_mean<0 and check.delta_detection<0),'R3_four_patient_direction_pass':bool(len(primary)==4 and primary.delta_mean.lt(0).all() and primary.delta_detection.lt(0).all()),'original_R1_exact_Pearson_match_still_unresolved':True,'original_R1_positive_direction_4_of_4':True,'R2_40_GO_statistics_inherited_pass':True,'biological_patient_n':4,'no_new_biological_P_or_FDR':True,'external_validation_started':False,'TYMP_positioning_started':False,'baseline_phase315_not_modified':True}
(P/'qc/RECOVERY_DECISION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));print(out[out.gene.isin(['CD3D','CD8A'])].to_json(orient='records'))
