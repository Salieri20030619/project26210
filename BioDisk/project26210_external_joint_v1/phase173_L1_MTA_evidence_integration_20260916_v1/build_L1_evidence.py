from pathlib import Path
import json,csv,hashlib,shutil
import pandas as pd,numpy as np
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=R/'phase173_L1_MTA_evidence_integration_20260916_v1'
P143=R/'phase143_sections25_27_candidate_report_remote_v1'
P160=R/'phase160_external_bridge_investigation_20260915_v1'
P162=R/'phase162_MTA_polyamine_part1_reassessment_20260915_v1'
P163=R/'phase163_polyamine_MTA_serial_receiver_20260915_v1'
P98=R.parent/'project26210_phase98_public_paired_metabolomics_v1'
sources=[];checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def track(p):
 p=Path(p);sources.append(dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size));return p
def read(p):
 return pd.read_csv(track(p),sep='\t')
def check(k,v):
 checks.append(dict(check=k,passed=bool(v)))
 if not v:raise RuntimeError(k)
E=json.loads(track(P143/'qc/EVIDENCE.json').read_text())
material=json.loads(track(P143/'report/report_materialized.json').read_text())
ledger=[]
for key in ['auc_pair','auc_cor','flux_pair','flux_cor','enzyme_f1','terminal_enzyme','f1_tymp']:
 a=pd.DataFrame(E[key]);a.to_csv(O/('L1_SOURCE_'+key+'.tsv'),sep='\t',index=False)
 for i,r in a.iterrows():
  r=r.to_dict()
  if key in ['auc_pair','flux_pair']:
   n=int(r['n_patients']);effect=float(r['median_patient_delta']);P=float(r['paired_wilcox_PValue']);F=float(r['FDR']);pos=round(float(r['fraction_patients_positive'])*n);endpoint=r['contrast']+'::'+r['pathway'];etype='median_paired_score_delta';family='four_features_per_contrast'
  elif key in ['auc_cor','flux_cor']:
   n=int(r['n_patients']);effect=float(r['median_patient_rho']);P=float(r['signed_rank_PValue']);F=float(r['FDR']);pos=round(float(r['fraction_positive'])*n);endpoint=r['predictor']+'::'+r['pathway'];etype='median_patient_rho';family='four_features_per_predictor'
  elif key=='enzyme_f1':
   n=int(r['n_patients']);effect=float(r['median_rho']);P=float(r['P']);F=float(r['FDR']);pos=int(r['positive_patients']);endpoint='F1::'+r['gene'];etype='median_patient_rho';family='original_12_genes'
  elif key=='terminal_enzyme':
   n=int(r['n']);effect=float(r['median_effect']);P=float(r['P']);F=float(r['FDR']);pos=int(r['positive']);endpoint=r['method']+'::'+r['gene'];etype='terminal_minus_early_log2CPM';family='original_39_gene_method_endpoints'
  else:
   continue
  tier='FDR_supported' if F<.05 else 'nominal_only' if P<.05 else 'direction_or_null'
  ledger.append(dict(evidence_id='P143_'+key+'_'+str(i+1),dataset_scope='original_five_scRNA_studies',compartment='malignant_epithelial',endpoint=endpoint,n_biological=n,positive=pos,effect_type=etype,effect=effect,P=P,FDR=F,FDR_family=family,evidence_tier=tier,source=str(P143/'qc/EVIDENCE.json'),source_json_key=key,source_row=i+1,independent_new_validation=False))
part=pd.DataFrame(ledger);part.to_csv(O/'L1_TYMP_METABOLIC_EVIDENCE.tsv',sep='\t',index=False)
# Confirm original report rendered rounded values match source enzyme table.
t=next(x for x in material if x.get('type')=='table' and '十二基因' in x.get('caption',''))
a=pd.DataFrame(E['enzyme_f1']).set_index('gene')
for row in t['rows']:
 check('report_enzyme_'+row[0],int(a.loc[row[0],'n_patients'])==int(row[1]) and abs(float(a.loc[row[0],'median_rho'])-float(row[2]))<=.0000051 and abs(float(a.loc[row[0],'P'])-float(row[4]))<=.0000051 and abs(float(a.loc[row[0],'FDR'])-float(row[5]))<=.0000051)
auc=read(P162/'FIRST_PART_AND_LEGACY_EVIDENCE/four_pathway_AUCell_paired_statistics_v1.tsv')
check('report_AUCell_source_matches',np.allclose(auc.FDR.astype(float),pd.DataFrame(E['auc_pair']).FDR.astype(float)))
contrast=read(P162/'PART1_GENE_CONTRAST_SUMMARY.tsv');contrast.to_csv(O/'L1_SUPPLEMENT_PART1_CONTRASTS.tsv',sep='\t',index=False)
pub160=read(P160/'META_BRIDGE_RESULTS.tsv')
pub160=pub160[pub160.candidate=='MET:5-methylthioadenosine (MTA)'].copy()
pub162=read(P162/'PUBLIC_META_AXIS_ASSOCIATIONS.tsv')
pub163=read(P163/'META_CONDITIONAL_EDGES.tsv')
old=read(P98/'analysis_v2/cross_cohort_meta_analysis_v2.tsv')
old.to_csv(O/'L1_LEGACY_MTA_FULL_CONTEXT.tsv',sep='\t',index=False)
pub162.to_csv(O/'L1_SOURCE_PUBLIC_AXIS_FULL.tsv',sep='\t',index=False)
pub163.to_csv(O/'L1_SOURCE_SERIAL_CONDITIONAL_FULL.tsv',sep='\t',index=False)
track(P160/'EXTERNAL_BRIDGE_AND_CLIENT_CONTEXT_CN.md');track(P162/'MTA_POLYAMINE_PART1_REASSESSMENT_CN.md');track(P163/'SERIAL_MTA_AND_RECEIVER_REVIEW_CN.md')
serial=read(P163/'COHORT_SERIAL_DECOMPOSITION.tsv');serial.to_csv(O/'L1_SERIAL_MODEL_CONTEXT.tsv',sep='\t',index=False)
public=[]
for i,r in pub160.iterrows():
 hist=old[old.endpoint==('A_B_TYMP_MTA_MARGINAL' if r.edge=='a' else 'B_C_MTA_SPP1_MARGINAL' if r.edge=='XY' else 'NO_MATCH')]
 historical=float(hist.iloc[0].fixed_q_across_endpoints) if r.adjustment=='raw' and len(hist) else np.nan
 public.append(dict(source_stage='Phase160',source=str(P160/'META_BRIDGE_RESULTS.tsv'),source_row=int(i)+1,candidate=r.candidate,outcome='MTA' if r.edge=='a' else 'tissue_SPP1',edge=r.edge,edge_definition={'a':'TYMP--MTA','XY':'MTA--tissue_SPP1','b':'MTA--tissue_SPP1_given_TYMP'}[r.edge],adjustment=r.adjustment,n=int(r.total_n),cohorts=int(r.cohorts),positive_cohorts=int(r.positive_cohorts),rho=float(r.random_rho),P=float(r.random_P),source_FDR=float(r.random_FDR),historical_fixed_FDR=historical,FDR_family='Phase160_MTA_context_single_endpoint_per_edge_adjustment',I2=float(r.I2),reused_same_patients=True))
for i,r in pub162.iterrows():
 if r.gene not in ['SRM','SMS','MTAP','ODC1','AMD1','MAT2A']:continue
 out='MTA' if 'MTA' in r.outcome else 'tissue_SPP1'
 edge={'a':'TYMP--'+r.gene,'XY':r.gene+'--'+out,'b':r.gene+'--'+out+'_given_TYMP','b_adjust_MTAP':r.gene+'--MTA_given_TYMP_MTAP'}[r.edge]
 public.append(dict(source_stage='Phase162',source=str(P162/'PUBLIC_META_AXIS_ASSOCIATIONS.tsv'),source_row=int(i)+1,candidate=r.gene,outcome=out,edge=r.edge,edge_definition=edge,adjustment=r.adjustment,n=int(r.total_n),cohorts=int(r.cohorts),positive_cohorts=int(r.positive_cohorts),rho=float(r.random_rho),P=float(r.random_P),source_FDR=float(r.random_FDR),historical_fixed_FDR=np.nan,FDR_family='Phase162_gene6_per_outcome_edge_adjustment',I2=float(r.I2),reused_same_patients=True))
pub=pd.DataFrame(public);pub.to_csv(O/'L1_PUBLIC_MTA_EVIDENCE.tsv',sep='\t',index=False)
# Read-only source rows and input identity; do not re-fit tests or recalculate FDR.
audit=[]
inventory=json.loads(track(P160/'INPUT_INVENTORY.json').read_text())
for x in inventory:
 cohort=x['cohort'];f=read(P162/('PUBLIC_PATIENTS_'+cohort+'.tsv'))
 feature='MET:5-methylthioadenosine (MTA)'
 n=len(f);ids=f.iloc[:,0]
 check('patient_id_'+cohort,ids.is_unique)
 check('MTA_observed_'+cohort,feature in f and f[feature].notna().all() and '5-methylthioadenosine (MTA)' in x['metabolite_hits'])
 audit.append(dict(cohort=cohort,patients=n,patient_id_unique=ids.is_unique,measured_feature='5-methylthioadenosine (MTA)',complete_MTA=int(f[feature].notna().sum()),measurement='observed_reference_tissue_metabolomics',original_met_file=x['met_file'],patient_values=str(P162/('PUBLIC_PATIENTS_'+cohort+'.tsv')),no_predicted_metabolite_used=True))
check('total_patients_255',sum(x['patients'] for x in audit)==255)
pd.DataFrame(audit).to_csv(O/'L1_PUBLIC_SAMPLE_MEASUREMENT_AUDIT.tsv',sep='\t',index=False)
track(P160/'PATIENT_ID_AUDIT.tsv')
rx=[]
for sid,gene,role,claim in [
('R-HSA-112265','TYMP','nucleoside_entry','Thymidine/deoxyuridine phosphorolysis; reaction does not establish direct MTA synthesis'),
('R-HSA-351222','AMD1','dcSAM_supply','SAM decarboxylation supplies aminopropyl donor; RNA is not substrate flux'),
('R-HSA-351215','SRM','MTA_generation','Spermidine synthesis produces MTA'),
('R-HSA-351210','SMS','MTA_generation','Spermine synthesis produces MTA'),
('R-HSA-1237160','MTAP','MTA_clearance_salvage','MTA phosphorolysis yields adenine and methylthioribose phosphate')]:
 file=track(O/'sources'/(sid+'.json'));j=json.loads(file.read_text())
 check('reaction_identity_'+sid,j['stId']==sid)
 rx.append(dict(reaction_id=sid,gene=gene,role=role,reaction_name=j['displayName'],claim=claim,evidence_type='curated_biochemical_reaction',project_inference=False,url='https://reactome.org/content/detail/'+sid,archive=str(file)))
pd.DataFrame(rx).to_csv(O/'L1_BIOCHEMICAL_LINKS.tsv',sep='\t',index=False)
# independent source-value equality checks, preserving all original tests
for r in pub.itertuples(index=False):
 origin=pub160 if r.source_stage=='Phase160' else pub162
 orig=origin.loc[r.source_row-1]
 check('public_source_'+r.source_stage+'_'+str(r.source_row),float(orig.random_rho)==r.rho and float(orig.random_P)==r.P and float(orig.random_FDR)==r.source_FDR and int(orig.total_n)==r.n)
for x in sources:check('hash_unchanged_'+Path(x['path']).name,sha(Path(x['path']))==x['sha256'])
pd.DataFrame(sources).to_csv(O/'SOURCE_MANIFEST.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'VALIDATION.tsv',sep='\t',index=False)
status=dict(stage='Phase173',step='L1-01_to_L1-04',source_evidence_rows=len(part),public_evidence_rows=len(pub),biochemical_sources=len(rx),public_patients=255,checks=len(checks),passed=sum(x['passed'] for x in checks),new_statistical_tests=0,new_figures=0,first_part_recomputed=False)
(O/'EVIDENCE_BUILD_STATUS.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status))
print('CANDIDATE SMS',pub162[(pub162.gene=='SMS')&(pub162.edge=='b_adjust_MTAP')&(pub162.adjustment=='myeloid')].to_string(index=False))
print('SERIAL columns',serial.columns.tolist())
