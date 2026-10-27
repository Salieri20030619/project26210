from pathlib import Path
import sys,json,hashlib,platform,datetime
import numpy as np
import pandas as pd
import scipy
from scipy.stats import wilcoxon,norm
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
S=P/'local_records_migration_20260908_150455/project'
O=P/'phase144_parts2_3_4_audit_remote_v1/numerical_v2'
(O/'evidence').mkdir(parents=True,exist_ok=True)
checks=[]; sources=[]; tables={}
def check(name,actual,expected):
 ok=bool(np.isclose(actual,expected,rtol=1e-7,atol=1e-10,equal_nan=True)) if isinstance(actual,(float,int,np.number)) else actual==expected
 checks.append(dict(check=name,actual=str(actual),expected=str(expected),status='PASS' if ok else 'REVIEW'))
def read(key,f):
 f=Path(f);raw=f.read_bytes(); sources.append(dict(key=key,path=str(f),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
 (O/'evidence'/f'{key}{f.suffix}').write_bytes(raw)
 if f.suffix=='.tsv':
  a=pd.read_csv(f,sep='\t',float_precision='round_trip');tables[key]=a;return a
 return raw.decode('utf-8-sig')
def bh(x):
 x=np.array(x,float);o=np.argsort(x);v=np.minimum.accumulate((x[o]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1];r=np.empty(len(x));r[o]=np.minimum(v,1);return r
def bhcheck(key,pcol,qcol,groups):
 a=tables[key]
 for gn,g in (a.groupby(groups,dropna=False) if groups else [('all',a)]):
  q=bh(g[pcol])
  for i,(_,r) in enumerate(g.iterrows()):check(f'{key}:BH:{gn}:{i}',q[i],r[qcol])
def patientcheck(key,patientkey,groups,effect,median,ncol,pcol,method='auto',correction=False):
 a=tables[key];b=tables[patientkey]
 if 'eligible' in b:b=b[b.eligible.astype(str).str.lower()=='true']
 for i,r in a.iterrows():
  g=b.copy()
  for k in groups:g=g[g[k]==r[k]]
  v=g[effect].dropna().to_numpy(float)
  check(f'{key}:{i}:n',len(v),r[ncol]);check(f'{key}:{i}:median',float(np.median(v)),r[median]);check(f'{key}:{i}:positive',int((v>0).sum()),r.n_positive)
  check(f'{key}:{i}:unique_patients',g.patient.nunique(),len(g))
  pv=wilcoxon(v,alternative='two-sided',method=method,correction=correction).pvalue
  check(f'{key}:{i}:P',float(pv),r[pcol])
base=S/'phase33_part2_y02_myeloid_monocle3_local_v1/analysis_local_v1'
for key,name in [('Y02','y02_monocle3_primary_summary_v1.tsv'),('Y02_patients','y02_monocle3_patient_effects_v1.tsv'),('Y02_states','y02_monocle3_state_summary_v1.tsv'),('Y02_concordance','y02_monocle3_dpt_concordance_summary_v1.tsv')]:read(key,base/name)
base=S/'phase34_part2_y03_myeloid_pyscenic_local_v1/analysis_local_v1'
for key,name in [('Y03','y03_candidate_endpoint_summary_v1.tsv'),('Y03_patients','y03_candidate_patient_effects_v1.tsv'),('Y03_studies','y03_candidate_study_stratified_summary_v1.tsv')]:read(key,base/name)
base=P/'phase35_part1_c01_c04_cell_communication_v1'
for key,name in [('C01','c01_cellchat_integrated_v1/c01_cellchat_endpoint_summary_v1.tsv'),('C01_patients','c01_cellchat_integrated_v1/c01_cellchat_patient_endpoint_effects_v1.tsv'),('C03','c03_nichenet_v1/c03_nichenet_ligand_activity_summary_v1.tsv'),('C03_patients','c03_nichenet_v1/c03_nichenet_patient_ligand_activity_v1.tsv'),('C04','c02_c04_analysis_v1/c04_receptor_tymp_f1_summary_v1.tsv'),('C02','c02_c04_analysis_v1/c02_receptor_specificity_summary_v1.tsv'),('LR_mapping','resource_audit_v2/cellchat_targeted_interactions_v1.tsv')]:read(key,base/name)
base=P/'phase36_part1_c05_tcga_joint_survival_v1'
for key,name in [('C05_logrank','analysis_v1/c05_joint_group_logrank_summary_v1.tsv'),('C05_cox','analysis_v1/c05_joint_group_cox_models_v1.tsv'),('C05_interaction','analysis_v1/c05_continuous_interaction_cox_v1.tsv'),('C05_cohort','analysis_v1/c05_tcga_analysis_cohort_v1.tsv'),('C05_survival','analysis_v1/c05_joint_group_survival_summary_v1.tsv'),('C05_registry','provenance/PHASE36_C05_ANALYSIS_REGISTRY_v1.md')]:read(key,base/name)
base=P/'phase109_targeted_gap_analysis_local_v1/analysis_v3'
for key,name in [('C109','cellchat_fine_state_specificity_summary_v1.tsv'),('C109_patients','cellchat_patient_fine_state_specificity_v1.tsv'),('A109','abundance_primary_results_v1.tsv'),('A109_patients','phase109_patient_summary_v1.tsv'),('A109_loso','abundance_leave_one_study_out_v1.tsv')]:read(key,base/name)
read('P109_registry',P/'phase109_targeted_gap_analysis_local_v1/PHASE109_PREREGISTRATION_v1.md')
base=P/'phase37_part4_s02_s04_spatial_metabolism_lr_v1'
for key,name in [('S02','analysis_v1/tables/s02_tjrcc12_inference_v1.tsv'),('S02_sections','analysis_v1/tables/s02_section_statistics_v1.tsv'),('S04','s04_reanalysis_v2/tables/s04_tjrcc12_inference_block_mean_v2.tsv'),('S04_coverage','lr_gene_coverage_v1/lr_pair_scope_summary_v1.tsv'),('S04_registry','provenance/PHASE37_S04_ANALYSIS_REGISTRY_v3.md')]:read(key,base/name)
read('S05',P/'phase38_part4_s05_boundary_sensitivity_v1/analysis_v1/tables/s05_draft_boundary_gradient_summary_v1.tsv')
read('S05_registry',P/'phase38_part4_s05_boundary_sensitivity_v1/provenance/PHASE38_S05_ANALYSIS_REGISTRY_v1.md')
read('S20',S/'phase20_spp1_vascular_immune_local_v1/server_results_v2_sync_v1/tables/primary_across_section_inference_v1.tsv')
base=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase70_cell2location_spatial_composition_v1/results_v1')
for key,name in [('S70_decisions','primary_decisions_v1.tsv'),('S70_meta','composition_adjusted_bridge_meta_v1.tsv'),('S70_summary','PHASE70_COMPOSITION_AWARE_SPATIAL_SUMMARY_v1.md')]:read(key,base/name)
read('S97',S/'phase97_spp1_vascular_spatial_bridge_local_v1/README.md')
read('scope_registry',P/'phase106_dataset_separation_audit_local_v1/DATASET_UNIVERSE_REGISTRY_v1.tsv')
read('old_checklist',S/'phase24_part1_analysis_completion_audit_local_v1/ANALYSIS_COMPLETION_AUDIT_v1.tsv')
patientcheck('Y02','Y02_patients',['study','endpoint'],'patient_effect','median_patient_effect','n_eligible_patients','wilcoxon_p_two_sided')
patientcheck('Y03','Y03_patients',['endpoint','tf'],'patient_effect','median_patient_effect','n_eligible_patients','wilcoxon_p_two_sided',method='exact')
bhcheck('Y02','wilcoxon_p_two_sided','BH_FDR_within_study_primary4',['study'])
bhcheck('Y03','wilcoxon_p_two_sided','BH_FDR_within_endpoint3TF',['endpoint'])
bhcheck('C01','wilcoxon_p_two_sided','BH_FDR_4_endpoints',[])
bhcheck('C03','wilcoxon_p_two_sided','BH_FDR_2_ligands',['family'])
bhcheck('C04','wilcoxon_p_two_sided','BH_FDR_within_outcome8features',['outcome'])
bhcheck('C109','wilcoxon_p','bh_q_4',[])
bhcheck('A109','permutation_p','bh_q_4',[])
bhcheck('S20','wilcoxon_p','wilcoxon_fdr_bh',[])
a=tables['C109'];b=tables['C109_patients']
for i,r in a.iterrows():
 v=b[r.endpoint].to_numpy(float)
 for label,x,y in [('n',len(v),r.n_patients),('median',np.median(v),r['median']),('positive',sum(v>0),r.n_positive),('P',wilcoxon(v,method='approx').pvalue,r.wilcoxon_p)]:check(f'C109:{i}:{label}',x,y)
 check(f'C109:{i}:LOO',all(np.median(np.delete(v,j))>0 for j in range(len(v))),bool(r.all_leave_one_patient_positive))
for key in ['C05_cox','C05_interaction']:
 a=tables[key]
 for i,r in a.iterrows():
  for label,x,y in [('HR',np.exp(r.beta),r.hazard_ratio),('z',r.beta/r.standard_error,r.z),('P',2*norm.sf(abs(r.z)),r.p_two_sided),('CIlo',np.exp(r.beta-norm.ppf(.975)*r.standard_error),r.ci95_low),('CIhi',np.exp(r.beta+norm.ppf(.975)*r.standard_error),r.ci95_high)]:check(f'{key}:{i}:{label}',x,y)
check('Y02:cell_count',int(tables['Y02_states'].n_cells.sum()),50761)
check('Y03:patient_count',tables['Y03_patients'].patient.nunique(),40)
check('Y03:BH_positive',int((tables['Y03'].BH_FDR_within_endpoint3TF<.05).sum()),9)
# Inventory all relevant original result families and rebuilt stages; metadata only for large objects.
roots=[P/x for x in ['phase33_part2_y02_myeloid_monocle3_v1','phase34_part2_y03_myeloid_pyscenic_v1','phase35_part1_c01_c04_cell_communication_v1','phase36_part1_c05_tcga_joint_survival_v1','phase37_part4_s02_s04_spatial_metabolism_lr_v1','phase38_part4_s05_boundary_sensitivity_v1','phase109_targeted_gap_analysis_local_v1']]
roots += [x for x in S.iterdir() if x.is_dir() and x.name.startswith(('phase8_','phase16_','phase17_','phase18_','phase19_','phase20_','phase21_','phase33_','phase34_','phase35_','phase36_','phase37_','phase38_','phase97_'))]
legacy=Path('/mnt/BioDisk/ludingyu/project26210/results')
roots += [legacy/x for x in ['part3','part3_3','part3_liana','part3_macro','part3_macro_subtype','part4_spatial']]
roots += [base.parent]
inv=[]
for root in roots:
 if root.exists():
  for f in root.rglob('*'):
   if f.is_file():
    st=f.stat();inv.append(dict(root=str(root),path=str(f),bytes=st.st_size,mtime_utc=datetime.datetime.fromtimestamp(st.st_mtime,datetime.timezone.utc).isoformat()))
pd.DataFrame(inv).to_csv(O/'FILE_INVENTORY.tsv',sep='\t',index=False)
pd.DataFrame(sources).to_csv(O/'SOURCE_SHA256.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'NUMERICAL_CHECKS.tsv',sep='\t',index=False)
status=dict(audit='EXISTING_RESULT_TABLE_AUDIT_NOT_FULL_MODEL_RERUN',python=sys.version,executable=sys.executable,pandas=pd.__version__,scipy=scipy.__version__,sources=len(sources),inventoried_files=len(inv),checks=len(checks),passed=sum(x['status']=='PASS' for x in checks),review=[x for x in checks if x['status']!='PASS'])
(O/'AUDIT_COMPUTATION_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
print(json.dumps(status,ensure_ascii=False))
print('SCHEMAS',json.dumps({k:list(tables[k].columns) for k in ['C01_patients','C03_patients','C05_logrank','C05_cohort','S02','S02_sections','S70_meta']},ensure_ascii=False))

from scipy.stats import chi2
patientcheck('C01','C01_patients',['endpoint'],'high_minus_low','median_high_minus_low','n_eligible_patients','wilcoxon_p_two_sided',method='exact')
for i,r in tables['C03'].iterrows():
 v=tables['C03_patients'].loc[tables['C03_patients'].ligand==r.ligand,r.endpoint].dropna().to_numpy(float)
 check(f'C03:{i}:median',np.median(v),r.median_patient_effect)
 check(f'C03:{i}:P',wilcoxon(v,method='approx',correction=True).pvalue,r.wilcoxon_p_two_sided)
 check(f'C03:{i}:n',len(v),r.n_eligible_patients)
bhcheck('C05_logrank','logrank_p_two_sided','BH_FDR_within_endpoint6axes',['endpoint'])
bhcheck('S02','wilcoxon_p_two_sided','BH_FDR_within_source4pathways',['source'])
bhcheck('S70_meta','wilcoxon_signed_rank_p','wilcoxon_bh_q',['preparation','analysis_family'])
for i,r in tables['C05_logrank'].iterrows():
 check(f'C05_logrank:{i}:P_from_chisq',chi2.sf(r.chisq,r.degrees_freedom),r.logrank_p_two_sided)
 c=tables['C05_cohort'];g=c[(c[r.endpoint+'_time']>0)&c[r.endpoint+'_event'].notna()]
 check(f'C05_logrank:{i}:n',len(g),r.n_patients)
 check(f'C05_logrank:{i}:events',g[r.endpoint+'_event'].sum(),r.n_events)
check('C05:unique_patients',tables['C05_cohort'].patient.nunique(),510)
# Sensitivity to SciPy version's changed auto handling of tied absolute ranks.
r=tables['Y03_patients'];sens=[]
for (endpoint,tf),g in r[r.eligible.astype(str).str.lower()=='true'].groupby(['endpoint','tf']):
 v=g.patient_effect.to_numpy(float)
 sens.append(dict(endpoint=endpoint,tf=tf,n=len(v),exact_p=wilcoxon(v,method='exact').pvalue,current_auto_p=wilcoxon(v,method='auto').pvalue))
sens=pd.DataFrame(sens);sens['current_auto_q']=sens.groupby('endpoint').current_auto_p.transform(lambda z:bh(z))
sens.to_csv(O/'Y03_SCIPY_AUTO_SENSITIVITY.tsv',sep='\t',index=False)
# Audit-only multiplicity annotation of secondary Cox effects, not a preregistered replacement.
for key in ['C05_cox','C05_interaction']:
 a=tables[key].copy();a['audit_only_BH_within_endpoint_model_term']=a.groupby(['endpoint','model','term']).p_two_sided.transform(lambda z:bh(z)) if 'model' in a else a.groupby(['endpoint','term']).p_two_sided.transform(lambda z:bh(z))
 a.to_csv(O/(key+'_AUDIT_BH.tsv'),sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'NUMERICAL_CHECKS.tsv',sep='\t',index=False)
status.update(checks=len(checks),passed=sum(x['status']=='PASS' for x in checks),review=[x for x in checks if x['status']!='PASS'])
(O/'AUDIT_COMPUTATION_STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
print('FINAL',json.dumps(status,ensure_ascii=False))
print('COX_MODELS',tables['C05_cox'].model.unique())
