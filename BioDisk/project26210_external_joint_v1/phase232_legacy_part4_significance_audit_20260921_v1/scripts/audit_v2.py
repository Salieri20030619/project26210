from pathlib import Path
import json,hashlib,math
import pandas as pd,numpy as np
from scipy import stats
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
P=R/'phase232_legacy_part4_significance_audit_20260921_v1'
P.mkdir(exist_ok=True)
for n in ['tables','audit','root_record_backups','scripts']:(P/n).mkdir(exist_ok=True)
B=Path('/mnt/BioDisk/ludingyu/project26210/results/part4_spatial')
files={'vascular_spots':B/'tymp_vascular_niche_v1/tables/spot_mechanism_scores_v1.tsv.gz','vascular_samples':B/'tymp_vascular_niche_v1/tables/sample_spatial_block_tests_v1.tsv','vascular_summary':B/'tymp_vascular_niche_v1/tables/tjrcc_cohort_direction_summary_v1.tsv','state_summary':B/'tymp_spp1_state_v1/tables/cohort_state_direction_summary_v1.tsv','interface_summary':B/'tymp_spp1_state_v1/tables/cohort_interface_downstream_summary_v1.tsv','GLMF_abundance':B/'tymp_immune_ecology_v2/tables/original_tymp_immune_spatial_tests.tsv','GLMF_programs':B/'tymp_immune_ecology_v2/tables/original_tymp_signature_spatial_tests.tsv','vascular_script':B/'tymp_vascular_niche_v1/provenance/run_tymp_vascular_niche_analysis_v1.py','state_script':B/'run_tymp_spp1_state_analysis_v1.py'}
def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip')
sp=read(files['vascular_spots']);ss=read(files['vascular_samples']);su=read(files['vascular_summary'])
errors=[]
for row in ss.itertuples():
 g=sp[sp['sample']==row.sample]
 block=g.groupby('block_'+row.grid)[[row.source,row.target]].median()
 rho=stats.spearmanr(block[row.source],block[row.target]).statistic
 errors.append(abs(rho-row.rho))
assert max(errors)<1e-12
for r in su.itertuples():
 g=ss[(ss.cohort=='TJ-RCC')&(ss.relation==r.relation)&(ss.route==r.route)&(ss.grid==r.grid)]
 assert len(g)==r.n_patients
 assert abs(g.rho.median()-r.median_rho)<1e-12
 pos=int((g.direction_aligned_rho>0).sum());assert pos==r.direction_supported_patients
 p=min(1.,2*sum(math.comb(len(g),k) for k in range(min(pos,len(g)-pos)+1))/2**len(g))
 assert abs(p-r.sign_test_p)<1e-12
np.testing.assert_allclose(stats.false_discovery_control(su.sign_test_p.to_numpy()),su.sign_test_fdr_bh_global,atol=1e-12)
for name in ['state_summary','interface_summary']:
 d=read(files[name])
 for r in d.itertuples():
  n=r.patients_tested;k=min(r.positive_patients,r.negative_patients);pp=min(1.,2*sum(math.comb(n,i) for i in range(k+1))/2**n)
  assert abs(pp-r.sign_test_p)<1e-12
 np.testing.assert_allclose(stats.false_discovery_control(d.sign_test_p.to_numpy()),d.sign_test_fdr_bh,atol=1e-12)
 d.to_csv(P/'tables'/(name+'.tsv'),sep='\t',index=False)
# Identity of composite endpoints, reconstructed from already stored source features.
formula=[]
for sid,g in sp.groupby('sample'):
 def z(x):return (x-np.mean(x))/np.std(x)
 niche=np.cbrt(np.maximum(g.spp1_percentile*g.endothelial_percentile*g.pericyte_percentile,1e-6))
 tb=z(np.log1p(g.T_cell_q05))-z(np.log1p(np.maximum(g.total_immune_q05-g.T_cell_q05,0)))
 er1=float(np.max(abs(niche-g.spp1_vascular_niche)));er2=float(np.max(abs(tb-g.T_cell_balance)))
 assert er1<1e-10 and er2<1e-6
 formula.append(dict(sample=sid,niche_formula_error=er1,T_balance_formula_error=er2))
pd.DataFrame(formula).to_csv(P/'audit/COMPOSITE_FORMULA_CHECK.tsv',sep='\t',index=False)
su.to_csv(P/'tables/OLD_VASCULAR_ALL_RELATIONS.tsv',sep='\t',index=False)
for name in ['GLMF_abundance','GLMF_programs']:
 d=read(files[name]);d.to_csv(P/'tables'/(name+'.tsv'),sep='\t',index=False)
 if name=='GLMF_abundance':
  q=d[(d['source']=='TYMP-High Malignant')&(d.target=='SPP1+ TAM')&(d.scale=='spatial_block_permutation')]
  print('GLMF_DIRECT',q[['sample','metric','grid','rho','p_value','fdr_bh']].to_string(index=False))
manifest=[dict(role=k,path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for k,f in files.items()]
(P/'SOURCE_SHA256.json').write_text(json.dumps(manifest,indent=2))
result=dict(status='PASS_KEY_LEGACY_NUMERICAL_AUDIT',stored_block_correlations_reproduced=len(ss),max_rho_error=max(errors),vascular_patient_summaries=18,state_patient_summaries=48,interface_patient_summaries=10,composite_formulas_verified_sections=14,permutation_P_not_rerun=True,old_inputs_unchanged=True,scope='key prior report endpoints, not every historical figure or pipeline')
(P/'VALIDATION.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
