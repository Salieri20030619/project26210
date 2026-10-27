from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd
from scipy.stats import norm
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase157_ordered_bridge_reassessment_20260915_v1';D=O/'01_public_TYMP_MTA';D.mkdir(exist_ok=False)
A=Path('/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1')
co=pd.read_csv(A/'analysis_v2/cohort_associations_v2.tsv',sep='\t');meta=pd.read_csv(A/'analysis_v2/cross_cohort_meta_analysis_v2.tsv',sep='\t');checks=[]
for ep in ['A_B_TYMP_MTA_MARGINAL','MTAP_B_MTAP_MTA_CONTEXT']:
 x=co[co.endpoint==ep];w=x.n-3;zz=np.arctanh(x.rho);z=(w*zz).sum()/w.sum();r=np.tanh(z);pv=2*norm.sf(abs(z)*np.sqrt(w.sum()));ref=meta[meta.endpoint==ep].iloc[0]
 checks.append(dict(endpoint=ep,pass_check=bool(np.isclose(r,ref.fixed_rho) and np.isclose(pv,ref.fixed_p)),rho=float(r),P=float(pv)))
assert all(x['pass_check'] for x in checks)
for f in [A/'analysis_v2/cohort_associations_v2.tsv',A/'analysis_v2/cross_cohort_meta_analysis_v2.tsv',A/'analysis_v2/paired_sample_counts_v2.tsv',A/'analysis_validation_v4/validation_status_v4.json',A/'analysis_phase99_v1/meta_tymp_mtap_associations_v1.tsv']:
 shutil.copy2(f,D/f.name)
(D/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE_PUBLIC_MEASURED_MTA_REASSESSMENT',checks=checks,paired_tumours=255,new_C4_available_in_these_samples=False,client_data_used=False,decision='Pooled TYMP-MTA not supported; preserve CPTAC nominal positive; old SPP1 endpoints not C4.'),indent=2))
B=Path('/mnt/BioDisk/wangmingyuan/project26210_phase77_mta_perturbation_signature_v1')
s=pd.read_csv(B/'training_v1/phase77_weighted_signatures_v1.tsv',sep='\t');direct=s[s.signature_id=='DIRECT_MTA_GSE249082']
print('STEP1',checks)
print('DIRECT_SIGNATURE',len(direct),direct.direction.value_counts().to_dict(),direct.groupby('direction').weight.agg(['min','max']).to_dict())
for f in [B/'03_run_phase77_singlecell_transfer_v1.R']:
 lines=f.read_text().splitlines()
 print('SCORE_CODE','\n'.join(line for line in lines if any(k in line.lower() for k in ['weight','score','zscore','scale(']))[:6000])
lr=pd.read_csv(P/'phase35_part1_c01_c04_cell_communication_v1/resource_lock_v1/cellchat_interactions_v1.tsv',sep='\t',low_memory=False)
print('LR_SAMPLE',lr.loc[lr.ligand.isin(['TGFB1','IL6','MIF','CSF1']),['ligand','ligand.symbol','receptor','receptor.symbol','annotation']].head(8).to_dict('records'))
print('RSCRIPT',[str(p) for p in [Path('/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'),Path('/mnt/BioDisk/ludingyu/project26210/.venv-pyscenic/bin/python')] if p.exists()])
