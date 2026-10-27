from pathlib import Path
import json,hashlib,datetime,shutil
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent;A=R/'phase223_spatial_cohort_reinclusion_20260920_v1';B=R/'phase222_spatial_S2_S5_execution_20260920_v1'
assert (P/'EXIT_STATUS').read_text().strip()=='0'
for sid,n in [('GLMF1',3199),('GLMF2',3880)]:
 d=pd.read_csv(P/'tables'/(sid+'_DISPLAY_INPUTS.tsv.gz'),sep='\t')
 old=pd.read_csv(A/'inputs'/(sid+'_FROZEN_TRANSFER_INPUTS.tsv.gz'),sep='\t')
 assert len(d)==n and d.spot_id.is_unique
 for c in old.columns:
  if old[c].dtype.kind in 'f':np.testing.assert_allclose(d[c],old[c],rtol=1e-13,atol=1e-13)
  else:assert d[c].equals(old[c]),c
 assert d.SPP1_share_display_eligible.notna().equals(d.source_eligible)
 assert not any('draft' in c.lower() or 'pathology' in c.lower() for c in d)
 assert np.isfinite(d.filter(like='display_program__').to_numpy()).all()
f=pd.read_csv(P/'tables/GLMF_CASE_EFFECTS.tsv',sep='\t')
orig=pd.read_csv(A/'results/ALL_14_SECTION_EFFECTS.tsv',sep='\t')
pd.testing.assert_frame_equal(f,orig[orig.dataset=='GLMF_GSE210041'].reset_index(drop=True))
e=pd.read_csv(P/'tables/EXISTING_PATIENT_EVIDENCE.tsv',sep='\t')
for dataset,source in [('TJRCC',A/'results/COHORT_SPECIFIC_SUMMARY.tsv'),('GSE175540',B/'PATIENT_RESULT_SUMMARY.tsv')]:
 o=pd.read_csv(source,sep='\t')
 if dataset=='TJRCC':o=o[o.dataset==dataset]
 else:o=o[o.model=='primary']
 for row in e[e.dataset==dataset].itertuples():
  base=o[o.endpoint==row.endpoint].iloc[0]
  for col in ['n_patients','positive','negative','median_effect','P','FDR']:np.testing.assert_allclose(getattr(row,col),base[col],equal_nan=True)
assert len(e)==9 and len(f)==8
before=json.loads((P/'audit/RDS_BEFORE.json').read_text())
for x in before:
 s=Path(x['path']).stat();assert (s.st_size,s.st_mtime_ns,s.st_ino)==(x['bytes'],x['mtime_ns'],x['inode'])
for x in json.loads((P/'audit/SOURCE_MANIFEST.json').read_text()):
 s=Path(x['path']).stat();assert (s.st_size,s.st_mtime_ns,s.st_ino)==(x['bytes'],x['mtime_ns'],x['inode'])
assert not any(f.suffix.lower() in ['.png','.pdf','.svg','.rds','.h5ad'] for f in P.rglob('*') if f.is_file())
result=dict(status='PASS',display_spots=7079,GLMF_sections=2,all_frozen_columns_preserved=True,case_effect_rows=8,patient_summary_rows=9,old_RDS_checked=len(before),old_RDS_changed=0,source_metadata_unchanged=True,new_figures=0,new_models=0,new_inferential_tests=0,original_images_visually_reviewed='lowres originals; not pathology diagnosis')
(P/'FINAL_VALIDATION.json').write_text(json.dumps(result,indent=2))
(P/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE_P1_P2_P3_PREPARATION',requires_content_review_before_new_figures=True,optional_14section_robustness_not_run=True,download_continues_independently=True),indent=2))
print(json.dumps(result))
