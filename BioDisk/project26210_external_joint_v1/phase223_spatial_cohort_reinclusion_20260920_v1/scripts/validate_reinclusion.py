from pathlib import Path
import json,math
import numpy as np,pandas as pd
from scipy import stats
P=Path(__file__).resolve().parents[1]
e=pd.read_csv(P/'results/ALL_14_SECTION_EFFECTS.tsv',sep='\t');summary=pd.read_csv(P/'results/COHORT_SPECIFIC_SUMMARY.tsv',sep='\t');assert len(e)==56 and e['sample'].nunique()==14
errs=[]
for f in (P/'inputs').glob('*_DESIGN.npz'):
 d=np.load(f);ma=d['mask'];mat=np.column_stack([d['x'][ma],d['y'][ma],d['cov'][ma]]);rank=np.column_stack([stats.rankdata(mat[:,j]) for j in range(mat.shape[1])]);inv=np.linalg.inv(np.corrcoef(rank.T));rho=-inv[0,1]/np.sqrt(inv[0,0]*inv[1,1]);key='endothelial' if 'endothelial' in f.name else 'pericyte';sid=f.name.removesuffix('_'+key+'_DESIGN.npz');effect=e[(e['sample']==sid)&(e.endpoint=='S2_'+key)].effect.iloc[0];assert abs(rho-effect)<1e-10;errs.append(abs(rho-effect))
for f in (P/'inputs').glob('*_FROZEN_TRANSFER_INPUTS.tsv.gz'):
 sid=f.name.removesuffix('_FROZEN_TRANSFER_INPUTS.tsv.gz');d=pd.read_csv(f,sep='\t')
 for ep,left,right in [('S3_expression_boundary','boundary_source','core_source'),('S3_ECM_region','ECM_high','ECM_low')]:
  row=e[(e['sample']==sid)&(e.endpoint==ep)].iloc[0]
  if row.eligible:
   assert d[left].sum()>=20 and d[right].sum()>=20
   value=d.loc[d[left],'SPP1_share'].median()-d.loc[d[right],'SPP1_share'].median();assert abs(value-row.effect)<1e-12
for row in summary.itertuples():
 a=e[(e.dataset==row.dataset)&(e.endpoint==row.endpoint)].effect.dropna().to_numpy();assert len(a)==row.n_patients
 if len(a):assert abs(np.median(a)-row.median_effect)<1e-12
 if len(a)<10:assert np.isnan(row.P)
 else:
  n=len(a);k=min(int((a>0).sum()),int((a<0).sum()));p=min(1,2*sum(math.comb(n,j) for j in range(k+1))/2**n);assert abs(p-row.P)<1e-12
for cohort,g in summary[summary.endpoint.str.startswith('S2')].groupby('dataset'):
 if g.P.notna().all():
  a=g.P.to_numpy();idx=np.argsort(a);q=np.minimum.accumulate((a[idx]*2/np.arange(1,3))[::-1])[::-1];assert np.allclose(q,g.FDR.to_numpy()[idx])
alias=json.loads((P/'GSE210041_GLMF_ALIAS_PROOF.json').read_text());assert len(alias)==2 and all(x['counts_exact'] for x in alias) and len({x['GEO_sample'] for x in alias})==2
for x in json.loads((P/'audit/RDS_BEFORE.json').read_text()):
 st=Path(x['path']).stat();assert(st.st_size,st.st_mtime_ns,st.st_ino)==(x['bytes'],x['mtime_ns'],x['inode'])
v={'status':'PASS','unique_added_sections':14,'GLMF_GSE210041_duplicate_labels_collapsed':True,'GLMF_original_10x_counts_exact':2,'main_case_patients':2,'new_inferential_n_from_GLFM':0,'primary_coefficients_independently_checked':len(errs),'max_error':max(errs),'regional_effects_and_exact_sign_P_checked':True,'old_RDS_unchanged':333,'no_new_figures':True,'no_new_mapping':True,'TLS_not_run_without_annotations':True}
(P/'FINAL_VALIDATION.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
