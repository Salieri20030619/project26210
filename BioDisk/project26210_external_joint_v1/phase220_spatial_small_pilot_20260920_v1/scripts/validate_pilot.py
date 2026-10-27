from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.optimize import lsq_linear
P=Path(__file__).resolve().parents[1];R=P.parent
assert (P/'EXIT_STATUS').read_text().strip()=='0'
r=pd.read_csv(P/'simulation/ALL_RECOVERY_RESULTS.tsv.gz',sep='\t');truth=pd.read_csv(P/'simulation/MIXTURE_TRUTH_AND_SEEDS.tsv.gz',sep='\t')
assert len(truth)==26*5*3*3 and len(r)==len(truth)*3 and truth.mixture_id.is_unique
assert r.groupby('mixture_id').size().eq(3).all() and r.within_macro_abs_error.notna().all()
assert np.isfinite(r[['estimate_C4','C4_abs_error','within_macro_abs_error']]).all().all()
assert np.allclose(abs(r.estimate_C4-r.true_C4),r.C4_abs_error)
assert all(abs(sum(json.loads(x).values())-1)<1e-10 for x in r.estimated_fractions)
meta=pd.read_csv(R/'phase219_spatial_reference_identifiability_20260920_v1/reference/PATIENT_STATE_METADATA.tsv',sep='\t')
for x in truth.itertuples():
 sources=json.loads(x.sources);fs=json.loads(x.true_fractions);assert abs(sum(fs.values())-1)<1e-10
 assert meta.iloc[x.C4_profile_row].patient==x.C4_patient and meta.iloc[x.C4_profile_row].study==x.study
 assert all(meta.iloc[i].study==x.study and meta.iloc[i].label==label for label,i in sources.items())
 assert abs(fs['C4']+fs['Other_mature_macrophage']-.3)<1e-12
checks=[]
for f in sorted((P/'spatial').glob('*/*_FIT_CACHE.npz')):
 z=np.load(f);A=z['A']*z['weights'][:,None];Y=z['y']*z['weights'];C=z['coef'];E=C@A.T-Y;grad=E@A
 scale=np.maximum(1,np.max(abs(Y@A),axis=1))
 dual=float(np.max(np.maximum(-grad,0)/scale[:,None]));complement=float(np.max(abs(grad*C)/np.maximum(1,np.sum(Y*Y,axis=1))[:,None]))
 assert (C>=0).all() and dual<1e-7 and complement<1e-7,(str(f),dual,complement)
 # A distinct bound-constrained solver validates one deterministic spot per map.
 j=len(Y)//2;fit=lsq_linear(A,Y[j],bounds=(0,np.inf),method='bvls',tol=1e-10,max_iter=500)
 rel=abs(np.linalg.norm(A@fit.x-Y[j])-np.linalg.norm(E[j]))/max(1,np.linalg.norm(Y[j]));assert rel<1e-7
 checks.append(dict(file=str(f.relative_to(P)),spots=len(Y),dual_violation=dual,complementarity=complement,alternate_solver_relative_objective_difference=rel))
old=json.loads((P/'audit/PREEXISTING_RDS_STAT.json').read_text())
for x in old:
 s=Path(x['path']).stat();assert (s.st_size,s.st_mtime_ns,s.st_ino)==(x['bytes'],x['mtime_ns'],x['inode'])
# Check inherited frozen identity files and source profile files were not changed by this pilot.
sourcehash=json.loads((R/'phase218_spatial_identity_alignment_20260920_v1/SOURCE_SHA256.json').read_text());assert all(hashlib.sha256(Path(f).read_bytes()).hexdigest()==h for f,h in sourcehash.items())
result={'status':'PASS','mixtures':len(truth),'fits':len(r),'unique_C4_donor_profiles':truth.C4_profile_row.nunique(),'heldout_truth_study_integrity':True,'two_section_three_panel_nnls_checks':checks,'preexisting_RDS_unchanged':len(old),'frozen_identity_hashes_unchanged':True,'validation_is_numerical_not_biological':True}
(P/'FINAL_VALIDATION.json').write_text(json.dumps(result,indent=2));print('VALIDATION_PASS',len(truth),len(r),len(checks),len(old))
