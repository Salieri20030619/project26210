
from pathlib import Path
import h5py,numpy as np,pandas as pd,json,hashlib
from anndata.io import read_elem
P=Path('/mnt/BioDisk/wangmingyuan')
O=P/'project26210_external_joint_v1/phase165_A_input_confirmation_20260915_v1'
S=P/'project26210_external_joint_v1_phase54_gse175540_spatial_v1'
F=P/'project26210_phase78_fine_state_spatial_v1/run_v1/sections_v1'
genes='TYMP ODC1 SRM SMS MTAP AMD1 MAT2A SLC16A3 SPP1 COL1A1 COL1A2 DCN LUM'.split()
rows=[]; cov=[]
for a in sorted((S/'run_v1/sections_v1').glob('*/processed_section_v1.h5ad')):
 gsm=a.parent.name
 with h5py.File(a) as f:
  obs=read_elem(f['obs']); var=read_elem(f['var']); xy=read_elem(f['obsm/spatial']); c=read_elem(f['layers/counts'])
  names=list(var.index); ab=pd.read_csv(F/gsm/'fine_state_abundance_v1.tsv.gz',sep='\t',index_col=0)
  selected=c[:,[names.index(g) for g in genes if g in names]].toarray()
  present=[g for g in genes if g in names]
  for g in genes:
   j=present.index(g) if g in present else None
   cov.append(dict(section=gsm,gene=g,present=j is not None,detected_spots=int((selected[:,j]>0).sum()) if j is not None else 0,total_counts=float(selected[:,j].sum()) if j is not None else 0))
  raw=list((S/'input_extract_v1').rglob(gsm+'*filtered*'))
  label=raw[0].name if raw else ''
  rows.append(dict(section=gsm,source_filename=label,preparation='FFPE' if len(names)==17943 else 'frozen',spots=len(obs),genes=len(names),coordinates_finite=bool(np.isfinite(xy).all()),coordinates_unique=len(np.unique(xy,axis=0))==len(obs),barcodes_unique=obs.index.is_unique,abundance_exact_order=list(obs.index)==list(ab.index),abundance_finite=bool(np.isfinite(ab.to_numpy()).all()),counts_nonnegative_integer=bool((c.data>=0).all() and np.equal(c.data,np.floor(c.data)).all()),total_counts_match=bool(np.allclose(np.asarray(c.sum(axis=1)).ravel(),obs.total_counts)),all_candidates_present=all(g in names for g in genes),input_path=str(a),abundance_path=str(F/gsm/'fine_state_abundance_v1.tsv.gz')))
d=pd.DataFrame(rows);d.to_csv(O/'SECTION_INPUT_QC.tsv',sep='\t',index=False)
pd.DataFrame(cov).to_csv(O/'CANDIDATE_COVERAGE.tsv',sep='\t',index=False)
pil=[]
for prep,x in d.groupby('preparation'):
 x=x.assign(distance=(x.spots-x.spots.median()).abs()).sort_values(['distance','section'])
 pil.append(x.iloc[0].to_dict())
pd.DataFrame(pil).to_csv(O/'PROVISIONAL_TECHNICAL_PILOT.tsv',sep='\t',index=False)
print(d.to_string(index=False,columns=list(d.columns)[:13]))
print('PILOT',[(x['section'],x['preparation']) for x in pil])
print('MIN_COVERAGE',pd.DataFrame(cov).groupby('gene').detected_spots.min().to_dict())

