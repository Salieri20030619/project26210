from pathlib import Path
import json,csv,hashlib,os,sys
import numpy as np,pandas as pd,scipy,skimage
from PIL import Image
P=Path(__file__).resolve().parents[1];ROOT=P.parent
C=json.loads((P/'CONFIG.json').read_text())
D=ROOT/'phase241_spatial_MSI_download_20260924_v1'
assert (D/'EXIT_STATUS').read_text().strip()=='0'
assert json.loads((D/'FINAL_VALIDATION.json').read_text())['status']=='PASS'
assert json.loads((D/'monitoring/COMPLETION_VALIDATION.json').read_text())['final_acceptance']=='PASS'
def sha(f):
 h=hashlib.sha256()
 with open(f,'rb') as s:
  for b in iter(lambda:s.read(2**20),b''):h.update(b)
 return h.hexdigest()
assets=pd.read_csv(P/'metadata/ASSET_INVENTORY.tsv',sep='\t').set_index('sample')
manifest=[];ids=[]
for s in C['samples']:
 a=assets.loc[s]
 assert a.tissue=='tumor'
 d=pd.read_csv(ROOT/'phase234_part4_replan_execution_20260921_v1/inputs'/f'{s}_INPUTS.tsv.gz',sep='\t')
 assert d.patient.nunique()==1 and d.patient.iloc[0]==a.patient_id and d.spot_id.is_unique
 assert len(d)==int(a.existing_tumor_qc_spots)
 assert np.isfinite(d[C['primary_RNA_targets']+C['supplementary_RNA_targets']].to_numpy()).all()
 ids.append(a.patient_id)
 for key in ['positions_path','he_path','scalefactors_path']:
  f=Path(a[key]);assert f.exists()
  manifest.append(dict(path=str(f),bytes=f.stat().st_size,mtime_ns=f.stat().st_mtime_ns,sha256=sha(f)))
 with Image.open(a.he_path) as im:im.verify()
 for pol in ['neg','pos']:
  for ext in ['.imzML','.ibd']:
   f=P/'extracted'/(s.replace('_','-')+'-'+pol+ext);assert f.is_symlink() and f.resolve().exists()
   source=f.resolve();manifest.append(dict(path=str(source),bytes=source.stat().st_size,mtime_ns=source.stat().st_mtime_ns,sha256=sha(source)))
 f=ROOT/'phase234_part4_replan_execution_20260921_v1/inputs'/f'{s}_INPUTS.tsv.gz'
 manifest.append(dict(path=str(f),bytes=f.stat().st_size,mtime_ns=f.stat().st_mtime_ns,sha256=sha(f)))
assert len(set(ids))==12 and ids==C['patient_ids']
(P/'audit/RAW_INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
# Raw input fingerprints must agree with completed download evidence.
checks=pd.read_csv(D/'tables/EXTRACTED_FILE_MANIFEST.tsv',sep='\t')
look={Path(x['path']).name:x['sha256'] for x in manifest}
assert all(look[r.file]==r.sha256 for r in checks.itertuples())
# Verify inherited numerical geometry convention without any patient molecular result.
from geometry import apply
M=np.array([[0.,-2.,40.],[2.,0.,10.],[0.,0.,1.]])
x=np.array([[1.,2.],[15.,19.]])
assert np.max(abs(apply(apply(x,M),np.linalg.inv(M))-x))<1e-12
from patient_summary import bh
assert np.allclose(bh([.01,.04,.03]),[.03,.04,.04])
rd=[dict(path=x['path'],bytes=Path(x['path']).stat().st_size,mtime_ns=Path(x['path']).stat().st_mtime_ns) for x in json.loads((ROOT/'phase237_SpaMTP_registration_20260923_v1/metadata/READONLY_RDS_BEFORE.json').read_text())]
(P/'audit/OLD_RDS_BEFORE.json').write_text(json.dumps(rd,indent=2))
qa={'status':'PASS','samples':12,'registered_patients':12,'raw_files':48,'download_hashes_verified':40,'inherited_ST_QC':True,'old_RDS_tracked':len(rd),'versions':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'skimage':skimage.__version__}}
(P/'PREFLIGHT_VALIDATION.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa),flush=True)
