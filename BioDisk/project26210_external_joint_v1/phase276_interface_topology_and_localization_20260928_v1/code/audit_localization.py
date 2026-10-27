from pathlib import Path
import os,json,datetime,hashlib
import numpy as np,pandas as pd,h5py
P=Path(__file__).resolve().parents[1];R=P.parent
roots=['/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial','/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/pathology_annotation','/mnt/BioDisk/wangmingyuan/project26210_pathology_boundary_ai_draft_v1','/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase54_gse175540_spatial_v1']
scans=[];candidate=[]
for root in roots:
 count=0;found=0
 for directory,dirs,files in os.walk(root):
  dirs[:]=[x for x in dirs if x not in ['.git','__pycache__']]
  for name in files:
   count+=1;low=name.lower();f=Path(directory)/name
   hit=low.endswith(('.tif','.tiff','.czi','.nd2','.svs','.mrxs','.qptiff','.geojson')) or any(x in low for x in ['multiplex','immunofluor','qupath','phenochart','cell_seg','cell_coord','pathology_polygons'])
   if not hit:continue
   found+=1;record=dict(path=str(f),bytes=f.stat().st_size,mtime_ns=f.stat().st_mtime_ns,status='CANDIDATE_NOT_VALIDATED',details='')
   if 'pathology_polygons' in low:
    a=json.loads(f.read_text());record.update(status=a.get('review_status','UNKNOWN'),details=a.get('annotation_source','UNKNOWN'),sha256=hashlib.sha256(f.read_bytes()).hexdigest())
   elif low.endswith(('.tif','.tiff','.qptiff')):
    try:
     import tifffile
     with tifffile.TiffFile(f) as tf:
      record['details']=json.dumps(dict(pages=len(tf.pages),series=[dict(shape=x.shape,axes=x.axes) for x in tf.series],ome=bool(tf.ome_metadata)))
      record['status']='IMAGE_METADATA_REQUIRES_ASSAY_PROVENANCE'
    except Exception as ex:record['details']=str(ex)
   candidate.append(record)
 scans.append(dict(root=root,files_scanned=count,candidates=found))
pd.DataFrame(scans).to_csv(P/'tables/LOCALIZATION_SCAN_SCOPE.tsv',sep='\t',index=False)
pd.DataFrame(candidate,columns=['path','bytes','mtime_ns','status','details','sha256']).to_csv(P/'tables/LOCALIZATION_CANDIDATES.tsv',sep='\t',index=False)
im=pd.read_csv(R/'phase229_spatial_F1_partition_eligibility_20260921_v1/tables/IMAGE_COORDINATE_SUPPORT.tsv',sep='\t')
rows=[]
for r in im.itertuples():
 path=r.image_path;status='';shape=''
 if '::' in path:
  filename,entry=path.split('::')
  with h5py.File(filename,'r') as f:
   assert entry in f;shape=str(f[entry].shape);assert len(f[entry].shape)==3 and f[entry].shape[2] in [3,4]
  status='H&E_RGB_SUPPORT_ONLY'
 else:
  assert Path(path).exists();status='H&E_ASSET_FROM_FROZEN_REGISTER'
 rows.append(dict(sample=r.sample,image_path=path,width=r.image_width,height=r.image_height,embedded_shape=shape,status=status,independent_tumor_boundary=False,SPP1_macrophage_CD8_multiplex_eligible=False))
pd.DataFrame(rows).to_csv(P/'tables/IMAGE_AND_LOCALIZATION_READINESS.tsv',sep='\t',index=False)
print('SCANS',scans);print('CANDIDATES',json.dumps(candidate,ensure_ascii=False,indent=2));print('HE_SECTIONS',len(rows))
(P/'qc/LOCALIZATION_AUDIT.json').write_text(json.dumps(dict(at=datetime.datetime.now().astimezone().isoformat(),scope=roots,HE_sections=len(rows),approved_pathology_boundaries=0,usable_SPP1_macrophage_CD8_multiplex_datasets=0,limitation='bounded local project/public-source inventory; not an assertion that no such public dataset exists worldwide',no_new_download=True,no_customer_wetlab_reanalysis=True),indent=2))
