from common import *
import h5py,re,gzip,io
from PIL import Image
from matplotlib.path import Path as PolyPath
ROOTS=[
'/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial',
'/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/pathology_annotation',
'/mnt/BioDisk/wangmingyuan/project26210_pathology_boundary_ai_draft_v1',
'/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase54_gse175540_spatial_v1']
def main():
 gate('boundary_audit_start')
 imfile=R/'phase229_spatial_F1_partition_eligibility_20260921_v1/tables/IMAGE_COORDINATE_SUPPORT.tsv'
 regfile=R/'phase234_part4_replan_execution_20260921_v1/tables/INPUT_SECTION_REGISTER.tsv'
 im=read(imfile);reg=read(regfile).set_index('sample');images=[];fields=[];manifest=[]
 for r in im.itertuples():
  path=r.image_path;embedded='::' in path
  if embedded:
   fn,key=path.split('::',1)
   with h5py.File(fn,'r') as f:
    node=f[key];shape=node.shape;assert len(shape)==3 and shape[-1] in [3,4];height,width=shape[:2]
   st=Path(fn).stat();identity=dict(path=fn,bytes=st.st_size,mtime_ns=st.st_mtime_ns)
  else:
   with Image.open(io.BytesIO(gzip.open(path,'rb').read()) if path.endswith('.gz') else path) as img:width,height=img.size;shape=(height,width,len(img.getbands()))
   st=Path(path).stat();identity=dict(path=path,bytes=st.st_size,mtime_ns=st.st_mtime_ns)
  d=read(R/'phase234_part4_replan_execution_20260921_v1/inputs'/f'{r.sample}_INPUTS.tsv.gz')
  inside=(d.x_hires_px>=0)&(d.y_hires_px>=0)&(d.x_hires_px<width)&(d.y_hires_px<height)
  image_match=width==r.image_width and height==r.image_height
  images.append(dict(sample=r.sample,dataset=reg.at[r.sample,'dataset'],patient=reg.at[r.sample,'patient'],image_path=path,width=width,height=height,shape=str(shape),registered_dimensions_match=image_match,spots=len(d),coordinate_finite=bool(np.isfinite(d[['x_hires_px','y_hires_px']]).all().all()),spots_inside=int(inside.sum()),fraction_inside=float(inside.mean()),spot_diameter_image_px=float(r.spot_diameter_fullres*r.scale),pixel_scale=float(r.scale),available_for_image_only_morphology=True,approved_independent_tumor_stroma_boundary=False,existing_multiplex_cell_localization=False))
  manifest.append(identity)
  with h5py.File(reg.at[r.sample,'input_path'],'r') as h:
   keys=list(h['obs'].keys())
   obskeys=[k for k in keys if any(x in k.lower() for x in ['annot','patholog','tumor','tumour','region','domain','cluster','leiden'])]
   unskeys=list(h.get('uns',{}).keys())
   fields.append(dict(sample=r.sample,obs_candidate_fields=';'.join(obskeys),uns_top_keys=';'.join(unskeys),interpretation='field_name_only_not_independent_pathology_provenance'))
 tab(images,'IMAGE_COORDINATE_AVAILABILITY.tsv');tab(fields,'SPATIAL_OBJECT_ANNOTATION_FIELDS.tsv')
 # Bounded inventory of annotation-like assets, retain all hits for classification.
 scans=[];candidates=[];polygons=[]
 for root in ROOTS:
  n=0;nh=0
  for dr,dirs,files in os.walk(root):
   dirs[:]=[x for x in dirs if x not in ['.git','__pycache__']]
   for name in files:
    n+=1;low=name.lower();f=Path(dr)/name
    hit=low.endswith(('.geojson','.czi','.nd2','.svs','.mrxs','.qptiff','.tif','.tiff')) or (any(x in low for x in ['patholog','annotation','annotated','tumor_mask','tumour_mask','stroma_mask','qupath','region_label','roi_label']) and low.endswith(('.json','.geojson','.csv','.tsv','.xml','.md','.txt')))
    if not hit:continue
    nh+=1;st=f.stat();rec=dict(path=str(f),bytes=st.st_size,mtime_ns=st.st_mtime_ns,kind='UNRESOLVED_ASSET_NOT_VALIDATED_BOUNDARY')
    if 'pathology_polygons' in low and low.endswith('.json'):
     q=json.loads(f.read_text());digest=sha(f);rec.update(kind='AI_TUMOR_NAT_POLYGONS',sha256=digest,review=q.get('review_status',''),annotation_source=q.get('annotation_source',''))
     sample=q['sample'];d=read(R/'phase234_part4_replan_execution_20260921_v1/inputs'/f'{sample}_INPUTS.tsv.gz');points=d[['x_hires_px','y_hires_px']].to_numpy();assigned=[]
     for j,poly in enumerate(q['polygons']):
      v=np.array(poly['vertices'],float);inside_poly=PolyPath(v).contains_points(points);assigned.append(inside_poly)
      area=.5*abs(np.dot(v[:,0],np.roll(v[:,1],-1))-np.dot(v[:,1],np.roll(v[:,0],-1)))
      valid_bbox=bool((v[:,0]>=0).all() and (v[:,0]<=q['image_width']).all() and (v[:,1]>=0).all() and (v[:,1]<=q['image_height']).all())
      polygons.append(dict(path=str(f),sample=sample,sha256=digest,polygon=j,label=poly['label'],vertices=len(v),area_px2=float(area),vertices_in_image=valid_bbox,finite_vertices=bool(np.isfinite(v).all()),spots_inside=int(inside_poly.sum()),coordinate_space=q['coordinate_space'],image_dimensions_match=bool(q['image_width']==next(z['width'] for z in images if z['sample']==sample) and q['image_height']==next(z['height'] for z in images if z['sample']==sample)),review_status=q.get('review_status',''),prior_RNA_reference='computational boundary' in q.get('annotation_source','').lower(),boundary_type='tumor_vs_NAT_not_intratumoral_tumor_stroma',eligible_as_independent_gold_standard=False))
     rec['spots_in_overlapping_polygons']=int((np.column_stack(assigned).sum(1)>1).sum())
    elif low.endswith(('.tif','.tiff','.qptiff')):
     try:
      import tifffile
      with tifffile.TiffFile(f) as t:rec.update(kind='RASTER_IMAGE_NOT_BOUNDARY',details=json.dumps(dict(pages=len(t.pages),series=[dict(shape=s.shape,axes=s.axes) for s in t.series],ome=bool(t.ome_metadata))))
     except Exception as ex:rec['details']=str(ex)
    elif low.endswith(('.csv','.tsv')):
     try:
      sep=',' if low.endswith('.csv') else '\t';head=pd.read_csv(f,sep=sep,nrows=0);rec['columns']=';'.join(head.columns.astype(str));rec['kind']='TABLE_PROVENANCE_REQUIRED'
     except Exception as ex:rec['details']=str(ex)
    elif low.endswith(('.md','.txt')):rec['kind']='PROVENANCE_DOCUMENT'
    candidates.append(rec)
  scans.append(dict(root=root,exists=Path(root).exists(),files_scanned=n,annotation_like_assets=nh))
 tab(scans,'BOUNDARY_SCAN_SCOPE.tsv');tab(candidates,'BOUNDARY_CANDIDATE_INVENTORY.tsv');poly=pd.DataFrame(polygons);tab(poly,'POLYGON_GEOMETRY_AND_PROVENANCE.tsv')
 # Canonical polygon identities; repeated copies are not independent annotations.
 unique=poly[['sample','sha256','review_status','prior_RNA_reference','boundary_type','eligible_as_independent_gold_standard']].drop_duplicates();tab(unique,'UNIQUE_BOUNDARY_ANNOTATIONS.tsv')
 assert len(images)==38 and len(unique)==2
 # Preserve exact provenance records used in audit, no mutation.
 for sample in ['GLMF1','GLMF2']:
  f=Path(ROOTS[2])/'annotations'/f'{sample}_pathology_polygons.json'
  (P/'resources'/f'{sample}_existing_polygon.json').write_bytes(f.read_bytes())
 note=Path(ROOTS[2])/'README_AI_ASSISTED_BOUNDARY_DRAFT.md';(P/'resources/EXISTING_BOUNDARY_PROVENANCE.md').write_bytes(note.read_bytes())
 manifest += [dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)) for f in [imfile,regfile,note]+list((Path(ROOTS[2])/'annotations').glob('*.json'))]
 dump(P/'qc/BOUNDARY_INPUT_MANIFEST.json',manifest)
 dump(P/'qc/BOUNDARY_AUDIT_COMPLETE.json',dict(status='COMPUTED',images=len(images),images_dimensions_match=sum(x['registered_dimensions_match'] for x in images),all_coordinates_finite=all(x['coordinate_finite'] for x in images),unique_preexisting_polygons=len(unique),independent_validated_tumor_stroma_boundaries=0,scan_roots=ROOTS,limitations='bounded existing project assets; not a global absence claim',boundary_mutation=False,new_images_generated=False))
 print('BOUNDARY_AUDIT_COMPLETE',len(images),len(candidates),len(unique),flush=True)
if __name__=='__main__':main()
