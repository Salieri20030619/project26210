from pathlib import Path
import json,hashlib,platform,sys,traceback,shutil
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse
from PIL import Image
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent
A=R/'phase223_spatial_cohort_reinclusion_20260920_v1'
B=R/'phase222_spatial_S2_S5_execution_20260920_v1'
def save(n,o): (P/n).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda v:v.item() if hasattr(v,'item') else str(v)))
def tab(n,d): d.to_csv(P/n,sep='\t',index=False)
def sha(f):
 h=hashlib.sha256()
 with Path(f).open('rb') as z:
  for b in iter(lambda:z.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def meta(f):
 f=Path(f);s=f.stat();return dict(path=str(f),bytes=s.st_size,mtime_ns=s.st_mtime_ns,inode=s.st_ino)
def run():
 sources=[];checks=[];summ=[];components=[];regions=[];ranges=[];defs=json.loads((R/'phase218_spatial_identity_alignment_20260920_v1/PROGRAM_DEFINITIONS.json').read_text())
 tracks=[meta(x['path']) for x in json.loads((A/'audit/RDS_BEFORE.json').read_text())]
 save('audit/RDS_BEFORE.json',tracks)
 selected={k:defs[k] for k in ['spp1_state_no_spp1','vascular_remodeling']}
 selected.update({k+'_no_LGALS3':[g for g in v if g!='LGALS3'] for k,v in list(selected.items())})
 save('inputs/DISPLAY_PROGRAM_DEFINITIONS.json',dict(programs=selected,normalization='log1p(raw_counts/spot_total*10000); gene-wise population z within all QC spots; unweighted gene mean',purpose='descriptive display transfer only; not exact reproduction of historical GLMF scores; no new correlations or P values'))
 for item in json.loads((A/'ASSET_ROLE_MANIFEST.json').read_text())[:2]:
  sid=item['sample'];rawdir=Path(item['raw_source']).parent/sid/'spatial'
  fp=A/'inputs'/(sid+'_FROZEN_TRANSFER_INPUTS.tsv.gz')
  d=pd.read_csv(fp,sep='\t')
  with h5py.File(item['input_path']) as h:
   obs=read_elem(h['obs']);var=read_elem(h['var']);X=read_elem(h['X']);xy=read_elem(h['obsm']['spatial'])
   embedded=h['uns/spatial'][sid]['images/hires'][:]
  assert np.array_equal(d.spot_id,obs.index)
  scale=json.loads((rawdir/'scalefactors_json.json').read_text());factor=scale['tissue_hires_scalef']
  pos=pd.read_csv(rawdir/'tissue_positions_list.csv',header=None,names=['spot_id','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('spot_id')
  aligned=pos.loc[obs.index]
  assert np.array_equal(xy,aligned[['pixel_col','pixel_row']].to_numpy())
  assert np.array_equal(aligned[['array_row','array_col']].to_numpy(),obs[['array_row','array_col']].to_numpy())
  assert np.all(aligned.in_tissue==1)
  im=np.array(Image.open(rawdir/'tissue_hires_image.png').convert('RGB'))
  ei=np.rint(embedded*255).astype(np.uint8) if embedded.max()<=1 else embedded.astype(np.uint8)
  delta=int(np.abs(im.astype(int)-ei.astype(int)).max());assert delta<=1
  height,width=im.shape[:2];p=xy*factor
  assert ((p[:,0]>=0)&(p[:,0]<width)&(p[:,1]>=0)&(p[:,1]<height)).all()
  q=pd.read_csv(item['abundance_path'],sep='\t',index_col=0)
  assert q.index.equals(obs.index)
  q.columns=[c.removeprefix('q05cell_abundance_w_sf_') for c in q.columns]
  q=q.div(q.sum(axis=1),axis=0)
  macro=q[['SPP1-state Macrophage','C1QC-state Macrophage','Other Macrophage']].sum(axis=1)
  np.testing.assert_allclose(d.SPP1_share,q['SPP1-state Macrophage']/macro,rtol=1e-12,atol=1e-12)
  d.insert(0,'sample',sid);d['x_hires_px']=p[:,0];d['y_hires_px']=p[:,1]
  d['array_row']=aligned.array_row.to_numpy();d['array_col']=aligned.array_col.to_numpy()
  d['SPP1_state_total_fraction']=q['SPP1-state Macrophage'].to_numpy()
  d['SPP1_share_display_eligible']=d.SPP1_share.where(d.source_eligible)
  # Descriptive score transfer with frozen gene sets; no outcome-driven selection.
  gs=sorted(set(sum(selected.values(),[])));gi=var.index.get_indexer(gs);assert (gi>=0).all()
  totals=np.asarray(X.sum(axis=1)).ravel();assert (totals>0).all()
  lx=np.log1p(X[:,gi].toarray().astype(float)/totals[:,None]*1e4)
  sd=lx.std(axis=0);Z=(lx-lx.mean(axis=0))/np.where(sd==0,1,sd)
  for key,genes in selected.items():d['display_program__'+key]=Z[:,[gs.index(g) for g in genes]].mean(axis=1)
  G=sparse.load_npz(A/'inputs'/(sid+'_adjacency.npz'))
  rows,cols=G.nonzero();length=np.linalg.norm(xy[rows]-xy[cols],axis=1)
  spacing=float(np.median(length));assert length.max()<=1.3*spacing+1e-6
  # Keep all components; do not select favorable tissue fragments.
  for comp,g in d.groupby('component'):
   components.append(dict(sample=sid,component=int(comp),spots=len(g),eligible_source_spots=int(g.source_eligible.sum()),x_min=g.x_hires_px.min(),x_max=g.x_hires_px.max(),y_min=g.y_hires_px.min(),y_max=g.y_hires_px.max()))
  path=Path('/mnt/BioDisk/wangmingyuan/project26210_pathology_boundary_ai_draft_v1')
  annotation=path/'annotations'/(sid+'_pathology_polygons.json')
  aj=json.loads(annotation.read_text());assert aj['review_status']=='not_pathologist_reviewed'
  for col in ['core_source','boundary_source','ECM_low','ECM_high']:
   regions.append(dict(sample=sid,region=col,definition='Phase223 frozen expression/composition definition',spots=int(d[col].sum()),pathology_confirmed=False))
  af=path/'tables'/(sid+'_spot_pathology_signed_distance.csv')
  ad=pd.read_csv(af).set_index('spot_id').loc[obs.index]
  np.testing.assert_allclose(ad[['x_image_px','y_image_px']].to_numpy(),p,atol=1e-5)
  tab('audit/'+sid+'_AI_DRAFT_LABELS_ONLY.tsv',pd.DataFrame(dict(spot_id=obs.index,draft_label=ad.pathology_label.to_numpy(),status='NOT_PATHOLOGIST_REVIEWED')))
  # No historical draft abundance columns imported: they use a different reference.
  tab('tables/'+sid+'_DISPLAY_INPUTS.tsv.gz',d)
  diam=scale['spot_diameter_fullres'];um_per_px=55/diam/factor
  out=dict(sample=sid,GEO_sample=item['GEO_sample'],patient=item['patient'],QC_spots=len(d),in_tissue_positions=int(pos.in_tissue.sum()),all_array_positions=len(pos),source_eligible=int(d.source_eligible.sum()),components=int(d.component.nunique()),image_width_px=width,image_height_px=height,image_array_max_delta=delta,spot_diameter_hires_px=diam*factor,um_per_hires_px=um_per_px,one_mm_hires_px=1000/um_per_px,scale_basis='nominal Visium 55um spot diameter; not independent scanner calibration',coordinate_orientation='x=column rightward; y=row downward; origin upper-left; no flip or rotation',TLS='NO_AUDITED_ANNOTATION',pathology='AI_DRAFT_ONLY',ai_draft_counts=ad.pathology_label.value_counts().to_dict())
  summ.append(out)
  checks.append(dict(sample=sid,barcode_exact=True,coordinates_exact=True,image_matches_embedded=True,all_spots_inside_image=True,frozen_share_exact=True,all_components_retained=True,program_genes_complete=True))
  for f in [fp,Path(item['input_path']),Path(item['posterior_path']),Path(item['abundance_path']),rawdir/'tissue_hires_image.png',rawdir/'tissue_lowres_image.png',rawdir/'tissue_positions_list.csv',rawdir/'scalefactors_json.json',annotation,af]:
   rec=meta(f)
   if f.suffix not in ['.h5ad']:rec['sha256']=sha(f)
   sources.append(rec)
  print(sid,out,flush=True)
 tab('tables/COMPONENT_COVERAGE.tsv',pd.DataFrame(components));tab('tables/REGION_BASIS.tsv',pd.DataFrame(regions))
 save('tables/IMAGE_COORDINATE_COVERAGE.json',summ);save('audit/DISPLAY_INPUT_CHECKS.json',checks)
 alltab=pd.concat([pd.read_csv(P/'tables'/(s+'_DISPLAY_INPUTS.tsv.gz'),sep='\t') for s in ['GLMF1','GLMF2']])
 for name in ['SPP1_share_display_eligible','SPP1_state_total_fraction','total_macro_fraction','endothelial_fraction','pericyte_fraction']+[c for c in alltab if c.startswith('display_program__')]:
  a=alltab[name].dropna();ranges.append(dict(variable=name,finite_spots=len(a),minimum=a.min(),maximum=a.max(),display_min=0 if not name.startswith('display_program__') else -float(a.abs().max()),display_max=1 if not name.startswith('display_program__') else float(a.abs().max()),rule='fixed 0-1 for fractions; symmetric pooled full range for within-section program z scores; no clipping'))
 tab('tables/PROPOSED_COLOR_LIMITS.tsv',pd.DataFrame(ranges))
 legacy=pd.read_csv(A/'results/COHORT_SPECIFIC_SUMMARY.tsv',sep='\t')
 legacy=legacy[legacy.dataset=='TJRCC'].copy();legacy['model']='primary';legacy['source_file']=str(A/'results/COHORT_SPECIFIC_SUMMARY.tsv')
 large=pd.read_csv(B/'PATIENT_RESULT_SUMMARY.tsv',sep='\t');large=large[large.model=='primary'].copy();large['dataset']='GSE175540';large['source_file']=str(B/'PATIENT_RESULT_SUMMARY.tsv')
 tab('tables/EXISTING_PATIENT_EVIDENCE.tsv',pd.concat([legacy,large],ignore_index=True))
 ef=pd.read_csv(A/'results/ALL_14_SECTION_EFFECTS.tsv',sep='\t');tab('tables/GLMF_CASE_EFFECTS.tsv',ef[ef.dataset=='GLMF_GSE210041'])
 shutil.copy2(R/'phase219_spatial_reference_identifiability_20260920_v1/program_check/PROGRAM_CHECK_SUMMARY.tsv',P/'tables/EXISTING_SHARED_GENE_CHECK.tsv')
 for f in [A/'results/ALL_14_SECTION_EFFECTS.tsv',A/'results/COHORT_SPECIFIC_SUMMARY.tsv',B/'PATIENT_RESULT_SUMMARY.tsv',R/'phase219_spatial_reference_identifiability_20260920_v1/program_check/PROGRAM_CHECK_SUMMARY.tsv']:
  sources.append(dict(**meta(f),sha256=sha(f)))
 save('audit/SOURCE_MANIFEST.json',sources)
 for x in tracks:assert meta(x['path'])==x
 for x in sources:
  m=meta(x['path']);assert all(m[k]==x[k] for k in m)
 save('audit/PRESERVATION.json',dict(old_RDS_checked=len(tracks),old_RDS_changed=0,source_files_changed=0))
 save('SOFTWARE.json',dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,h5py=h5py.__version__))
 save('STATUS.json',dict(status='DISPLAY_DATA_PREPARED',new_figures=0,new_inferential_tests=0,new_mapping=0))
 (P/'EXIT_STATUS').write_text('0\n')
try:run()
except Exception:
 (P/'ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');raise
