from pathlib import Path
from collections import deque
import json,hashlib,gzip,io,traceback,sys
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from anndata.io import read_elem
from PIL import Image
P=Path(__file__).resolve().parents[1];R=P.parent
A=R/'phase223_spatial_cohort_reinclusion_20260920_v1';S=R/'phase221_spatial_S1_identity_eligibility_freeze_20260920_v1';B=R/'phase222_spatial_S2_S5_execution_20260920_v1'
def save(n,o):(P/n).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(n,d):d.to_csv(P/n,sep='\t',index=False)
def meta(f):
 f=Path(f);s=f.stat();return dict(path=str(f),bytes=s.st_size,mtime_ns=s.st_mtime_ns,inode=s.st_ino)
def run():
 old=[meta(x['path']) for x in json.loads((A/'audit/RDS_BEFORE.json').read_text())];save('audit/RDS_BEFORE.json',old)
 records=json.loads((A/'ASSET_ROLE_MANIFEST.json').read_text())
 pat=pd.read_csv(S/'FROZEN_SECTION_PATIENT_MAP.tsv',sep='\t').set_index('section')
 avail=pd.read_csv(S/'S3_S4_AVAILABILITY.tsv',sep='\t').set_index('section')
 inventory=json.loads((R/'phase217_spatial_asset_progress_inventory_20260920_v1/ASSET_INVENTORY.json').read_text())
 for sid,v in pat.iterrows():
  records.append(dict(sample=sid,dataset='GSE175540',patient=v.patient_id,input_path=str(R.parent/'project26210_external_joint_v1_phase54_gse175540_spatial_v1/run_v1/sections_v1'/sid/'processed_section_v1.h5ad')))
 defs=json.loads((R/'phase227_paper_spatial_replication_plan_20260921_v1/evidence/AUTHOR_PROGRAM_GENESETS.json').read_text())
 programs={k:defs[k] for k in ['POSTNCAF','APOEMph']}
 programs.update(ECM4=['COL1A1','COL1A2','DCN','LUM'],CD8_expression=['CD8A','CD8B'],Tcell_expression=['CD3D','CD3E','TRAC'])
 save('PROGRAM_DEFINITIONS.json',dict(programs=programs,method='log1p raw counts normalized to 10000; per-gene within-section population z; equal mean; not inferred cell fractions',technical_gate='all listed genes present; at least20 spots detect >=2 program genes; nonzero variance',gene_sets_not_new_selected=True))
 summaries=[];gene_qc=[];assets=[];sources=[];region_summary=[]
 for rec in records:
  sid=rec['sample'];legacy=rec['dataset']!='GSE175540'
  fp=(A/'inputs'/(sid+'_FROZEN_TRANSFER_INPUTS.tsv.gz')) if legacy else S/'frozen_inputs'/(sid+'_SPOT_INPUTS.tsv.gz')
  gp=(A/'inputs'/(sid+'_adjacency.npz')) if legacy else S/'frozen_inputs'/(sid+'_adjacency.npz')
  d=pd.read_csv(fp,sep='\t').rename(columns={'x_pixels':'x','y_pixels':'y','SPP1_share_within_mature_macro':'SPP1_share','neighbor_macro_fraction':'neighbor_macro','neighbor_log1p_UMI':'neighbor_depth'})
  d.insert(0,'sample',sid);d.insert(1,'patient',rec['patient']);d.insert(2,'dataset',rec['dataset'])
  G=sparse.load_npz(gp);deg=np.asarray(G.sum(axis=1)).ravel()
  assert G.shape==(len(d),len(d)) and (G!=G.T).nnz==0 and deg.max()<=6
  with h5py.File(rec['input_path']) as h:
   obs=read_elem(h['obs']);var=read_elem(h['var']);xy=read_elem(h['obsm']['spatial'])
   assert np.array_equal(obs.index,d.spot_id)
   np.testing.assert_allclose(xy,d[['x','y']],rtol=0,atol=1e-8)
   node=h['layers/counts'] if 'counts' in h.get('layers',{}) else h['X'];X=read_elem(node).astype(float)
   assert np.isfinite(X.data).all() and (X.data>=0).all() and np.allclose(X.data,np.round(X.data))
   raw_layer='layers/counts' if 'counts' in h.get('layers',{}) else 'X'
   totals=np.asarray(X.sum(axis=1)).ravel();assert (totals>0).all()
   if legacy:
    spatial=h['uns/spatial'];lib=list(spatial)[0];shape=spatial[lib]['images/hires'].shape
    scal=read_elem(spatial[lib]['scalefactors']);scale=float(scal['tissue_hires_scalef']);diam=float(scal['spot_diameter_fullres'])
    image_path=rec['input_path']+'::uns/spatial/'+lib+'/images/hires'
   else:
    matched={x['kind']:x['path'] for x in inventory if sid in x['path'] and x['kind'] in ['hires_image','scalefactors','coordinates']}
    with gzip.open(matched['scalefactors'],'rt') as f:scal=json.load(f)
    scale=float(scal['tissue_hires_scalef']);diam=float(scal['spot_diameter_fullres'])
    with gzip.open(matched['hires_image'],'rb') as f:
     im=Image.open(f);w,hh=im.size;shape=(hh,w,3)
    image_path=matched['hires_image'];sources.extend(meta(v) for v in matched.values())
   inside=(xy[:,0]*scale>=0)&(xy[:,0]*scale<shape[1])&(xy[:,1]*scale>=0)&(xy[:,1]*scale<shape[0])
   # Record image ranges rather than silently dropping any spots.
   assets.append(dict(sample=sid,image_path=image_path,image_width=shape[1],image_height=shape[0],scale=scale,spot_diameter_fullres=diam,inside_image=int(inside.sum()),QC_spots=len(d),all_inside=bool(inside.all()),verification='coordinate range and image header; not new pathology review'))
  tumor=d.malignant_fraction.to_numpy()>=.2;tn=G@tumor.astype(float)
  interface=(tumor&(tn<deg))|((~tumor)&(tn>0))
  core=tumor&(tn==deg)&(deg>0)
  source=d.source_eligible.to_numpy(bool)
  np.testing.assert_array_equal(source,(d.total_macro_fraction>=.1)&(deg>=3))
  if legacy:
   np.testing.assert_array_equal(d.core_source,core&source)
   np.testing.assert_array_equal(d.boundary_source,interface&source)
   ecml=d.ECM_low.to_numpy(bool);ecmh=d.ECM_high.to_numpy(bool)
  else:
   rp=B/'cache'/(sid+'_REGION_MEMBERS.tsv.gz');rg=pd.read_csv(rp,sep='\t');assert np.array_equal(rg.spot_id,d.spot_id)
   np.testing.assert_array_equal(rg.tumor_expression_region,tumor)
   np.testing.assert_array_equal(rg.core_source,core&source);np.testing.assert_array_equal(rg.interface_source,interface&source)
   ecml=rg.ECM_low_source.to_numpy(bool);ecmh=rg.ECM_high_source.to_numpy(bool);sources.append(meta(rp))
   for c in ['TLS_positive','TLS_explicit_negative','TLS_peri_source','TLS_far_source']:d[c]=rg[c].to_numpy()
  d['degree']=deg;d['tumor_expression_region']=tumor;d['expression_interface']=interface
  d['core_source']=core&source;d['interface_source']=interface&source
  d['partition']=np.where(interface,'expression_interface',np.where(core,'malignant_expression_core',np.where(tumor,'isolated_malignant','nonmalignant_background')))
  d['ECM_low_source']=ecml;d['ECM_high_source']=ecmh
  # Graph-hop distance does not cross gaps or invent physical pathology borders.
  hops=np.full(len(d),np.inf);queue=deque(np.flatnonzero(interface).tolist())
  for i in queue:hops[i]=0
  while queue:
   i=queue.popleft()
   for j in G.indices[G.indptr[i]:G.indptr[i+1]]:
    if not np.isfinite(hops[j]):hops[j]=hops[i]+1;queue.append(int(j))
  d['expression_interface_distance_hops']=np.where(np.isfinite(hops),hops,np.nan)
  d['distance_component_reachable']=np.isfinite(hops)
  d['x_hires_px']=xy[:,0]*scale;d['y_hires_px']=xy[:,1]*scale
  # Reuse composition partitions; never construct boundaries from these target programs.
  pgates={}
  for key,genes in programs.items():
   ix=var.index.get_indexer(genes);present=ix>=0
   if present.all():
    values=X[:,ix].toarray();z=np.log1p(values/totals[:,None]*1e4)
    sd=z.std(axis=0);z=(z-z.mean(axis=0))/np.where(sd==0,1,sd)
    score=z.mean(axis=1);det=(values>0).sum(axis=1);passgate=int((det>=2).sum())>=20 and score.std()>0
    d['program_'+key]=score;d['detected_genes_'+key]=det
   else:
    score=np.full(len(d),np.nan);passgate=False;det=np.zeros(len(d))
    d['program_'+key]=score;d['detected_genes_'+key]=np.nan
   pgates[key]=bool(passgate)
   gene_qc.append(dict(sample=sid,program=key,listed_genes=len(genes),present_genes=int(present.sum()),missing=';'.join(np.array(genes)[~present]),spots_detect_ge2=int((det>=2).sum()),technical_program_pass=bool(passgate),not_cell_identity=True))
  for reg,g in d.groupby('partition'):
   region_summary.append(dict(sample=sid,dataset=rec['dataset'],patient=rec['patient'],partition=reg,spots=len(g),source_spots=int(g.source_eligible.sum())))
  summary=dict(sample=sid,dataset=rec['dataset'],patient=rec['patient'],spots=len(d),source_spots=int(source.sum()),core_source=int((core&source).sum()),interface_source=int((interface&source).sum()),ECM_low_source=int(ecml.sum()),ECM_high_source=int(ecmh.sum()),expression_boundary_eligible=int((core&source).sum())>=20 and int((interface&source).sum())>=20,ECM_region_eligible=int(ecml.sum())>=20 and int(ecmh.sum())>=20,components=int(d.component.nunique()),no_interface_path_spots=int((~np.isfinite(hops)).sum()),histological_boundary_qualified=False,CAF_abundance_identity=False,CD8_abundance_identity=False,old_SPP1_not_C4=True,counts_layer=raw_layer,all_spots_in_image=bool(inside.all()),program_POSTNCAF_eligible=pgates['POSTNCAF'],program_APOEMph_eligible=pgates['APOEMph'],program_ECM4_eligible=pgates['ECM4'],program_CD8_expression_eligible=pgates['CD8_expression'],TLS_eligible=False)
  if not legacy:summary.update(TLS_eligible=int(d.TLS_peri_source.sum())>=20 and int(d.TLS_far_source.sum())>=20,TLS_annotation_type=avail.loc[sid,'TLS_annotation_type'])
  else:summary.update(TLS_annotation_type='NO_AUDITED_ANNOTATION')
  summaries.append(summary);sources.extend([meta(fp),meta(gp),meta(rec['input_path'])])
  tab('spots/'+sid+'_F1_SPATIAL_INPUTS.tsv.gz',d)
  print('F1',sid,len(d),summary['expression_boundary_eligible'],pgates,flush=True)
 summary=pd.DataFrame(summaries);assert len(summary)==38 and summary.spots.sum()==114640
 tab('tables/SECTION_ELIGIBILITY.tsv',summary);tab('tables/PROGRAM_TECHNICAL_QC.tsv',pd.DataFrame(gene_qc));tab('tables/PARTITION_COVERAGE.tsv',pd.DataFrame(region_summary));tab('tables/IMAGE_COORDINATE_SUPPORT.tsv',pd.DataFrame(assets))
 endpointcols=['expression_boundary_eligible','ECM_region_eligible','TLS_eligible','program_POSTNCAF_eligible','program_APOEMph_eligible','program_ECM4_eligible','program_CD8_expression_eligible','histological_boundary_qualified','CAF_abundance_identity','CD8_abundance_identity']
 rows=[]
 for dataset,g in summary.groupby('dataset'):
  for ep in endpointcols:
   gg=g[g[ep]];rows.append(dict(dataset=dataset,endpoint=ep,sections=len(gg),patient_units=gg.patient.nunique()))
 tab('tables/COHORT_ENDPOINT_ELIGIBILITY.tsv',pd.DataFrame(rows))
 save('audit/SOURCE_MANIFEST.json',sources)
 for x in old:assert meta(x['path'])==x
 for x in sources:assert meta(x['path'])==x
 save('audit/PRESERVATION.json',dict(RDS_checked=len(old),RDS_changed=0,input_metadata_changed=0))
 save('SOFTWARE.json',dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,h5py=h5py.__version__))
 save('STATUS.json',dict(status='F1_COMPUTE_COMPLETE',sections=38,spots=114640,new_inference=False,new_figures=0))
 (P/'EXIT_STATUS').write_text('0\n')
try:run()
except Exception:
 (P/'ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');raise
