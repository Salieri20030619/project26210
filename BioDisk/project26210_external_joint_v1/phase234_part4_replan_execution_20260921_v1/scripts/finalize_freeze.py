from common import *
import h5py, shutil, traceback
from anndata.io import read_elem
from scipy.sparse.csgraph import dijkstra
def main():
 assert (P/'A0_COMPLETE.json').exists()
 assert not (P/'A1_COMPLETE.json').exists()
 # Preserve pre-association diameter diagnostic, then resolve against independently known array geometry.
 shutil.copy2(P/'tables/SCALE_ELIGIBILITY.tsv',P/'audit/INITIAL_DIAMETER_SCALE_DIAGNOSTIC.tsv')
 regs=pd.read_csv(P/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t');checks=[];scales=[]
 for ri,r in regs.iterrows():
  sid=r['sample'];d=pd.read_csv(P/'inputs'/(sid+'_INPUTS.tsv.gz'),sep='\t')
  with h5py.File(r.input_path) as h:obs=read_elem(h['obs'])
  G=sparse.load_npz(P/'graphs'/(sid+'_ring.npz')).tocsr();rr,cc=G.nonzero()
  ar=obs.array_row.to_numpy(float);ac=obs.array_col.to_numpy(float)
  dr=np.abs(ar[rr]-ar[cc]);dc=np.abs(ac[rr]-ac[cc])
  lattice=bool((((dr==0)&(dc==2))|((dr==1)&(dc==1))).all())
  uv=np.column_stack([ac*50,ar*(50*np.sqrt(3))])
  xy=d[['x','y']].to_numpy();design=np.column_stack([np.ones(len(d)),uv])
  beta=np.linalg.lstsq(design,xy,rcond=None)[0];pred=design@beta
  sv=np.linalg.svd(beta[1:],compute_uv=False);factor=1/np.sqrt(sv.prod())
  rmse=float(np.sqrt(np.mean((pred-xy)**2))*factor)
  ratio=float(sv.max()/sv.min());physical=lattice and rmse<5 and ratio<1.1
  median_px=float(np.median(np.linalg.norm(xy[rr]-xy[cc],axis=1)))
  checks.append(dict(sample=sid,array_lattice_edges_verified=lattice,image_array_affine_rmse_um=rmse,image_array_anisotropy_ratio=ratio,image_pixel_to_um=factor,calibrated_median_edge_um=median_px*factor,diameter_based_median_edge_um=r.median_edge_um,method='Visium array_row/array_col: horizontal pitch100um, row pitch sqrt(3)*50um; all first-ring edges100um',physical_eligible=physical))
  regs.loc[ri,'physical_eligible']=physical;regs.loc[ri,'pixel_to_um']=factor;regs.loc[ri,'median_edge_um']=100.
  d['physical_x_um']=uv[:,0];d['physical_y_um']=uv[:,1]
  if physical:
   WG=G.astype(float);WG.data[:]=100.
   rows=[];cols=[];vals=[]
   for start in range(0,len(d),256):
    dst=dijkstra(WG,directed=False,indices=np.arange(start,min(start+256,len(d))),limit=800)
    a,b=np.where(np.isfinite(dst)&(dst>0));rows.extend((a+start).tolist());cols.extend(b.tolist());vals.extend(dst[a,b].tolist())
   D=sparse.csr_matrix((vals,(rows,cols)),shape=G.shape)
   assert (D-D.T).nnz==0 and not D.diagonal().any()
   assert np.array_equal(d.component.to_numpy()[np.array(rows)],d.component.to_numpy()[np.array(cols)])
   sparse.save_npz(P/'graphs'/(sid+'_distance800.npz'),D)
   ids=np.flatnonzero(d.expression_interface.to_numpy(bool))
   dd=dijkstra(WG,directed=False,indices=ids,min_only=True) if len(ids) else np.full(len(d),np.inf)
   d['expression_distance_um']=np.where(np.isfinite(dd),dd,np.nan)
   for radius in [800,400,200]:
    H=D.copy();H.data=(H.data<=radius).astype(float);H.eliminate_zeros();nn=np.asarray(H.sum(axis=1)).ravel()
    scales.append(dict(sample=sid,scale=str(radius),physical_eligible=True,sources_ge3=int((nn>=3).sum()),median_neighbors=float(np.median(nn)),min_neighbors=int(nn.min()),max_neighbors=int(nn.max()),median_section_fraction=float(np.median(nn/max(len(d)-1,1))),max_section_fraction=float(nn.max()/max(len(d)-1,1)),median_edge_um=100))
  else:
   for radius in [800,400,200]:scales.append(dict(sample=sid,scale=str(radius),physical_eligible=False,reason='ARRAY_IMAGE_GEOMETRY_FAILED'))
  deg=np.asarray(G.sum(axis=1)).ravel();scales.append(dict(sample=sid,scale='ring',physical_eligible=True,sources_ge3=int((deg>=3).sum()),median_neighbors=float(np.median(deg)),median_edge_um=100))
  a=pd.read_csv(r.abundance_path,sep='\t',index_col=0)
  if r.dataset=='GSE175540':a=a[[c for c in a if c.startswith('q05__')]];a.columns=[c.removeprefix('q05__') for c in a]
  else:a.columns=[c.removeprefix('q05cell_abundance_w_sf_') for c in a]
  assert np.array_equal(a.index,d.spot_id)
  for k,c in [('endothelial_q05','Endothelial'),('pericyte_q05','SMC/Pericyte')]:d[k]=a[c].to_numpy()
  macrocols=[c for c in a if any(k in c for k in ['SPP1','C1QC']) or c=='Other Macrophage']
  d['total_macro_q05']=a[macrocols].sum(axis=1).to_numpy()
  tab('inputs/'+sid+'_INPUTS.tsv.gz',d)
  print('A0_SCALE',sid,physical,round(rmse,3),round(ratio,4),flush=True)
 ep=pd.read_csv(P/'tables/FROZEN_ENDPOINTS.tsv',sep='\t')
 for x in ['V1_SPP1_q05','V2_SPP1_fraction','V3_SPP1_macro_share','V4_SPP1_program']:
  for y in ['endothelial_q05','pericyte_q05']:
   ep.loc[len(ep)]=['SPP1_vascular',x+'__'+y,x,y,'parallel_absolute_target']
 for y in ['endothelial_q05','pericyte_q05','POSTNCAF','ECM4']:
  ep.loc[len(ep)]=['comparison_partners','total_macro_q05__'+y,'total_macro_q05',y,'comparison']
 tab('tables/FROZEN_ENDPOINTS.tsv',ep)
 tab('tables/INPUT_SECTION_REGISTER.tsv',regs);tab('tables/SCALE_ELIGIBILITY.tsv',pd.DataFrame(scales));tab('audit/PHYSICAL_CALIBRATION.tsv',pd.DataFrame(checks))
 rules=json.loads((P/'IMPLEMENTATION_FREEZE_BEFORE_ASSOCIATION.json').read_text())
 save('audit/INITIAL_IMPLEMENTATION_FREEZE.json',rules)
 rules.update(physical_resolution_time=stamp(),radius_distance='shortest physical path along verified Visium array lattice; all edges100um; image-array affine compatibility verified (RMSE<5um,anisotropy<1.1); diameter metadata disagreement retained, not used to define radii',vascular_targets='both q05 and total13-normalized Endothelial/SMC parallel for all V1-V4; four targets included in vascular family',adjustment_target_list='endothelial_fraction,pericyte_fraction,POSTNCAF,ECM4,cytotoxicity,tcell_dysfunction,immunosuppression_core,invasion_emt; additional q05 vessel targets use target q05 and same fixed covariates',common_source_qualification='separate Q0/Q1/Q2 own masks; commonQ2 mask identical by construction, patient matching separately')
 save('IMPLEMENTATION_FREEZE_BEFORE_ASSOCIATION.json',rules)
 save('A0_FINAL_FREEZE.json',dict(time=stamp(),sections=38,spot_count=114640,physical_sections=int(regs.physical_eligible.sum()),endpoint_count=len(ep),lr_candidates=len(pd.read_csv(P/'tables/FROZEN_LR_CANDIDATES.tsv',sep='\t')),association_started=False,parameters_sha256=sha(P/'IMPLEMENTATION_FREEZE_BEFORE_ASSOCIATION.json'),endpoints_sha256=sha(P/'tables/FROZEN_ENDPOINTS.tsv')))
if __name__=='__main__':
 try:main()
 except Exception:(P/'ERROR_SCALE.txt').write_text(traceback.format_exc());raise
