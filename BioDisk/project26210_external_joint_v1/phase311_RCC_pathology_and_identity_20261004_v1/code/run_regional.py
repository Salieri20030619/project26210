from pathlib import Path
import json,datetime
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import rankdata
P=Path(__file__).resolve().parents[1];old=P.parent/'phase310_RCC_mapped_spatial_inputs_20261004_v1'
a=ad.read_h5ad(P/'resources/proseg_region6_identity.h5ad');C=a.layers['maxpost_counts'].toarray();gn=list(a.var_names);ids=a.obs_names.astype(str)
ann=json.loads((P/'resources/CELL_CLUSTER_ANNOTATION.json').read_text())
groups=a.obs.identity_cluster.astype(str).map({k:v['label'] for k,v in ann['clusters'].items()}).fillna('Unknown').to_numpy()
markers=json.loads((P/'resources/SEGMENTATION_IDENTITY_PROTOCOL.json').read_text())['identity_plan']['markers']
def dc(gs):
 js=[gn.index(g) for g in gs if g in gn];return (C[:,js]>0).sum(axis=1) if js else np.zeros(len(a),int)
t=dc(markers['T']);my=dc(markers['myeloid']);fb=dc(markers['fibroblast']);en=dc(markers['endothelial']);ep=dc(markers['epithelial']);bline=dc(markers['B'])
macro=(groups=='Myeloid')&((dc(['C1QA','C1QB','C1QC'])>=1)|((dc(['CD68'])>=1)&(dc(['CD163'])>=1)))&(t<2)&(ep<2)&(en<2)
fib=(groups=='Fibroblast')&(fb>=2)&(dc(['COL1A1','COL1A2','COL3A1'])>=1)&(dc(['MYH11','ACTG2'])==0)&(t<2)&(my<2)&(ep<2)&(en<2)
cd8=(groups=='T')&(t>=2)&(dc(['CD8A','CD8B'])>=1)&(my<2)&(ep<2)&(en<2)
macro=macro&(bline<2);fib=fib&(bline<2);cd8=cd8&(bline<2)
cell=a.obs.copy();cell['original_cell_id']=ids;cell['reviewed_cluster_lineage']=groups;cell['mature_myeloid']=macro;cell['strict_fibroblast']=fib;cell['strict_CD8']=cd8
oc=pd.read_csv(old/'tables/PILOT_CELL_METADATA.tsv',sep='\t',dtype={'cell_id':str}).set_index('cell_id')
og=pd.read_csv(old/'tables/PILOT_FEATURE_AUDIT.tsv',sep='\t')
assert set(ids).issubset(set(oc.index))
alignment=oc.loc[ids];nucxy=alignment[['nucleus_x_proxy_um','nucleus_y_proxy_um']].to_numpy()
common=alignment.nuclear_RNA.ge(5).to_numpy() & np.isfinite(nucxy).all(axis=1)
cell['paired_geometry_eligible']=common
cell.to_csv(P/'tables/CELL_IDENTITY_FROZEN.tsv.gz',sep='\t',index=False,compression='gzip')
counts={'QC_cells':len(a),'same_ID_geometry_cells':int(common.sum()),'myeloid_all':int(macro.sum()),'fibroblast_all':int(fib.sum()),'CD8_all':int(cd8.sum()),'myeloid_common':int((macro&common).sum()),'fibroblast_common':int((fib&common).sum()),'CD8_common':int((cd8&common).sum())}
(P/'qc/IDENTITY_MEMBERSHIP_FROZEN.json').write_text(json.dumps({'time':datetime.datetime.now().astimezone().isoformat(),'counts':counts,'annotation':ann,'SPP1_not_used_for_identity':True},indent=2))
# All subsequent descriptive effects use one exact common set, no cell-type refitting.
idx=np.flatnonzero(common);centers_mask=macro[common];F=fib[common];T=cd8[common];pid=ids[common]
A=sparse.load_npz(old/'resources/pilot_assigned_counts_qv20.npz')[alignment.matrix_row.to_numpy(int)]
N=sparse.load_npz(old/'resources/pilot_nuclear_counts_qv20.npz')[alignment.matrix_row.to_numpy(int)]
spidx=gn.index('SPP1');old_sp=int(og.loc[og.gene.eq('SPP1'),'matrix_col'].iloc[0])
assays={'Proseg':(C[common,spidx],C[common].sum(axis=1)),'OriginalAssigned':(A[common,old_sp].toarray().ravel(),np.asarray(A[common].sum(axis=1)).ravel()),'OriginalNuclear':(N[common,old_sp].toarray().ravel(),np.asarray(N[common].sum(axis=1)).ravel())}
geometries={'Proseg_centroid':a.obsm['spatial'][common],'Original_nuclear_centroid_proxy':nucxy[common]}
rows=[];contrasts=[];detailed=[]
def corr(x,y,controls=None):
 if len(x)<10 or np.unique(x).size<2 or np.unique(y).size<2:return np.nan,'insufficient_variation_or_centers',None
 X=rankdata(x);Y=rankdata(y);rank=1
 if controls is not None:
  Z=np.column_stack([np.ones(len(x))]+[rankdata(v) for v in controls.T]);rank=int(np.linalg.matrix_rank(Z))
  if rank<Z.shape[1]:return np.nan,'rank_deficient_controls',rank
  X=X-Z@np.linalg.lstsq(Z,X,rcond=None)[0];Y=Y-Z@np.linalg.lstsq(Z,Y,rcond=None)[0]
 X=X-X.mean();Y=Y-Y.mean();den=np.linalg.norm(X)*np.linalg.norm(Y)
 if den<1e-10:return np.nan,'zero_residual_variance',rank
 return float(X@Y/den),'',rank
for geom,xy in geometries.items():
 assert np.isfinite(xy).all()
 tree=cKDTree(xy);ftree=cKDTree(xy[F]) if F.any() else None;ttree=cKDTree(xy[T]) if T.any() else None
 for radius in [50,100]:
  edge=np.minimum.reduce([xy[:,0]-xy[:,0].min(),xy[:,0].max()-xy[:,0],xy[:,1]-xy[:,1].min(),xy[:,1].max()-xy[:,1]])
  ix=np.flatnonzero(centers_mask&(edge>=radius));cen=xy[ix];neighbor=tree.query_ball_point(cen,radius,return_length=True)-1
  fn=ftree.query_ball_point(cen,radius,return_length=True) if ftree else np.zeros(len(ix),int)
  ff=np.divide(fn,neighbor,out=np.full(len(ix),np.nan),where=neighbor>0)
  dist=ttree.query(cen)[0] if ttree else np.full(len(ix),np.nan)
  good=np.isfinite(ff)&np.isfinite(dist);ix=ix[good];ff=ff[good];dist=dist[good];neighbor=neighbor[good];cen=cen[good]
  structural_reason=''
  if len(ix)<100:structural_reason+='myeloid_centers_lt100;'
  if int(F.sum())<100:structural_reason+='fibroblasts_lt100;'
  if int(T.sum())<50:structural_reason+='CD8_lt50;'
  for assay,(spc,total) in assays.items():
   expr=np.log1p(spc[ix]/np.maximum(total[ix],1)*1e4);pos=spc[ix]>0
   spp_reason=structural_reason
   if int(pos.sum())<30:spp_reason+='SPP1_detected_lt30;'
   if int((~pos).sum())<30:spp_reason+='SPP1_undetected_lt30;'
   controls=np.column_stack([total[ix],neighbor,cen])
   for endpoint,x,y in [('SPP1_vs_fibroblast_fraction',expr,ff),('SPP1_vs_nearest_CD8_distance',expr,dist),('fibroblast_fraction_vs_nearest_CD8_distance',ff,dist)]:
    reason=spp_reason if endpoint.startswith('SPP1') else structural_reason
    for model in ['unadjusted','depth_density_xy']:
     rho=np.nan;rank=None
     if not reason:rho,reason,rank=corr(x,y,controls if model!='unadjusted' else None)
     rows.append({'region':'tumor1_region6','geometry':geom,'assay':assay,'radius_um':radius,'endpoint':endpoint,'model':model,'centers':len(ix),'SPP1_detected':int(pos.sum()),'SPP1_undetected':int((~pos).sum()),'fibroblasts':int(F.sum()),'CD8_cells':int(T.sum()),'rho':rho,'P':None,'FDR':None,'design_rank':rank,'not_assessable_reason':reason,'unit':'one unpaired named region; descriptive cells are not biological replicates'})
   contrasts.append({'geometry':geom,'assay':assay,'radius_um':radius,'positive_cells':int(pos.sum()),'negative_cells':int((~pos).sum()),'median_fibroblast_fraction_positive':float(np.median(ff[pos])) if pos.any() else None,'median_fibroblast_fraction_negative':float(np.median(ff[~pos])) if (~pos).any() else None,'median_CD8_distance_positive':float(np.median(dist[pos])) if pos.any() else None,'median_CD8_distance_negative':float(np.median(dist[~pos])) if (~pos).any() else None,'group_comparison_eligible':not bool(spp_reason),'reason':spp_reason})
   for k,j in enumerate(ix):
    detailed.append({'geometry':geom,'assay':assay,'radius_um':radius,'original_cell_id':pid[j],'SPP1_count':float(spc[j]),'SPP1_normalized':float(expr[k]),'library_counts':float(total[j]),'all_neighbor_count':int(neighbor[k]),'fibroblast_fraction':float(ff[k]),'nearest_CD8_distance_um':float(dist[k]),'x':float(cen[k,0]),'y':float(cen[k,1])})
pd.DataFrame(rows).to_csv(P/'tables/REGIONAL_ALL_SPATIAL_EFFECTS.tsv',sep='\t',index=False)
pd.DataFrame(contrasts).to_csv(P/'tables/REGIONAL_ALL_GROUP_CONTRASTS.tsv',sep='\t',index=False)
pd.DataFrame(detailed,columns=['geometry','assay','radius_um','original_cell_id','SPP1_count','SPP1_normalized','library_counts','all_neighbor_count','fibroblast_fraction','nearest_CD8_distance_um','x','y']).to_csv(P/'tables/REGIONAL_MODEL_INPUTS.tsv.gz',sep='\t',index=False,compression='gzip')
summary={'status':'REGIONAL_DESCRIPTIVE_ANALYSIS_COMPLETE','identity_counts':counts,'models':len(rows),'assessable':int(pd.DataFrame(rows).rho.notna().sum()),'tumor_boundary_endpoint':'not_assessable: no independent tumor compartment mapping; not replaced by EPCAM','biological_n':'one region only; donor mapping unresolved; no inferential P/FDR'}
(P/'qc/REGIONAL_ANALYSIS_COMPLETE.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2));print(pd.DataFrame(rows).query("geometry=='Proseg_centroid' and assay=='Proseg'").to_string(index=False),flush=True)
