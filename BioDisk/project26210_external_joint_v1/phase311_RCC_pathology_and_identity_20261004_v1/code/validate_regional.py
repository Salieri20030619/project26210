from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd,anndata as ad
from sklearn.neighbors import BallTree
from scipy.stats import rankdata
from scipy import sparse
P=Path(__file__).resolve().parents[1];O=P.parent/'phase310_RCC_mapped_spatial_inputs_20261004_v1'
assert (P/'logs/run_regional.exit').read_text().strip()=='0'
checks=[]
def ck(n,b):
 checks.append({'check':n,'pass':bool(b)})
 if not b:raise AssertionError(n)
effects=pd.read_csv(P/'tables/REGIONAL_ALL_SPATIAL_EFFECTS.tsv',sep='\t')
ck('72_planned_models',len(effects)==72)
ck('no_inferential_P_FDR',effects[['P','FDR']].isna().all().all())
frozen=pd.read_csv(P/'tables/CELL_IDENTITY_FROZEN.tsv.gz',sep='\t',dtype={'original_cell_id':str})
a=ad.read_h5ad(P/'resources/proseg_region6_identity.h5ad')
ck('frozen_id_order',np.array_equal(frozen.original_cell_id,a.obs_names.astype(str)))
ck('mutually_exclusive_targets',(frozen[['mature_myeloid','strict_fibroblast','strict_CD8']].sum(axis=1)<=1).all())
common=frozen.paired_geometry_eligible.to_numpy();f=frozen.loc[common].reset_index(drop=True)
old=pd.read_csv(O/'tables/PILOT_CELL_METADATA.tsv',sep='\t',dtype={'cell_id':str}).set_index('cell_id').loc[f.original_cell_id]
geos={'Proseg_centroid':a.obsm['spatial'][common],'Original_nuclear_centroid_proxy':old[['nucleus_x_proxy_um','nucleus_y_proxy_um']].to_numpy()}
inputs=pd.read_csv(P/'tables/REGIONAL_MODEL_INPUTS.tsv.gz',sep='\t',dtype={'original_cell_id':str})
ck('model_input_ID_qualification',set(inputs.original_cell_id).issubset(set(f.loc[f.mature_myeloid,'original_cell_id'])))
C=a.layers['maxpost_counts'].toarray()[common];g=list(a.var_names).index('SPP1')
features=pd.read_csv(O/'tables/PILOT_FEATURE_AUDIT.tsv',sep='\t');g0=int(features.loc[features.gene.eq('SPP1'),'matrix_col'].iloc[0])
A=sparse.load_npz(O/'resources/pilot_assigned_counts_qv20.npz')[old.matrix_row.to_numpy(int)]
N=sparse.load_npz(O/'resources/pilot_nuclear_counts_qv20.npz')[old.matrix_row.to_numpy(int)]
assays={'Proseg':(C[:,g],C.sum(axis=1)),'OriginalAssigned':(A[:,g0].toarray().ravel(),np.asarray(A.sum(axis=1)).ravel()),'OriginalNuclear':(N[:,g0].toarray().ravel(),np.asarray(N.sum(axis=1)).ravel())}
imap={c:i for i,c in enumerate(f.original_cell_id)}
maxgeo=maxrho=0.
for geo,xy in geos.items():
 at=BallTree(xy);ft=BallTree(xy[f.strict_fibroblast]) if f.strict_fibroblast.any() else None;tt=BallTree(xy[f.strict_CD8]) if f.strict_CD8.any() else None
 for rad in [50,100]:
  base=inputs[(inputs.geometry==geo)&(inputs.radius_um==rad)&(inputs.assay=='Proseg')]
  ix=np.array([imap[x] for x in base.original_cell_id],int)
  if not len(ix):continue
  centers=xy[ix];nc=at.query_radius(centers,r=rad,count_only=True)-1
  fc=ft.query_radius(centers,r=rad,count_only=True) if ft else np.zeros(len(ix),int)
  nd=tt.query(centers,k=1)[0].ravel() if tt else np.repeat(np.nan,len(ix))
  ff=fc/nc;err=max(float(np.max(np.abs(ff-base.fibroblast_fraction))),float(np.max(np.abs(nd-base.nearest_CD8_distance_um))));maxgeo=max(maxgeo,err)
  ck('BallTree_geometry_'+geo+str(rad),err<1e-8)
  # Direct all-pairs Euclidean check on first 8 centers.
  for j in range(min(8,len(ix))):
   all_d=np.sqrt(((xy-centers[j])**2).sum(axis=1))
   ck('bruteforce_neighborhood_'+geo+str(rad)+'_'+str(j),int((all_d<=rad).sum())-1==int(nc[j]) and int((all_d[f.strict_fibroblast]<=rad).sum())==int(fc[j]))
   ck('bruteforce_CD8_'+geo+str(rad)+'_'+str(j),np.isclose(all_d[f.strict_CD8].min(),nd[j],atol=1e-9))
  for assay,(sp,lib) in assays.items():
   d=inputs[(inputs.geometry==geo)&(inputs.radius_um==rad)&(inputs.assay==assay)]
   ck('common_centers_'+geo+str(rad)+assay,np.array_equal(d.original_cell_id,base.original_cell_id))
   ck('SPP1_source_counts_'+geo+str(rad)+assay,np.array_equal(d.SPP1_count,sp[ix]))
   ck('library_source_counts_'+geo+str(rad)+assay,np.array_equal(d.library_counts,lib[ix]))
   sx=np.log1p(sp[ix]/np.maximum(lib[ix],1)*1e4)
   Z=np.column_stack([np.ones(len(ix)),rankdata(lib[ix]),rankdata(nc),rankdata(centers[:,0]),rankdata(centers[:,1])])
   Q,_=np.linalg.qr(Z)
   for endpoint,x,y in [('SPP1_vs_fibroblast_fraction',sx,ff),('SPP1_vs_nearest_CD8_distance',sx,nd),('fibroblast_fraction_vs_nearest_CD8_distance',ff,nd)]:
    for model in ['unadjusted','depth_density_xy']:
     row=effects[(effects.geometry==geo)&(effects.assay==assay)&(effects.radius_um==rad)&(effects.endpoint==endpoint)&(effects.model==model)].iloc[0]
     if pd.isna(row.rho):
      ck('NA_reason_'+geo+str(rad)+assay+endpoint+model,pd.notna(row.not_assessable_reason) and bool(row.not_assessable_reason));continue
     rx,ry=rankdata(x),rankdata(y)
     if model!='unadjusted':rx=rx-Q@(Q.T@rx);ry=ry-Q@(Q.T@ry)
     corr=float(np.corrcoef(rx,ry)[0,1]);err=abs(corr-row.rho);maxrho=max(maxrho,err)
     ck('effect_QR_recompute_'+geo+str(rad)+assay+endpoint+model,err<1e-8)
ck('all_unassessable_have_reason',effects.loc[effects.rho.isna(),'not_assessable_reason'].notna().all())
rec={'status':'PASS','time':datetime.datetime.now().astimezone().isoformat(),'named_checks':len(checks),'checks':checks,'max_geometry_error':maxgeo,'max_rho_error':maxrho,'model_count':len(effects),'no_population_inference':True,'method':'independent BallTree and exact Euclidean distances; raw input counts; QR rank-residual correlations; all endpoint/qualification rows'}
(P/'qc/REGIONAL_INDEPENDENT_VALIDATION.json').write_text(json.dumps(rec,indent=2));print(json.dumps({k:v for k,v in rec.items() if k!='checks'},indent=2))
