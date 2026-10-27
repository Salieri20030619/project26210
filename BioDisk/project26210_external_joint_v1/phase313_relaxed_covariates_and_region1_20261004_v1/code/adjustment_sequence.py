from pathlib import Path
import json,datetime
import numpy as np,pandas as pd
from scipy.stats import rankdata
P=Path(__file__).resolve().parents[1];Q=P.parent/'phase312_relaxed_entry_spatial_20261004_v1'
d=pd.read_csv(Q/'tables/ALL_MODEL_INPUTS.tsv.gz',sep='\t',dtype={'original_cell_id':str})
d=d[(d.condition=='L2_cluster_broad')&(d.center_set=='condition_members')];rows=[];qa=[]
old=pd.read_csv(Q/'tables/ALL_SPATIAL_EFFECTS.tsv',sep='\t')
old=old[(old.condition=='L2_cluster_broad')&(old.center_set=='condition_members')]
models={'unadjusted':[],'depth_only':['library_counts'],'depth_plus_density':['library_counts','all_neighbor_count'],'depth_density_xy':['library_counts','all_neighbor_count','x','y']}
for (geo,assay,rad),b in d.groupby(['geometry','assay','radius_um'],sort=False):
 for ep,x,y in [('SPP1_vs_fibroblast_fraction',b.SPP1_normalized,b.fibroblast_fraction),('SPP1_vs_nearest_CD8_distance',b.SPP1_normalized,b.nearest_CD8_distance_um),('fibroblast_fraction_vs_nearest_CD8_distance',b.fibroblast_fraction,b.nearest_CD8_distance_um)]:
  for model,covs in models.items():
   X=rankdata(x);Y=rankdata(y);Z=np.column_stack([np.ones(len(b))]+[rankdata(b[c]) for c in covs]);rk=np.linalg.matrix_rank(Z);reason=''
   if rk<Z.shape[1] or len(np.unique(x))<2 or len(np.unique(y))<2:rho=np.nan;reason='rank_or_variance'
   else:
    X=X-Z@np.linalg.lstsq(Z,X,rcond=None)[0];Y=Y-Z@np.linalg.lstsq(Z,Y,rcond=None)[0];rho=float(X@Y/(np.linalg.norm(X)*np.linalg.norm(Y)))
    U,s,_=np.linalg.svd(Z,full_matrices=False);U=U[:,:rk];x2=rankdata(x);y2=rankdata(y);x2-=U@(U.T@x2);y2-=U@(U.T@y2);rho2=float(np.corrcoef(x2,y2)[0,1]);assert abs(rho-rho2)<1e-10
   if model in ['unadjusted','depth_density_xy']:
    ref=old[(old.geometry==geo)&(old.assay==assay)&(old.radius_um==rad)&(old.endpoint==ep)&(old.model==model)].iloc[0];assert abs(rho-ref.rho)<1e-10
   rows.append({'geometry':geo,'assay':assay,'radius_um':rad,'endpoint':ep,'model':model,'centers':len(b),'rho':rho,'P':None,'FDR':None,'controls':';'.join(covs),'design_rank':int(rk),'reason':reason})
out=pd.DataFrame(rows);out.to_csv(P/'tables/A_FIXED_ADJUSTMENT_SEQUENCE.tsv',sep='\t',index=False)
rec={'status':'PASS','time':datetime.datetime.now().astimezone().isoformat(),'models':len(out),'old_endpoint_model_rows_matched':72,'all_new_effects_independent_SVD_verified':True,'no_P_FDR':True,'interpretation':'fixed-sequence change point, not causal/confounder identification'}
(P/'qc/A_ADJUSTMENT_VALIDATION.json').write_text(json.dumps(rec,indent=2))
print(json.dumps(rec,indent=2));print(out[(out.geometry=='Proseg_centroid')&(out.assay=='Proseg')&(out.radius_um==50)][['endpoint','model','rho']].to_string(index=False),flush=True)
