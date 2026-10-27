from pathlib import Path
import json,datetime,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import rankdata
from sklearn.neighbors import BallTree
P=Path(__file__).resolve().parents[1];Q=P.parent/'phase311_RCC_pathology_and_identity_20261004_v1';O=P.parent/'phase310_RCC_mapped_spatial_inputs_20261004_v1'
checks=[]
def ck(n,b):
 checks.append({'check':n,'pass':bool(b)})
 if not b:raise AssertionError(n)
ck('analysis_exit0',(P/'logs/analyze_relaxed.exit').read_text().strip()=='0')
a=ad.read_h5ad(Q/'resources/proseg_region6_identity.h5ad');C=a.layers['maxpost_counts'].toarray();gn=list(a.var_names)
m=pd.read_csv(P/'tables/ALL_CONDITION_MEMBERS.tsv.gz',sep='\t',dtype={'original_cell_id':str})
strict=pd.read_csv(Q/'tables/CELL_IDENTITY_FROZEN.tsv.gz',sep='\t',dtype={'original_cell_id':str})
ck('complete_cell_ID_order',np.array_equal(m.original_cell_id,a.obs_names.astype(str)))
clusters=a.obs.identity_cluster.astype(str).to_numpy()
def marks(gs):return sum((C[:,gn.index(g)]>0).astype(int) for g in gs if g in gn)
cd=(clusters=='3')&(marks(['CD8A','CD8B'])>0)
expected={'L2_cluster_broad':[(clusters=='1'),(clusters=='0'),cd],'L1_single_marker_supported':[(clusters=='1')&(marks(['CD68','C1QA','C1QB','C1QC','CD163','AIF1','LYZ','FCGR3A'])>0),(clusters=='0')&(marks(['COL1A1','COL1A2','COL3A1','LUM','PDGFRA','DPT','FBLN1'])>0),cd],'L0_strict_members_relaxed_counts':[strict.mature_myeloid.to_numpy(),strict.strict_fibroblast.to_numpy(),strict.strict_CD8.to_numpy()]}
for name,arr in expected.items():
 for role,b in zip(['myeloid','fibroblast','CD8'],arr):ck('membership_'+name+role,np.array_equal(m[name+'_'+role],b))
old=pd.read_csv(O/'tables/PILOT_CELL_METADATA.tsv',sep='\t',dtype={'cell_id':str}).set_index('cell_id').loc[a.obs_names]
common=(old.nuclear_RNA>=5).to_numpy()&np.isfinite(old[['nucleus_x_proxy_um','nucleus_y_proxy_um']]).all(axis=1).to_numpy()
ck('common_geometry_mask',np.array_equal(common,m.common_geometry))
f=m.loc[common].reset_index(drop=True);idmap={x:i for i,x in enumerate(f.original_cell_id)}
geo={'Proseg_centroid':a.obsm['spatial'][common],'Original_nuclear_centroid_proxy':old.loc[common,['nucleus_x_proxy_um','nucleus_y_proxy_um']].to_numpy()}
og=pd.read_csv(O/'tables/PILOT_FEATURE_AUDIT.tsv',sep='\t');sj=int(og.loc[og.gene.eq('SPP1'),'matrix_col'].iloc[0])
A=sparse.load_npz(O/'resources/pilot_assigned_counts_qv20.npz')[old.matrix_row.to_numpy(int)][common]
N=sparse.load_npz(O/'resources/pilot_nuclear_counts_qv20.npz')[old.matrix_row.to_numpy(int)][common]
assays={'Proseg':(C[common,gn.index('SPP1')],C[common].sum(axis=1)),'OriginalAssigned':(A[:,sj].toarray().ravel(),np.asarray(A.sum(axis=1)).ravel()),'OriginalNuclear':(N[:,sj].toarray().ravel(),np.asarray(N.sum(axis=1)).ravel())}
d=pd.read_csv(P/'tables/ALL_SPATIAL_EFFECTS.tsv',sep='\t');inp=pd.read_csv(P/'tables/ALL_MODEL_INPUTS.tsv.gz',sep='\t',dtype={'original_cell_id':str})
ck('432_models_with216_per_center_set',len(d)==432 and (d.groupby('center_set').size()==216).all())
ck('full_factorial_no_duplicate_keys',not d.duplicated(['condition','center_set','geometry','assay','radius_um','endpoint','model']).any())
ck('no_P_FDR',d[['P','FDR']].isna().all().all())
maxgeo=maxrho=0.;ninputs=0
for geom,xy in geo.items():
 alltree=BallTree(xy);edge=np.minimum.reduce([xy[:,0]-xy[:,0].min(),xy[:,0].max()-xy[:,0],xy[:,1]-xy[:,1].min(),xy[:,1].max()-xy[:,1]])
 for name in expected:
  fm=f[name+'_fibroblast'].to_numpy(bool);tm=f[name+'_CD8'].to_numpy(bool);ft=BallTree(xy[fm]);tt=BallTree(xy[tm])
  for centerset in ['condition_members','matched_strict_centers']:
   mc=f[name+'_myeloid'].to_numpy(bool) if centerset=='condition_members' else f['L0_strict_members_relaxed_counts_myeloid'].to_numpy(bool)
   for radius in [50,100]:
    tag=geom+name+centerset+str(radius);base=inp[(inp.geometry==geom)&(inp.condition==name)&(inp.center_set==centerset)&(inp.radius_um==radius)&(inp.assay=='Proseg')]
    expected_ix=np.flatnonzero(mc&(edge>=radius));pos=xy[expected_ix]
    nc=alltree.query_radius(pos,r=radius,count_only=True)-1;fc=ft.query_radius(pos,r=radius,count_only=True);nd=tt.query(pos,k=1)[0].ravel()
    good=(nc>0)&np.isfinite(nd);ix=expected_ix[good];nc=nc[good];fc=fc[good];nd=nd[good];pos=pos[good];frac=fc/nc
    ck('all_centers_'+tag,np.array_equal(base.original_cell_id.to_numpy(),f.original_cell_id.to_numpy()[ix]))
    err=max(float(np.max(np.abs(frac-base.fibroblast_fraction))),float(np.max(np.abs(nd-base.nearest_CD8_distance_um))));maxgeo=max(maxgeo,err)
    ck('independent_BallTree_'+tag,err<1e-8 and np.array_equal(nc,base.all_neighbor_count) and np.array_equal(fc,base.fibroblast_neighbor_count))
    for j in range(min(3,len(ix))):
     distances=np.sqrt(((xy-pos[j])**2).sum(axis=1))
     ck('bruteforce_'+tag+str(j),int((distances<=radius).sum())-1==int(nc[j]) and int((distances[fm]<=radius).sum())==int(fc[j]) and np.isclose(distances[tm].min(),nd[j],atol=1e-9))
    for assay,(sp,lib) in assays.items():
     sub=inp[(inp.geometry==geom)&(inp.condition==name)&(inp.center_set==centerset)&(inp.radius_um==radius)&(inp.assay==assay)];ninputs+=len(sub)
     ck('assay_counts_'+tag+assay,np.array_equal(sub.original_cell_id,base.original_cell_id) and np.array_equal(sub.SPP1_count,sp[ix]) and np.array_equal(sub.library_counts,lib[ix]))
     expr=np.log1p(sp[ix]/np.maximum(lib[ix],1)*1e4);ck('expression_'+tag+assay,np.allclose(expr,sub.SPP1_normalized,atol=1e-12))
     Z=np.column_stack([np.ones(len(ix)),rankdata(lib[ix]),rankdata(nc),rankdata(pos[:,0]),rankdata(pos[:,1])]);U,s,_=np.linalg.svd(Z,full_matrices=False);rk=np.linalg.matrix_rank(Z);U=U[:,:rk]
     for ep,x,y in [('SPP1_vs_fibroblast_fraction',expr,frac),('SPP1_vs_nearest_CD8_distance',expr,nd),('fibroblast_fraction_vs_nearest_CD8_distance',frac,nd)]:
      for model in ['unadjusted','depth_density_xy']:
       r=d[(d.geometry==geom)&(d.condition==name)&(d.center_set==centerset)&(d.radius_um==radius)&(d.assay==assay)&(d.endpoint==ep)&(d.model==model)].iloc[0]
       ck('coverage_'+tag+assay+ep+model,int(r.centers)==len(ix) and int(r.SPP1_detected)==int((sp[ix]>0).sum()))
       if pd.isna(r.rho):ck('NA_reason_'+tag+assay+ep+model,pd.notna(r.not_assessable_reason));continue
       rx,ry=rankdata(x),rankdata(y)
       if model!='unadjusted':rx=rx-U@(U.T@rx);ry=ry-U@(U.T@ry)
       corr=float(np.corrcoef(rx,ry)[0,1]);error=abs(corr-r.rho);maxrho=max(maxrho,error);ck('SVD_effect_'+tag+assay+ep+model,error<1e-8)
ck('all_input_rows_used',ninputs==len(inp))
historical=pd.read_csv(Q/'tables/REGIONAL_ALL_SPATIAL_EFFECTS.tsv',sep='\t')
j=d[(d.condition=='L0_strict_members_relaxed_counts')&(d.center_set=='condition_members')].merge(historical,on=['geometry','assay','radius_um','endpoint','model'],suffixes=('_new','_old'))
ck('strict_old56_values_preserved',j.rho_old.notna().sum()==56 and (j.rho_new-j.rho_old).abs().max()<1e-12)
ck('strict16_gate_releases',int((j.rho_old.isna()&j.rho_new.notna()).sum())==16)
direction=pd.read_csv(P/'tables/COMPLETE_DIRECTIONAL_MATCH.tsv',sep='\t');eps=['SPP1_vs_fibroblast_fraction','SPP1_vs_nearest_CD8_distance','fibroblast_fraction_vs_nearest_CD8_distance']
ck('directional_match_not_significance',np.array_equal(direction.all_three_positive,direction[eps].notna().all(axis=1)&(direction[eps]>0).all(axis=1)))
rec={'status':'PASS','time':datetime.datetime.now().astimezone().isoformat(),'named_checks':len(checks),'checks':checks,'models_verified':len(d),'model_input_rows_verified':ninputs,'max_geometry_error':maxgeo,'max_rho_error':maxrho,'historical_strict56_unchanged':True,'sample_unit':'one named region; no inferential statistics','method':'frozen cluster IDs and raw matrices; independent BallTree + exact Euclidean checks; SVD residual correlation; coverage/IDs/complete factor family'}
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(rec,indent=2));print(json.dumps({k:v for k,v in rec.items() if k!='checks'},indent=2),flush=True)
