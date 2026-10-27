from pathlib import Path
import json,datetime,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import rankdata
P=Path(__file__).resolve().parents[1];Q=P.parent/'phase311_RCC_pathology_and_identity_20261004_v1';O=P.parent/'phase310_RCC_mapped_spatial_inputs_20261004_v1'
protocol=json.loads((P/'resources/RELAXED_EXPLORATION_PROTOCOL.json').read_text())
a=ad.read_h5ad(Q/'resources/proseg_region6_identity.h5ad');C=a.layers['maxpost_counts'].toarray();gn=list(a.var_names)
strict=pd.read_csv(Q/'tables/CELL_IDENTITY_FROZEN.tsv.gz',sep='\t',dtype={'original_cell_id':str})
assert np.array_equal(strict.original_cell_id.to_numpy(),a.obs_names.astype(str))
old=pd.read_csv(O/'tables/PILOT_CELL_METADATA.tsv',sep='\t',dtype={'cell_id':str}).set_index('cell_id').loc[a.obs_names]
common=old.nuclear_RNA.ge(5).to_numpy()&np.isfinite(old[['nucleus_x_proxy_um','nucleus_y_proxy_um']]).all(axis=1).to_numpy()
assert np.array_equal(common,strict.paired_geometry_eligible.to_numpy())
lineage=strict.reviewed_cluster_lineage.to_numpy()
def detect(gs):
 js=[gn.index(g) for g in gs if g in gn];return (C[:,js]>0).sum(axis=1) if js else np.zeros(len(C),int)
mygenes=['CD68','C1QA','C1QB','C1QC','CD163','AIF1','LYZ','FCGR3A'];fbgenes=['COL1A1','COL1A2','COL3A1','LUM','PDGFRA','DPT','FBLN1']
loose_cd8=(lineage=='T')&(detect(['CD8A','CD8B'])>=1)
conditions={'L2_cluster_broad':((lineage=='Myeloid'),(lineage=='Fibroblast'),loose_cd8),'L1_single_marker_supported':((lineage=='Myeloid')&(detect(mygenes)>=1),(lineage=='Fibroblast')&(detect(fbgenes)>=1),loose_cd8),'L0_strict_members_relaxed_counts':tuple(strict[x].to_numpy(bool) for x in ['mature_myeloid','strict_fibroblast','strict_CD8'])}
for i in range(3):
 assert np.all(~conditions['L0_strict_members_relaxed_counts'][i]|conditions['L1_single_marker_supported'][i])
 assert np.all(~conditions['L1_single_marker_supported'][i]|conditions['L2_cluster_broad'][i])
members=strict[['original_cell_id','reviewed_cluster_lineage']].copy();members['common_geometry']=common
cover=[];flags=[]
spidx=gn.index('SPP1')
for name,(m,f,t) in conditions.items():
 for role,x in zip(['myeloid','fibroblast','CD8'],[m,f,t]):members[name+'_'+role]=x
 cover.append({'condition':name,'QC_myeloid':int(m.sum()),'QC_fibroblast':int(f.sum()),'QC_CD8':int(t.sum()),'common_myeloid':int((m&common).sum()),'common_fibroblast':int((f&common).sum()),'common_CD8':int((t&common).sum()),'QC_myeloid_SPP1_detected':int((m&(C[:,spidx]>0)).sum()),'common_myeloid_SPP1_detected':int((m&common&(C[:,spidx]>0)).sum())})
 for role,x in zip(['myeloid','fibroblast','CD8'],[m,f,t]):
  for flag,v in {'MYH11_or_ACTG2':detect(['MYH11','ACTG2'])>0,'two_T_markers':detect(['CD3D','CD3E','CD3G','TRAC','CD247'])>=2,'two_myeloid_markers':detect(mygenes)>=2,'two_B_markers':detect(['MS4A1','CD79A','CD79B','CD19'])>=2,'two_epithelial_markers':detect(['EPCAM','KRT8','KRT7','KRT18','KRT19','CDH1'])>=2}.items():
   flags.append({'condition':name,'role':role,'flag':flag,'population':int(x.sum()),'flag_cells':int((x&v).sum()),'fraction':float(v[x].mean()) if x.any() else None})
members.to_csv(P/'tables/ALL_CONDITION_MEMBERS.tsv.gz',sep='\t',index=False,compression='gzip')
pd.DataFrame(cover).to_csv(P/'tables/COVERAGE_BY_CONDITION.tsv',sep='\t',index=False)
pd.DataFrame(flags).to_csv(P/'tables/LINEAGE_FLAGS_BY_CONDITION.tsv',sep='\t',index=False)
og=pd.read_csv(O/'tables/PILOT_FEATURE_AUDIT.tsv',sep='\t');s0=int(og.loc[og.gene.eq('SPP1'),'matrix_col'].iloc[0])
A=sparse.load_npz(O/'resources/pilot_assigned_counts_qv20.npz')[old.matrix_row.to_numpy(int)]
N=sparse.load_npz(O/'resources/pilot_nuclear_counts_qv20.npz')[old.matrix_row.to_numpy(int)]
assays={'Proseg':(C[common,spidx],C[common].sum(axis=1)),'OriginalAssigned':(A[common,s0].toarray().ravel(),np.asarray(A[common].sum(axis=1)).ravel()),'OriginalNuclear':(N[common,s0].toarray().ravel(),np.asarray(N[common].sum(axis=1)).ravel())}
geos={'Proseg_centroid':a.obsm['spatial'][common],'Original_nuclear_centroid_proxy':old.loc[common,['nucleus_x_proxy_um','nucleus_y_proxy_um']].to_numpy()}
ids=a.obs_names.astype(str).to_numpy()[common];strict_centers=conditions['L0_strict_members_relaxed_counts'][0][common]
results=[];inputs=[];contrasts=[]
def correlation(x,y,Z=None):
 if len(x)<30 or np.unique(x).size<2 or np.unique(y).size<2:return np.nan,'insufficient_variation_or_centers',None
 x=rankdata(x);y=rankdata(y);rk=1
 if Z is not None:
  Z=np.column_stack([np.ones(len(x))]+[rankdata(z) for z in Z.T]);rk=int(np.linalg.matrix_rank(Z))
  if rk<Z.shape[1]:return np.nan,'rank_deficient_controls',rk
  x=x-Z@np.linalg.lstsq(Z,x,rcond=None)[0];y=y-Z@np.linalg.lstsq(Z,y,rcond=None)[0]
 x=x-x.mean();y=y-y.mean();d=np.linalg.norm(x)*np.linalg.norm(y)
 if d<1e-10:return np.nan,'zero_residual_variance',rk
 return float(x@y/d),'',rk
for geo,xy in geos.items():
 total_tree=cKDTree(xy)
 edge=np.minimum.reduce([xy[:,0]-xy[:,0].min(),xy[:,0].max()-xy[:,0],xy[:,1]-xy[:,1].min(),xy[:,1].max()-xy[:,1]])
 for name,(mf,ff,tf) in conditions.items():
  m=mf[common];f=ff[common];t=tf[common];ft=cKDTree(xy[f]);tt=cKDTree(xy[t])
  for center_set in ['condition_members','matched_strict_centers']:
   centers=m if center_set=='condition_members' else strict_centers
   for rad in [50,100]:
    ix=np.flatnonzero(centers&(edge>=rad));pos=xy[ix]
    all_n=total_tree.query_ball_point(pos,rad,return_length=True)-1
    fib_n=ft.query_ball_point(pos,rad,return_length=True);cd8_dist=tt.query(pos)[0]
    valid=(all_n>0)&np.isfinite(cd8_dist);ix=ix[valid];pos=pos[valid];all_n=all_n[valid];fib_n=fib_n[valid];cd8_dist=cd8_dist[valid];fibfrac=fib_n/all_n
    structural=''
    if len(ix)<30:structural+='centers_lt30;'
    if int(f.sum())<20:structural+='fibroblasts_lt20;'
    if int(t.sum())<20:structural+='CD8_lt20;'
    for assay,(sp,lib) in assays.items():
     expr=np.log1p(sp[ix]/np.maximum(lib[ix],1)*1e4);det=sp[ix]>0
     sig_reason=structural
     if det.sum()<10:sig_reason+='SPP1_detected_lt10;'
     if (~det).sum()<10:sig_reason+='SPP1_undetected_lt10;'
     z=np.column_stack([lib[ix],all_n,pos])
     base={'condition':name,'center_set':center_set,'geometry':geo,'assay':assay,'radius_um':rad,'centers':len(ix),'SPP1_detected':int(det.sum()),'SPP1_undetected':int((~det).sum()),'fibroblast_reference_cells':int(f.sum()),'CD8_reference_cells':int(t.sum())}
     for endpoint,x,y in [('SPP1_vs_fibroblast_fraction',expr,fibfrac),('SPP1_vs_nearest_CD8_distance',expr,cd8_dist),('fibroblast_fraction_vs_nearest_CD8_distance',fibfrac,cd8_dist)]:
      reason=sig_reason if endpoint.startswith('SPP1') else structural
      for model in ['unadjusted','depth_density_xy']:
       if reason:rho,why,rk=np.nan,reason,None
       else:rho,why,rk=correlation(x,y,z if model!='unadjusted' else None)
       results.append({**base,'endpoint':endpoint,'model':model,'rho':rho,'P':None,'FDR':None,'design_rank':rk,'not_assessable_reason':why,'evidence':'post-hoc descriptive expansion; one named region'})
     contrasts.append({**base,'group_comparison_eligible':not bool(sig_reason),'median_fibro_fraction_SPP1_detected':float(np.median(fibfrac[det])) if det.any() else None,'median_fibro_fraction_SPP1_undetected':float(np.median(fibfrac[~det])) if (~det).any() else None,'median_CD8_distance_SPP1_detected':float(np.median(cd8_dist[det])) if det.any() else None,'median_CD8_distance_SPP1_undetected':float(np.median(cd8_dist[~det])) if (~det).any() else None})
     for k,j in enumerate(ix):inputs.append({'condition':name,'center_set':center_set,'geometry':geo,'assay':assay,'radius_um':rad,'original_cell_id':ids[j],'SPP1_count':float(sp[j]),'SPP1_normalized':float(expr[k]),'library_counts':float(lib[j]),'all_neighbor_count':int(all_n[k]),'fibroblast_neighbor_count':int(fib_n[k]),'fibroblast_fraction':float(fibfrac[k]),'nearest_CD8_distance_um':float(cd8_dist[k]),'x':float(pos[k,0]),'y':float(pos[k,1])})
d=pd.DataFrame(results);d.to_csv(P/'tables/ALL_SPATIAL_EFFECTS.tsv',sep='\t',index=False)
pd.DataFrame(inputs).to_csv(P/'tables/ALL_MODEL_INPUTS.tsv.gz',sep='\t',index=False,compression='gzip')
pd.DataFrame(contrasts).to_csv(P/'tables/ALL_GROUP_CONTRASTS.tsv',sep='\t',index=False)
keys=['condition','center_set','geometry','assay','radius_um','model']
wide=d.pivot(index=keys,columns='endpoint',values='rho').reset_index();endpoints=protocol['endpoints']
wide['all_three_assessable']=wide[endpoints].notna().all(axis=1)
wide['all_three_positive']=wide['all_three_assessable']&(wide[endpoints]>0).all(axis=1)
wide['positive_edges']=(wide[endpoints]>0).sum(axis=1)
wide.to_csv(P/'tables/COMPLETE_DIRECTIONAL_MATCH.tsv',sep='\t',index=False)
primary=d[(d.center_set=='condition_members')&(d.geometry=='Proseg_centroid')&(d.assay=='Proseg')&(d.radius_um==50)]
primary.to_csv(P/'tables/PRIMARY_LOOSE_TO_STRICT_RESULTS.tsv',sep='\t',index=False)
hist=pd.read_csv(Q/'tables/REGIONAL_ALL_SPATIAL_EFFECTS.tsv',sep='\t')
j=d[(d.condition=='L0_strict_members_relaxed_counts')&(d.center_set=='condition_members')].merge(hist,on=['geometry','assay','radius_um','endpoint','model'],suffixes=('_relaxedgate','_historical'))
j.to_csv(P/'tables/STRICT_MEMBER_GATE_CHANGE_COMPARISON.tsv',sep='\t',index=False)
rec={'status':'ANALYSIS_COMPLETE','time':datetime.datetime.now().astimezone().isoformat(),'models_primary_and_sensitivities':len(d[d.center_set=='condition_members']),'matched_center_models':len(d[d.center_set=='matched_strict_centers']),'matched_L0_duplicates_reference':72,'assessable':int(d.rho.notna().sum()),'conditions':cover,'newly_assessable_strict_member_models':int((j.rho_historical.isna()&j.rho_relaxedgate.notna()).sum()),'max_historical_effect_change':float((j.rho_historical-j.rho_relaxedgate).abs().max()),'no_P_FDR':True,'no_new_segmentation':True}
(P/'qc/ANALYSIS_COMPLETE.json').write_text(json.dumps(rec,indent=2))
manifest={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Q/'resources/proseg_region6_identity.h5ad',Q/'tables/CELL_IDENTITY_FROZEN.tsv.gz',P/'resources/RELAXED_EXPLORATION_PROTOCOL.json']}
(P/'qc/INPUT_SHA256.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(rec,indent=2));print(primary[['condition','endpoint','model','centers','SPP1_detected','rho','not_assessable_reason']].to_string(index=False),flush=True)
