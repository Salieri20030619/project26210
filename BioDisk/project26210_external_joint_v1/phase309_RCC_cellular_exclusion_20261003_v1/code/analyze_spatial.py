from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
import json,datetime,hashlib
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata
P=Path(__file__).resolve().parents[1]
proto=json.loads((P/'resources/SPATIAL_PROTOCOL_FROZEN.json').read_text())
a=pd.read_csv(P/'tables/CELL_IDENTITY_CANDIDATES.tsv.gz',sep='\t').set_index('cell_id')
xy=a[['x_centroid','y_centroid']].to_numpy();q=a.QC.to_numpy();cd8=a.CD8.to_numpy();fib=a.Fibroblast_candidate.to_numpy();tum=a.Tumor.to_numpy();ca9=a.Tumor_CA9.to_numpy()
assert cd8.sum()>=100 and fib.sum()>=100 and tum.sum()>=1000
end=q&(a.Endothelial_marker_support>=2).to_numpy()&(a.T_marker_support<2).to_numpy()&(a.Myeloid_marker_support<2).to_numpy()&(a.Epithelial_marker_support<2).to_numpy()
mural=q&(a.Mural_marker_support>=2).to_numpy()&(a.Endothelial_marker_support<2).to_numpy()&(a.T_marker_support<2).to_numpy()&(a.Myeloid_marker_support<2).to_numpy()&(a.Epithelial_marker_support<2).to_numpy()
out=a.loc[cd8,['x_centroid','y_centroid','RNA_counts','author_group']].copy();pos=xy[cd8]
def neighbor_counts(mask,radius):return cKDTree(xy[mask]).query_ball_point(pos,radius,return_length=True,workers=4) if mask.any() else np.zeros(len(pos),dtype=int)
for label,mask in [('QC',q),('Fibroblast',fib),('Tumor',tum),('Tumor_CA9',ca9),('Endothelial',end),('Mural',mural)]:
 for radius in ([50,100] if label=='Fibroblast' else [100]):out[f'{label}_count_{radius}um']=neighbor_counts(mask,radius)
for label,mask in [('Tumor',tum),('Tumor_CA9',ca9)]:
 d,j=cKDTree(xy[mask]).query(pos,workers=4);out[f'{label}_distance_um']=d;out[f'{label}_nearest_cell']=a.index.to_numpy()[mask][j]
out['QC_count_100um']-=1
origin=xy[q].min(axis=0);tiles=np.floor((pos-origin)/500).astype(int);out['tile_x']=tiles[:,0];out['tile_y']=tiles[:,1];out['tile_500um']=[f'{i}_{j}' for i,j in tiles]
def residual_rho(data,x,y,controls):
 n=len(data);xx=rankdata(data[x].to_numpy());yy=rankdata(data[y].to_numpy());columns=[np.ones(n)];used=[]
 for c in controls:
  z=rankdata(data[c].to_numpy());sd=z.std()
  if sd>0:columns.append((z-z.mean())/sd);used.append(c)
 design=np.column_stack(columns);rank=np.linalg.matrix_rank(design)
 if n-rank<10 or np.ptp(xx)==0 or np.ptp(yy)==0:return np.nan,'NOT_ENOUGH_RESIDUAL_VARIATION',used,rank
 if rank!=design.shape[1]:return np.nan,'RANK_DEFICIENT',used,rank
 rx=xx-design@np.linalg.lstsq(design,xx,rcond=None)[0];ry=yy-design@np.linalg.lstsq(design,yy,rcond=None)[0]
 if np.linalg.norm(rx)<1e-9 or np.linalg.norm(ry)<1e-9:return np.nan,'ZERO_RESIDUAL',used,rank
 return float(rx@ry/np.linalg.norm(rx)/np.linalg.norm(ry)),'',used,rank
rows=[];blockrows=[];group_rows=[]
for target in ['Tumor','Tumor_CA9']:
 distance=target+'_distance_um';d=out[out[distance]<=200].copy()
 for radius in [50,100]:
  exposure=f'Fibroblast_count_{radius}um';cov=[f'{target}_count_100um','QC_count_100um','Endothelial_count_100um','Mural_count_100um','RNA_counts','x_centroid','y_centroid']
  for model,controls in [('unadjusted',[]),('density_vascular_depth_xy',cov)]:
   rho,reason,used,rank=residual_rho(d,exposure,distance,controls)
   rows.append({'donors':1,'target':target,'fibro_radius_um':radius,'model':model,'CD8_cells':len(d),'rho':rho,'reason':reason,'control_variables':';'.join(used),'design_rank':rank,'P':np.nan,'FDR':np.nan,'unit':'within-one-section descriptive association; cells are not independent samples'})
  for tile,g in d.groupby('tile_500um'):
   if len(g)<30:rho=np.nan;reason='CD8_LT30';used=[];rank=0
   else:rho,reason,used,rank=residual_rho(g,exposure,distance,cov)
   blockrows.append({'target':target,'fibro_radius_um':radius,'tile':tile,'CD8_cells':len(g),'rho':rho,'reason':reason,'design_rank':rank,'P':np.nan,'FDR':np.nan})
  lo,hi=np.quantile(d[exposure],[.25,.75]);eligible=hi>lo
  for label,mask in [('low',d[exposure]<=lo),('high',d[exposure]>=hi)]:
   g=d[mask] if eligible else d.iloc[:0]
   group_rows.append({'target':target,'fibro_radius_um':radius,'group':label,'q25':lo,'q75':hi,'eligible':eligible,'CD8_cells':len(g),'median_cancer_distance_um':g[distance].median(),'within20um_fraction':(g[distance]<=20).mean() if len(g) else np.nan,'within50um_fraction':(g[distance]<=50).mean() if len(g) else np.nan,'median_local_cancer_count_100um':g[f'{target}_count_100um'].median(),'median_local_QC_count_100um':g.QC_count_100um.median(),'meaning':'unadjusted illustrative group contrast; not matched treatment effect'})
out.to_csv(P/'tables/CD8_SPATIAL_READOUTS.tsv.gz',sep='\t')
pd.DataFrame(rows).to_csv(P/'tables/SPATIAL_ASSOCIATIONS_ALL.tsv',sep='\t',index=False)
blocks=pd.DataFrame(blockrows);blocks.to_csv(P/'tables/SPATIAL_BLOCK_EFFECTS_ALL.tsv',sep='\t',index=False)
pd.DataFrame(group_rows).to_csv(P/'tables/FIBRO_NEIGHBORHOOD_CONTRASTS.tsv',sep='\t',index=False)
summary=[]
for (target,radius),g in blocks.groupby(['target','fibro_radius_um']):
 v=g.rho.dropna();summary.append({'target':target,'fibro_radius_um':radius,'eligible_blocks':len(v),'positive_blocks':int((v>0).sum()),'negative_blocks':int((v<0).sum()),'median_rho':v.median(),'donors':1,'P':np.nan,'FDR':np.nan,'unit':'technical spatial partitions, not biological replication'})
pd.DataFrame(summary).to_csv(P/'tables/SPATIAL_BLOCK_SUMMARY.tsv',sep='\t',index=False)
des=[]
for target in ['Tumor','Tumor_CA9']:
 for cohort,mask in [('all_QC_CD8',np.ones(len(out),dtype=bool)),('within200um_cancer',out[target+'_distance_um']<=200)]:
  v=out.loc[mask,target+'_distance_um'];des.append({'target':target,'cohort':cohort,'CD8_cells':len(v),'median_distance_um':v.median(),'q25_um':v.quantile(.25),'q75_um':v.quantile(.75),'within20um_fraction':(v<=20).mean(),'within50um_fraction':(v<=50).mean(),'donors':1})
pd.DataFrame(des).to_csv(P/'tables/CD8_DISTANCE_DESCRIPTIVE.tsv',sep='\t',index=False)
(P/'qc/SPATIAL_ANALYSIS_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','time':datetime.datetime.now().astimezone().isoformat(),'all_CD8':len(out),'donors':1,'association_models':len(rows),'block_models':len(blockrows),'protocol_sha256':hashlib.sha256((P/'resources/SPATIAL_PROTOCOL_FROZEN.json').read_bytes()).hexdigest(),'identity_sha256':hashlib.sha256((P/'tables/CELL_IDENTITY_CANDIDATES.tsv.gz').read_bytes()).hexdigest()},indent=2))
print('ASSOCIATIONS\n',pd.DataFrame(rows).to_string(index=False),flush=True);print('CONTRASTS\n',pd.DataFrame(group_rows).to_string(index=False),flush=True);print('BLOCK_SUMMARY\n',pd.DataFrame(summary).to_string(index=False),flush=True)
