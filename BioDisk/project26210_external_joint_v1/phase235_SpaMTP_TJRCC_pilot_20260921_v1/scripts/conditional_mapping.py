from pathlib import Path
import json,datetime
import numpy as np,pandas as pd
from scipy import ndimage,stats
from scipy.spatial import cKDTree
P=Path(__file__).resolve().parents[1];R=P.parent
B=Path('/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial/TJ_RCC_figshare_24599295/extracted/ST_visium')
def normparams(mask,size=128):
 yy,xx=np.where(mask);y0,y1=yy.min(),yy.max();x0,x1=xx.min(),xx.max()
 h,w=y1-y0+1,x1-x0+1;scale=min((size-16)/h,(size-16)/w)
 nh,nw=max(1,int(round(h*scale))),max(1,int(round(w*scale)))
 return dict(x0=int(x0),y0=int(y0),sx=nw/w,sy=nh/h,offsetx=(size-nw)//2,offsety=(size-nh)//2)
def normxy(xy,p):
 out=xy.astype(float).copy();out[:,0]=(xy[:,0]-p['x0']+.5)*p['sx']-.5+p['offsetx'];out[:,1]=(xy[:,1]-p['y0']+.5)*p['sy']-.5+p['offsety'];return out
def orient(xy,flip,k,dy,dx):
 out=xy.copy()
 if flip:out[:,0]=127-out[:,0]
 for _ in range(k):out=np.column_stack([out[:,1],127-out[:,0]])
 out+=np.array([dx,dy]);return out
def main():
 rules=dict(time=datetime.datetime.now().astimezone().isoformat(),scope='conditional exploratory association under plausible shape-only mappings;not validated alignment',candidate_registration='best translation within each dihedral orientation;keep orientations whose dice>=best-0.03;no molecular objective',peak_candidates='all frozen5ppm matches;no significance selection',SM_quantity='pixel intensity/TIC * feature_count;nearest nonzero SM pixel within1.5 normalized SM pitch',RNA_targets=['TYMP_RNA','polyamine_mta_context','V1_SPP1_q05','V3_SPP1_macro_share','V4_SPP1_program'],summary='existing frozen4x4 spatial quantile block arithmetic means;>=10 mapped ST spots/block and >=8 blocks;>=50%QC spots mapped,>=50 mapped spots;Spearman effect only',statistics='two cases;no P,FDR,spot bootstrap or counting mappings as repeats',registration_interval='minimum/maximum effect across retained orientations, not a confidence interval',new_figures=0)
 (P/'metadata/CONDITIONAL_ASSOCIATION_FREEZE.json').write_text(json.dumps(rules,indent=2))
 poses=pd.read_csv(P/'tables/PAIRING_ORIENTATION_DIAGNOSTIC.tsv',sep='\t')
 cand=pd.read_csv(P/'tables/TARGETED_MASS_CANDIDATES_PASS5PPM.tsv',sep='\t')
 effects=[];coverage=[];blocks=[];mapping_tables=[];numchecks=[]
 for sample in ['R114_T','Y12_T']:
  pos=pd.read_csv(B/sample/'outs/spatial/tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','px_row','px_col'])
  xy=pos[['px_col','px_row']].to_numpy();spacing=np.median(cKDTree(xy).query(xy,k=2)[0][:,1])
  low=xy.min(axis=0)-spacing;hi=xy.max(axis=0)+spacing
  gy,gx=np.meshgrid(np.linspace(low[1],hi[1],256),np.linspace(low[0],hi[0],256),indexing='ij')
  dist,ix=cKDTree(xy).query(np.column_stack([gx.ravel(),gy.ravel()]))
  maskST=((pos.in_tissue.to_numpy()[ix]==1)&(dist<.66*spacing)).reshape(256,256)
  maskST=ndimage.binary_closing(maskST,iterations=2);pn=normparams(maskST)
  d=pd.read_csv(R/'phase234_part4_replan_execution_20260921_v1/inputs'/(sample+'_INPUTS.tsv.gz'),sep='\t')
  points=d[['x','y']].to_numpy();STnorm=normxy((points-low)/(hi-low)*255,pn)
  for pol in ['neg','pos']:
   label=sample.replace('_','-')+'-'+pol;z=np.load(P/'results/matrices'/(label+'.npz'));X=z['intensity'];tic=z['tic'];xySM=z['xy'];positive=tic>0
   mask=np.zeros((xySM[:,1].max(),xySM[:,0].max()),bool);mask[xySM[:,1]-1,xySM[:,0]-1]=positive
   pm=normparams(mask);SMnorm=normxy(xySM[positive]-1,pm)
   pp=poses[(poses['sample']==sample)&(poses.polarity==pol)]
   best=pp.sort_values('dice',ascending=False).groupby(['flip','rotation90'],as_index=False).first()
   best=best[best.dice>=pp.dice.max()-.03].copy()
   pitch=np.median(cKDTree(SMnorm).query(SMnorm,k=2)[0][:,1])
   for r in best.itertuples():
    transformed=orient(SMnorm,r.flip,int(r.rotation90),r.translation_y,r.translation_x)
    distance,nearest=cKDTree(transformed).query(STnorm);ok=distance<=1.5*pitch
    identity=f'{sample}_{pol}_f{int(r.flip)}_r{r.rotation90}_y{r.translation_y}_x{r.translation_x}'
    coverage.append(dict(sample=sample,polarity=pol,mapping=identity,dice=r.dice,flip=r.flip,rotation90=r.rotation90,translation_y=r.translation_y,translation_x=r.translation_x,mapped_spots=int(ok.sum()),QC_spots=len(d),mapped_fraction=float(ok.mean()),molecular_evidence_used_to_select=False,mapping_validated=False))
    idx=np.flatnonzero(positive)[nearest]
    mp=pd.DataFrame(dict(sample=sample,polarity=pol,mapping=identity,spot_id=d.spot_id,SM_pixel_index=idx,distance_normalized=distance,SM_pitch_normalized=pitch,included=ok))
    mapping_tables.append(mp)
    for c in cand[(cand['sample']==sample)&(cand.polarity==pol)].itertuples():
     norm=X[:,c.feature_index].astype(float)/np.where(tic>0,tic,1)*X.shape[1]
     frame=d.loc[ok,['block_4']+rules['RNA_targets']].copy();frame['SM_candidate_TIC']=norm[idx[ok]]
     grouped=frame.groupby('block_4').mean();n=frame.groupby('block_4').size();grouped=grouped[n>=10]
     for b,values in grouped.iterrows():
      blocks.append(dict(sample=sample,polarity=pol,mapping=identity,molecule=c.molecule,mz=c.observed_mz,block=int(b),spots=int(n.loc[b]),**values.to_dict()))
     for target in rules['RNA_targets']:
      x=grouped.SM_candidate_TIC.to_numpy();y=grouped[target].to_numpy()
      qualified=ok.mean()>=.5 and ok.sum()>=50 and len(grouped)>=8 and len(np.unique(x))>=3 and len(np.unique(y))>=3
      effect=float(stats.spearmanr(x,y).statistic) if qualified else np.nan
      if qualified:
       xr=stats.rankdata(x);yr=stats.rankdata(y);xr-=xr.mean();yr-=yr.mean()
       check=xr@yr/np.sqrt((xr@xr)*(yr@yr));assert abs(check-effect)<1e-12
       numchecks.append(abs(check-effect))
      effects.append(dict(sample=sample,polarity=pol,mapping=identity,molecule=c.molecule,mz=c.observed_mz,RNA_target=target,block_n=len(grouped),mapped_spots=int(ok.sum()),mapped_fraction=float(ok.mean()),rho=effect,reason='' if qualified else 'coverage_or_variation_not_eligible',registration='conditional_shape_only',metabolite_identity='mass_candidate'))
  print('CONDITIONAL',sample,flush=True)
 ee=pd.DataFrame(effects);ee.to_csv(P/'tables/CONDITIONAL_MAPPING_EFFECTS.tsv',sep='\t',index=False)
 pd.DataFrame(coverage).to_csv(P/'tables/CONDITIONAL_MAPPING_COVERAGE.tsv',sep='\t',index=False)
 pd.DataFrame(blocks).to_csv(P/'tables/CONDITIONAL_REGIONAL_VALUES.tsv',sep='\t',index=False)
 pd.concat(mapping_tables,ignore_index=True).to_csv(P/'results/CONDITIONAL_SPOT_PIXEL_MAPS.tsv.gz',sep='\t',index=False)
 summary=[]
 for (sample,pol,molecule,mz,target),g in ee.groupby(['sample','polarity','molecule','mz','RNA_target']):
  v=g.rho.dropna().to_numpy()
  summary.append(dict(sample=sample,polarity=pol,molecule=molecule,mz=mz,RNA_target=target,retained_mappings=len(g),estimable_mappings=len(v),rho_min=float(v.min()) if len(v) else np.nan,rho_max=float(v.max()) if len(v) else np.nan,positive_mappings=int((v>0).sum()),negative_mappings=int((v<0).sum()),direction_stable=(len(v)==len(g) and ((v>0).all() or (v<0).all())),interval_type='mapping-range not CI',evidence='conditional exploratory case;not validated registration'))
 su=pd.DataFrame(summary);su.to_csv(P/'tables/CONDITIONAL_DIRECTION_STABILITY.tsv',sep='\t',index=False)
 (P/'CONDITIONAL_ASSOCIATION_VALIDATION.json').write_text(json.dumps(dict(status='PASS',mapping_hypotheses=len(coverage),section_mapping_effects=len(ee),estimable=int(ee.rho.notna().sum()),independent_rank_checks=len(numchecks),max_error=max(numchecks) if numchecks else None,biological_patients=2,new_P_values=0,new_figures=0,validated_registration=False),indent=2))
 print(su[su.molecule.isin(['Spermidine','Spermine'])].to_string(index=False),flush=True)
if __name__=='__main__':main()
