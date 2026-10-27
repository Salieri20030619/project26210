from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree
from skimage.transform import resize
P=Path(__file__).resolve().parents[1];R=P.parent
B=Path('/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial/TJ_RCC_figshare_24599295/extracted/ST_visium')
def normalize(mask,size=128):
 yy,xx=np.where(mask)
 crop=mask[yy.min():yy.max()+1,xx.min():xx.max()+1]
 scale=min((size-16)/crop.shape[0],(size-16)/crop.shape[1])
 shape=tuple(max(1,int(round(v*scale))) for v in crop.shape)
 scaled=resize(crop.astype(float),shape,order=0,preserve_range=True,anti_aliasing=False)>0.5
 out=np.zeros((size,size),bool);y=(size-shape[0])//2;x=(size-shape[1])//2;out[y:y+shape[0],x:x+shape[1]]=scaled
 return out
def main():
 # Image/coordinate feasibility only; never choose registration using TYMP/SPP1/metabolite associations.
 plan=dict(method='coarse morphology-only orientation diagnostic, positive TIC mask vs independent Visium in_tissue field;8 dihedral orientations, translations-8..8 grid step2 after isotropic bbox fit',acceptance='NO molecular association accepted from silhouette alone; independent landmarks or transform metadata needed',mass_or_RNA_target_used_for_alignment=False,physical_SM_pixel_size_unverified=True,new_figures=0)
 (P/'metadata/PAIRING_DIAGNOSTIC_PLAN.json').write_text(json.dumps(plan,indent=2))
 allrows=[];pairrows=[];sources=[]
 for sample in ['R114_T','Y12_T']:
  raw=B/sample/'outs/spatial'
  posfile=raw/'tissue_positions_list.csv'
  pos=pd.read_csv(posfile,header=None,names=['barcode','in_tissue','array_row','array_col','px_row','px_col'])
  xy=pos[['px_col','px_row']].to_numpy();spacing=float(np.median(cKDTree(xy).query(xy,k=2)[0][:,1]))
  low=xy.min(axis=0)-spacing;hi=xy.max(axis=0)+spacing;gridx=np.linspace(low[0],hi[0],256);gridy=np.linspace(low[1],hi[1],256)
  gy,gx=np.meshgrid(gridy,gridx,indexing='ij');query=np.column_stack([gx.ravel(),gy.ravel()])
  dist,ix=cKDTree(xy).query(query);st=((pos.in_tissue.to_numpy()[ix]==1)&(dist<.66*spacing)).reshape(256,256)
  st=ndimage.binary_closing(st,iterations=2);target=normalize(st)
  frozen=R/'phase234_part4_replan_execution_20260921_v1/inputs'/(sample+'_INPUTS.tsv.gz')
  d=pd.read_csv(frozen,sep='\t',usecols=['spot_id'])
  assert d.spot_id.isin(pos.barcode).all()
  for f in [posfile,raw/'scalefactors_json.json',raw/'tissue_hires_image.png',frozen]:
   ss=f.stat();sources.append(dict(path=str(f),bytes=ss.st_size,mtime_ns=ss.st_mtime_ns,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
  for pol in ['neg','pos']:
   z=np.load(P/'results/matrices'/(sample.replace('_','-')+'-'+pol+'.npz'));coord=z['xy'];tic=z['tic'];h=int(coord[:,1].max());w=int(coord[:,0].max())
   mask=np.zeros((h,w),bool);mask[coord[:,1]-1,coord[:,0]-1]=tic>0;source=normalize(mask)
   opts=[]
   for flip in [False,True]:
    for k in range(4):
     img=np.rot90(np.fliplr(source) if flip else source,k)
     for dy in range(-8,9,2):
      for dx in range(-8,9,2):
       moved=ndimage.shift(img.astype(np.uint8),(dy,dx),order=0,mode='constant')>0
       dice=2*np.count_nonzero(moved&target)/(moved.sum()+target.sum())
       opts.append(dict(sample=sample,polarity=pol,flip=flip,rotation90=k,translation_y=dy,translation_x=dx,dice=float(dice)))
   opt=pd.DataFrame(opts);allrows.extend(opts);best=opt.sort_values('dice',ascending=False).iloc[0]
   perorientation=opt.groupby(['flip','rotation90']).dice.max().sort_values(ascending=False)
   pairrows.append(dict(sample=sample,polarity=pol,ST_QC_spots=len(d),ST_original_in_tissue=int(pos.in_tissue.sum()),SM_nonzero_pixels=int(mask.sum()),name_match=True,ST_image_available=True,SM_image_available=False,SM_grid_width=w,SM_grid_height=h,SM_coordinate_unit='unverified grid',best_mask_dice=float(best.dice),orientation_margin=float(perorientation.iloc[0]-perorientation.iloc[1]),best_flip=bool(best.flip),best_rotation90=int(best.rotation90),mapping_accepted=False,reason='No supplied SM-to-ST transform/independent anatomical landmarks; shape fit alone is not verified registration'))
  print('PAIRING',sample,flush=True)
 pd.DataFrame(allrows).to_csv(P/'tables/PAIRING_ORIENTATION_DIAGNOSTIC.tsv',sep='\t',index=False)
 pd.DataFrame(pairrows).to_csv(P/'tables/PAIRING_FEASIBILITY.tsv',sep='\t',index=False)
 (P/'audit/READONLY_ST_SOURCES.json').write_text(json.dumps(sources,indent=2))
 print(pd.DataFrame(pairrows).to_string(index=False),flush=True)
if __name__=='__main__':main()
