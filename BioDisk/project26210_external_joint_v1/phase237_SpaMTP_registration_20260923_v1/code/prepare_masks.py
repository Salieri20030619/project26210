from geometry import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=[];manifest=json.loads((P/'metadata/INPUT_MANIFEST.json').read_text());saved=[]
fig,ax=plt.subplots(2,4,figsize=(14,7))
for i,sample in enumerate(['R114_T','Y12_T']):
 d=ST/sample/'outs/spatial';rgb=np.asarray(Image.open(d/'tissue_hires_image.png').convert('RGB'))/255;sf=json.loads((d/'scalefactors_json.json').read_text());pos=pd.read_csv(d/'tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','px_row','px_col']);xy=pos[['px_col','px_row']].to_numpy()*sf['tissue_hires_scalef'];pitch=np.median(cKDTree(xy).query(xy,k=2)[0][:,1]);pts=xy[pos.in_tissue.to_numpy()==1]
 lo=np.floor(pts.min(0)-pitch).astype(int);hi=np.ceil(pts.max(0)+pitch).astype(int);lo=np.maximum(lo,0);hi=np.minimum(hi,np.array(rgb.shape[:2][::-1])-1)
 roi=np.zeros(rgb.shape[:2],bool);roi[lo[1]:hi[1]+1,lo[0]:hi[0]+1]=True
 od=-np.log(np.clip(rgb,.01,1)).mean(2);threshold=float(filters.threshold_otsu(od[roi]));masks={}
 for name,mult in [('HE',1.),('HE_low',.85),('HE_high',1.15)]:
  mask=(od>threshold*mult)&roi;mask=morphology.binary_closing(mask,morphology.disk(max(1,int(round(pitch*.2)))));mask=ndimage.binary_fill_holes(largest(mask));masks[name]=mask
 T=normalize_matrix(masks['HE']);ny,nx=np.mgrid[:N,:N];original=apply(np.c_[nx.ravel(),ny.ravel()],np.linalg.inv(T));dist,idx=cKDTree(xy).query(original);provider=((pos.in_tissue.to_numpy()[idx]==1)&(dist<.66*pitch)).reshape(N,N);provider=ndimage.binary_fill_holes(largest(morphology.binary_closing(provider,morphology.disk(1))))
 target={name:warp(mask,T)>.5 for name,mask in masks.items()};target['provider']=provider
 np.savez_compressed(P/'inputs'/(sample+'_ST.npz'),HE=target['HE'],HE_low=target['HE_low'],HE_high=target['HE_high'],provider=provider,T=T,rgb=rgb.astype(np.float32),hires_scale=sf['tissue_hires_scalef'],pitch=pitch,xy_hires=xy,barcode=pos.barcode.to_numpy(dtype=str),in_tissue=pos.in_tissue.to_numpy(),optical_density=od.astype(np.float32))
 for j,name in enumerate(['HE','HE_low','HE_high','provider']):ax[i,j].imshow(target[name],cmap='gray');ax[i,j].set_title(sample+' '+name);ax[i,j].axis('off')
 for pol in ['neg','pos']:
  label=sample.replace('_','-')+'-'+pol;z=np.load(OLD/'results/matrices'/(label+'.npz'));xy=z['xy']-1;raw=np.zeros((xy[:,1].max()+1,xy[:,0].max()+1),bool);raw[xy[:,1],xy[:,0]]=z['tic']>0;S=normalize_matrix(raw);source=warp(raw,S)>.5
  erode=ndimage.binary_erosion(raw,iterations=1);dilate=ndimage.binary_dilation(raw,iterations=1)
  np.savez_compressed(P/'inputs'/(label+'_SM.npz'),source=source,erode=warp(erode,S)>.5,dilate=warp(dilate,S)>.5,S=S,raw_mask=raw,xy=xy,tic=z['tic'])
  rows.append(dict(sample=sample,polarity=pol,label=label,SM_mask_pixels=int(raw.sum()),ST_HE_mask_pixels=int(masks['HE'].sum()),ST_original_tissue_spots=int(pos.in_tissue.sum()),HE_otsu_threshold=threshold,ST_spot_pitch_hires=pitch,HE_provider_dice=2*np.count_nonzero(target['HE']&provider)/(target['HE'].sum()+provider.sum()),normalization='isotropic'))
for q in (OLD/'results/SpaMTP_objects').glob('*.rds'):saved.append(dict(path=str(q),bytes=q.stat().st_size,mtime_ns=q.stat().st_mtime_ns))
(P/'metadata/READONLY_RDS_BEFORE.json').write_text(json.dumps(saved,indent=2));pd.DataFrame(rows).to_csv(P/'tables/INPUT_GEOMETRY_QC.tsv',sep='\t',index=False)
fig.tight_layout();fig.savefig(P/'diagnostics/mask_definitions.png',dpi=140);plt.close(fig);print(pd.DataFrame(rows).to_string(index=False))
