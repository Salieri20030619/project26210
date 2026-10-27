from common import *
import h5py,gzip,io
from PIL import Image
from skimage.color import rgb2hed,rgb2hsv,rgb2gray
from skimage.filters import sobel,laplace
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import adjusted_rand_score
def decode(path):
 if '::' in path:
  f,k=path.split('::',1)
  with h5py.File(f,'r') as h:z=h[k][:]
 else:
  with Image.open(io.BytesIO(gzip.open(path,'rb').read()) if path.endswith('.gz') else path) as im:z=np.asarray(im.convert('RGB'))
 z=z[:,:,:3]
 if z.dtype!=np.uint8:
  if np.nanmax(z)<=1:z=np.clip(z*255,0,255)
  z=z.astype(np.uint8)
 return z
def patches(z,size):
 rgb=z.astype(float)/255;hsv=rgb2hsv(rgb);hed=rgb2hed(rgb);gray=rgb2gray(rgb)
 tissue=(gray<.95)&(hsv[:,:,1]>.05);edges=sobel(gray);sharp=np.abs(laplace(gray))
 rows=[]
 for y in range(0,len(z),size):
  for x in range(0,z.shape[1],size):
   s=np.s_[y:min(y+size,len(z)),x:min(x+size,z.shape[1])]
   m=tissue[s];frac=m.mean()
   if m.size<size*size/2 or frac<.5:continue
   vec=[float(frac)]
   for arr in [hed[:,:,0],hed[:,:,1],gray,hsv[:,:,1]]:
    v=arr[s][m];vec += [float(v.mean()),float(np.std(v)),float(np.quantile(v,.9))]
   vec += [float(edges[s][m].mean()),float(sharp[s][m].mean())]
   rows.append([x,y]+vec)
 d=pd.DataFrame(rows,columns=['x','y','tissue_fraction','H_mean','H_sd','H_q90','E_mean','E_sd','E_q90','gray_mean','gray_sd','gray_q90','sat_mean','sat_sd','sat_q90','edge_mean','laplace_abs_mean'])
 if len(d)<30:return d,tissue,None
 cols=d.columns[2:];scaler=StandardScaler();v=scaler.fit_transform(d[cols])
 labels=KMeans(n_clusters=3,random_state=287,n_init=20).fit_predict(v)
 # IDs ordered by H optical density; no biological or molecular labels.
 rank=d.assign(label=labels).groupby('label').H_mean.mean().sort_values().index.tolist()
 mapping={v:i for i,v in enumerate(rank)};d['morphology_class']=[mapping[x] for x in labels]
 grid=np.full(((len(z)+size-1)//size,(z.shape[1]+size-1)//size),-1,dtype=np.int16)
 for r in d.itertuples():grid[r.y//size,r.x//size]=r.morphology_class
 return d,tissue,grid
def main():
 gate('image_morphology_start')
 assets=read(R/'phase286_CAF_identity_boundary_audit_20260929_v1/tables/IMAGE_COORDINATE_AVAILABILITY.tsv')
 qc=[];spots=[];centroids=[];files=[]
 for r in assets.itertuples():
  z=decode(r.image_path);assert z.shape[:2]==(r.height,r.width)
  d=read(R/'phase234_part4_replan_execution_20260921_v1/inputs'/f'{r.sample}_INPUTS.tsv.gz',usecols=['x_hires_px','y_hires_px','spot_id'])
  a,t,grid=patches(z,64);b,_,grid2=patches(z,96)
  for size,g in [(64,a),(96,b)]:
   g.to_csv(P/'tables'/f'{r.sample}_IMAGE_PATCHES_{size}.tsv.gz',sep='\t',index=False)
   if 'morphology_class' in g:
    for cl,sub in g.groupby('morphology_class'):
     rec=dict(sample=r.sample,patch_size=size,morphology_class=cl,patches=len(sub))
     rec.update(sub.drop(columns=['x','y','morphology_class']).mean().to_dict());centroids.append(rec)
  x=d.x_hires_px.to_numpy();y=d.y_hires_px.to_numpy()
  la=np.full(len(d),-1);lb=la.copy()
  if grid is not None:la=grid[(y//64).astype(int),(x//64).astype(int)]
  if grid2 is not None:lb=grid2[(y//96).astype(int),(x//96).astype(int)]
  good=(la>=0)&(lb>=0);ari=adjusted_rand_score(la[good],lb[good]) if good.sum()>=30 else np.nan
  d['sample']=r.sample;d['morphology_class_64']=la;d['morphology_class_96']=lb
  d['center_tissue_mask']=t[np.clip(y.astype(int),0,len(t)-1),np.clip(x.astype(int),0,t.shape[1]-1)]
  spots.append(d)
  unique,num=np.unique(la[la>=0],return_counts=True)
  qc.append(dict(sample=r.sample,dataset=r.dataset,patient=r.patient,width=r.width,height=r.height,decoded=True,image_tissue_fraction=float(t.mean()),patches_64=len(a),patches_96=len(b),spots=len(d),assigned_spots_64=int((la>=0).sum()),spot_centers_on_tissue=int(d.center_tissue_mask.sum()),classes_spot_counts=json.dumps(dict(zip(unique.astype(str),num.tolist()))),cross_scale_common_spots=int(good.sum()),cross_scale_ARI=float(ari),RNA_used=False,independent_image_clustering=True,validated_tumor_stroma_boundary=False,interpretation='unsupervised_stain_texture_classes_not_histological_identity'))
  fn=r.image_path.split('::')[0];st=Path(fn).stat();files.append(dict(sample=r.sample,path=r.image_path,bytes=st.st_size,mtime_ns=st.st_mtime_ns))
  np.savez_compressed(P/'inputs'/f'{r.sample}_IMAGE_REGIONS.npz',grid64=grid,grid96=grid2)
  # Unmodified raw image exports solely for image inspection, not analytical figures.
  if r.sample in ['GLMF1','GLMF2','H46_T','J38_T','GSM5924030','GSM5924031']:
   Image.fromarray(z).save(P/'resources'/f'{r.sample}_SOURCE_IMAGE.png')
  print('IMAGE_DONE',r.sample,len(a),round(float(ari),3),flush=True)
 tab(qc,'IMAGE_MORPHOLOGY_QC.tsv');tab(centroids,'IMAGE_CLASS_FEATURES.tsv')
 pd.concat(spots).to_csv(P/'tables/IMAGE_ONLY_SPOT_CLASSES.tsv.gz',sep='\t',index=False)
 dump(P/'qc/IMAGE_INPUT_MANIFEST.json',files)
 dump(P/'qc/IMAGE_MORPHOLOGY_COMPLETE.json',dict(status='COMPUTED',images=len(qc),spots=sum(r['spots'] for r in qc),RNA_features_used=False,biological_tumor_stroma_labels_assigned=False,new_analytical_figures=0))
 print('IMAGE_MORPHOLOGY_COMPLETE',flush=True)
if __name__=='__main__':main()
