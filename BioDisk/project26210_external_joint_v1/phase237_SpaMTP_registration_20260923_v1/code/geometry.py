from pathlib import Path
import json,hashlib,datetime,math
import numpy as np,pandas as pd
from scipy import ndimage,optimize
from scipy.spatial import cKDTree
from skimage import filters,morphology,measure,transform,color
from PIL import Image
P=Path(__file__).resolve().parents[1];OLD=P.parent/'phase235_SpaMTP_TJRCC_pilot_20260921_v1'
ST=Path('/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial/TJ_RCC_figshare_24599295/extracted/ST_visium')
N=256

def sha(q):return hashlib.sha256(q.read_bytes()).hexdigest()
def apply(x,M):return np.asarray(x)@M[:2,:2].T+M[:2,2]
def normalize_matrix(mask):
 y,x=np.where(mask);s=224/max(x.max()-x.min()+1,y.max()-y.min()+1);cx=(x.max()+x.min())/2;cy=(y.max()+y.min())/2
 return np.array([[s,0,127.5-s*cx],[0,s,127.5-s*cy],[0,0,1.]])
def warp(mask,M,shape=(N,N),order=0):return transform.warp(mask.astype(float),inverse_map=np.linalg.inv(M),output_shape=shape,order=order,preserve_range=True,mode='constant',cval=0)
def largest(mask):
 l,n=ndimage.label(mask)
 if not n:raise ValueError('empty mask')
 sizes=np.bincount(l.ravel());sizes[0]=0;return l==sizes.argmax()
def boundary(mask,maxn=650):
 cs=measure.find_contours(mask.astype(float),.5);q=np.concatenate([c[:,::-1] for c in cs if len(c)>10]);idx=np.linspace(0,len(q)-1,min(maxn,len(q))).astype(int);return q[idx]
def sdt(mask):return ndimage.distance_transform_edt(~mask)-ndimage.distance_transform_edt(mask)
def geometry(source,target):
 sp=boundary(source);tp=boundary(target);return dict(source=source,target=target,sp=sp,tp=tp,sd=sdt(source),td=sdt(target),sc=np.array(ndimage.center_of_mass(source))[::-1],tc=np.array(ndimage.center_of_mass(target))[::-1],scale=np.sqrt(target.sum()/source.sum()))
def matrix(par,g,flip,affine=False):
 angle,logs,dx,dy=par[:4];a=math.radians(angle);rot=np.array([[math.cos(a),-math.sin(a)],[math.sin(a),math.cos(a)]])
 stretch=np.eye(2)
 if affine:
  aniso,shear=par[4:];stretch=np.array([[math.exp(aniso),shear],[0,math.exp(-aniso)]])
 L=rot@stretch@np.diag([-1 if flip else 1,1])*g['scale']*math.exp(logs);t=g['tc']+np.array([dx,dy])-L@g['sc'];M=np.eye(3);M[:2,:2]=L;M[:2,2]=t;return M

def distance_values(M,g,sp=None,tp=None):
 sp=g['sp'] if sp is None else sp;tp=g['tp'] if tp is None else tp
 a=apply(sp,M);b=apply(tp,np.linalg.inv(M));scale=np.sqrt(abs(np.linalg.det(M[:2,:2])))
 da=np.abs(ndimage.map_coordinates(g['td'],[a[:,1],a[:,0]],order=1,mode='constant',cval=256));db=np.abs(ndimage.map_coordinates(g['sd'],[b[:,1],b[:,0]],order=1,mode='constant',cval=256))*scale
 return da,db

def fit(g,flip,k,affine=False,initial=None,train=None):
 center=k*90;bounds=[(center-45,center+45),(-.25,.25),(-32,32),(-32,32)]+([(-.12,.12),(-.12,.12)] if affine else [])
 sp,tp=(g['sp'],g['tp']) if train is None else train
 def fun(v):
  M=matrix(v,g,flip,affine);da,db=distance_values(M,g,sp,tp);return (np.mean(np.minimum(da,64)**1.4)+np.mean(np.minimum(db,64)**1.4))/2
 starts=[initial] if initial is not None else [np.array([center+offset,0,0,0]+([0,0] if affine else []),float) for offset in [-20,0,20]]
 res=[]
 for x in starts:
  rr=optimize.minimize(fun,x,method='Powell',bounds=bounds,options=dict(maxiter=100,xtol=.015,ftol=1e-5));res.append(rr)
 rr=min(res,key=lambda r:r.fun);M=matrix(rr.x,g,flip,affine);w=warp(g['source'],M)>.5;inter=np.count_nonzero(w&g['target']);dice=2*inter/(w.sum()+g['target'].sum());da,db=distance_values(M,g)
 return dict(flip=bool(flip),rotation_seed=k,angle_degrees=float(rr.x[0]),parameters=rr.x.tolist(),matrix=M.tolist(),objective=float(rr.fun),optimizer_success=bool(rr.success),dice=float(dice),surface_mean=float((da.mean()+db.mean())/2),surface_p95=float(np.percentile(np.r_[da,db],95)),determinant=float(np.linalg.det(M[:2,:2])),scale=float(np.sqrt(abs(np.linalg.det(M[:2,:2])))),model='affine' if affine else 'similarity')
