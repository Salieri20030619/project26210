from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,to_rgba
from matplotlib.patches import FancyBboxPatch,Rectangle
from matplotlib.collections import LineCollection
from scipy.ndimage import gaussian_filter
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase192_sections28_212_figures_20260917_v1';q=r/'part2_rebuild_v1'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':10,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
purple='#993F80';teal='#408E91';gold='#BD9753'
cm=LinearSegmentedColormap.from_list('expression',['#759DBF','#FAF9F5',purple])
seq=LinearSegmentedColormap.from_list('density',['#F6F3F7','#BC99BC','#72316C'])
reg=json.loads((d/'inputs/FEATURE_REGISTRY.json').read_text())
m=pd.read_csv(d/'inputs/MACRO_DISPLAY_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str})
pal=['#538B9E','#8BA672',purple,'#AA7998','#BD9753','#9B929F','#50978D','#C88773','#7892BE','#B3A877','#77768E','#B684AC']
def save(fig,name):
 fig.savefig(r/'phase198_fig2_contrast_corr_revision_20260918_v5/panels/C4_no_nested_letters.pdf',dpi=400)
 pass
 plt.close(fig);print('SAVED',name,flush=True)
def title(fig,k,s):fig.suptitle(k+'  |  '+s,x=.035,y=.985,ha='left',fontsize=12,fontweight='bold')
def zrow(a):
 a=np.asarray(a,float);assert np.isfinite(a).all()
 return (a-a.mean(1,keepdims=True))/np.maximum(a.std(1,keepdims=True),1e-8)
def density(x,y,weights=None,bins=180,bounds=None):
 if bounds is None:
  xr=np.ptp(x)*.05;yr=np.ptp(y)*.05;bounds=[[min(x)-xr,max(x)+xr],[min(y)-yr,max(y)+yr]]
 h,xe,ye=np.histogram2d(x,y,bins=bins,range=bounds,weights=weights)
 return gaussian_filter(h.T,2),[(xe[0]+xe[1])/2,(xe[-2]+xe[-1])/2,(ye[0]+ye[1])/2,(ye[-2]+ye[-1])/2],bounds
x=m.UMAP1.values;y=m.UMAP2.values
base,ext,bounds=density(x,y)
def outline(ax,xx=x,yy=y,mask=None,b=bounds,color=purple):
 mask=(m.candidate=='4').values if mask is None else mask
 if mask.sum()<5:return
 h,e,_=density(xx[mask],yy[mask],bounds=b)
 ax.contour(h,levels=[h.max()*.12],extent=e,colors=[color],linewidths=.8,zorder=4)
def clean(ax):ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect('equal');ax.spines[['top','right','left','bottom']].set_visible(False)
fig=plt.figure(figsize=(7.48,5.2));title(fig,'2.9','A distinct SPP1-enriched macrophage state')
gs=fig.add_gridspec(2,2,height_ratios=[3.2,1],top=.87,bottom=.06,wspace=.12)
for j in [0,1]:
 ax=fig.add_subplot(gs[0,j]);ax.contourf(base,levels=[base.max()*.006,base.max()],extent=ext,colors=['#E7E8E8'])
 if j==0:
  for c,col in zip(reg['clusters'],pal):
   mask=(m.candidate==c).values;h,e,_=density(x[mask],y[mask],bounds=bounds)
   ax.contourf(h,levels=[h.max()*.16,h.max()*.48,h.max()+1e-5],extent=e,colors=[to_rgba(col,.40),to_rgba(col,.85)])
   # Label at density mode instead of outlier-prone mean.
   iy,ix=np.unravel_index(h.argmax(),h.shape)
   ax.text(e[0]+ix/(h.shape[1]-1)*(e[1]-e[0]),e[2]+iy/(h.shape[0]-1)*(e[3]-e[2]),'C'+c,fontsize=7,weight='bold',ha='center',bbox=dict(fc='white',ec='none',alpha=.8,pad=.6))
 else:
  mask=(m.candidate=='4').values;h,e,_=density(x[mask],y[mask],bounds=bounds)
  ax.contourf(h,levels=np.linspace(h.max()*.03,h.max()+1e-5,8),extent=e,cmap=seq);outline(ax)
 ax.set_title(['Frozen macrophage states','C4 localization'][j],loc='left',fontsize=9);clean(ax)
ax=fig.add_subplot(gs[1,:]);ax.axis('off')
for k,(c,col) in enumerate(zip(reg['clusters'],pal)):
 xx=(k%3)/3;yy=.94-(k//3)*.25
 ax.add_patch(Rectangle((xx,yy-.07),.017,.08,transform=ax.transAxes,fc=col,ec='none'))
 ax.text(xx+.025,yy-.03,'C'+c+'  '+reg['labels'][c],transform=ax.transAxes,va='center',fontsize=6.8)
save(fig,'02_macrophage_density_atlas')
