from pathlib import Path
import json
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase202_candidate_labels_comparison_20260918_v6'
base=(r/'phase192_sections28_212_figures_20260917_v1/scripts/render_figures.py').read_text().split('# 01 actual biochemical relations only')[0]
exec(compile(base,'frozen_density_setup','exec'))
out=r/'phase202_candidate_labels_comparison_20260918_v6'
loccolors={'4':'#993F80','10':'#146A98','15':'#397348'}
locstyles={'4':'--','10':'-.','15':':'}
centers={};grids={}
for c in ['4','10','15']:
 hh,ee,_=density(x[m.candidate==c],y[m.candidate==c],bounds=bounds)
 iy,ix=np.unravel_index(hh.argmax(),hh.shape)
 centers[c]=(ee[0]+ix/(hh.shape[1]-1)*(ee[1]-ee[0]),ee[2]+iy/(hh.shape[0]-1)*(ee[3]-ee[2]))
 grids[c]=(hh,ee)
def labels(ax,outline_all=True):
 for c in ['4','10','15']:
  hh,ee=grids[c]
  if outline_all:ax.contour(hh,levels=[hh.max()*.12],extent=ee,colors=[loccolors[c]],linewidths=1.2,linestyles=locstyles[c],zorder=6)
  xx,yy=centers[c]
  ax.text(xx,yy,'C'+c,color=loccolors[c],fontsize=8,weight='bold',ha='center',va='center',zorder=7,bbox=dict(fc='white',ec=loccolors[c],lw=.45,alpha=.92,pad=1.2))
# Fig2D: exact original atlas and C4 filled density, adding two candidate boundaries.
fig=plt.figure(figsize=(7.48,5.2));title(fig,'2.9','SPP1-expressing macrophage candidates')
gs=fig.add_gridspec(2,2,height_ratios=[3.2,1],top=.87,bottom=.06,wspace=.12)
for j in [0,1]:
 ax=fig.add_subplot(gs[0,j]);ax.contourf(base,levels=[base.max()*.006,base.max()],extent=ext,colors=['#E7E8E8'])
 if j==0:
  for c,col in zip(reg['clusters'],pal):
   mask=(m.candidate==c).values;hh,ee,_=density(x[mask],y[mask],bounds=bounds)
   ax.contourf(hh,levels=[hh.max()*.16,hh.max()*.48,hh.max()+1e-5],extent=ee,colors=[to_rgba(col,.40),to_rgba(col,.85)])
   iy,ix=np.unravel_index(hh.argmax(),hh.shape)
   xx=ee[0]+ix/(hh.shape[1]-1)*(ee[1]-ee[0]);yy=ee[2]+iy/(hh.shape[0]-1)*(ee[3]-ee[2])
   ax.text(xx,yy,'C'+c,fontsize=7,weight='bold',ha='center',color=loccolors.get(c,'#222222'),bbox=dict(fc='white',ec=loccolors.get(c,'none'),lw=.45,alpha=.88,pad=.8))
 else:
  hh,ee=grids['4'];ax.contourf(hh,levels=np.linspace(hh.max()*.03,hh.max()+1e-5,8),extent=ee,cmap=seq)
  labels(ax)
 ax.set_title(['Frozen macrophage states','C4 / C10 / C15 locations'][j],loc='left',fontsize=9);clean(ax)
ax=fig.add_subplot(gs[1,:]);ax.axis('off')
for k,(c,col) in enumerate(zip(reg['clusters'],pal)):
 xx=(k%3)/3;yy=.94-(k//3)*.25
 ax.add_patch(Rectangle((xx,yy-.07),.017,.08,transform=ax.transAxes,fc=col,ec='none'))
 ax.text(xx+.025,yy-.03,'C'+c+'  '+reg['labels'][c],transform=ax.transAxes,va='center',fontsize=6.8)
fig.savefig(out/'panels/candidate_locations.pdf');plt.close(fig)
# Fig3B: reconstruct the frozen expression surfaces exactly, add consistent candidate locations.
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
cmap=LinearSegmentedColormap.from_list('contourcontrast',['#DBEAF1','#92BBC9','#F8F5E9','#E7B47C','#B95534'])
fig,axs=plt.subplots(2,2,figsize=(7.8,6.65));title(fig,'2.11','Functional expression with candidate locations')
fields=['expr:SPP1','score:Cholesterol metabolism','score:Lysosome','score:Antigen processing and presentation']
audit=[];oldmat=np.load(r/'phase193_section211_contours_20260917_v2/tables/UNCHANGED_LOCAL_MEANS.npz')
for ax,f in zip(axs.ravel(),fields):
 hh,ee,_=density(x,y,weights=m[f].values,bounds=bounds)
 v=np.ma.masked_where(base<base.max()*.008,np.divide(hh,base,out=np.zeros_like(hh),where=base>0))
 assert np.allclose(v.filled(np.nan),oldmat[f],equal_nan=True,atol=1e-12)
 lo,hi=np.percentile(v.compressed(),[2,98]);levels=np.linspace(lo,hi,8)
 xx=np.linspace(ee[0],ee[1],v.shape[1]);yy=np.linspace(ee[2],ee[3],v.shape[0])
 im=ax.contourf(xx,yy,v,levels=levels,cmap=cmap,extend='both',antialiased=True)
 ax.contour(xx,yy,v,levels=levels[1:-1],colors='#5C6566',linewidths=.32,alpha=.65)
 labels(ax);clean(ax)
 ax.set_title(f.split(':')[1].replace('Antigen processing and presentation','Antigen presentation'),loc='left',fontsize=9)
 cb=fig.colorbar(im,ax=ax,fraction=.033,pad=.025,shrink=.79,ticks=levels[::2])
 cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'));cb.ax.tick_params(labelsize=6)
 cb.ax.set_title('Local mean\nlog1p(CP10K)',fontsize=6,pad=7)
 audit.append({'feature':f,'original_values_identical':True,'display_limits':[lo,hi]})
handles=[Line2D([0],[0],color=loccolors[c],lw=1.3,ls=locstyles[c],label='C'+c+' location') for c in ['4','10','15']]
handles.append(Line2D([0],[0],color='#5C6566',lw=.65,label='Equal-expression contour'))
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.035),ncol=2,frameon=False,fontsize=8)
fig.text(.5,.013,'Colored outlines mark cluster locations, not expression maxima.',ha='center',fontsize=7,color='#555555')
fig.subplots_adjust(top=.88,left=.025,right=.94,bottom=.18,wspace=.19,hspace=.24)
fig.savefig(out/'panels/functional_candidate_locations.pdf');plt.close(fig)
(out/'qc/LOCATOR_VALIDATION.json').write_text(json.dumps({'frozen_coordinate_sources':str(r/'phase192_sections28_212_figures_20260917_v1/inputs/MACRO_DISPLAY_METADATA.tsv.gz'),'boundary_density_fraction':.12,'label_positions':centers,'colors':loccolors,'line_styles':locstyles,'expression_fields':audit,'new_inference':False},indent=2))
print('CANDIDATE_LOCATORS_READY')
