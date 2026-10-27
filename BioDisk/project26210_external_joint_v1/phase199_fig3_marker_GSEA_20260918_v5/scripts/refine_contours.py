from pathlib import Path
import json,subprocess,hashlib
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
old=r/'phase192_sections28_212_figures_20260917_v1'
exec((old/'scripts/render_figures.py').read_text().split('# 01 actual')[0])
d=r/'phase199_fig3_marker_GSEA_20260918_v5'
for sub in ['figures','qc','tables','scripts']:(d/sub).mkdir(parents=True,exist_ok=True)
from matplotlib.lines import Line2D
cmap=LinearSegmentedColormap.from_list('contourcontrast',['#DBEAF1','#92BBC9','#F8F5E9','#E7B47C','#B95534'])
fig,axs=plt.subplots(2,2,figsize=(7.8,6.65))
title(fig,'2.11','Functional expression contours across macrophage states')
fields=['expr:SPP1','score:Cholesterol metabolism','score:Lysosome','score:Antigen processing and presentation']
audit=[];matrices={}
c4h,ce,_=density(x[m.candidate=='4'],y[m.candidate=='4'],bounds=bounds)
for ax,f,letter in zip(axs.ravel(),fields,'ABCD'):
 h,e,_=density(x,y,weights=m[f].values,bounds=bounds)
 v=np.ma.masked_where(base<base.max()*.008,np.divide(h,base,out=np.zeros_like(h),where=base>0))
 lo,hi=np.percentile(v.compressed(),[2,98]);levels=np.linspace(lo,hi,8)
 xx=np.linspace(e[0],e[1],v.shape[1]);yy=np.linspace(e[2],e[3],v.shape[0])
 im=ax.contourf(xx,yy,v,levels=levels,cmap=cmap,extend='both',antialiased=True)
 ax.contour(xx,yy,v,levels=levels[1:-1],colors='#5C6566',linewidths=.32,alpha=.65)
 ax.contour(c4h,levels=[c4h.max()*.12],extent=ce,colors=[purple],linewidths=1.5,linestyles='--',zorder=5)
 clean(ax)
 label=f.split(':')[1].replace('Antigen processing and presentation','Antigen presentation')
 ax.set_title(label,loc='left',fontsize=9)
 cb=fig.colorbar(im,ax=ax,fraction=.033,pad=.025,shrink=.79,ticks=levels[::2])
 from matplotlib.ticker import FormatStrFormatter
 cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'));cb.ax.tick_params(labelsize=6)
 cb.ax.set_title('Local mean\nlog1p(CP10K)',fontsize=6,pad=7)
 audit.append({'feature':f,'lower_display_p2':lo,'upper_display_p98':hi,'levels':levels.tolist(),'visible_min':float(v.min()),'visible_max':float(v.max())})
 matrices[f]=v.filled(np.nan)
fig.legend(handles=[Line2D([0],[0],color=purple,lw=1.5,ls='--',label='C4 boundary'),Line2D([0],[0],color='#5C6566',lw=.65,label='Equal-expression contour')],loc='lower center',bbox_to_anchor=(.5,.055),ncol=2,frameon=False,fontsize=8)
fig.text(.5,.027,'Cool to warm: low to high within each panel; contour intervals are equally spaced.',ha='center',fontsize=7,color='#555555')
fig.subplots_adjust(top=.88,left=.025,right=.94,bottom=.14,wspace=.19,hspace=.24)

fig.savefig(d/'panels/contours_no_internal_letters.pdf')
plt.close(fig)
oldmat=np.load(r/'phase193_section211_contours_20260917_v2/tables/UNCHANGED_LOCAL_MEANS.npz')
assert set(matrices)==set(oldmat.files)
for key,mat in matrices.items():assert np.allclose(mat,oldmat[key],equal_nan=True,atol=1e-12)
(d/'qc/CONTOUR_UNCHANGED_CHECK.json').write_text(json.dumps({'features':4,'local_means_unchanged':True,'only_nested_letters_removed':True}))
print('CONTOURS_VERIFIED')
