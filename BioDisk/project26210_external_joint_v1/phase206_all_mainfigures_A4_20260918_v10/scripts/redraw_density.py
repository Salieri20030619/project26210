from pathlib import Path
import json
r=Path.cwd();dest=r/'phase206_all_mainfigures_A4_20260918_v10'
s=(r/'phase192_sections28_212_figures_20260917_v1/scripts/render_figures.py').read_text().split('# 01 actual biochemical relations only')[0];exec(compile(s,'frozen_density','exec'))
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
plt.rcParams.update({'font.size':7,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6})
colors={'4':'#993F80','10':'#146A98','15':'#397348'};styles={'4':'--','10':'-.','15':':'}
def mark(ax):
 for c in colors:
  hh,ee,_=density(x[m.candidate==c],y[m.candidate==c],bounds=bounds);ax.contour(hh,levels=[hh.max()*.12],extent=ee,colors=[colors[c]],linestyles=styles[c],linewidths=.65);iy,ix=np.unravel_index(hh.argmax(),hh.shape);xx=ee[0]+ix/(hh.shape[1]-1)*(ee[1]-ee[0]);yy=ee[2]+iy/(hh.shape[0]-1)*(ee[3]-ee[2]);ax.text(xx,yy,'C'+c,ha='center',va='center',fontsize=6,fontweight='bold',color=colors[c],bbox=dict(fc='white',ec='none',alpha=.85,pad=.4))
fig,ax=plt.subplots(figsize=(45/25.4,49/25.4));fig.subplots_adjust(left=.01,right=.99,top=.99,bottom=.1);ax.contourf(base,levels=[base.max()*.006,base.max()],extent=ext,colors=['#E7E8E8'])
for c,col in zip(reg['clusters'],pal):
 mask=(m.candidate==c).values;hh,ee,_=density(x[mask],y[mask],bounds=bounds);ax.contourf(hh,levels=[hh.max()*.16,hh.max()*.48,hh.max()+1e-5],extent=ee,colors=[to_rgba(col,.40),to_rgba(col,.85)])
 if c not in colors:
  iy,ix=np.unravel_index(hh.argmax(),hh.shape);ax.text(ee[0]+ix/(hh.shape[1]-1)*(ee[1]-ee[0]),ee[2]+iy/(hh.shape[0]-1)*(ee[3]-ee[2]),'C'+c,ha='center',fontsize=6,color='#343434',bbox=dict(fc='white',ec='none',alpha=.7,pad=.1))
mark(ax);clean(ax);fig.savefig(dest/'panels/2D.pdf');plt.close(fig)
cmap=LinearSegmentedColormap.from_list('contours',['#DBEAF1','#92BBC9','#F8F5E9','#E7B47C','#B95534']);fig=plt.figure(figsize=(94/25.4,115/25.4));old=np.load(r/'phase193_section211_contours_20260917_v2/tables/UNCHANGED_LOCAL_MEANS.npz');fields=['expr:SPP1','score:Cholesterol metabolism','score:Lysosome','score:Antigen processing and presentation']
for i,f in enumerate(fields):
 row,col=divmod(i,2);ax=fig.add_axes([.015+col*.50,.55-row*.46,.40,.37]);hh,ee,_=density(x,y,weights=m[f].values,bounds=bounds);v=np.ma.masked_where(base<base.max()*.008,np.divide(hh,base,out=np.zeros_like(hh),where=base>0));assert np.allclose(v.filled(np.nan),old[f],equal_nan=True,atol=1e-12);lo,hi=np.percentile(v.compressed(),[2,98]);lev=np.linspace(lo,hi,8);im=ax.contourf(np.linspace(ee[0],ee[1],v.shape[1]),np.linspace(ee[2],ee[3],v.shape[0]),v,levels=lev,cmap=cmap,extend='both');ax.contour(np.linspace(ee[0],ee[1],v.shape[1]),np.linspace(ee[2],ee[3],v.shape[0]),v,levels=lev[1:-1],colors='#5C6566',linewidths=.25,alpha=.65);mark(ax);clean(ax);ax.set_title(['SPP1','Cholesterol','Lysosome','Antigen presentation'][i],fontsize=7,pad=4)
 ca=fig.add_axes([.425+col*.5,.61-row*.46,.014,.23]);cb=fig.colorbar(im,cax=ca,ticks=[lo,hi]);cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'));cb.ax.tick_params(labelsize=6,length=1,pad=1)
fig.legend(handles=[Line2D([],[],color=colors[c],ls=styles[c],lw=.8,label='C'+c) for c in colors],loc='lower center',bbox_to_anchor=(.5,.015),ncol=3,fontsize=6,frameon=False,handlelength=1.5,columnspacing=1);fig.savefig(dest/'panels/3B.pdf');plt.close(fig)
# Frozen trajectory expression curves, only layout changes.
dy=np.load(r/'phase192_sections28_212_figures_20260917_v1/inputs/DYNAMIC_DISPLAY.npz');groups={'Monocyte':['FCN1','VCAN','S100A8','S100A9','LYZ'],'SPP1':['SPP1'],'Lipid':['APOE','LPL','CYP27A1','ABCA1','ABCG1','NPC1','SOAT1'],'Lysosome':['CTSB','CTSD','CTSL','LAMP1','LAMP2','LIPA','ATP6V0D1'],'Antigen':['HLA-DRA','HLA-DPA1','HLA-DPB1','CD74']}
fig,axs=plt.subplots(2,3,figsize=(93/25.4,54/25.4));fig.subplots_adjust(left=.12,right=.98,bottom=.14,top=.89,wspace=.65,hspace=.8)
for ax,(label,genes) in zip(axs.ravel(),groups.items()):
 inds=[list(dy['genes']).index(g) for g in genes]
 for mode,col,ls in [('primary','#993F80','-'),('sensitivity','#398B92','--')]:ax.plot(dy['grid'],dy[mode][inds].mean(0),color=col,ls=ls,lw=.85)
 ax.set_title(label,fontsize=7,pad=3);ax.set_xticks([0,1],['0','1']);ax.yaxis.set_major_locator(MaxNLocator(2));ax.tick_params(labelsize=6,length=2,pad=1)
axs.ravel()[-1].axis('off');axs.ravel()[-1].legend(handles=[Line2D([],[],color='#993F80',label='Primary'),Line2D([],[],color='#398B92',ls='--',label='Sensitivity')],fontsize=6,loc='center',frameon=False,handlelength=1)
fig.text(.49,.01,'Pseudotime rank',ha='center',fontsize=6);fig.savefig(dest/'panels/3E.pdf');plt.close(fig)
print('DENSITY_READY')
