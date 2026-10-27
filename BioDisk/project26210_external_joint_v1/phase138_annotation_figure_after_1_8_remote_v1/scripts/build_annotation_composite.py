import io,sys,json,zipfile
import numpy as np,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,to_rgba
from scipy.cluster.hierarchy import linkage,dendrogram
from scipy.ndimage import gaussian_filter,binary_dilation,label as cc_label
from matplotlib import transforms
p='/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad'
a=ad.read_h5ad(p,backed='r');xy=np.asarray(a.obsm['X_umap']);labs=a.obs.cell_type_broad_v1.astype(str).to_numpy()
panels={'Dendritic cell':['CD1C','HLA-DRA','FCER1A','CD74'],'Macrophage':['APOC1','C1QA','C1QB'],'Monocyte':['S100A9','S100A8'],'Mast cell':['TPSAB1','KIT'],'B cell':['MS4A1','CD79A','CD79B'],'NK cell':['NKG7','KLRB1'],'T cell':['CD3D','CD3E','CD3G'],'Plasma cell':['MZB1','JCHAIN'],'Epithelial cell':['EPCAM','KRT8','KRT17','KRT18','KRT19'],'Endothelial':['VWF','PECAM1'],'Pericyte/SMC':['ACTA2','MYH11','ABCC9','PDGFRB','RGS5'],'Fibroblast':['COL1A1','DCN','COL3A1'],'Erythrocyte':['HBB','HBA1']}
types=list(panels)+['Epithelial/Myeloid mixed'];genes=[g for gs in panels.values() for g in gs];assert len(genes)==len(set(genes));assert set(genes)<=set(a.var_names)
x=a[:,genes].X.toarray();assert np.isfinite(x).all() and x.min()>=0 and len(labs)==258120 and set(labs)==set(types)
means=np.array([x[labs==t].mean(axis=0) for t in types]);frac=np.array([(x[labs==t]>0).mean(axis=0) for t in types]);norm=(means-means.min(axis=0))/(np.ptp(means,axis=0)+1e-12)
Z=linkage(norm,method='average',metric='euclidean',optimal_ordering=True);tree=dendrogram(Z,no_plot=True);rows=tree['leaves'];rt=[types[i] for i in rows]
gene_order=[];blocks=[]
for t in rt:
 if t in panels:
  start=len(gene_order);gene_order.extend(genes.index(g) for g in panels[t]);blocks.append((t,start,len(gene_order)-1))
ng=len(genes);nr=len(types);vv=norm[np.ix_(rows,gene_order)];ff=frac[np.ix_(rows,gene_order)];gn=[genes[i] for i in gene_order]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.linewidth':.6})
cmap=LinearSegmentedColormap.from_list('reference',['#EDF3F4','#A6BAC9','#866093','#4B1856'])
groups={'Myeloid':['Macrophage','Monocyte','Dendritic cell','Mast cell'],'Lymphoid':['T cell','NK cell','B cell','Plasma cell'],'Stromal':['Endothelial','Pericyte/SMC','Fibroblast'],'Epithelial':['Epithelial cell']}
colors={'Macrophage':'#DD8992','Monocyte':'#C74753','Dendritic cell':'#9F3342','Mast cell':'#E8B6BE','T cell':'#85518C','NK cell':'#AE7DAF','B cell':'#547E89','Plasma cell':'#C39ABF','Endothelial':'#69B392','Pericyte/SMC':'#CCA051','Fibroblast':'#AA7D41','Epithelial cell':'#8CB963','Erythrocyte':'#AF7764','Epithelial/Myeloid mixed':'#AAAAAA'}

from pathlib import Path
import hashlib
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
w=Path(__file__).resolve().parents[1]
markers=['CLEC9A','C1QA','S100A8','TPSAB1','MS4A1','NKG7','CD3D','MZB1','EPCAM','PECAM1','RGS5','COL1A1','HBB']
marker_types=['Dendritic cell','Macrophage','Monocyte','Mast cell','B cell','NK cell','T cell','Plasma cell','Epithelial cell','Endothelial','Pericyte/SMC','Fibroblast','Erythrocyte']
marker_colors=['#7AAFA0','#D76A73','#E1A078','#AF7694','#8D719D','#A5A76B','#81729D','#BC8CA6','#D4AA62','#6AAE99','#7DA2BF','#A3BC78','#C57B65']
mx=a[:,markers].X.toarray()
scales=[float(np.quantile(mx[:,i][mx[:,i]>0],.99))for i in range(13)]
a.file.close()
np.savez_compressed(w/'source_data/marker_dotplot_values.npz',fraction=ff,colour=vv,genes=np.array(gn),cell_types=np.array(rt),tree=Z)
W,H=190,270
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.45})
fig=plt.figure(figsize=(W/25.4,H/25.4))
def axmm(x,y,ww,hh):return fig.add_axes([x/W,1-(y+hh)/H,ww/W,hh/H])
def text(x,y,t,**kw):
 args=dict(fontsize=6.5,ha='left',va='top');args.update(kw);return fig.text(x/W,1-y/H,t,**args)
def title(letter,y,t):
 text(4,y,letter,fontsize=11,fontweight='bold');text(10,y+.4,t,fontsize=8.5,fontweight='bold')
title('A',3,'Single-cell identities and major compartments')
ax=axmm(5,11,57,50)
rng=np.random.default_rng(26210);ix=rng.permutation(len(labs));cols=np.array([to_rgba(colors[t])for t in labs])
ax.scatter(xy[ix,0],xy[ix,1],c=cols[ix],s=.05,linewidths=0,rasterized=True,alpha=.78)
xmin,xmax=xy[:,0].min()-1.5,xy[:,0].max()+1.5
ymin,ymax=xy[:,1].min()-1.5,xy[:,1].max()+1.5
edgesx=np.linspace(xmin,xmax,260);edgesy=np.linspace(ymin,ymax,330)
xc=(edgesx[:-1]+edgesx[1:])/2;yc=(edgesy[:-1]+edgesy[1:])/2
contour_info={}
for group,members in groups.items():
 ids=np.isin(labs,members);hist=np.histogram2d(xy[ids,0],xy[ids,1],bins=[edgesx,edgesy])[0].T
 smooth=gaussian_filter(hist,2.3);sv=np.sort(smooth.ravel())[::-1]
 thr=sv[np.searchsorted(np.cumsum(sv),smooth.sum()*.975)];mask=smooth>=thr
 cc,ncc=cc_label(mask);kept=np.zeros_like(mask)
 for k in range(1,ncc+1):
  component=cc==k
  if smooth[component].sum()>=max(25,ids.sum()*.001):kept|=component
 soft=gaussian_filter(binary_dilation(kept,iterations=3).astype(float),3.0)
 ax.contour(xc,yc,soft,levels=[.35],colors=['#343434'],linewidths=.45,linestyles=[(0,(2.1,2.1))],zorder=5)
 contour_info[group]={'cells':int(ids.sum()),'threshold':float(thr)}
ax.set_xlim(xmin,xmax);ax.set_ylim(ymin,ymax);ax.set_aspect('equal');ax.axis('off')
for end in [(.22,.04),(.04,.19)]:
 ax.annotate('',xy=end,xytext=(.04,.04),xycoords='axes fraction',arrowprops={'arrowstyle':'-|>','lw':.55,'color':'#333333'})
ax.text(.13,0,'UMAP1',transform=ax.transAxes,fontsize=6.5,ha='center')
ax.text(0,.12,'UMAP2',transform=ax.transAxes,fontsize=6.5,rotation=90,va='center')
text(72,11,'258,120 cells',fontsize=10,fontweight='bold')
text(72,17,'75 libraries | 40 patients | 14 cell identities',fontsize=7)
headers={'Myeloid':'#D98B95','Lymphoid':'#85518C','Stromal':'#69B392','Epithelial':'#8CB963'}
for j,group in enumerate(['Myeloid','Lymphoid','Stromal']):
 left=72+j*38
 fig.add_artist(Rectangle((left/W,1-28/H),33/W,4.5/H,color=headers[group],transform=fig.transFigure))
 text(left+16.5,24.1,group,fontsize=7,fontweight='bold',ha='center',color='white')
 for i,ct in enumerate(groups[group]):
  yy=32+i*5
  fig.add_artist(Line2D([left/W],[1-yy/H],marker='o',markersize=3,c=colors[ct],ls='',transform=fig.transFigure))
  text(left+3,yy-1.2,ct,fontsize=6.5)
 if group=='Stromal':
  fig.add_artist(Rectangle((left/W,1-53/H),33/W,4.5/H,color=headers['Epithelial'],transform=fig.transFigure))
  text(left+16.5,49.1,'Epithelial',fontsize=7,fontweight='bold',ha='center',color='white')
  fig.add_artist(Line2D([left/W],[1-57/H],marker='o',markersize=3,c=colors['Epithelial cell'],ls='',transform=fig.transFigure))
  text(left+3,55.8,'Epithelial cell',fontsize=6.5)
text(72,62,'Other: erythrocyte (brown); epithelial/myeloid mixed (grey).',fontsize=6.5,color='#555555')

# B retains the original 38-marker matrix, gene blocks and original tree.
title('B',67,'Marker expression across cell identities')
ax=axmm(39,95,130,35)
dx=10*np.arange(ng)+5;dy=10*np.arange(nr)+5;xx,yy=np.meshgrid(dx,dy)
ax.scatter(xx.ravel(),yy.ravel(),s=ff.ravel()*25+.08,c=vv.ravel(),cmap=cmap,vmin=0,vmax=1,edgecolors='#5A5A5A',linewidths=.2)
ax.set_xlim(0,ng*10);ax.set_ylim(nr*10,0)
ax.set_xticks(dx,gn,rotation=90,fontsize=6.5);ax.set_yticks(dy,rt,fontsize=6.5);ax.tick_params(length=1.5,pad=1)
tr=transforms.blended_transform_factory(ax.transData,ax.transAxes)
for ct,start,end in blocks:
 ax.plot([start*10+1,start*10+1,(end+1)*10-1,(end+1)*10-1],[1.002,1.04,1.04,1.002],transform=tr,color='#222222',lw=.5,clip_on=False)
 ax.text((start+end+1)*5,1.065,ct,rotation=90,transform=tr,ha='center',va='bottom',fontsize=6.5)
ta=axmm(170,95,13,35)
dendrogram(Z,orientation='right',ax=ta,no_labels=True,color_threshold=0,above_threshold_color='#454545')
ta.set_ylim(ax.get_ylim());ta.axis('off')
for c in ta.collections:c.set_linewidth(.45)
ca=axmm(46,142,35,1.8);ca.imshow(np.linspace(0,1,256)[None,:],cmap=cmap,aspect='auto')
ca.set_xticks([0,127,255],['0','0.5','1'],fontsize=6.5);ca.set_yticks([]);ca.tick_params(length=1.5,pad=1)
text(45,148,'Mean expression (gene-scaled)',fontsize=6.5)
text(109,142,'Expressing cells (%)',fontsize=6.5)
for j,f in enumerate([.2,.4,.6,.8,1]):
 xx=111+j*13
 fig.add_artist(Line2D([xx/W],[1-147/H],marker='o',markersize=np.sqrt(f*25+.08),mfc='#888888',mec='#444444',mew=.2,ls='',transform=fig.transFigure))
 text(xx+2,146,str(int(f*100)),fontsize=6.5)

# C retains the 13 original expression panels; no embedded identity overview.
title('C',154,'Marker-specific expression in the same UMAP')
lims=[xy[:,0].min()-1,xy[:,0].max()+1,xy[:,1].min()-1,xy[:,1].max()+1]
positions=[(j//5,j%5)for j in range(10)]+[(2,1),(2,2),(2,3)]
for j,(row,col) in enumerate(positions):
 left=4+col*37.5;top=162+row*35
 text(left+16.25,top,markers[j],ha='center',fontsize=7.7,fontstyle='italic',fontweight='bold')
 text(left+16.25,top+3.6,marker_types[j],ha='center',fontsize=6.5,color='#555555')
 ax=axmm(left+1,top+7,31,24.5)
 v=mx[:,j];pos=np.flatnonzero(v>0);pos=pos[np.argsort(v[pos],kind='stable')]
 cm=LinearSegmentedColormap.from_list('expr_'+str(j),['#E3E3E3',marker_colors[j]])
 ax.scatter(xy[:,0],xy[:,1],s=.045,c='#DDDDDD',linewidths=0,rasterized=True)
 ax.scatter(xy[pos,0],xy[pos,1],c=v[pos],vmin=0,vmax=scales[j],cmap=cm,s=.085,linewidths=0,rasterized=True)
 ax.set_xlim(lims[:2]);ax.set_ylim(lims[2:]);ax.set_aspect('equal');ax.axis('off')
 ca=ax.inset_axes([.65,.005,.29,.018]);ca.imshow(np.linspace(0,1,128)[None,:],cmap=cm,aspect='auto')
 ca.set_xticks([0,127],['0',f'{scales[j]:.1f}'],fontsize=6.5);ca.set_yticks([]);ca.tick_params(length=0,pad=1);ca.spines[:].set_visible(False)

from matplotlib.patches import FancyArrowPatch
text(7,238,'Shared axes',fontsize=6.5,color='#555555')
for end in [(27,255),(12,244)]:
 fig.add_artist(FancyArrowPatch((12/W,1-255/H),(end[0]/W,1-end[1]/H),transform=fig.transFigure,arrowstyle='-|>',mutation_scale=7,lw=.6,color='#333333'))
text(20,257,'UMAP1',ha='center',fontsize=6.5)
text(8,249,'UMAP2',ha='center',va='center',rotation=90,fontsize=6.5)

fig.canvas.draw()
renderer=fig.canvas.get_renderer()
inactive=set(id(t)for ax in fig.axes if not ax.axison for t in ax.get_xticklabels()+ax.get_yticklabels())
outside=[]
for textart in fig.findobj(match=Text):
 if textart.get_visible() and textart.get_text().strip() and id(textart)not in inactive:
  box=textart.get_window_extent(renderer)
  if box.x0<0 or box.y0<0 or box.x1>fig.bbox.width or box.y1>fig.bbox.height:outside.append({'text':textart.get_text(),'bounds':[float(v)for v in box.bounds]})
for ext in ['pdf','png']:
 path=w/'figures'/('Figure1_annotation_composite_v1.'+ext);fig.savefig(path,dpi=600)
 print('SAVED',path.name,path.stat().st_size,flush=True)
plt.close(fig)
meta={'status':'DRAWN','width_mm':W,'height_mm':H,'cells':len(labs),'cell_type_counts':{ct:int((labs==ct).sum())for ct in types},
 'B_marker_count':len(genes),'B_gene_order':gn,'B_celltype_order':rt,'C_markers':markers,'C_colour_maxima':scales,'A_contours':contour_info,
 'same_UMAP_coordinates':True,'source_coordinate_sha256':hashlib.sha256(xy.tobytes()).hexdigest(),'C_expression_sha256':hashlib.sha256(mx.tobytes()).hexdigest(),
 'B_fraction_sha256':hashlib.sha256(ff.tobytes()).hexdigest(),'B_colour_sha256':hashlib.sha256(vv.tobytes()).hexdigest(),'outside_page_text':outside,'new_inferential_tests':False}
(w/'qc/FIGURE_METADATA.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta),flush=True)
