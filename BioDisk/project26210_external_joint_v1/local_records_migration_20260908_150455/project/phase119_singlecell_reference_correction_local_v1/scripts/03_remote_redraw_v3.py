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
buf=io.BytesIO();pack=zipfile.ZipFile(buf,'w',zipfile.ZIP_DEFLATED)
def save(fig,name):
 for ext in ['png','pdf']:
  b=io.BytesIO();fig.savefig(b,format=ext,dpi=300,bbox_inches='tight',pad_inches=.14);pack.writestr(name+'.'+ext,b.getvalue())
 plt.close(fig)
def legends(fig,rect,size=55):
 ax=fig.add_axes(rect);ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1)
 ax.text(.5,.97,'Fraction of cells\nin group (%)',ha='center',va='top',fontsize=9)
 for i,f in enumerate([.2,.4,.6,.8,1]):
  xx=.1+i*.2;ax.scatter([xx],[.65],s=f*size,facecolors='#888888',edgecolors='#444444',linewidths=.3);ax.text(xx,.50,str(int(f*100)),ha='center',fontsize=8)
 ax.text(.5,.29,'Mean expression\nin group',ha='center',fontsize=9)
 ca=ax.inset_axes([.03,.10,.94,.07]);ca.imshow(np.linspace(0,1,256)[None,:],cmap=cmap,aspect='auto');ca.set_xticks([0,127,255],['0.0','0.5','1.0'],fontsize=8);ca.set_yticks([]);ca.tick_params(length=2)
def dot_horizontal():
 fig=plt.figure(figsize=(15,6.4));ax=fig.add_axes([.14,.18,.67,.55]);dx=10*np.arange(ng)+5;dy=10*np.arange(nr)+5;xx,yy=np.meshgrid(dx,dy)
 ax.scatter(xx.ravel(),yy.ravel(),s=ff.ravel()*68+.08,c=vv.ravel(),cmap=cmap,vmin=0,vmax=1,edgecolors='#5A5A5A',linewidths=.25)
 ax.set_xlim(0,ng*10);ax.set_ylim(nr*10,0);ax.set_xticks(dx,gn,rotation=90,fontsize=8);ax.set_yticks(dy,rt,fontsize=9);ax.tick_params(length=2)
 tr=transforms.blended_transform_factory(ax.transData,ax.transAxes)
 for t,s,e in blocks:
  ax.plot([s*10+1,s*10+1,(e+1)*10-1,(e+1)*10-1],[1.002,1.035,1.035,1.002],transform=tr,color='#222222',lw=.8,clip_on=False)
  ax.text((s+e+1)*5,1.055,t,rotation=90,transform=tr,ha='center',va='bottom',fontsize=8)
 ta=fig.add_axes([.812,.18,.065,.55]);dendrogram(Z,orientation='right',ax=ta,no_labels=True,color_threshold=0,above_threshold_color='#454545');ta.set_ylim(ax.get_ylim());ta.axis('off'); [c.set_linewidth(.65) for c in ta.collections]
 legends(fig,[.886,.29,.106,.29],68);fig.text(.14,.035,'Colour scaled within each gene; tree based on average-linkage clustering of cell-type marker profiles.',fontsize=8,color='#555555')
 save(fig,'01_marker_dotplot_tree_horizontal_v3')
def dot_vertical():
 fig=plt.figure(figsize=(10.4,13));ax=fig.add_axes([.30,.16,.51,.69]);dx=10*np.arange(nr)+5;dy=10*np.arange(ng)+5;xx,yy=np.meshgrid(dx,dy)
 ax.scatter(xx.ravel(),yy.ravel(),s=ff.T.ravel()*70+.08,c=vv.T.ravel(),cmap=cmap,vmin=0,vmax=1,edgecolors='#5A5A5A',linewidths=.25)
 ax.set_xlim(0,nr*10);ax.set_ylim(ng*10,0);ax.set_xticks(dx,rt,rotation=90,fontsize=8);ax.set_yticks(dy,gn,fontsize=8);ax.tick_params(length=2)
 ta=fig.add_axes([.30,.853,.51,.10]);dendrogram(Z,orientation='top',ax=ta,no_labels=True,color_threshold=0,above_threshold_color='#454545');ta.set_xlim(ax.get_xlim());ta.axis('off'); [c.set_linewidth(.65) for c in ta.collections]
 tr=transforms.blended_transform_factory(ax.transAxes,ax.transData)
 for t,s,e in blocks:
  ax.plot([-.16,-.19,-.19,-.16],[s*10+1,s*10+1,(e+1)*10-1,(e+1)*10-1],transform=tr,color='#222222',lw=.75,clip_on=False)
  ax.text(-.215,(s+e+1)*5,t.replace(' cell','\ncell').replace('Pericyte/SMC','Pericyte/\nSMC'),transform=tr,ha='right',va='center',fontsize=7)
 legends(fig,[.85,.38,.14,.16],70);save(fig,'02_marker_dotplot_tree_vertical_v3')
groups={'Myeloid':['Macrophage','Monocyte','Dendritic cell','Mast cell'],'Lymphoid':['T cell','NK cell','B cell','Plasma cell'],'Stromal':['Endothelial','Pericyte/SMC','Fibroblast'],'Epithelial':['Epithelial cell']}
colors={'Macrophage':'#DD8992','Monocyte':'#C74753','Dendritic cell':'#9F3342','Mast cell':'#E8B6BE','T cell':'#85518C','NK cell':'#AE7DAF','B cell':'#547E89','Plasma cell':'#C39ABF','Endothelial':'#69B392','Pericyte/SMC':'#CCA051','Fibroblast':'#AA7D41','Epithelial cell':'#8CB963','Erythrocyte':'#AF7764','Epithelial/Myeloid mixed':'#AAAAAA'}
fig=plt.figure(figsize=(7.3,9));ax=fig.add_axes([.08,.31,.84,.65]);rng=np.random.default_rng(26210);ix=rng.permutation(len(labs));cols=np.array([to_rgba(colors[t]) for t in labs]);ax.scatter(xy[ix,0],xy[ix,1],c=cols[ix],s=.28,linewidths=0,rasterized=True,alpha=.78)
xmin,xmax=xy[:,0].min()-1.5,xy[:,0].max()+1.5;ymin,ymax=xy[:,1].min()-1.5,xy[:,1].max()+1.5
edgesx=np.linspace(xmin,xmax,260);edgesy=np.linspace(ymin,ymax,330);xc=(edgesx[:-1]+edgesx[1:])/2;yc=(edgesy[:-1]+edgesy[1:])/2
contour_info={}
for group,members in groups.items():
 ids=np.isin(labs,members);h=np.histogram2d(xy[ids,0],xy[ids,1],bins=[edgesx,edgesy])[0].T;smooth=gaussian_filter(h,2.3);sv=np.sort(smooth.ravel())[::-1];thr=sv[np.searchsorted(np.cumsum(sv),smooth.sum()*.975)];mask=smooth>=thr
 cc,ncc=cc_label(mask);kept=np.zeros_like(mask)
 for k in range(1,ncc+1):
  component=cc==k
  if smooth[component].sum()>=max(25,ids.sum()*.001):kept|=component
 expanded=binary_dilation(kept,iterations=3);soft=gaussian_filter(expanded.astype(float),3.0)
 ax.contour(xc,yc,soft,levels=[.35],colors=['#343434'],linewidths=.8,linestyles=[(0,(2.1,2.1))],zorder=5)
 contour_info[group]={'n_cells':int(ids.sum()),'method':'97.5% smoothed density support; component minimum max(25,0.1% group mass); 3-bin outline padding','members':members}
ax.set_xlim(xmin,xmax);ax.set_ylim(ymin,ymax);ax.set_aspect('equal');ax.axis('off');ax.text(0,1.02,'A',transform=ax.transAxes,fontsize=17,fontweight='bold');ax.text(.10,1.02,'N = 258,120',transform=ax.transAxes,fontstyle='italic',fontsize=10)
# Grouped legends with coloured header bands, following the reference.
headers={'Myeloid':'#D98B95','Lymphoid':'#85518C','Stromal':'#69B392','Epithelial':'#8CB963','Other':'#999999'}
for col,group in enumerate(['Myeloid','Lymphoid','Stromal']):
 la=fig.add_axes([.09+col*.30,.085,.27,.20]);la.set_xlim(0,1);la.set_ylim(0,1);la.axis('off');la.add_patch(plt.Rectangle((0,.84),.94,.16,color=headers[group]));la.text(.47,.92,group,ha='center',va='center',color='white',fontsize=10,fontstyle='italic',fontweight='bold')
 for j,t in enumerate(groups[group]):la.scatter(.04,.72-j*.15,s=19,color=colors[t]);la.text(.12,.72-j*.15,t,va='center',fontsize=8)
 if group=='Stromal':
  la.add_patch(plt.Rectangle((0,.05),.94,.15,color=headers['Epithelial']));la.text(.47,.125,'Epithelial',ha='center',va='center',color='white',fontsize=9,fontweight='bold');la.scatter(.04,-.035,s=19,color=colors['Epithelial cell'],clip_on=False);la.text(.12,-.035,'Epithelial cell',va='center',fontsize=8)
fig.text(.10,.035,'Other: erythrocyte (brown); epithelial/myeloid mixed (grey).',fontsize=8,color='#666666');save(fig,'03_UMAP_major_classes_dashed_v3')
dot_horizontal();dot_vertical()
meta={'status':'PASS_SOURCE_CHECKS','cells':len(labs),'cell_type_counts':{t:int((labs==t).sum()) for t in types},'genes':genes,'tree_row_order_bottom_to_top':rt,'tree_method':'average linkage; Euclidean distance of gene-scaled group means; optimal leaf ordering','dot_size':'fraction expressing X>0; area proportional','dot_colour':'gene-wise min-max group mean log1p expression','umap_contours':contour_info,'notes':'Same frozen UMAP and broad annotation; no cell identifiers or expression matrices exported','versions':{'numpy':np.__version__,'matplotlib':matplotlib.__version__,'anndata':ad.__version__}}
pack.writestr('aggregate_validation_v3.json',json.dumps(meta,indent=2));pack.close();sys.stdout.buffer.write(buf.getvalue())
