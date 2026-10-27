import io,sys,json,zipfile
import numpy as np,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,to_rgba
from matplotlib.lines import Line2D
p='/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad'
a=ad.read_h5ad(p,backed='r')
types=['Dendritic cell','Macrophage','Monocyte','Mast cell','B cell','NK cell','T cell','Plasma cell','Epithelial cell','Endothelial','Pericyte/SMC','Fibroblast','Erythrocyte','Epithelial/Myeloid mixed']
genes=['CLEC9A','C1QA','S100A8','TPSAB1','MS4A1','NKG7','CD3D','MZB1','EPCAM','PECAM1','RGS5','COL1A1','HBB']
colors=['#7AAFA0','#D76A73','#E1A078','#AF7694','#8D719D','#A5A76B','#81729D','#BC8CA6','#D4AA62','#6AAE99','#7DA2BF','#A3BC78','#C57B65','#A7A7A7']
labels=a.obs['cell_type_broad_v1'].astype(str).to_numpy();xy=np.asarray(a.obsm['X_umap']);x=a[:,genes].X.toarray()
assert len(labels)==258120 and set(labels)==set(types) and np.isfinite(xy).all() and np.isfinite(x).all() and x.min()>=0
counts=np.array([(labels==t).sum() for t in types]);assert counts.sum()==258120
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'pdf.fonttype':42,'axes.linewidth':.6})
buffer=io.BytesIO();pack=zipfile.ZipFile(buffer,'w',zipfile.ZIP_DEFLATED)
lims=[xy[:,0].min()-1,xy[:,0].max()+1,xy[:,1].min()-1,xy[:,1].max()+1]
def save(fig,name):
 for ext in ['png','pdf']:
  b=io.BytesIO();fig.savefig(b,format=ext,dpi=300,bbox_inches='tight',pad_inches=.12);pack.writestr(name+'.'+ext,b.getvalue())
 plt.close(fig)
def frame(ax,arrow=False):
 ax.set_xlim(lims[:2]);ax.set_ylim(lims[2:]);ax.set_aspect('equal');ax.axis('off')
 if arrow:
  for end in [(.22,.04),(.04,.19)]:ax.annotate('',xy=end,xytext=(.04,.04),xycoords='axes fraction',arrowprops={'arrowstyle':'-|>','lw':.75,'color':'#333333'})
  ax.text(.13,0,'UMAP 1',transform=ax.transAxes,fontsize=6,ha='center');ax.text(0,.12,'UMAP 2',transform=ax.transAxes,fontsize=6,rotation=90,va='center')
def landscape(ax):
 order=np.random.default_rng(26210).permutation(len(labels));col=np.array([to_rgba(colors[types.index(t)]) for t in labels]);ax.scatter(xy[order,0],xy[order,1],c=col[order],s=.28,linewidths=0,rasterized=True,alpha=.85);frame(ax,True)
scales=[]
for j in range(13):scales.append(float(np.quantile(x[:,j][x[:,j]>0],.99)))
def feature(ax,j,arrow=False):
 v=x[:,j];ix=np.flatnonzero(v>0);ix=ix[np.argsort(v[ix],kind='stable')];cmap=LinearSegmentedColormap.from_list('expr',['#E3E3E3',colors[j]])
 ax.scatter(xy[:,0],xy[:,1],s=.18,c='#DDDDDD',linewidths=0,rasterized=True)
 ax.scatter(xy[ix,0],xy[ix,1],c=v[ix],vmin=0,vmax=scales[j],cmap=cmap,s=.35,linewidths=0,rasterized=True)
 frame(ax,arrow);ax.set_title(genes[j],fontstyle='italic',fontweight='bold',fontsize=12,pad=18);ax.text(.5,1.01,types[j],transform=ax.transAxes,ha='center',fontsize=8,color='#555555')
 ca=ax.inset_axes([.64,.005,.29,.015]);ca.imshow(np.linspace(0,1,128)[None,:],cmap=cmap,aspect='auto');ca.set_xticks([0,127],['0',f'{scales[j]:.1f}'],fontsize=6);ca.set_yticks([]);ca.tick_params(length=0,pad=2);ca.spines[:].set_visible(False)
fig,ax=plt.subplots(figsize=(7,7.8));landscape(ax);ax.set_title('Cell identities',loc='left',fontweight='bold');ax.text(1,1.02,'258,120 cells',transform=ax.transAxes,ha='right',fontsize=10)
handles=[Line2D([],[],marker='o',ls='',markersize=5,color=c,label=t) for t,c in zip(types,colors)]
ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.04),ncol=2,frameon=False,fontsize=9,handletextpad=.4,columnspacing=1.5);save(fig,'03_cell_identity_UMAP_v1')
for start,end,name in [(0,7,'04_marker_UMAP_immune_v1'),(7,13,'05_marker_UMAP_stromal_epithelial_v1')]:
 fig,axes=plt.subplots(1,end-start,figsize=(2.5*(end-start),4.2))
 for j,ax in zip(range(start,end),np.ravel(axes)):feature(ax,j,j==start)
 fig.subplots_adjust(wspace=.15,bottom=.14,top=.83);fig.text(.02,.025,'Log1p expression. Gene-specific colour limits: 99th percentile of positive expression.',fontsize=9,color='#555555');save(fig,name)
fig,axes=plt.subplots(3,5,figsize=(14,12.3));axes=axes.ravel();landscape(axes[0]);axes[0].set_title('Cell identities',fontweight='bold',pad=18)
for j in range(13):feature(axes[j+1],j,j in [0,4,9])
axes[-1].axis('off');axes[-1].legend(handles=handles,loc='center',frameon=False,fontsize=9,title='Cell type',title_fontsize=11)
fig.subplots_adjust(hspace=.32,wspace=.20,bottom=.06,top=.95);fig.text(.10,.017,'258,120 cells | Identical UMAP coordinates in every panel | Log1p expression with gene-specific colour scales',fontsize=10,color='#555555');save(fig,'00_singlecell_marker_overview_v1')
for j,g in enumerate(genes):
 fig,ax=plt.subplots(figsize=(3.6,4.5));feature(ax,j,True);save(fig,'Marker_'+g+'_v1')
meta={'status':'PASS_SOURCE_CHECKS','cells':len(labels),'counts':dict(zip(types,counts.tolist())),'genes':genes,'gene_colour_upper':scales,'source':p,'source_shape':list(a.shape),'n_studies':int(a.obs.study.nunique()),'n_patient_keys':len(a.obs[['study','curated_patient_id']].drop_duplicates()),'n_sample_keys':len(a.obs[['study','sample']].drop_duplicates()),'versions':{'numpy':np.__version__,'matplotlib':matplotlib.__version__,'anndata':ad.__version__},'notes':'Frozen coordinates and annotations. Every cell included. Markers plotted across all cell identities. No cell identifiers or per-cell matrices transmitted.'}
pack.writestr('aggregate_validation.json',json.dumps(meta,indent=2));pack.close();sys.stdout.buffer.write(buffer.getvalue())
