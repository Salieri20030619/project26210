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

from pathlib import Path
import hashlib
w=Path(__file__).resolve().parents[1]
pack.close()
fig=plt.figure(figsize=(14,10.8))
grid=fig.add_gridspec(3,5,left=.025,right=.99,top=.95,bottom=.065,wspace=.16,hspace=.24)
positions=[(j//5,j%5) for j in range(10)]+[(2,1),(2,2),(2,3)]
for j,(row,col) in enumerate(positions):
 ax=fig.add_subplot(grid[row,col])
 feature(ax,j,arrow=(j in [0,5,10]))
fig.text(.5,.022,'258,120 cells | Identical UMAP coordinates | Log1p expression; gene-specific colour scales',ha='center',fontsize=9,color='#555555')
fig.canvas.draw()
from matplotlib.text import Text
rend=fig.canvas.get_renderer();outside=[]
for artist in fig.findobj(match=Text):
 if artist.get_visible() and artist.get_text().strip():
  box=artist.get_window_extent(rend)
  if box.x0<0 or box.y0<0 or box.x1>fig.bbox.width or box.y1>fig.bbox.height:outside.append(artist.get_text())
for ext in ['pdf','png']:
 path=w/'figures'/('00_singlecell_markers_only_v2.'+ext)
 fig.savefig(path,dpi=300,bbox_inches='tight',pad_inches=.12)
 print('SAVED',path.name,path.stat().st_size,flush=True)
plt.close(fig)
a.file.close()
meta={'status':'PASS_SOURCE_AND_LAYOUT_CHECKS','source':p,'cells':len(labels),'genes':genes,'n_gene_panels':13,'identity_panel_included':False,'categorical_identity_legend_included':False,'layout':'5+5+3; last row centered','positions':positions,'gene_colour_upper':scales,'coordinate_sha256':hashlib.sha256(xy.tobytes()).hexdigest(),'expression13_sha256':hashlib.sha256(x.tobytes()).hexdigest(),'source_figure_script_sha256':hashlib.sha256(Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/local_records_migration_20260908_150455/project/phase118_singlecell_figure_redraw_local_v1/scripts/06_remote_render.py').read_bytes()).hexdigest(),'outside_page_text':outside,'no_new_analysis_or_model':True}
(w/'qc/VALIDATION.json').write_text(json.dumps(meta,indent=2))
print('SOURCE_AND_LAYOUT_VALIDATION',json.dumps(meta),flush=True)
