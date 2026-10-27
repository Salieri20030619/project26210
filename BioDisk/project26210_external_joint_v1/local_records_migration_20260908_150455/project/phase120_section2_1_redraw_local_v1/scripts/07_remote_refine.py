from pathlib import Path
import io,sys,json,tempfile,subprocess,zipfile
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter
root=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1');b=io.BytesIO();z=zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.linewidth':.65})
def save(fig,name):
 for ext in ['png','pdf']:
  out=io.BytesIO();fig.savefig(out,format=ext,dpi=300,bbox_inches='tight',pad_inches=.12);z.writestr(name+'.'+ext,out.getvalue())
 plt.close(fig)
# Render the continuous inferCNV checkpoint directly on the server.
r_code=r'''
o<-readRDS(commandArgs(TRUE)[1]);out<-commandArgs(TRUE)[2];m<-log2(o@expr.data);ref<-unlist(o@reference_grouped_cell_indices);obs<-unlist(o@observation_grouped_cell_indices);ord<-c(ref,obs[order(colMeans(abs(m[,obs,drop=FALSE])))]);m<-m[,ord,drop=FALSE]
f<-file(file.path(out,'matrix.bin'),'wb');writeBin(as.numeric(m),f,size=4,endian='little');close(f)
write.table(o@gene_order$chr,file.path(out,'chromosomes.txt'),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(c(nrow(m),ncol(m),length(ref)),file.path(out,'shape.txt'),row.names=FALSE,col.names=FALSE)
'''
with tempfile.TemporaryDirectory(prefix='26210_21_render_') as temp:
 script=Path(temp)/'render.R';script.write_text(r_code)
 rr=subprocess.run(['/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript','--vanilla',str(script),str(root/'03_malignancy_cnv/infercnv_v2/patients/GSE156632__RCC2/infercnv_run/preliminary.infercnv_obj'),temp],capture_output=True,timeout=180)
 if rr.returncode:raise RuntimeError(rr.stderr.decode(errors='replace'))
 ng,nc,nref=np.loadtxt(Path(temp)/'shape.txt',dtype=int);matrix=np.fromfile(Path(temp)/'matrix.bin',dtype='<f4').reshape((ng,nc),order='F').T;chrom=np.loadtxt(Path(temp)/'chromosomes.txt',dtype=str)
 starts=np.r_[0,np.flatnonzero(chrom[1:]!=chrom[:-1])+1];ends=np.r_[starts[1:],ng];centers=(starts+ends-1)/2
 fig=plt.figure(figsize=(11,7));ax=fig.add_axes([.055,.14,.76,.76]);im=ax.imshow(matrix,aspect='auto',interpolation='nearest',vmin=-.5,vmax=.5,cmap=LinearSegmentedColormap.from_list('cnv',['#287BA5','#FFFFFF','#C83E50']),rasterized=True)
 for val in ends[:-1]:ax.axvline(val-.5,color='#444444',lw=.55)
 ax.axhline(nref-.5,color='#444444',lw=.9);ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('Chromosomes',labelpad=8)
 for ii,center in enumerate(centers):ax.text(center,1.015+(0.035 if ii in [13,20] else 0),chrom[starts[ii]].replace('chr',''),transform=ax.get_xaxis_transform(),ha='center',fontsize=8)
 ax.text(1.015,1-nref/(2*nc),'Reference\n200 cells',transform=ax.transAxes,va='center',fontsize=10);ax.text(1.015,(nc-nref)/(2*nc),'Tumour epithelial\nquery\n701 cells',transform=ax.transAxes,va='center',fontsize=10)
 fig.text(.055,.97,'GSE156632  RCC2',fontsize=13,fontweight='bold');ca=fig.add_axes([.84,.19,.135,.025]);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-.5,0,.5]);ca.set_title('log2 relative signal',fontsize=8,pad=6);save(fig,'01_CNV_chromosome_heatmap_v2')
a=ad.read_h5ad(root/'02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad',backed='r')
epi=a.obs.cell_type_broad_v1.astype(str).eq('Epithelial cell').to_numpy();idx=np.flatnonzero(epi);xy=np.asarray(a.obsm['X_umap'])[idx]
calls=pd.read_csv(root/'03_malignancy_cnv/malignancy_cnv_summary_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz',sep='\t').set_index('exact_cell_id')
ids=a.obs_names[idx].astype(str);status=np.array(['Normal-tissue epithelial' if i not in calls.index else calls.loc[i,'malignancy_status_v1'] for i in ids]);assert len(ids)==62652
names={'Normal-tissue epithelial':'Normal-tissue epithelial','malignant_consensus':'Malignant consensus','malignant_infercnv_only':'Malignant inferCNV only','non_malignant_consensus':'Non-malignant consensus','method_discordant':'Method-discordant','indeterminate':'Indeterminate'}
cols={'Normal-tissue epithelial':'#97B8CD','malignant_consensus':'#B779A8','malignant_infercnv_only':'#DBB56C','non_malignant_consensus':'#6BAF9E','method_discordant':'#BCBCBC','indeterminate':'#E2E2E2'}
def axesoff(ax):
 ax.set_aspect('equal');ax.set_xlim(xy[:,0].min()-1,xy[:,0].max()+1);ax.set_ylim(xy[:,1].min()-1,xy[:,1].max()+1);ax.axis('off')
fig,ax=plt.subplots(figsize=(6.7,7));
for s in names:
 ii=status==s;ax.scatter(xy[ii,0],xy[ii,1],s=.75,color=cols[s],linewidths=0,rasterized=True,label=f'{names[s]} ({ii.sum():,})')
axesoff(ax);ax.set_title('Epithelial CNV classification',loc='left',fontweight='bold');ax.legend(loc='upper center',bbox_to_anchor=(.5,-.01),frameon=False,fontsize=8,ncol=2,markerscale=4);plt.close(fig)
genes=['NDUFA4L2','CA9','SLC17A3','MET','TYMP'];expr=a[idx,genes].X.toarray();assert np.isfinite(expr).all()
fig,axs=plt.subplots(1,5,figsize=(14,4));cmap=LinearSegmentedColormap.from_list('expr',['#E6E6E6','#6C367D'])
for j,(g,ax) in enumerate(zip(genes,axs)):
 v=expr[:,j];o=np.argsort(v,kind='stable');hi=np.quantile(v[v>0],.99);sc=ax.scatter(xy[o,0],xy[o,1],c=v[o],cmap=cmap,s=.4,vmin=0,vmax=hi,linewidths=0,rasterized=True);axesoff(ax);ax.set_title(g,fontstyle='italic',fontweight='bold');fig.colorbar(sc,ax=ax,orientation='horizontal',shrink=.55,pad=.025,ticks=[0,hi])
fig.subplots_adjust(wspace=.12);plt.close(fig)

# Reference H/J: black-background density displays on the frozen epithelial UMAP.
xe=np.linspace(xy[:,0].min()-1,xy[:,0].max()+1,190);ye=np.linspace(xy[:,1].min()-1,xy[:,1].max()+1,230)
def density(weights):
 h=np.histogram2d(xy[:,0],xy[:,1],bins=[xe,ye],weights=weights)[0].T;h=gaussian_filter(h,2.5);return h/max(h.sum(),1e-12)/((xe[1]-xe[0])*(ye[1]-ye[0]))
def dark(ax,h):
 im=ax.imshow(h,extent=[xe[0],xe[-1],ye[0],ye[-1]],origin='lower',cmap='magma',aspect='equal',interpolation='bilinear',vmin=0);ax.set_facecolor('#000004');ax.set_xticks([]);ax.set_yticks([]);return im
fig,axs=plt.subplots(1,2,figsize=(7.3,6.3));masks=[np.isin(status,['malignant_consensus','malignant_infercnv_only']),status=='non_malignant_consensus']
for ax,mask,title in zip(axs,masks,['CNV-supported malignant','Non-malignant consensus']):
 im=dark(ax,density(mask.astype(float)));ax.scatter(xy[mask,0],xy[mask,1],s=.10,color='white',alpha=.65,linewidths=0,rasterized=True);ax.set_title(title+'\n'+format(mask.sum(),',')+' cells',fontsize=10);cb=fig.colorbar(im,ax=ax,orientation='horizontal',shrink=.8,pad=.025);cb.set_label('Within-group cell density',fontsize=8)
fig.subplots_adjust(wspace=.12);save(fig,'07_CNV_group_density_reference_H_v1')
base=density(np.ones(len(xy)));fig,axs=plt.subplots(1,5,figsize=(14,5.3))
for j,(ax,g) in enumerate(zip(axs,genes)):
 dark(ax,base);v=expr[:,j];order=np.argsort(v,kind='stable');vmax=np.quantile(v[v>0],.99);sc=ax.scatter(xy[order,0],xy[order,1],c=v[order],cmap='gray',s=.16,vmin=0,vmax=vmax,linewidths=0,rasterized=True);ax.set_title(g,fontstyle='italic',fontweight='bold');cb=fig.colorbar(sc,ax=ax,orientation='horizontal',shrink=.75,pad=.02,ticks=[0,vmax]);cb.set_label('Log1p expression',fontsize=7)
fig.subplots_adjust(wspace=.16);fig.text(.08,.02,'Common epithelial cell-density background; expression is overlaid at the original cell positions.',fontsize=9);save(fig,'08_marker_density_reference_J_v1')
# Source statistics only, no new inferential tests.
statdir=root/'04_subclones_tymp/cnv_clone_tymp_cnmf_association_v1'
stats=pd.read_csv(statdir/'cohort_factor_contrast_summary_v1.tsv',sep='\t');print(stats.columns.tolist(),file=sys.stderr)
z.writestr('formal_aggregate_statistics_v2.csv',stats.to_csv(index=False))
meta={'cells_epithelial':len(ids),'status_counts':{names[s]:int((status==s).sum()) for s in names},'cnv_heatmap':'GSE156632 RCC2 preliminary inferCNV checkpoint; continuous log2 relative signal; 7029 genes; 200 reference + 701 tumour query; fixed color clip +/-0.5','markers':genes,'scope':'CNV classification and expression displays; no reassignment or new inferential test','versions':{'anndata':ad.__version__,'matplotlib':matplotlib.__version__}}
z.writestr('remote_aggregate_validation_v2.json',json.dumps(meta,indent=2));z.close();sys.stdout.buffer.write(b.getvalue())
