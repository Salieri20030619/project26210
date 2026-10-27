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
o<-readRDS(commandArgs(TRUE)[1]);out<-commandArgs(TRUE)[2]
m<-log2(o@expr.data);ref<-unlist(o@reference_grouped_cell_indices);obs<-unlist(o@observation_grouped_cell_indices)
ord<-c(ref,obs[order(colMeans(abs(m[,obs,drop=FALSE])))]);m<-m[,rev(ord),drop=FALSE]
chr<-as.character(o@gene_order$chr);rr<-rle(chr);en<-cumsum(rr$lengths);st<-c(1,head(en,-1)+1);cent<-(st+en)/2
pal<-colorRampPalette(c('#287BA5','#FFFFFF','#C83E50'))(201)
draw<-function(){
 par(mar=c(4,2,3,7),xpd=NA,las=1,family='sans')
 image(seq_len(nrow(m)),seq_len(ncol(m)),m,col=pal,zlim=c(-.5,.5),useRaster=TRUE,axes=FALSE,xlab='',ylab='')
 abline(v=head(en,-1)+.5,col='#333333',lwd=.55);abline(h=length(obs)+.5,col='#333333',lwd=1)
 axis(3,at=cent,labels=sub('chr','',rr$values),tick=FALSE,cex.axis=.70,line=-.4)
 mtext('Chromosomes',side=1,line=1,cex=.9);mtext('GSE156632  RCC2',side=3,line=1.2,adj=0,font=2,cex=1)
 text(nrow(m)*1.01,length(obs)+length(ref)/2,paste0('Reference\n',length(ref),' cells'),adj=0,cex=.78)
 text(nrow(m)*1.01,length(obs)/2,paste0('Tumour epithelial\nquery\n',length(obs),' cells'),adj=0,cex=.78)
 box(lwd=.7)
 par(fig=c(.81,.95,.13,.18),new=TRUE,mar=c(1,0,0,0));image(seq(-.5,.5,length.out=201),1,matrix(seq(-.5,.5,length.out=201),ncol=1),col=pal,axes=FALSE,xlab='',ylab='');axis(1,at=c(-.5,0,.5),labels=c('-0.5','0','0.5'),cex.axis=.6,tck=-.1);mtext('log2 relative signal',side=3,line=.1,cex=.65)
}
png(file.path(out,'01_CNV_chromosome_heatmap_v1.png'),width=3000,height=2300,res=300);draw();dev.off()
pdf(file.path(out,'01_CNV_chromosome_heatmap_v1.pdf'),width=10,height=7.67,useDingbats=FALSE);draw();dev.off()
'''
with tempfile.TemporaryDirectory(prefix='26210_21_render_') as temp:
 script=Path(temp)/'render.R';script.write_text(r_code)
 rr=subprocess.run(['/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript','--vanilla',str(script),str(root/'03_malignancy_cnv/infercnv_v2/patients/GSE156632__RCC2/infercnv_run/preliminary.infercnv_obj'),temp],capture_output=True,timeout=180)
 if rr.returncode:raise RuntimeError(rr.stderr.decode(errors='replace'))
 for f in Path(temp).glob('01_CNV*'):z.writestr(f.name,f.read_bytes())
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
axesoff(ax);ax.set_title('Epithelial CNV classification',loc='left',fontweight='bold');ax.legend(loc='upper center',bbox_to_anchor=(.5,-.01),frameon=False,fontsize=8,ncol=2,markerscale=4);save(fig,'02_epithelial_CNV_status_UMAP_v1')
genes=['NDUFA4L2','CA9','SLC17A3','MET','TYMP'];expr=a[idx,genes].X.toarray();assert np.isfinite(expr).all()
fig,axs=plt.subplots(1,5,figsize=(14,4));cmap=LinearSegmentedColormap.from_list('expr',['#E6E6E6','#6C367D'])
for j,(g,ax) in enumerate(zip(genes,axs)):
 v=expr[:,j];o=np.argsort(v,kind='stable');hi=np.quantile(v[v>0],.99);sc=ax.scatter(xy[o,0],xy[o,1],c=v[o],cmap=cmap,s=.4,vmin=0,vmax=hi,linewidths=0,rasterized=True);axesoff(ax);ax.set_title(g,fontstyle='italic',fontweight='bold');fig.colorbar(sc,ax=ax,orientation='horizontal',shrink=.55,pad=.025,ticks=[0,hi])
fig.subplots_adjust(wspace=.12);save(fig,'03_epithelial_marker_UMAP_v1')
# Source statistics only, no new inferential tests.
statdir=root/'04_subclones_tymp/cnv_clone_tymp_cnmf_association_v1'
stats=pd.read_csv(statdir/'cohort_factor_contrast_summary_v1.tsv',sep='\t');print(stats.columns.tolist(),file=sys.stderr)
z.writestr('formal_aggregate_statistics.csv',stats.to_csv(index=False))
meta={'cells_epithelial':len(ids),'status_counts':{names[s]:int((status==s).sum()) for s in names},'cnv_heatmap':'GSE156632 RCC2 preliminary inferCNV checkpoint; continuous log2 relative signal; 7029 genes; 200 reference + 701 tumour query; fixed color clip +/-0.5','markers':genes,'scope':'CNV classification and expression displays; no reassignment or new inferential test','versions':{'anndata':ad.__version__,'matplotlib':matplotlib.__version__}}
z.writestr('remote_aggregate_validation.json',json.dumps(meta,indent=2));z.close();sys.stdout.buffer.write(b.getvalue())
