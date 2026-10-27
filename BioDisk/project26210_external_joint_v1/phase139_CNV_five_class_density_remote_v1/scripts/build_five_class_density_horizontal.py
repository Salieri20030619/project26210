from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
w=Path(__file__).resolve().parents[1]
for d in ['figures','qc','source_data']: (w/d).mkdir(exist_ok=True)
root=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
a=ad.read_h5ad(root/'02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad',backed='r')
idx=np.flatnonzero(a.obs.cell_type_broad_v1.astype(str).eq('Epithelial cell').to_numpy())
xy=np.asarray(a.obsm['X_umap'])[idx]
ids=a.obs_names[idx].astype(str)
calls=pd.read_csv(root/'03_malignancy_cnv/malignancy_cnv_summary_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz',sep='\t').set_index('exact_cell_id')
status=calls.malignancy_status_v1.reindex(ids).fillna('outside_tumour_calls').to_numpy()
groups=[('malignant_consensus','Malignant consensus',10844),('malignant_infercnv_only','Malignant inferCNV only',1106),('non_malignant_consensus','Non-malignant consensus',498),('method_discordant','Method-discordant',22216),('indeterminate','Indeterminate',4006)]
assert len(xy)==62652 and sum(n for _,_,n in groups)==len(calls)==38670
xe=np.linspace(xy[:,0].min()-1,xy[:,0].max()+1,190);ye=np.linspace(xy[:,1].min()-1,xy[:,1].max()+1,230)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'axes.linewidth':.5})
fig=plt.figure(figsize=(280/25.4,125/25.4))
fig.text(.025,.96,'CNV classification of tumour-tissue epithelial cells',fontsize=11,weight='bold')

rows=[]
for j,(key,title,n) in enumerate(groups):
 mask=status==key
 assert mask.sum()==n
 # Original 189 x 229 histogram grid and Gaussian smoothing, sigma=2.5 bins.
 h=np.histogram2d(xy[:,0],xy[:,1],bins=[xe,ye],weights=mask.astype(float))[0].T
 h=gaussian_filter(h,2.5);h=h/h.sum()/((xe[1]-xe[0])*(ye[1]-ye[0]))
 x=.025+j*.195;y=.23
 ax=fig.add_axes([x,y,.17,.595])
 im=ax.imshow(h,extent=[xe[0],xe[-1],ye[0],ye[-1]],origin='lower',cmap='magma',aspect='equal',interpolation='bilinear',vmin=0)
 ax.scatter(xy[mask,0],xy[mask,1],s=.09,color='white',alpha=.65,linewidths=0,rasterized=True)
 ax.set_facecolor('#000004');ax.set_xticks([]);ax.set_yticks([])
 ax.set_xlabel('UMAP1',labelpad=2,fontsize=7);ax.set_ylabel('UMAP2',labelpad=2,fontsize=7)
 title=title.replace('Malignant consensus','Malignant\nconsensus').replace('Malignant inferCNV only','Malignant\ninferCNV only').replace('Non-malignant consensus','Non-malignant\nconsensus')
 ax.set_title(chr(65+j)+'  '+title+'\n'+format(n,',')+' cells ('+format(n/38670*100,'.2f')+'%)',fontsize=8,pad=6)
 ca=fig.add_axes([x+.012,.125,.146,.014])
 cb=fig.colorbar(im,cax=ca,orientation='horizontal')
 cb.ax.tick_params(labelsize=6,length=2,pad=1)
 cb.locator=matplotlib.ticker.MaxNLocator(nbins=3);cb.update_ticks()
 cb.set_label('Within-group cell density',fontsize=6.5,labelpad=2)
 rows.append(dict(status=key,cells=n,percent=n/38670*100,density_max=float(h.max()),density_integral=float(h.sum()*(xe[1]-xe[0])*(ye[1]-ye[0]))))
name='10_CNV_five_class_density_reference_H_horizontal_v3'
fig.savefig(w/'figures'/(name+'.pdf'),dpi=600)
fig.savefig(w/'figures'/(name+'.png'),dpi=600)
pd.DataFrame(rows).to_csv(w/'source_data/class_counts_and_density.tsv',sep='\t',index=False)
meta=dict(cells_total=38670,epithelial_coordinate_scope=62652,normal_tissue_excluded=23982,width_mm=280,height_mm=125,classes=rows,counts_complete=True,original_coordinates=True,histogram_bins=[189,229],gaussian_sigma_bins=2.5,normalization='within_group_integral_1',new_inferential_tests=False)
(w/'qc/VALIDATION.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta))

