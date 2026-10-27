from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.stats import rankdata
w=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase143_sections25_27_candidate_report_remote_v1')
for d in ['scripts','figures','qc','report','fonts']: (w/d).mkdir(parents=True,exist_ok=True)
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory/trajectory_v1')
a=ad.read_h5ad(b/'objects/part1_malignant_trajectory_all_methods_v1.h5ad',backed='r')
g=pd.read_csv(b/'source_data/monocle3_principal_graph_segments_v1.tsv',sep='\t')
xy=np.asarray(a.obsm['X_umap']);t=a.obs.Monocle3_pseudotime_v1.to_numpy(float);valid=np.isfinite(t)
r=np.full(len(t),np.nan);r[valid]=(rankdata(t[valid])-1)/max(valid.sum()-1,1)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(190/25.4,150/25.4))
ax.scatter(xy[:,0],xy[:,1],s=1.5,c='#dddddd',linewidths=0,rasterized=True)
o=np.argsort(np.nan_to_num(r,nan=-1));o=o[np.isfinite(r[o])]
sc=ax.scatter(xy[o,0],xy[o,1],c=r[o],s=2.2,cmap='viridis',vmin=0,vmax=1,linewidths=0,rasterized=True)
segments=np.stack([g[['x','y']].to_numpy(),g[['xend','yend']].to_numpy()],axis=1)
ax.add_collection(LineCollection(segments,colors='white',linewidths=2.0,zorder=3))
ax.add_collection(LineCollection(segments,colors='#242424',linewidths=.7,zorder=4))
ax.set_aspect('equal');ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
ax.set_title('Monocle3 trajectory on the original UMAP',loc='left',fontweight='bold')
cb=fig.colorbar(sc,ax=ax,shrink=.6,pad=.025);cb.set_label('Global pseudotime rank');cb.set_ticks([0,.5,1],labels=['Early','0.5','Late'])
ax.spines[['top','right']].set_visible(False)
fig.subplots_adjust(left=.1,right=.89,top=.9,bottom=.16)
fig.text(.1,.075,f'{len(xy):,} malignant cells | {len(g)} original principal-graph segments',fontsize=8)
fig.text(.1,.04,'Lines: frozen Monocle3 principal graph. Color: global ordering; not longitudinal cell motion.',fontsize=7)
for ext in ['pdf','png']:fig.savefig(w/'figures'/('14_Monocle3_original_UMAP_trajectory_v1.'+ext),dpi=450)
a.obs[['Monocle3_pseudotime_v1','Monocle3_pseudotime_patient_v1']].to_csv(w/'qc/trajectory_values.tsv.gz',sep='\t')
(w/'qc/TRAJECTORY_VALIDATION.json').write_text(json.dumps(dict(cells=len(xy),finite_pseudotime=int(valid.sum()),segments=len(g),models_rerun=False,original_graph=True,color='global_pseudotime_rank',heatmap_color_uses_different_within_patient_rank=True),indent=2))
print('TRAJECTORY_DONE',valid.sum(),len(g))
