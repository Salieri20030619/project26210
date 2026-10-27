from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter1d
from scipy.stats import spearmanr
w=Path(__file__).resolve().parents[1];prev=w.parent/'phase123_section2_2_reference_redraw_remote_v1';d=np.load(prev/'data/plot_values_server_only.npz');cells=d['cells'];genes=d['genes'];x=d['expression'];assert x.shape==(11950,36)
u=pd.read_csv('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/05_cnmf/full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz',sep='\t').set_index('exact_cell_id').loc[cells];f1=u.F1.to_numpy()
# Shared F1 ordering and tie-preserving bins: never independently sort a gene.
binid=np.minimum(79,np.floor((pd.Series(f1).rank(method='average').to_numpy()-1)*80/len(f1)).astype(int));labels=np.sort(np.unique(binid))[::-1];means=np.array([x[binid==b].mean(axis=0) for b in labels]);counts=np.array([(binid==b).sum() for b in labels]);scores=np.array([f1[binid==b].mean() for b in labels]);assert counts.sum()==11950 and np.all(np.diff(scores)<=0)
smooth=gaussian_filter1d(means*counts[:,None],sigma=1.5,axis=0,mode='nearest')/gaussian_filter1d(counts.astype(float),sigma=1.5,mode='nearest')[:,None]
z=(smooth-smooth.mean(axis=0))/np.where(smooth.std(axis=0)>0,smooth.std(axis=0),1)
cor=np.array([spearmanr(means[:,j],scores).statistic for j in range(len(genes))]);roworder=np.argsort(-np.nan_to_num(cor),kind='stable')
# Linear interpolation only along the common horizontal trend axis.
grid=np.linspace(0,len(labels)-1,800);display=np.array([np.interp(grid,np.arange(len(labels)),z[:,j]) for j in roworder]);scorebar=np.interp(grid,np.arange(len(labels)),scores)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.linewidth':.7})
cmap=LinearSegmentedColormap.from_list('reference_E',['#1C6791','#F8F4ED','#B7232B'])
fig=plt.figure(figsize=(5.5,8.4));ax=fig.add_axes([.16,.08,.48,.78]);im=ax.imshow(display,aspect='auto',cmap=cmap,vmin=-2,vmax=2,interpolation='nearest');ax.set_yticks(range(len(genes)),genes[roworder],fontsize=8.5);ax.yaxis.tick_right();ax.tick_params(axis='y',length=0,pad=4);ax.set_xticks([0,799],['F1-high','F1-low']);ax.tick_params(axis='x',length=0,pad=6);ax.set_xlabel('Common F1-ranked bins',labelpad=7)
for spine in ax.spines.values():spine.set_linewidth(.8)
top=fig.add_axes([.16,.872,.48,.030]);top.imshow(scorebar[None,:],aspect='auto',cmap=cmap,vmin=0,vmax=1,interpolation='nearest');top.set_xticks([]);top.set_yticks([]);top.set_title('F1 usage',fontsize=8,loc='left',pad=4)
ca=fig.add_axes([.16,.95,.32,.016]);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-2,0,2]);ca.set_title('Smoothed RNA trend (z-score)',fontsize=8,pad=5);cb.ax.tick_params(length=2,labelsize=8);fig.text(.02,.98,'E',fontsize=18,fontweight='bold',va='top')
for ext in ['png','pdf']:fig.savefig(w/'figures'/('E_F1_ordered_smooth_gradient_v2.'+ext),dpi=350,bbox_inches='tight',pad_inches=.1)
plt.close(fig)
np.savez_compressed(w/'data/E_binning_and_smoothing_server_only.npz',bin_labels=labels,counts=counts,F1_mean=scores,raw_bin_mean=means,smoothed_mean=smooth,z_score=z,genes=genes,row_order=roworder)
qc={'status':'PASS','cells':len(cells),'genes':len(genes),'requested_bins':80,'occupied_tie_preserving_bins':len(labels),'smallest_bin':int(counts.min()),'largest_bin':int(counts.max()),'all_genes_share_same_cell_and_bin_order':True,'smoothing':'count-weighted Gaussian sigma 1.5 bins; horizontal linear display interpolation; no vertical interpolation','row_order':'descending descriptive Spearman relation of raw bin means with mean F1; gene membership unchanged','metric':'RNA expression trend, not VIPER inferred activity','prior_inferential_results_unchanged':True}
(w/'qc/E_VALIDATION.json').write_text(json.dumps(qc,indent=2));print(json.dumps(qc))
