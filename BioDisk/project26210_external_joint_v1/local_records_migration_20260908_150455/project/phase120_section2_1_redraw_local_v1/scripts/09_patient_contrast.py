from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).resolve().parents[2];d=Path('C:/Users/31937/Downloads/26210/02_2.1恶性上皮亚克隆_v1');s=r/'phase32_part1_m03_cnv_clone_trees_local_v1/input_snapshot_local_v1'
c=pd.read_csv(s/'formal_cell_cnv_clone_assignments_v1.tsv.gz',sep='\t');u=pd.read_csv(s/'consensus_usage_rank7_v1.tsv.gz',sep='\t');c=c.merge(u,on='exact_cell_id',validate='one_to_one');c=c[c.clone_cells>=20]
a=c.groupby(['infercnv_unit','formal_inferred_cnv_clone_v1']).agg(F1=('F1','median'),TYMP=('TYMP_log1p_CP10K_v1','median')).reset_index();pts=[]
for pt,g in a.groupby('infercnv_unit'):
 if len(g)<2:continue
 g=g.sort_values('formal_inferred_cnv_clone_v1');hi=g.F1.idxmax();pts.append((pt,g.loc[hi,'TYMP'],g.drop(hi).TYMP.mean()))
t=pd.DataFrame(pts,columns=['patient','F1_high_clone','other_clones']);t['delta']=t.F1_high_clone-t.other_clones
stats=pd.read_csv(d/'核查记录/formal_aggregate_statistics.csv');f1=stats[(stats.factor=='F1')&np.isclose(stats.distance_threshold,.25)].iloc[0]
assert len(t)==14 and np.isclose(t.delta.median(),f1.median_patient_tymp_delta,atol=1e-9)
plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'axes.linewidth':.7});fig,axs=plt.subplots(1,2,figsize=(10.5,5.3),gridspec_kw={'width_ratios':[1,1.5]})
for row in t.itertuples():axs[0].plot([0,1],[row.other_clones,row.F1_high_clone],'-o',c='#B57AA3' if row.delta>0 else '#8CA9B9',alpha=.8,lw=.9,ms=4)
axs[0].set_xticks([0,1],['Other primary\nclones','F1-high\nprimary clone']);axs[0].set_ylabel('TYMP median log1p CP10K');axs[0].set_title('Patient-matched clone contrast',loc='left',fontweight='bold');axs[0].text(.03,.97,f'n = 14 patients\nMedian delta = {f1.median_patient_tymp_delta:.3f}\nP = {f1.wilcoxon_two_sided_p:.4f}\nFDR = {f1.wilcoxon_BH_FDR_within_threshold:.4f}',transform=axs[0].transAxes,va='top',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.8})
t=t.sort_values('delta');axs[1].hlines(range(len(t)),0,t.delta,color=['#B57AA3' if v>0 else '#8CA9B9' for v in t.delta],lw=1.4);axs[1].scatter(t.delta,range(len(t)),c=['#B57AA3' if v>0 else '#8CA9B9' for v in t.delta],s=27);axs[1].set_yticks(range(len(t)),[p.replace('GSE','').replace('__',' / ') for p in t.patient],fontsize=8);axs[1].axvline(0,c='#999999',ls='--',lw=.7);axs[1].set_xlabel('Within-patient TYMP difference');axs[1].set_title('9/14 patients positive',loc='left',fontweight='bold')
for ax in axs:ax.spines[['top','right']].set_visible(False)
fig.tight_layout(w_pad=3)
for ext in ['png','pdf']:fig.savefig(d/('09_patient_level_TYMP_clone_contrast_v1.'+ext),dpi=300,bbox_inches='tight',pad_inches=.12)
plt.close(fig)
print('Patient plot verified against frozen formal median',t.delta.median())
