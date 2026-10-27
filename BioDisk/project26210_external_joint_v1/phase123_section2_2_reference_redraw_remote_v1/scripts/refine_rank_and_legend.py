from pathlib import Path
import json,hashlib,sys,textwrap
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle,Patch
from scipy.cluster.hierarchy import linkage,dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
w=Path(__file__).resolve().parents[1];out=w/'figures';qc=w/'qc';base=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1');r=base/'05_cnmf'
archive=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/local_records_migration_20260908_150455/project')
source_files=[]
def read(path):
 path=Path(path);source_files.append(path);return pd.read_csv(path,sep='\t')
u=read(r/'full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz').set_index('exact_cell_id');load=read(r/'full_ranks4_8_v1/consensus_loading_rank7_v1.tsv.gz').set_index('gene');m=read(r/'input_v1/cnmf_cell_metadata_v1.tsv.gz').set_index('exact_cell_id').loc[u.index]
factors=[f'F{i}' for i in range(1,8)];assert len(u)==11950 and len(load)==2500 and u.index.is_unique
assert np.allclose(u.sum(axis=1),1) and np.isfinite(u.values).all() and (u.values>=0).all()
m['factor']=u.idxmax(axis=1);m['max_usage']=u.max(axis=1)
rank=read(r/'full_ranks4_8_v1/rank_qc_v1.tsv');evidence=read(r/'tymp_rank_evidence_v1/factor_evidence_v1.tsv');pcorr=read(r/'tymp_rank_evidence_v1/factor_patient_stratified_correlations_v1.tsv.gz')
select=read(r/'tymp_study_sensitivity_v1/selected_rank_v1.tsv').iloc[0];assert select.selected_rank==7 and select.selected_TYMP_factor=='F1'
phase57=archive/'phase57_part1_by_cohort_signal_audit_local_v1/analysis_local_v1';patients=read(phase57/'patient_effects_long_v1.tsv');synth=read(phase57/'cross_cohort_synthesis_v1.tsv')
tymp_id=next(x for x in synth.signal_id if x=='F1__TYMP_RNA' or ('TYMP' in x and 'strict' not in x.lower() and 'broad' not in x.lower() and 'axis' not in x.lower()))
pt=patients[patients.signal_id==tymp_id].copy();st=synth[synth.signal_id==tymp_id].iloc[0];assert len(pt)==20 and (pt.effect>0).sum()==17
colors=['#BA7868','#7295BE','#79A697','#C9AD63','#A58BB4','#8DB6C3','#C78FAD'];fd=dict(zip(factors,colors))
studies=['GSE156632','GSE159115','GSE178481','GSE207493','GSE242299'];sd=dict(zip(studies,['#7296BC','#E3AD66','#76A893','#9F86B4','#C47F85']))
blue=LinearSegmentedColormap.from_list('reference_blue',['#FFFFFF','#D8DBE9','#8A98BA','#315B91'])
div=LinearSegmentedColormap.from_list('rna_div',['#3779A3','#F8F7F3','#BA514C'])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':10,'axes.labelsize':9,'axes.linewidth':.65,'pdf.fonttype':42,'savefig.facecolor':'white','xtick.major.width':.65,'ytick.major.width':.65})
figlist=[]
def save(fig,name):
 for ext in ['png','pdf']:
  f=out/(name+'.'+ext);fig.savefig(f,dpi=320,bbox_inches='tight',pad_inches=.10);figlist.append(f)
 plt.close(fig)
def clean(ax):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3)
def letter(fig,s):fig.text(.015,.985,s,fontsize=16,fontweight='bold',va='top')

runs=read(r/'full_ranks4_8_v1/run_qc_v1.tsv')
fig,axs=plt.subplots(1,2,figsize=(8.8,3.4));a=axs[0]
a.plot(rank['rank'],rank.mean_component_silhouette,'o-',c='#7295BE',lw=1.3,label='Mean silhouette');a.plot(rank['rank'],rank.min_run_fraction,'s-',c='#BA7868',lw=1.3,label='Minimum run coverage');a.axvspan(6.84,7.16,color='#EEE3C8',zorder=0);a.set(xlabel='Number of programs (k)',ylabel='Stability',xticks=rank['rank'],ylim=(0,1.05));a.legend(frameon=False,fontsize=7,loc='lower left');a.set_title('Twenty starts per rank',loc='left');clean(a)
a=axs[1];rng=np.random.default_rng(26210)
for k in range(4,9):
 vals=runs.loc[runs['rank']==k,'mse'].to_numpy();assert len(vals)==20;col='#D8B76E' if k==7 else '#8FAABA'
 a.boxplot(vals,positions=[k],widths=.43,showfliers=False,patch_artist=True,boxprops={'facecolor':col,'edgecolor':'#555555','linewidth':.7},medianprops={'color':'#333333','linewidth':1},whiskerprops={'color':'#777777','linewidth':.7},capprops={'color':'#777777','linewidth':.7})
 a.scatter(k+rng.uniform(-.16,.16,len(vals)),vals,c='#333333',s=7,alpha=.4,zorder=3,linewidths=0)
a.set(xlabel='Number of programs (k)',ylabel='Reconstruction MSE',xticks=list(range(4,9)));a.set_title('Frozen selection: k = 7, candidate F1',loc='left');clean(a);fig.tight_layout(pad=2);save(fig,'01_cNMF_rank_selection_v1')
# Program loading enrichment: exact original FDR, independent of clone DE.
ora=read(r/'cnv_clone_program_function_v1/rank7_program_pathway_ORA_v1.tsv');chosen=[]
for f in factors:
 for term in ora[(ora.analysis==f)&(ora.FDR<.05)].sort_values(['FDR','pathway']).pathway.head(2):
  if term not in chosen:chosen.append(term)
short=lambda s:s.split('::')[-1].replace('HALLMARK_','').replace('REACTOME_','').replace('_',' ').capitalize()
fig,ax=plt.subplots(figsize=(8.1,max(4.4,len(chosen)*.30+1.4)))
for j,f in enumerate(factors):
 for i,term in enumerate(chosen):
  rr=ora[(ora.analysis==f)&(ora.pathway==term)]
  if rr.empty:continue
  row=rr.iloc[0]
  if row.FDR>=.05:continue
  ax.scatter(j,i,s=15+row.overlap_size*5,c=[-np.log10(max(row.FDR,1e-300))],cmap=blue,vmin=0,vmax=max(5,-np.log10(ora.loc[ora.FDR>0,'FDR'].min())),edgecolors='#52647D',linewidths=.4)
ax.set_xticks(range(7),factors);ax.set_yticks(range(len(chosen)),[textwrap.fill(short(t),37) for t in chosen],fontsize=7.5);ax.invert_yaxis();ax.set_xlim(-.6,6.6);ax.set_xlabel('cNMF program');ax.set_title('Functions enriched among top-loading genes',loc='left',pad=12);clean(ax);ax.grid(axis='y',c='#ECECEC',lw=.5);ax.set_axisbelow(True)
sm=plt.cm.ScalarMappable(norm=plt.Normalize(0,max(5,-np.log10(ora.loc[ora.FDR>0,'FDR'].min()))),cmap=blue);fig.colorbar(sm,ax=ax,fraction=.035,pad=.05,label='−log10(FDR)');fig.subplots_adjust(left=.39,right=.91,bottom=.14,top=.89);fig.text(.39,.035,'Dot area: number of overlapping genes; only FDR < 0.05 shown.',fontsize=7);save(fig,'08_program_function_enrichment_v1')
