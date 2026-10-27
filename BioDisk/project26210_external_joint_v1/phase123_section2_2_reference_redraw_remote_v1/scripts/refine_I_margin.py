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
# Select program genes by loading only; do not label these as differential tests.
genes=[];gene_factor=[]
for f in factors:
 choices=[g for g in load[f].sort_values(ascending=False,kind='stable').index if g not in genes][:5]
 genes.extend(choices);gene_factor.extend([f]*len(choices))
i_genes=genes.copy();i_gf=gene_factor.copy()
if 'TYMP' not in genes:genes.append('TYMP');gene_factor.append('TYMP')
# I: reference blue, exact gene contribution matrix and source factor groups.
vals=load.loc[i_genes,factors];norm=vals.div(vals.max(axis=1),axis=0);assert np.isfinite(norm).all().all()
fig=plt.figure(figsize=(5.3,9));ax=fig.add_axes([.27,.075,.39,.84]);im=ax.imshow(norm,cmap=blue,vmin=0,vmax=1,aspect='auto',interpolation='nearest');ax.set_xticks(range(7),factors,fontsize=9);ax.set_yticks(range(len(i_genes)),i_genes,fontsize=8.5,fontstyle='italic');ax.yaxis.tick_right();ax.tick_params(length=0,pad=4);ax.set_xlabel('cNMF program',labelpad=8)
ax.set_xticks(np.arange(-.5,7,1),minor=True);ax.set_yticks(np.arange(-.5,len(i_genes),1),minor=True);ax.grid(which='minor',color='#555555',lw=.45);ax.tick_params(which='minor',length=0)
for j,f in enumerate(factors):
 lo=j*5-.5;hi=(j+1)*5-.5;ax.plot([-.78,-.78],[lo,hi],clip_on=False,c=fd[f],lw=2.2);ax.text(-1,lo+2.5,f+' genes',ha='right',va='center',fontsize=8,color='#444444')
ax.set_xlim(-.5,6.5)
ca=fig.add_axes([.27,.95,.39,.015]);fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[0,.5,1]);ca.set_title('Relative gene contribution',fontsize=8);letter(fig,'I');save(fig,'05_I_gene_contribution_reference_blue_v2')
norm.to_csv(w/'qc/I_gene_contribution_aggregate.tsv',sep='\t')
