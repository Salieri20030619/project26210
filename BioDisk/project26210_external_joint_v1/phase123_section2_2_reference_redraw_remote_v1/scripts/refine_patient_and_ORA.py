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
# Formal F1-TYMP audit and frozen F1-CNV endpoint shown side by side.
cnv=pcorr[(pcorr['rank']==7)&(pcorr.factor=='F1')&(pcorr.outcome=='infercnv_burden_over_threshold_v1')].dropna(subset=['rho']).copy();cnv['study']=cnv.curated_patient_id.str.split('__').str[0];ce=evidence[(evidence['rank']==7)&(evidence.factor=='F1')].iloc[0]
fig,axs=plt.subplots(1,2,figsize=(10.1,6.2))
for ax,d,value,title in [(axs[0],pt,'effect','F1–TYMP RNA'),(axs[1],cnv,'rho','F1–CNV burden')]:
 d=d.sort_values(['study','curated_patient_id'],kind='stable').reset_index(drop=True);ys=np.arange(len(d));ax.axvline(0,c='#AAAAAA',ls='--',lw=.8);ax.hlines(ys,0,d[value],color=[sd[s] for s in d.study],lw=1.2);ax.scatter(d[value],ys,c=[sd[s] for s in d.study],s=28,zorder=3,edgecolors='white',linewidths=.45);ax.set_yticks(ys,[x.split('__',1)[1] for x in d.curated_patient_id],fontsize=7);ax.invert_yaxis();ax.set_xlim(-1,1);ax.set_xlabel('Within-patient Spearman ρ');ax.set_title(title,loc='left',fontweight='bold');clean(ax)
axs[0].text(0,-.115,f'n = {len(pt)} patients; {(pt.effect>0).sum()}/{len(pt)} positive\nMedian ρ = {st.median_effect:.3f}; P = {st.wilcoxon_p_two_sided:.4f}; FDR = {st.wilcoxon_fdr_bh_audit18:.4f}',transform=axs[0].transAxes,fontsize=8,va='top')
axs[1].text(0,-.115,f'n = {len(cnv)} patients; {(cnv.rho>0).sum()}/{len(cnv)} positive\nMedian ρ = {ce.median_patient_rho_CNV:.3f}; P = {ce.wilcoxon_greater_p_CNV:.3f}; FDR = {ce.wilcoxon_greater_fdr_CNV:.3f}',transform=axs[1].transAxes,fontsize=8,va='top')
fig.subplots_adjust(left=.15,right=.985,top=.94,bottom=.22,wspace=.59);save(fig,'07_F1_TYMP_and_CNV_patient_evidence_v1')

# Display all 17 existing significant up-list pathways; explicitly state top-250 selection.
from scipy.stats import hypergeom
rora=read(r/'cnv_clone_program_function_v1/F1_high_clone_DE_pathway_ORA_v1.tsv');de=read(r/'cnv_clone_program_function_v1/F1_high_clone_vs_other_paired_pseudobulk_DE_v1.tsv.gz')
checks=[]
for direction,d in rora.groupby('direction',sort=False):
 p=hypergeom.sf(d.overlap_size.to_numpy()-1,len(de),d.pathway_size.to_numpy(),d.selected_size.to_numpy());order=np.argsort(p);adj=np.empty(len(p));adj[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);err=float(np.max(np.abs(adj-d.FDR.to_numpy())));assert err<1e-10,err
 checks.append({'direction':direction,'tested_pathways':len(d),'FDR_significant':int((d.FDR<.05).sum()),'min_FDR':float(d.FDR.min()),'independent_hypergeometric_BH_max_error':err})
up=rora[(rora.direction=='up')&(rora.FDR<.05)].sort_values('FDR');assert len(up)==17
fig,ax=plt.subplots(figsize=(9.0,7.5));xx=-np.log10(up.FDR.to_numpy());ys=np.arange(len(up));ax.hlines(ys,0,xx,color='#D8DEEA',lw=2);ax.scatter(xx,ys,s=20+up.overlap_size.to_numpy()*6,c=xx,cmap=blue,vmin=0,vmax=4,edgecolors='#36547B',linewidths=.55,zorder=3);ax.set_yticks(ys,[textwrap.fill(s.split('::')[-1],45) for s in up.pathway],fontsize=8);ax.invert_yaxis();ax.set_xlim(0,max(xx)+.4);ax.axvline(-np.log10(.05),ls='--',c='#999999',lw=.7);ax.set_xlabel('−log10(FDR)');ax.set_title('F1-high clone: enrichment of the up-ranked gene set',loc='left',pad=14,fontsize=10);clean(ax)
fig.subplots_adjust(left=.51,right=.97,top=.92,bottom=.13);fig.text(.04,.045,'Top 250 genes ranked by edgeR P within each direction; 14 paired patients.\nUp: 17/201 pathways FDR < 0.05. Down: 0/201. Dot area reflects gene overlap.',fontsize=8,color='#444444');save(fig,'10_F1_high_clone_up_pathway_enrichment_v1')
(w/'qc/ORA_INDEPENDENT_RECHECK.json').write_text(json.dumps({'status':'PASS','tested_gene_universe':len(de),'selection_rule':'top 250 ranked by edgeR P separately within positive/negative logFC, not an FDR-filtered gene list','families':checks},indent=2))
old=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/local_records_migration_20260908_150455/project/phase116_sections_2_1_2_2_figure_package_local_v1/deliverable_v1/05_关键统计表/cnmf_F1_high_clone_DE_pathway_ORA_v1.tsv')
assert hashlib.sha256(old.read_bytes()).hexdigest()==hashlib.sha256((r/'cnv_clone_program_function_v1/F1_high_clone_DE_pathway_ORA_v1.tsv').read_bytes()).hexdigest()
(w/'qc/HISTORICAL_SUMMARY_CORRECTION.md').write_text('旧摘要“F1-high clone pathway ORA全部FDR=1”不正确。原始服务器表与Phase116历史打包表SHA256完全一致，均为上调17/201条FDR<0.05、下调0/201条。分别按原方向家族独立复算超几何P和BH FDR通过。基因集是各方向按edgeR P排序的前250基因，不是只纳入FDR显著基因。此次没有重跑DE或富集，也没有更换表；纠正的是历史摘要遗漏上调分支的表述。',encoding='utf8')
val=json.loads((w/'qc/VALIDATION.json').read_text());val['CNV_patient_n']=len(cnv);val['CNV_positive_patients']=int((cnv.rho>0).sum());val['clone_DE_ORA_up_FDR_significant']=17;val['clone_DE_ORA_down_FDR_significant']=0;val['historical_summary_corrected']=True;val['figure_files']=20;(w/'qc/VALIDATION.json').write_text(json.dumps(val,indent=2));print(json.dumps(checks),flush=True)
