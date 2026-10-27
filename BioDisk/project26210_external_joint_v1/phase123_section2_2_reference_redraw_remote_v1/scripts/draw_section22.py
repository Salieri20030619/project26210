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
# Fixed rank selection. Reuse original runs, no new model fit.
fig,axs=plt.subplots(1,2,figsize=(8.8,3.4));a=axs[0]
a.plot(rank['rank'],rank.mean_component_silhouette,'o-',c='#7295BE',lw=1.3,label='Mean silhouette');a.plot(rank['rank'],rank.min_run_fraction,'s-',c='#BA7868',lw=1.3,label='Minimum run coverage');a.axvspan(6.84,7.16,color='#EEE3C8',zorder=0);a.set(xlabel='Number of programs (k)',ylabel='Stability',xticks=rank['rank'],ylim=(0,1.05));a.legend(frameon=False,fontsize=7,loc='lower left');a.set_title('Twenty starts per rank',loc='left');clean(a)
a=axs[1];a.errorbar(rank['rank'],rank.median_run_mse,yerr=rank.iqr_run_mse/2,fmt='o-',lw=1.1,c='#79A697',capsize=3);a.axvspan(6.84,7.16,color='#EEE3C8',zorder=0);a.set(xlabel='Number of programs (k)',ylabel='Reconstruction MSE',xticks=rank['rank']);a.set_title('Frozen selection: k = 7, candidate F1',loc='left');clean(a);fig.tight_layout(pad=2);save(fig,'01_cNMF_rank_selection_v1')
# Select program genes by loading only; do not label these as differential tests.
genes=[];gene_factor=[]
for f in factors:
 choices=[g for g in load[f].sort_values(ascending=False,kind='stable').index if g not in genes][:5]
 genes.extend(choices);gene_factor.extend([f]*len(choices))
i_genes=genes.copy();i_gf=gene_factor.copy()
if 'TYMP' not in genes:genes.append('TYMP');gene_factor.append('TYMP')
hpath=base/'02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad';a=ad.read_h5ad(hpath,backed='r');idx=a.obs_names.get_indexer(u.index);assert (idx>=0).all()
xy=np.asarray(a.obsm['X_umap'])[idx];ix=np.argsort(idx);xx=a[idx[ix],genes].X.toarray();expr=np.empty_like(xx);expr[ix]=xx;assert np.isfinite(expr).all();a.file.close()
np.savez_compressed(w/'data'/'plot_values_server_only.npz',xy=xy,expression=expr,genes=np.array(genes),cells=np.array(u.index,dtype=str));print('RNA_AND_COORDINATES_READY',flush=True)
# Frozen UMAP: program labels, F1 usage, TYMP RNA.
fig,axs=plt.subplots(1,3,figsize=(10.8,4.1))
for f in factors:
 ii=m.factor==f;axs[0].scatter(xy[ii,0],xy[ii,1],s=2.1,c=fd[f],linewidths=0,rasterized=True,label=f)
axs[0].set_title('Dominant cNMF program',loc='left');axs[0].legend(frameon=False,ncol=4,markerscale=3,loc='upper center',bbox_to_anchor=(.5,-.015),fontsize=7,handletextpad=.15,columnspacing=.65)
for ax,v,title,cmap in [(axs[1],u.F1.to_numpy(),'F1 usage',blue),(axs[2],expr[:,genes.index('TYMP')],'TYMP RNA',div)]:
 order=np.argsort(v,kind='stable');vmax=float(np.quantile(v[v>0],.99)) if (v>0).any() else 1
 cmap=blue if title=='F1 usage' else LinearSegmentedColormap.from_list('tymp',['#E9E9E9','#A45654'])
 im=ax.scatter(xy[order,0],xy[order,1],c=v[order],cmap=cmap,s=2.1,vmin=0,vmax=vmax,linewidths=0,rasterized=True);ax.set_title(title,loc='left');cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.038,pad=.045,shrink=.65);cb.set_label('Normalized usage' if title=='F1 usage' else 'Log1p expression',fontsize=7)
for ax in axs:ax.set_aspect('equal');ax.axis('off')
fig.subplots_adjust(left=.02,right=.985,top=.90,bottom=.19,wspace=.22);save(fig,'02_cNMF_F1_TYMP_UMAP_v1')
# E: all retained cells, ordered by existing dominant program and usage.
order=m.assign(_factor=pd.Categorical(m.factor,categories=factors,ordered=True)).sort_values(['_factor','max_usage'],ascending=[True,False],kind='stable').index
positions=u.index.get_indexer(order);z=(expr-expr.mean(axis=0))/np.where(expr.std(axis=0)>0,expr.std(axis=0),1)
fig=plt.figure(figsize=(11,7.8));ax=fig.add_axes([.10,.09,.70,.77]);im=ax.imshow(z[positions].T,aspect='auto',interpolation='nearest',cmap=div,vmin=-2,vmax=2,rasterized=True)
ax.set_yticks(range(len(genes)),genes,fontsize=7.3);ax.yaxis.tick_right();ax.tick_params(axis='y',length=0,pad=4);ax.set_xticks([]);ax.set_xlabel('11,950 malignant cells ordered by dominant program',labelpad=8)
ordered=m.loc[order];bounds=np.r_[0,ordered.factor.value_counts().reindex(factors,fill_value=0).cumsum().to_numpy()];centers=(bounds[:-1]+bounds[1:]-1)/2
for b in bounds[1:-1]:ax.axvline(b-.5,c='white',lw=.65)
top=fig.add_axes([.10,.88,.70,.019]);top.imshow(np.array([factors.index(x) for x in ordered.factor])[None,:],aspect='auto',cmap=ListedColormap(colors),vmin=0,vmax=6,interpolation='nearest');top.set_xticks(centers,factors,fontsize=8);top.xaxis.tick_top();top.tick_params(length=0);top.set_yticks([])
ann=fig.add_axes([.10,.859,.70,.015]);ann.imshow(np.array([studies.index(x) for x in ordered.study])[None,:],aspect='auto',cmap=ListedColormap(list(sd.values())),vmin=0,vmax=4,interpolation='nearest');ann.axis('off')
ca=fig.add_axes([.10,.94,.19,.018]);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-2,0,2]);ca.set_title('RNA expression (gene z-score)',fontsize=8,pad=4);cb.ax.tick_params(labelsize=7,length=2)
fig.legend(handles=[Patch(facecolor=sd[s],label=s) for s in studies],loc='upper right',bbox_to_anchor=(.985,.975),frameon=False,ncol=3,fontsize=7,handlelength=.9,columnspacing=.9);letter(fig,'E');save(fig,'03_E_program_gene_expression_heatmap_v1')
# H: actual 20 x 7 run components, not synthetic correlation blocks.
components=read(w/'data/rank7_run_component_loadings.tsv.gz').set_index('gene');cm=read(w/'data/rank7_component_labels.tsv').set_index('component_id').loc[components.columns];assert components.shape==(2500,140)
corr=np.corrcoef(components.values.T);corr=np.clip(corr,-1,1);np.fill_diagonal(corr,1);tree=linkage(squareform(np.maximum(0,1-corr),checks=False),method='average',optimal_ordering=True)
fig=plt.figure(figsize=(7.2,7.4));ta=fig.add_axes([.14,.78,.72,.13]);dd=dendrogram(tree,ax=ta,no_labels=True,color_threshold=0,above_threshold_color='#AAAAAA',link_color_func=lambda k:'#9C9C9C');ta.axis('off');leaves=dd['leaves'];groups=cm.consensus_cluster.to_numpy()[leaves].astype(int)
ax=fig.add_axes([.14,.12,.72,.65]);im=ax.imshow(corr[np.ix_(leaves,leaves)],cmap=div,vmin=-1,vmax=1,interpolation='nearest',aspect='auto');ax.set_xticks([]);ax.set_yticks([])
starts=np.r_[0,np.flatnonzero(groups[1:]!=groups[:-1])+1];ends=np.r_[starts[1:],140]
for lo,hi in zip(starts,ends):
 ax.add_patch(Rectangle((lo-.5,lo-.5),hi-lo,hi-lo,fill=False,edgecolor='#333333',lw=.8))
 if hi-lo>=5:ax.text((lo+hi-1)/2,142.5,'F'+str(groups[lo]),ha='center',va='top',fontsize=8)
ax.set_xlabel('Original consensus labels; 140 run components',labelpad=19)
for rect,arr in [([.14,.771,.72,.009],groups[None,:]-1),([.124,.12,.011,.65],(groups-1)[:,None])]:
 aa=fig.add_axes(rect);aa.imshow(arr,cmap=ListedColormap(colors),vmin=0,vmax=6,aspect='auto',interpolation='nearest');aa.axis('off')
ca=fig.add_axes([.60,.94,.24,.018]);fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-1,0,1]);ca.set_title('Loading-profile Pearson r',fontsize=8);letter(fig,'H');save(fig,'04_H_consensus_component_correlation_v1')
pd.DataFrame(corr,index=components.columns,columns=components.columns).to_csv(w/'data/H_component_correlation.tsv.gz',sep='\t')
# I: reference blue, exact gene contribution matrix and source factor groups.
vals=load.loc[i_genes,factors];norm=vals.div(vals.max(axis=1),axis=0);assert np.isfinite(norm).all().all()
fig=plt.figure(figsize=(5.3,9));ax=fig.add_axes([.27,.075,.39,.84]);im=ax.imshow(norm,cmap=blue,vmin=0,vmax=1,aspect='auto',interpolation='nearest');ax.set_xticks(range(7),factors,fontsize=9);ax.set_yticks(range(len(i_genes)),i_genes,fontsize=8.5,fontstyle='italic');ax.yaxis.tick_right();ax.tick_params(length=0,pad=4);ax.set_xlabel('cNMF program',labelpad=8)
ax.set_xticks(np.arange(-.5,7,1),minor=True);ax.set_yticks(np.arange(-.5,len(i_genes),1),minor=True);ax.grid(which='minor',color='#555555',lw=.45);ax.tick_params(which='minor',length=0)
for j,f in enumerate(factors):
 lo=j*5-.5;hi=(j+1)*5-.5;ax.plot([-.78,-.78],[lo,hi],clip_on=False,c=fd[f],lw=2.2);ax.text(-1,lo+2.5,f+' genes',ha='right',va='center',fontsize=8,color='#444444')
ca=fig.add_axes([.27,.95,.39,.015]);fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[0,.5,1]);ca.set_title('Relative gene contribution',fontsize=8);letter(fig,'I');save(fig,'05_I_gene_contribution_reference_blue_v1')
norm.to_csv(w/'qc/I_gene_contribution_aggregate.tsv',sep='\t')
# J: counts and study composition; retain concentration of F1 in GSE207493.
ct=pd.crosstab(m.factor,m.study).reindex(index=factors,columns=studies,fill_value=0);frac=ct.div(ct.sum(axis=1),axis=0);assert int(ct.values.sum())==11950
f1frac=float(frac.loc['F1','GSE207493']);assert abs(f1frac-float(select.candidate_dominant_study_fraction))<1e-10
fig,axs=plt.subplots(1,2,figsize=(9.4,3.8));bottom=np.zeros(7)
for s in studies:axs[0].bar(factors,ct[s],bottom=bottom,color=sd[s],width=.73);bottom+=ct[s].to_numpy()
axs[0].set_ylabel('Cells');axs[0].set_title('Program-dominant cell counts',loc='left');clean(axs[0]);bottom=np.zeros(7)
for s in studies:axs[1].bar(factors,frac[s],bottom=bottom,color=sd[s],width=.73,label=s);bottom+=frac[s].to_numpy()
axs[1].set(ylim=(0,1.04),ylabel='Study fraction');axs[1].set_title('F1: '+format(f1frac*100,'.2f')+'% GSE207493',loc='left');clean(axs[1]);fig.legend(handles=[Patch(facecolor=sd[s],label=s) for s in studies],loc='lower center',ncol=5,frameon=False,fontsize=7,bbox_to_anchor=(.51,.005));fig.subplots_adjust(left=.08,right=.985,top=.88,bottom=.19,wspace=.30);save(fig,'06_J_program_study_composition_v1');ct.to_csv(qc/'program_study_cell_counts.tsv',sep='\t')
# Formal F1-TYMP audit and frozen F1-CNV endpoint shown side by side.
cnv=pcorr[(pcorr['rank']==7)&(pcorr.factor=='F1')&(pcorr.outcome=='infercnv_burden_over_threshold_v1')].copy();cnv['study']=cnv.curated_patient_id.str.split('__').str[0];ce=evidence[(evidence['rank']==7)&(evidence.factor=='F1')].iloc[0]
fig,axs=plt.subplots(1,2,figsize=(10.1,6.2))
for ax,d,value,title in [(axs[0],pt,'effect','F1–TYMP RNA'),(axs[1],cnv,'rho','F1–CNV burden')]:
 d=d.sort_values(['study','curated_patient_id'],kind='stable').reset_index(drop=True);ys=np.arange(len(d));ax.axvline(0,c='#AAAAAA',ls='--',lw=.8);ax.hlines(ys,0,d[value],color=[sd[s] for s in d.study],lw=1.2);ax.scatter(d[value],ys,c=[sd[s] for s in d.study],s=28,zorder=3,edgecolors='white',linewidths=.45);ax.set_yticks(ys,[x.split('__',1)[1] for x in d.curated_patient_id],fontsize=7);ax.invert_yaxis();ax.set_xlim(-1,1);ax.set_xlabel('Within-patient Spearman ρ');ax.set_title(title,loc='left',fontweight='bold');clean(ax)
axs[0].text(0,-.115,f'n = {len(pt)} patients; {(pt.effect>0).sum()}/{len(pt)} positive\nMedian ρ = {st.median_effect:.3f}; P = {st.wilcoxon_p_two_sided:.4f}; FDR = {st.wilcoxon_fdr_bh_audit18:.4f}',transform=axs[0].transAxes,fontsize=8,va='top')
axs[1].text(0,-.115,f'n = {len(cnv)} patients; {(cnv.rho>0).sum()}/{len(cnv)} positive\nMedian ρ = {ce.median_patient_rho_CNV:.3f}; P = {ce.wilcoxon_greater_p_CNV:.3f}; FDR = {ce.wilcoxon_greater_fdr_CNV:.3f}',transform=axs[1].transAxes,fontsize=8,va='top')
fig.subplots_adjust(left=.15,right=.985,top=.94,bottom=.22,wspace=.59);save(fig,'07_F1_TYMP_and_CNV_patient_evidence_v1')
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
sm=plt.cm.ScalarMappable(norm=plt.Normalize(0,max(5,-np.log10(ora.loc[ora.FDR>0,'FDR'].min()))),cmap=blue);fig.colorbar(sm,ax=ax,fraction=.035,pad=.05,label='−log10(FDR)');fig.subplots_adjust(left=.39,right=.91,bottom=.14,top=.89);save(fig,'08_program_function_enrichment_v1')
# Existing paired clone DE: show genuine differential evidence rather than calling loadings DE.
de=read(r/'cnv_clone_program_function_v1/F1_high_clone_vs_other_paired_pseudobulk_DE_v1.tsv.gz');deltas=read(r/'cnv_clone_program_function_v1/top_DE_patient_logCPM_deltas_v1.tsv.gz');deora=read(r/'cnv_clone_program_function_v1/F1_high_clone_DE_pathway_ORA_v1.tsv')
sel=de[de.gene.isin(deltas.gene)].sort_values(['edgeR_FDR','gene']).head(24).gene.tolist();hm=deltas[deltas.gene.isin(sel)].pivot(index='gene',columns='curated_patient_id',values='logCPM_delta').reindex(sel)
fig,(ax,ha)=plt.subplots(1,2,figsize=(11.2,6.4),gridspec_kw={'width_ratios':[.95,1.1]});sig=de.edgeR_FDR<.05
ax.scatter(de.edgeR_logFC,-np.log10(de.edgeR_FDR.clip(lower=1e-300)),s=4,c=np.where(sig,'#BA7868','#C6CCD4'),alpha=.65,linewidths=0,rasterized=True);ax.axhline(-np.log10(.05),c='#777777',ls='--',lw=.7);ax.axvline(0,c='#BBBBBB',lw=.6);ax.set(xlabel='F1-high clone vs other: log2 fold change',ylabel='−log10(edgeR FDR)');ax.set_title('Patient-blocked differential expression',loc='left');clean(ax)
im=ha.imshow(hm,cmap=div,vmin=-3,vmax=3,aspect='auto',interpolation='nearest');ha.set_yticks(range(len(sel)),sel,fontsize=7);ha.yaxis.tick_right();ha.set_xticks(range(hm.shape[1]),[x.split('__',1)[1] for x in hm.columns],rotation=90,fontsize=6.5);ha.tick_params(length=0);ha.set_title('Paired expression difference by patient',loc='left');cb=fig.colorbar(im,ax=ha,orientation='horizontal',fraction=.035,pad=.19);cb.set_label('F1-high clone − other: logCPM',fontsize=7)
fig.subplots_adjust(left=.08,right=.90,top=.91,bottom=.21,wspace=.32);save(fig,'09_F1_high_clone_differential_evidence_v1')
# Save complete evidence tables on server. No recomputed inferential P/FDR.
summary={'status':'PASS_DATA_AND_RENDER_CHECKS','cells':len(u),'genes_input':len(load),'rank':7,'runs':20,'components':components.shape[1],'no_cnmf_refit':True,'umap_coordinates_unchanged':True,'selected_genes':genes,'TYMP_in_loading_matrix':'TYMP' in load.index,'F1_dominant_cells':int(ct.loc['F1'].sum()),'F1_GSE207493_fraction':f1frac,'F1_TYMP_patient_n':len(pt),'F1_TYMP_positive':int((pt.effect>0).sum()),'F1_TYMP_median_rho':float(st.median_effect),'F1_TYMP_P':float(st.wilcoxon_p_two_sided),'F1_TYMP_FDR_audit18':float(st.wilcoxon_fdr_bh_audit18),'CNV_patient_n':len(cnv),'CNV_median_rho':float(ce.median_patient_rho_CNV),'CNV_FDR_rank_screen':float(ce.wilcoxon_greater_fdr_CNV),'clone_DE_patients':hm.shape[1],'clone_DE_ORA_FDR_significant':int((deora.FDR<.05).sum()),'figure_files':len(figlist),'versions':{'python':sys.version,'matplotlib':matplotlib.__version__,'numpy':np.__version__,'pandas':pd.__version__,'anndata':ad.__version__}}
(qc/'VALIDATION.json').write_text(json.dumps(summary,indent=2));pd.DataFrame(gene_factor,index=genes,columns=['selection_program']).to_csv(qc/'GENE_SELECTION.tsv',sep='\t')
pt.to_csv(qc/'F1_TYMP_patients.tsv',sep='\t',index=False);pd.DataFrame([st]).to_csv(qc/'F1_TYMP_formal_summary.tsv',sep='\t',index=False);pd.DataFrame([ce]).to_csv(qc/'F1_CNV_rank_screen_summary.tsv',sep='\t',index=False)
source_manifest=[]
for f in sorted(set(source_files)):
 source_manifest.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size})
(qc/'SOURCE_MANIFEST.json').write_text(json.dumps(source_manifest,indent=2));print(json.dumps(summary),flush=True)
