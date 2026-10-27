from pathlib import Path
import json,sys
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from matplotlib.ticker import MaxNLocator
w=Path(__file__).resolve().parents[1];q=w/'qc'
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory')
def read(p):return pd.read_csv(p,sep='\t')
def indexed(p):return pd.read_csv(p,sep='\t',index_col=0)
assert json.loads((q/'INDEPENDENT_VALIDATION.json').read_text())['status']=='PASS'
m=indexed(b/'trajectory_v1/input/cell_metadata_v1.tsv.gz');u=indexed(b/'trajectory_v1/input/scvi_umap_v1.tsv.gz').loc[m.index]
a=indexed(b/'metabolism_v1/source_data/four_pathway_AUCell_cell_scores_v1.tsv.gz').loc[m.index]
f=indexed(b/'metabolism_v1/source_data/scFEA_four_axis_cell_scores_v1.tsv.gz').loc[m.index]
tr=indexed(b/'trajectory_v1/source_data/trajectory_all_methods_cell_scores_v1.tsv.gz').loc[m.index]
e=indexed(w/'data/ENZYME_LOG1P_CP10K.tsv.gz').loc[m.index]
curves=read(w/'data/PATIENT_FITTED_CURVES.tsv.gz');states=read(w/'data/PATIENT_STATE_COMPOSITION.tsv.gz')
ts=read(b/'trajectory_v1/tables/trajectory_all_methods_correlation_summary_v1.tsv')
auc=['AUC_Pyrimidine','AUC_CysMet','AUC_ArgPro','AUC_Glycerophospholipid'];flux=['Pyrimidine_flux','CysMet_flux','ArgPro_flux','Choline_proxy']
genes=['TYMP','TK1','SLC29A1','SLC29A2','ODC1','SRM','SMS','MTAP','CHKA','PCYT1A','PCYT2','CEPT1']
methods=['Monocle3','Slingshot','scTour'];tcols=[x+'_pseudotime_patient_v1' for x in methods]
labels=dict(zip(auc,['Pyrimidine','Cysteine / methionine','Arginine / proline','Glycerophospholipid']))
labels.update(dict(zip(flux,['Pyrimidine flux proxy','Cys / Met flux proxy','Arg / Pro flux proxy','Choline proxy (M21 / M83)'])))
labels.update({g:g for g in genes});labels.update({'F1':'F1 score','TYMP_log1p_CP10K_v1':'TYMP RNA','infercnv_burden_over_threshold_v1':'CNV burden','progression_anchor_v1':'Root / orientation anchor'})
statecolors={'Base':'#8CA8B6','Hypoxia':'#CB8B4B','Cycle':'#8766AB','TYMP_state':'#C94C64'}
blue='#367F9C';red='#C85269';gray='#B2BAC2'
seq=LinearSegmentedColormap.from_list('activity',['#E9EDF1','#B5C9D3','#608BA5','#3A477D','#702F6B'])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':11,'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':300})
written=[]
def fmt(x):return f'{x:.2g}' if x<.001 else f'{x:.3f}'
def save(fig,section,name,title,foot):
 fig.suptitle(title,x=.035,ha='left',y=.985,fontsize=15,fontweight='bold',color='#243442')
 fig.text(.035,.017,foot,ha='left',va='bottom',fontsize=8,color='#4E5963')
 fig.tight_layout(rect=[.015,.065,.99,.94],h_pad=2,w_pad=2)
 dest=w/('figures_'+section)
 for ext in ['png','pdf']:fig.savefig(dest/(name+'.'+ext),facecolor='white',dpi=300)
 written.append({'section':section,'name':name,'title':title,'footnote':foot});plt.close(fig);print('SAVED',section,name,flush=True)
def umap(ax,v,title,label,cmap=seq,limits=None,diverge=False):
 v=np.asarray(v,float);order=np.argsort(np.nan_to_num(v,nan=-np.inf));kw={}
 if limits is not None:kw.update(vmin=limits[0],vmax=limits[1])
 if diverge:
  lim=max(abs(np.nanmin(v)),abs(np.nanmax(v)));kw={'norm':TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim)};cmap='RdBu_r'
 artist=ax.scatter(u.iloc[order,0],u.iloc[order,1],c=v[order],s=2.1,cmap=cmap,rasterized=True,linewidths=0,**kw)
 ax.set_title(title,loc='left',pad=9);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
 ax.set_xlim(u.iloc[:,0].min()-.4,u.iloc[:,0].max()+.4);ax.set_ylim(u.iloc[:,1].min()-.4,u.iloc[:,1].max()+.4)
 cb=ax.figure.colorbar(artist,ax=ax,fraction=.045,pad=.025,shrink=.72);cb.set_label(label,fontsize=8);cb.ax.tick_params(labelsize=7);cb.locator=MaxNLocator(4);cb.update_ticks()
 return artist
# 2.6 A: exact AUCell values on the frozen malignant UMAP.
fig,axs=plt.subplots(2,3,figsize=(13,8.5))
for ax,c,title in zip(axs.flat,['F1','TYMP_log1p_CP10K_v1']+auc,['A  F1 program','B  TYMP RNA']+[chr(67+i)+'  '+labels[x] for i,x in enumerate(auc)]):
 umap(ax,m[c],title,'F1 usage' if c=='F1' else 'log1p CP10K' if c.startswith('TYMP') else 'AUCell score')
save(fig,'2_6','01_A_four_pathway_UMAP_v1','Metabolic pathway activity across malignant epithelial cells','11,950 malignant cells | 26 patients | frozen UMAP | color = observed score; pathway-specific full color ranges')
def paired(key,cols,number):
 prefix='four_pathway_AUCell' if key=='AUCell' else 'scFEA_four_axis'
 ds=read(b/('metabolism_v1/tables/'+prefix+'_patient_deltas_v1.tsv'));ss=read(b/('metabolism_v1/tables/'+prefix+'_paired_statistics_v1.tsv'))
 fig,axs=plt.subplots(2,4,figsize=(15,8),sharey='col')
 contrasts=['TYMP-positive_vs_negative','F1-high_clone_vs_other']
 for row,contrast in enumerate(contrasts):
  for j,c in enumerate(cols):
   ax=axs[row,j];d=ds[(ds.pathway==c)&(ds.contrast==contrast)].sort_values('curated_patient_id');r=ss[(ss.pathway==c)&(ss.contrast==contrast)].iloc[0]
   for rr in d.itertuples():ax.plot([0,1],[rr.low_score,rr.high_score],color=gray,lw=.7,alpha=.7,zorder=1)
   ax.scatter(np.zeros(len(d)),d.low_score,color=blue,s=20,zorder=2);ax.scatter(np.ones(len(d)),d.high_score,color=red,s=20,zorder=2)
   ax.scatter([0,1],[d.low_score.median(),d.high_score.median()],marker='_',s=420,color='#202832',linewidths=2.5,zorder=3)
   ax.set_xticks([0,1],['TYMP-','TYMP+'] if row==0 else ['Other clones','F1-high clone']);ax.set_xlim(-.35,1.35)
   ax.set_title(labels[c],loc='left',fontsize=10,pad=40)
   ax.text(0,1.025,f'n={int(r.n_patients)} | positive {int(round(r.n_patients*r.fraction_patients_positive))}/{int(r.n_patients)}\nmedian delta={r.median_patient_delta:.3g}; P={fmt(r.paired_wilcox_PValue)}; q={fmt(r.FDR)}',transform=ax.transAxes,fontsize=8,va='bottom')
   if j==0:ax.set_ylabel('Patient median AUCell' if key=='AUCell' else 'Patient median module-z proxy')
 save(fig,'2_6',f'{number:02d}_patient_paired_{key}_v1',('Three pathways increase in TYMP+ cells and F1-high clones' if key=='AUCell' else 'scFEA predicts changes in transcriptome-derived flux proxies'),'Each line = one patient; >=20 cells per arm | two-sided paired Wilcoxon | BH across four pathways within each contrast\n'+('TYMP+ = raw UMI >0; TYMP- = raw UMI =0.' if key=='AUCell' else 'scFEA axes are module-z aggregates; choline is a partial proxy, and MTA is not an explicit model output.'))
paired('AUCell',auc,2)
def corrplot(key,cols,number):
 prefix='four_pathway_AUCell' if key=='AUCell' else 'scFEA_four_axis'
 ss=read(b/('metabolism_v1/tables/'+prefix+'_correlation_summary_v1.tsv'))
 predictors=['F1','TYMP_log1p_CP10K_v1'];mat=np.zeros((4,2));fig,ax=plt.subplots(figsize=(9,7))
 for i,c in enumerate(cols):
  for j,p in enumerate(predictors):
   r=ss[(ss.pathway==c)&(ss.predictor==p)].iloc[0];mat[i,j]=r.median_patient_rho
 im=ax.imshow(mat,cmap='RdBu_r',vmin=-.3,vmax=.3,aspect='auto')
 for i,c in enumerate(cols):
  for j,p in enumerate(predictors):
   r=ss[(ss.pathway==c)&(ss.predictor==p)].iloc[0]
   ax.text(j,i,f'rho = {r.median_patient_rho:+.3f}\nP = {fmt(r.signed_rank_PValue)}; q = {fmt(r.FDR)}\npositive {int(round(r.n_patients*r.fraction_positive))}/{int(r.n_patients)}',ha='center',va='center',fontsize=9,color='white' if abs(mat[i,j])>.2 else '#243442')
 ax.set_xticks([0,1],['F1 score','TYMP RNA']);ax.set_yticks(range(4),[labels[c] for c in cols]);ax.tick_params(length=0)
 fig.colorbar(im,ax=ax,fraction=.035,pad=.03,label='Median patient Spearman rho')
 save(fig,'2_6',f'{number:02d}_patient_correlations_{key}_v1',('F1 is positively associated with three metabolic pathways' if key=='AUCell' else 'F1 is positively associated with four scFEA proxy axes'),'20 independent patients | >=30 cells per patient | two-sided signed-rank test on patient rho\nBH across four pathways separately for each predictor; color limits fixed at +/-0.3.')
corrplot('AUCell',auc,3)
# Enzyme DotPlot: no gene-wise scaling; patients contribute equally within each quartile.
dot=read(q/'ENZYME_DOTPLOT_PATIENT_PROFILES.tsv');dg=dot.groupby(['gene','quartile']).agg(mean_log1p=('mean_log1p','mean'),fraction=('fraction_detected','mean')).reset_index()
fig,ax=plt.subplots(figsize=(13.5,5.2));vmax=dg.mean_log1p.max()
for j,gene in enumerate(genes):
 for k in range(1,5):
  r=dg[(dg.gene==gene)&(dg.quartile==k)].iloc[0]
  sc=ax.scatter(j,k-1,s=r.fraction*500,c=[r.mean_log1p],cmap=seq,vmin=0,vmax=vmax,edgecolors='#5A6772',linewidths=.35)
ax.set_xticks(range(12),genes,rotation=35,ha='right');ax.set_yticks(range(4),['Q1 - low F1','Q2','Q3','Q4 - high F1']);ax.set_xlim(-.65,11.65);ax.set_ylim(-.65,3.65)
for cut in [3.5,7.5]:ax.axvline(cut,color='#DCE1E6',lw=1)
for center,label in [(1.5,'Thymidine / pyrimidine'),(5.5,'Polyamine / MTA'),(9.5,'Kennedy pathway')]:ax.text(center,3.9,label,ha='center',fontsize=10)
fig.colorbar(sc,ax=ax,pad=.025,fraction=.035,label='Patient-equal mean log1p CP10K')
handles=[ax.scatter([],[],s=v*500,color='#CBD4DC',edgecolors='#5A6772',linewidths=.35,label=f'{v*100:.0f}%') for v in [.1,.3,.6,1]]
ax.legend(handles=handles,title='Detected cells',loc='center left',bbox_to_anchor=(1.24,.5),frameon=False,labelspacing=1.6)
save(fig,'2_6','04_enzyme_F1_quartile_dotplot_v1','Enzyme and transporter expression across the F1 continuum',f'{dot.patient.nunique()} patients | within-patient F1 rank quartiles; tied scores stay together; >=5 cells per quartile\nDot area = patient-equal detected fraction; color = unscaled mean log1p RNA; RNA expression does not measure enzyme flux.')
gs=read(q/'ENZYME_F1_STATISTICS.tsv').set_index('gene').loc[genes];gp=read(q/'ENZYME_F1_PATIENT_EFFECTS.tsv')
fig,ax=plt.subplots(figsize=(11,7.5));rng=np.random.default_rng(26210)
for i,g in enumerate(genes):
 d=gp.loc[gp.gene==g,'rho'].dropna().values;r=gs.loc[g]
 ax.scatter(d,i+rng.uniform(-.13,.13,len(d)),s=18,color=blue,alpha=.6)
 ax.scatter([r.median_rho],[i],marker='D',s=44,color=red,zorder=5)
 ax.text(1.03,i,f'{r.median_rho:+.3f}   {int(r.positive_patients)}/{int(r.n_patients)}   P={fmt(r.P)}   q={fmt(r.FDR)}',transform=ax.get_yaxis_transform(),fontsize=8,va='center')
ax.axvline(0,color=gray,lw=1);ax.set_yticks(range(12),genes);ax.invert_yaxis();ax.set_xlim(-.75,.85);ax.set_xlabel('Within-patient Spearman rho with F1')
ax.text(1.03,1.025,'Median rho | positive / n | patient-level P and q',transform=ax.transAxes,fontsize=8)
fig.subplots_adjust(right=.6)
save(fig,'2_6','05_enzyme_F1_patient_associations_v1','Patient-level evidence for enzyme and transporter co-expression','Dots = independent patients; diamond = median rho | >=30 cells per patient | two-sided signed-rank; BH family = 12 prespecified genes')
fig,axs=plt.subplots(2,2,figsize=(10,8.5))
for ax,c in zip(axs.flat,flux):umap(ax,f[c],labels[c],'Mean module z score',diverge=True)
save(fig,'2_6','06_scFEA_proxy_UMAP_v1','scFEA flux proxies vary across malignant epithelial cells','11,950 malignant cells | original scFEA outputs; no model refitting | symmetric full color range around zero per axis\nCholine proxy is partial; MTA is not modeled explicitly. These values are not measured metabolite concentrations.')
paired('scFEA',flux,7);corrplot('scFEA',flux,8)
# Complete exploratory screen: display all pathways, not only nominal positives.
ks=read(q/'ALL85_TYMP_PATIENT_STATISTICS.tsv').sort_values(['FDR','pathway']).reset_index(drop=True)
fig,axs=plt.subplots(1,2,figsize=(16,16))
for ax,part in zip(axs,[ks.iloc[:43],ks.iloc[43:]]):
 y=np.arange(len(part));colors=np.where(part.FDR<.05,red,gray)
 ax.barh(y,part.median_delta,color=colors,height=.7)
 ax.set_yticks(y,[s.replace(' / ','/').replace(' metabolism',' metab.') for s in part.pathway],fontsize=7)
 ax.invert_yaxis();ax.axvline(0,color='#5D6974',lw=.7);ax.set_xlabel('Median paired AUCell delta (TYMP+ minus TYMP-)')
 for i,r in enumerate(part.itertuples()):ax.text(1.02,i,f'{int(r.positive_patients)}/{int(r.n_patients)}  q={fmt(r.FDR)}',transform=ax.get_yaxis_transform(),fontsize=7,va='center')
 ax.text(1.02,1.02,'positive / n; FDR',transform=ax.transAxes,fontsize=7)
lim=max(abs(ks.median_delta.min()),abs(ks.median_delta.max()))*1.1
for ax in axs:ax.set_xlim(-lim,lim)
save(fig,'2_6','09_S_all85_pathways_TYMP_screen_v1',f'Exploratory TYMP comparison: {(ks.FDR<.05).sum()} of 85 pathways pass FDR < 0.05','14 paired patients; >=20 cells per arm | all 85 KEGG pathways displayed | red = BH q<0.05 within the 85-pathway family\nThis exploratory family is separate from the four prespecified pathways; full P, q and patient deltas are retained in server statistics.')
# 2.7: five distinct methods in a frozen embedding.
fig,axs=plt.subplots(2,3,figsize=(14,9))
ax=axs[0,0]
for st,col in statecolors.items():
 ix=tr.malignant_state_v1==st;ax.scatter(u.loc[ix].iloc[:,0],u.loc[ix].iloc[:,1],s=2.1,color=col,label=st.replace('_',' '),linewidths=0,rasterized=True)
ax.set_title('A  Malignant expression states',loc='left');ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.legend(frameon=False,markerscale=4,fontsize=8,loc='lower left')
for ax,method,col,letter in zip([axs[0,1],axs[0,2],axs[1,0]],methods,tcols,['B','C','D']):
 umap(ax,tr[col],letter+'  '+method,'Patient-scaled pseudotime',limits=(0,1))
 if method=='Monocle3':
  seg=read(b/'trajectory_v1/source_data/monocle3_principal_graph_segments_v1.tsv')
  lines=np.stack([seg[['x','y']].values,seg[['xend','yend']].values],axis=1)
  ax.add_collection(LineCollection(lines,colors='#444A51',linewidths=.55,alpha=.7,zorder=4))
umap(axs[1,1],tr.CytoTRACE2_Score,'E  CytoTRACE2','Potency score (higher = more potent)',limits=(0,1))
umap(axs[1,2],tr.CellRank_max_fate_probability_v1,'F  CellRank','Maximum fate probability',limits=(0,1))
save(fig,'2_7','01_five_method_trajectory_overview_v1','Five complementary views of malignant epithelial state ordering','11,950 cells | 26 patients | frozen states, UMAP and models | Monocle3 overlay = original learned principal graph\nPseudotime, potency and fate probability have distinct meanings; cross-sectional ordering does not establish clinical progression.')
fig,axs=plt.subplots(1,3,figsize=(13,4.8))
for ax,st in zip(axs,['Cycle','Hypoxia','TYMP_state']):umap(ax,tr['CellRank_fate_'+st+'_v1'],st.replace('_',' ')+' fate','Fate probability',limits=(0,1))
save(fig,'2_7','02_CellRank_three_fates_v1','CellRank assigns distinct terminal-state probabilities','Original kernel: 0.8 pseudotime + 0.2 connectivity | probabilities sum to 1 in each cell\nTerminal states were defined from expression states; TYMP-state fate is not an independent validation of TYMP-driven progression.')
vec=read(b/'trajectory_v1/source_data/cellrank_pseudotime_transition_vectors_v1.tsv.gz').set_index('exact_cell_id').loc[m.index]
assert np.allclose(vec[['UMAP1','UMAP2']],u,atol=1e-6)
fig,ax=plt.subplots(figsize=(9,7.5));umap(ax,tr.CellRank_fate_TYMP_state_v1,'Pseudotime-informed transition field','TYMP-state fate probability',limits=(0,1))
xe=np.linspace(u.iloc[:,0].min(),u.iloc[:,0].max()+1e-6,29);ye=np.linspace(u.iloc[:,1].min(),u.iloc[:,1].max()+1e-6,23)
vg=vec.copy();vg['xb']=np.digitize(vg.UMAP1,xe);vg['yb']=np.digitize(vg.UMAP2,ye)
vg=vg.groupby(['xb','yb']).agg(x=('UMAP1','mean'),y=('UMAP2','mean'),dx=('dUMAP1','mean'),dy=('dUMAP2','mean'),n=('UMAP1','size'));vg=vg[vg.n>=10]
vg.to_csv(q/'CELLRANK_DISPLAY_VECTOR_BINS.tsv',sep='\t')
ax.quiver(vg.x,vg.y,vg.dx,vg.dy,color='#293642',angles='xy',scale_units='xy',scale=.4,width=.003,headwidth=3.2,zorder=5)
save(fig,'2_7','03_CellRank_transition_vectors_v1','CellRank transition directions on the malignant UMAP','Arrows = mean original transition vector in 28 x 22 spatial bins with >=10 cells; common display scale\nPseudotime-informed transitions; no spliced/unspliced layers, RNA velocity or lineage tracing.')
def curvefig(features,name,title,kind):
 fig,axs=plt.subplots(len(features),3,figsize=(13,2.45*len(features)+1.4),squeeze=False,sharex=True,sharey='row')
 for i,feature in enumerate(features):
  for j,method in enumerate(methods):
   ax=axs[i,j];sub=curves[(curves.feature==feature)&(curves.method==method)]
   wide=sub.pivot(index='rank',columns='patient',values='fitted');x=wide.index.values;v=wide.values
   med=np.median(v,axis=1);lo,hi=np.quantile(v,[.25,.75],axis=1)
   ax.fill_between(x,lo,hi,color=blue,alpha=.18,linewidth=0);ax.plot(x,med,color=blue,lw=2)
   ax.set_xlim(0,1);ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(4))
   if j==0:ax.set_ylabel(labels[feature]+'\n'+('AUCell score' if kind=='auc' else 'Module-z proxy' if kind=='flux' else 'log1p CP10K'))
   if i==0:ax.set_title(method,loc='left',fontweight='bold')
   if i==len(features)-1:ax.set_xlabel('Within-patient pseudotime rank')
   if kind in ['auc','flux']:
    rr=ts[(ts.target==feature)&(ts.method==method+'_pseudotime_patient_v1')].iloc[0]
    ax.text(.03,.96,f'rho={rr.median_patient_rho:+.3f}; P={fmt(rr.signed_rank_PValue)}\nq={fmt(rr.FDR)}; positive {int(round(rr.fraction_positive*rr.n_patients))}/{int(rr.n_patients)}',transform=ax.transAxes,va='top',fontsize=8,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':2})
 save(fig,'2_7',name,title,'20 patients | line = median of patient-specific fits; band = patient IQR, not confidence interval | same y scale per feature across methods\n20 patient-rank bins; fixed Gaussian bandwidth 0.075; '+('P/q use original unsmoothed patient correlations (BH: 12 targets per method).' if kind!='rna' else 'RNA trends are descriptive and do not establish metabolic flux or protein activity.'))
curvefig(auc,'04_four_pathway_pseudotime_curves_v1','Pathway trajectories show method-dependent directions','auc')
curvefig(flux,'05_scFEA_pseudotime_curves_v1','scFEA proxy trajectories differ across pseudotime methods','flux')
curvefig(genes[:4],'06_thymidine_axis_RNA_curves_v1','Thymidine-axis RNA expression along pseudotime','rna')
curvefig(genes[4:8],'07_polyamine_MTA_axis_RNA_curves_v1','Polyamine / MTA-axis RNA expression along pseudotime','rna')
curvefig(genes[8:],'08_Kennedy_axis_RNA_curves_v1','Kennedy-pathway RNA expression along pseudotime','rna')
targets=['progression_anchor_v1','F1','TYMP_log1p_CP10K_v1','infercnv_burden_over_threshold_v1']+auc+flux
meth=tcols+['CytoTRACE2_Score'];mat=np.array([[ts[(ts.target==target)&(ts.method==method)].iloc[0].median_patient_rho for method in meth] for target in targets])
fig,ax=plt.subplots(figsize=(11,10));im=ax.imshow(mat,cmap='RdBu_r',vmin=-.35,vmax=.35,aspect='auto')
for i,target in enumerate(targets):
 for j,method in enumerate(meth):
  r=ts[(ts.target==target)&(ts.method==method)].iloc[0]
  ax.text(j,i,f'{r.median_patient_rho:+.3f}\nP={fmt(r.signed_rank_PValue)}; q={fmt(r.FDR)}',ha='center',va='center',fontsize=7.5,color='white' if abs(r.median_patient_rho)>.25 else '#243442')
ax.set_xticks(range(4),methods+['CytoTRACE2\n(potency)']);ax.set_yticks(range(12),[labels[x] for x in targets]);ax.tick_params(length=0)
fig.colorbar(im,ax=ax,fraction=.035,pad=.025,label='Median patient Spearman rho')
save(fig,'2_7','09_patient_method_consistency_v1','Patient-level associations vary across trajectory methods','20 independent patients; CNV correlations use 19 estimable patients | all 48 original endpoints shown; BH across 12 targets per method\nAnchor was used to orient/root models, so anchor association is a diagnostic, not independent validation; CytoTRACE2 measures potency.')
fig,axs=plt.subplots(1,3,figsize=(13,5))
for ax,method in zip(axs,methods):
 sub=states[states.method==method].groupby(['rank','state']).fraction.mean().unstack().loc[:,list(statecolors)]
 assert np.allclose(sub.sum(axis=1),1)
 ax.stackplot(sub.index,sub.values.T,colors=list(statecolors.values()),labels=[x.replace('_',' ') for x in statecolors],alpha=.95)
 ax.set_title(method,loc='left');ax.set_ylim(0,1);ax.set_xlim(0,1);ax.set_xlabel('Within-patient pseudotime rank')
axs[0].set_ylabel('Patient-equal state fraction');axs[-1].legend(frameon=False,loc='center left',bbox_to_anchor=(1.02,.5),fontsize=8)
save(fig,'2_7','10_state_composition_along_pseudotime_v1','Malignant state composition follows method-specific orderings','20 patients contribute equally | original Base / Hypoxia / Cycle / TYMP-state labels | fixed within-patient display smoothing\nFractions sum to 1 at each rank; expression-state ordering does not establish a low-to-high clinical malignancy timeline.')
(q/'FIGURE_REGISTRY.json').write_text(json.dumps(written,indent=2))
(q/'SOFTWARE_VERSIONS.json').write_text(json.dumps({'python':sys.version,'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},indent=2))
assert len(written)==19
print('DRAWING_COMPLETE',len(written),flush=True)
