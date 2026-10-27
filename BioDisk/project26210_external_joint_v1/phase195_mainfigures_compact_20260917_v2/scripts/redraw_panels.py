from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase195_mainfigures_compact_20260917_v2'
for s in ['panels','figures','qc','tables']:(d/s).mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
cm=LinearSegmentedColormap.from_list('div',['#568AA7','#FAF8F0','#A5417B'])
def save(fig,n):
 fig.savefig(d/'panels'/f'{n}.pdf',bbox_inches='tight',pad_inches=.045,dpi=400);plt.close(fig);print(n,flush=True)
# Regenerate four protein contours from original KDE grids, with no new fit.
p=r/'phase131_sections2_5_2_6_reference_correction_remote_v1'
stats=pd.read_csv(r/'phase128_section2_5_metaviper_redraw_remote_v1/qc/m12_f1_activity_patient_summary_v1.tsv',sep='\t')
fig,axs=plt.subplots(2,2,figsize=(6.4,5.0),layout='constrained');axs=axs.ravel()
pc=LinearSegmentedColormap.from_list('protein',['#FFFFFF','#D4C8E8','#7D63B1','#2C155C'])
for ax,g in zip(axs,['CASP9','CHMP1A','PRSS8','PXN']):
 z=np.load(p/'data'/f'F1_{g}_DENSITY_GRID.npz')
 lev=[.025,.06,.12,.22,.36,.52,.70,.86,1.000001]
 im=ax.contourf(z['x'],z['y'],z['relative_density'],levels=lev,cmap=pc,vmin=0,vmax=1)
 ax.contour(z['x'],z['y'],z['relative_density'],levels=lev[:-1],colors='#695279',linewidths=.35)
 ax.set_xticks([0,.5,1]);ax.set_xlabel('F1 usage');ax.set_title(g,fontsize=10,weight='bold')
 ax.axvline(0,color='#999999',ls='--',lw=.4);ax.axvline(1,color='#999999',ls='--',lw=.4)
 ax.tick_params(labelsize=7)
axs[0].set_ylabel('metaVIPER NES')
fig.colorbar(im,ax=axs,shrink=.7,fraction=.016,ticks=[.025,1],label='Relative density')
save(fig,'protein_compact')
stats.to_csv(d/'tables/ORIGINAL_PROTEIN_STATISTICS.tsv',sep='\t',index=False)
# Functional heatmap from frozen patient-balanced values.
v=pd.read_csv(r/'phase192_sections28_212_figures_20260917_v1/tables/PROGRAM_PATIENT_EQUAL_SCORES.tsv',sep='\t',index_col=0)
z=(v.values-v.values.mean(1,keepdims=True))/v.values.std(1,keepdims=True)
fig,ax=plt.subplots(figsize=(5.4,4.3),layout='constrained');im=ax.imshow(z,aspect='auto',cmap=cm,vmin=-2.5,vmax=2.5)
labels=[x.replace('Antigen processing and presentation','Antigen presentation').replace('Interferon Gamma Response','IFN-gamma response') for x in v.index]
ax.set_yticks(range(10),labels,fontsize=7);ax.set_xticks(range(12),v.columns,fontsize=7);ax.tick_params(length=0)
ax.add_patch(Rectangle((1.5,-.5),1,10,fill=False,ec='#993F80',lw=1.3))
fig.colorbar(im,ax=ax,fraction=.025,pad=.02,label='Row z-score',ticks=[-2,0,2]);save(fig,'program_compact')
# New biological gene-level panel using already frozen 26-patient contrasts.
z=np.load(r/'part2_rebuild_v1/04_markers_annotation/run_v4_macrophage_only/paired_patient_effects/cluster_4.npz',allow_pickle=True)
groups={'Identity':['SPP1','GPNMB','CCL18'],'Lipid':['APOE','LPL','CYP27A1','ABCA1','ABCG1','NPC1','SOAT1'],'Lysosome':['CTSB','CTSD','CTSL','LAMP1','LAMP2','LIPA'],'Antigen':['HLA-DRA','HLA-DPA1','HLA-DPB1','CD74'],'Backbone':['C1QA','C1QB','C1QC','TYROBP']}
genes=sum(groups.values(),[]);inds=[list(z['genes']).index(g) for g in genes]
val=z['effects'][:,inds].T;assert val.shape==(24,26) and np.isfinite(val).all()
fig,ax=plt.subplots(figsize=(7,3.8),layout='constrained');im=ax.imshow(val,aspect='auto',cmap=cm,vmin=-3,vmax=3)
ax.set_yticks(range(24),genes,fontsize=6.8,fontstyle='italic');ax.set_xticks(range(26),[f'P{i+1:02d}' for i in range(26)],rotation=90,fontsize=6);ax.tick_params(length=0)
start=0
for group,gs in groups.items():
 ax.axhline(start-.5,color='white',lw=1);start+=len(gs)
studies=z['study']
for i in range(1,26):
 if studies[i]!=studies[i-1]:ax.axvline(i-.5,color='white',lw=1.6)
fig.colorbar(im,ax=ax,fraction=.025,pad=.02,label='C4 - other macrophages\nDelta log2(CPM + 1)',extend='both',ticks=[-3,0,3])
save(fig,'functional_gene_pairs')
pd.DataFrame(val,index=genes,columns=z['patient']).to_csv(d/'tables/FUNCTIONAL_GENE_PAIRS.tsv',sep='\t')
pd.DataFrame({'display':[f'P{i+1:02d}' for i in range(26)],'patient':z['patient'],'study':studies}).to_csv(d/'tables/GENE_PATIENT_KEY.tsv',sep='\t',index=False)
# Contribution: all 16 pre-existing same-library patients, no new eligibility filter.
co=pd.read_csv(r/'phase177_L2_myeloid_state_integration_20260916_v1/source_tables/c4_contribution.tsv',sep='\t')
assert len(co)==16
fig,ax=plt.subplots(figsize=(6.6,1.65),layout='constrained')
vv=co[['C4_fraction','C4_SPP1_RNA_share']].values.T
im=ax.imshow(vv,aspect='auto',cmap='YlGnBu',vmin=0,vmax=1)
ax.set_yticks([0,1],['C4 cell fraction','C4 share of SPP1 RNA'],fontsize=8)
ax.set_xticks(range(16),[f'S{i+1:02d}' for i in range(16)],rotation=90,fontsize=6);ax.tick_params(length=0)
for i in range(2):
 for j in range(16):ax.text(j,i,f'{vv[i,j]*100:.1f}',ha='center',va='center',fontsize=6,color='white' if vv[i,j]>.5 else '#252525')
cb=fig.colorbar(im,ax=ax,fraction=.023,pad=.02,ticks=[0,.5,1],label='Share (%)');cb.ax.set_yticklabels(['0','50','100'])
save(fig,'spp1_contribution')
co.insert(0,'display',[f'S{i+1:02d}' for i in range(16)]);co.to_csv(d/'tables/SPP1_CONTRIBUTION.tsv',sep='\t',index=False)
# Regenerate dynamics with tighter geometry; same 35 genes, same matrices, same joint z.
dy=np.load(r/'phase192_sections28_212_figures_20260917_v1/inputs/DYNAMIC_DISPLAY.npz')
join=np.c_[dy['primary'],dy['sensitivity']];zz=(join-join.mean(1,keepdims=True))/np.maximum(join.std(1,keepdims=True),1e-8)
fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
for i,ax in enumerate(axs):
 im=ax.imshow(zz[:,i*40:(i+1)*40],aspect='auto',cmap=cm,vmin=-2.5,vmax=2.5)
 ax.set_yticks(range(35),dy['genes'] if i==0 else ['']*35,fontsize=6.8,fontstyle='italic')
 ax.set_xticks([0,19.5,39],['Early','Intermediate','Late']);ax.tick_params(length=0)
 ax.set_title(['Primary | 14 patients','Sensitivity | 24 patients'][i],fontsize=10)
 ax.set_xlabel('Within-partition pseudotime rank')
 for b in [4.5,7.5,14.5,21.5,25.5,29.5]:ax.axhline(b,color='white',lw=.8)
fig.colorbar(im,ax=axs,fraction=.015,pad=.015,label='Joint row z-score',ticks=[-2,0,2])
save(fig,'dynamic_compact')
(d/'tables/NEW_PANEL_METHODS.json').write_text(json.dumps({'gene_patients':26,'contribution_patients':16,'new_inference':False,'gene_color_limit':[-3,3],'dynamic_values_unchanged':True,'density_KDE_refit':False},indent=2))
