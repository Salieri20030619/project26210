from pathlib import Path
import json,ast
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize,PowerNorm
from matplotlib.cm import ScalarMappable
from matplotlib.patches import FancyBboxPatch,Rectangle
from matplotlib.lines import Line2D
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase198_fig2_contrast_corr_revision_20260918_v5'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
purple='#993F80';teal='#398B92';cm=LinearSegmentedColormap.from_list('effects',['#568AA7','#FAF8F0','#A5417B'])
def panel(fig,n):
 fig.savefig(d/'panels'/f'{n}.pdf',bbox_inches='tight',pad_inches=.06,dpi=400)
 plt.close(fig);print(n,flush=True)
# Top20 preserves original descriptive dot semantics.
o=r/'phase131_sections2_5_2_6_reference_correction_remote_v1'
ds=pd.read_csv(o/'data/PATHWAY_DOTPLOT_SUMMARY.tsv',sep='\t')
sel=json.loads((d/'tables/TOP20_SELECTION.json').read_text());paths=sel['paths']
v=json.loads((o/'qc/SECTION26_VALIDATION.json').read_text());factors=['F'+str(i) for i in range(1,8)];npats=[v['paired_patient_counts_per_program'][f] for f in factors]
row_stats=ds[ds.pathway.isin(paths)].groupby('pathway').mean_AUC.agg(['mean','std'])
ds['row_z']=(ds.mean_AUC-ds.pathway.map(row_stats['mean']))/ds.pathway.map(row_stats['std'])
ds[ds.pathway.isin(paths)].to_csv(d/'tables/TOP20_DISPLAY_VALUES.tsv',sep='\t',index=False)
zlim=float(abs(ds.loc[ds.pathway.isin(paths),'row_z']).max())
fig=plt.figure(figsize=(8.7,11.0))
vmax=ds[ds.pathway.isin(paths)].mean_AUC.max()
for side,group in enumerate(['TYMP_negative','TYMP_positive']):
 ax=fig.add_axes([.43+side*.215,.115,.20,.82])
 sub=ds[ds.group==group].set_index(['pathway','factor'])
 for i,path in enumerate(paths):
  for j,f in enumerate(factors):
   if (path,f) not in sub.index:ax.plot(j,i,marker='x',ms=3,color='#BBBBBB');continue
   z=sub.loc[(path,f)]
   sc=ax.scatter(j,i,s=z.fraction*50,c=[z.row_z],cmap='RdBu_r',vmin=-zlim,vmax=zlim,edgecolors='#444444',lw=.3)
 ax.set(xlim=(-.6,6.6),ylim=(19.6,-.6))
 ax.set_xticks(range(7),[f'{f}\nn={n}' for f,n in zip(factors,npats)],fontsize=6.5)
 ax.set_yticks(range(20),paths if side==0 else ['']*20,fontsize=7)
 ax.grid(color='#E2E2E2',lw=.4);ax.set_axisbelow(True);ax.tick_params(length=0)
 ax.set_title(['TYMP-','TYMP+'][side],fontsize=9)
ca=fig.add_axes([.91,.52,.015,.30]);cb=fig.colorbar(sc,cax=ca);cb.set_ticks([-2,-1,0,1,2]);cb.ax.tick_params(labelsize=6);cb.set_label('Within-pathway z-score of mean AUC',fontsize=7)
la=fig.add_axes([.895,.12,.10,.29]);la.axis('off')
la.legend(handles=[Line2D([],[],ls='',marker='o',ms=np.sqrt(50*x),mfc='#A2B7C7',mec='#444444',mew=.4,label=f'{x*100:.0f}%') for x in [.1,.3,.6,1]],title='Cells with\nAUC > 0',frameon=False,fontsize=6,title_fontsize=7,loc='upper left')
fig.text(.43,.035,'cNMF program | paired patients',fontsize=7)
fig.text(.43,.005,'Colors compare groups within each pathway | Original AUC values retained',fontsize=7)
panel(fig,'top20_pathways')
