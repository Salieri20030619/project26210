from pathlib import Path
import json,ast
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import FancyBboxPatch,Rectangle
from matplotlib.lines import Line2D
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase197_top20_network_dynamics_20260917_v4'
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
fig=plt.figure(figsize=(8.7,8.8))
vmax=ds[ds.pathway.isin(paths)].mean_AUC.max()
for side,group in enumerate(['TYMP_negative','TYMP_positive']):
 ax=fig.add_axes([.43+side*.215,.115,.20,.82])
 sub=ds[ds.group==group].set_index(['pathway','factor'])
 for i,path in enumerate(paths):
  for j,f in enumerate(factors):
   if (path,f) not in sub.index:ax.plot(j,i,marker='x',ms=3,color='#BBBBBB');continue
   z=sub.loc[(path,f)]
   sc=ax.scatter(j,i,s=z.fraction*50,c=[z.mean_AUC],cmap='viridis',vmin=0,vmax=vmax,edgecolors='#444444',lw=.3)
 ax.set(xlim=(-.6,6.6),ylim=(19.6,-.6))
 ax.set_xticks(range(7),[f'{f}\nn={n}' for f,n in zip(factors,npats)],fontsize=6.5)
 ax.set_yticks(range(20),paths if side==0 else ['']*20,fontsize=7)
 ax.grid(color='#E2E2E2',lw=.4);ax.set_axisbelow(True);ax.tick_params(length=0)
 ax.set_title(['TYMP-','TYMP+'][side],fontsize=9)
ca=fig.add_axes([.91,.52,.015,.30]);cb=fig.colorbar(sc,cax=ca);cb.ax.tick_params(labelsize=6);cb.set_label('Patient-equal mean AUC',fontsize=7)
la=fig.add_axes([.895,.12,.10,.29]);la.axis('off')
la.legend(handles=[Line2D([],[],ls='',marker='o',ms=np.sqrt(50*x),mfc='#A2B7C7',mec='#444444',mew=.4,label=f'{x*100:.0f}%') for x in [.1,.3,.6,1]],title='Cells with\nAUC > 0',frameon=False,fontsize=6,title_fontsize=7,loc='upper left')
fig.text(.43,.035,'cNMF program | paired patients',fontsize=7)
panel(fig,'top20_pathways')
# Membership network: original 24 genes and four KEGG gene sets.
z=np.load(r/'part2_rebuild_v1/04_markers_annotation/run_v4_macrophage_only/paired_patient_effects/cluster_4.npz',allow_pickle=True)
allgenes=['SPP1','GPNMB','CCL18','APOE','LPL','CYP27A1','ABCA1','ABCG1','NPC1','SOAT1','CTSB','CTSD','CTSL','LAMP1','LAMP2','LIPA','HLA-DRA','HLA-DPA1','HLA-DPB1','CD74','C1QA','C1QB','C1QC','TYROBP']
median={g:float(np.median(z['effects'][:,list(z['genes']).index(g)])) for g in allgenes}
pathnames=['Cholesterol metabolism','Lysosome','Phagosome','Antigen processing and presentation']
gmt=pd.read_csv(r/'phase177_L2_myeloid_state_integration_20260916_v1/source_tables/program_genes.tsv',sep='\t')
sets={p:set(gmt.loc[gmt.program==p,'gene_symbols'].iloc[0].split(';')) for p in pathnames}
display=allgenes[3:20];edges=[(p,g) for p in pathnames for g in display if g in sets[p]]
fig,ax=plt.subplots(figsize=(7.1,6.0));fig.subplots_adjust(left=.02,right=.95,top=.99,bottom=.20)
ax.set(xlim=(-.85,2.95),ylim=(-.045,1.045));ax.axis('off')
py=dict(zip(pathnames,[.87,.62,.37,.12]));gy=dict(zip(display,np.linspace(.98,.02,len(display))))
for p,g in edges:
 ax.plot([.4,1.95],[py[p],gy[g]],color='#A2ADB2',alpha=.60,lw=.75,zorder=1)
for p in pathnames:
 label=p.replace('Cholesterol metabolism','Cholesterol\nmetabolism').replace('Antigen processing and presentation','Antigen\npresentation')
 ax.add_patch(FancyBboxPatch((-.73,py[p]-.057),1.12,.114,boxstyle='round,pad=.012',fc='#E7EEEE',ec=teal,lw=.8,zorder=2))
 ax.text(-.17,py[p],label,ha='center',va='center',fontsize=9,color='#294D52')
norm=Normalize(-3,3)
for g in display:
 ax.plot(1.96,gy[g],marker='o',ms=10,mfc=cm(norm(median[g])),mec='#707070',mew=.5,zorder=3)
 ax.text(2.08,gy[g],g,va='center',fontsize=9,fontstyle='italic')
ca=fig.add_axes([.89,.30,.020,.50]);cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cm),cax=ca,extend='both',ticks=[-3,0,3]);cb.set_label('Median paired delta log2(CPM + 1)',fontsize=8)
# Seven state/backbone genes kept without invented pathway membership edges.
for k,g in enumerate(allgenes[:3]+allgenes[20:]):
 xx=.09+k*.119
 fig.add_artist(Line2D([xx],[.092],marker='o',ms=10,mfc=cm(norm(median[g])),mec='#777777',mew=.5,ls='',transform=fig.transFigure))
 fig.text(xx,.052,g,fontsize=8,ha='center',fontstyle='italic')
fig.text(.03,.147,'State identity and macrophage backbone',fontsize=8,weight='bold')
fig.text(.03,.012,'Edges: KEGG gene-set membership | Colors: C4 versus other macrophages | 26 patients',fontsize=7)
panel(fig,'gene_function_network')
pd.DataFrame(edges,columns=['pathway','gene']).to_csv(d/'tables/GENE_PATHWAY_EDGES.tsv',sep='\t',index=False)
pd.DataFrame({'gene':allgenes,'median_delta':[median[g] for g in allgenes]}).to_csv(d/'tables/NETWORK_GENE_VALUES.tsv',sep='\t',index=False)
# Source-share annotation: two compact side-by-side composition strips.
co=pd.read_csv(r/'phase177_L2_myeloid_state_integration_20260916_v1/source_tables/c4_contribution.tsv',sep='\t')
fr=co[['C4_fraction','C4_SPP1_RNA_share']].median().values
fig,axs=plt.subplots(1,2,figsize=(12.5,1.65))
fig.subplots_adjust(left=.01,right=.99,bottom=.22,top=.94,wspace=.12)
for ax,value,label in zip(axs,fr,['Macrophage cell abundance','Captured SPP1 RNA contribution']):
 ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 ax.text(.02,.86,label,fontsize=10,weight='bold')
 ax.add_patch(Rectangle((.02,.33),.76,.27,fc='#E1E5E5',ec='none'))
 ax.add_patch(Rectangle((.02,.33),.76*value,.27,fc=purple,ec='none'))
 ax.text(.82,.465,f'{value:.1%}',ha='left',va='center',fontsize=18,weight='bold',color=purple)
fig.text(.5,.03,'Purple: C4 | Gray: other macrophages | Patient medians, n = 16',ha='center',fontsize=8,color='#555555')
panel(fig,'C4_source_share')
# Curves replace the entire dynamic heatmap; no new curve fitting.
dy=np.load(r/'phase192_sections28_212_figures_20260917_v1/inputs/DYNAMIC_DISPLAY.npz')
groups={'Monocyte genes':['FCN1','VCAN','S100A8','S100A9','LYZ'],'SPP1':['SPP1'],'Lipid-handling genes':['APOE','LPL','CYP27A1','ABCA1','ABCG1','NPC1','SOAT1'],'Lysosomal genes':['CTSB','CTSD','CTSL','LAMP1','LAMP2','LIPA','ATP6V0D1'],'Antigen-presentation genes':['HLA-DRA','HLA-DPA1','HLA-DPB1','CD74']}
fig,axs=plt.subplots(2,3,figsize=(6.8,5.2),layout='constrained');axs=axs.ravel();axs[-1].axis('off');records=[]
for ax,(label,genes) in zip(axs,groups.items()):
 idx=[list(dy['genes']).index(g) for g in genes]
 for mode,col,ls in [('primary',purple,'-'),('sensitivity',teal,'--')]:
  vals=dy[mode][idx].mean(0)
  ax.plot(dy['grid'],vals,color=col,lw=1.7,ls=ls)
  records.extend({'panel':label,'mode':mode,'rank':float(t),'mean_logCP10K':float(v)} for t,v in zip(dy['grid'],vals))
 ax.set_title(label.replace('-','-\n') if len(label)>19 else label,fontsize=9)
 ax.set_xticks([0,.5,1],['Early','Mid','Late']);ax.tick_params(labelsize=7);ax.grid(axis='y',color='#EEEEEE',lw=.5)
 ax.set_xlabel('Pseudotime rank',fontsize=8)
axs[0].set_ylabel('Mean log1p(CP10K)',fontsize=8)
axs[-1].legend(handles=[Line2D([],[],color=purple,lw=1.7,label='Primary: 14 patients'),Line2D([],[],color=teal,ls='--',lw=1.7,label='Sensitivity: 24 patients')],loc='center',ncol=1,frameon=False,fontsize=8)
panel(fig,'expression_dynamics_curves')
pd.DataFrame(records).to_csv(d/'tables/DYNAMIC_CURVES.tsv',sep='\t',index=False)
(d/'tables/DYNAMIC_GENE_GROUPS.json').write_text(json.dumps(groups,indent=2))

# Short, wide tissue panel: unchanged 255 samples and row values.
bulk=pd.read_csv(r/'phase192_sections28_212_figures_20260917_v1/tables/TISSUE_MTA_ORDERED_PATIENTS.tsv',sep='\t',index_col=0)
features=['MET:5-methylthioadenosine (MTA)','SPP1','LST1','AIF1','CSF1R','CD68','FCER1G','SMS','MTAP']
fig=plt.figure(figsize=(12.5,3.2));start=0
for i,(co,dat) in enumerate(bulk.groupby('cohort',sort=False)):
 ax=fig.add_axes([.07+.86*start/255,.25,.86*len(dat)/255-.006,.62])
 vv=dat[features].values.T.copy();vv[0]=np.log2(vv[0]+1)
 zz=(vv-vv.mean(1,keepdims=True))/vv.std(1,keepdims=True)
 im=ax.imshow(zz,aspect='auto',cmap=cm,vmin=-2.5,vmax=2.5,interpolation='nearest')
 ax.set_title(co.replace('_',' ')+f' | n = {len(dat)}',fontsize=9)
 ax.set_yticks(range(9),['MTA']+features[1:] if i==0 else ['']*9,fontsize=7)
 ax.set_xticks([]);ax.tick_params(length=0)
 ax.add_patch(Rectangle((-.5,.5),len(dat),1,fill=False,ec=purple,lw=.9))
 start+=len(dat)
fig.text(.5,.14,'Patients ordered by measured MTA within each cohort',ha='center',fontsize=8)
ca=fig.add_axes([.951,.30,.011,.48]);fig.colorbar(im,cax=ca,label='Within-cohort row z-score',ticks=[-2,0,2]);ca.tick_params(labelsize=7)
panel(fig,'tissue_compact')
