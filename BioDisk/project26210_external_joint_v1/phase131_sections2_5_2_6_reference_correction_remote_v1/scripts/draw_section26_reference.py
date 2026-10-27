from pathlib import Path
import json,subprocess,sys
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
w=Path(__file__).resolve().parents[1];q=w/'qc';out=w/'figures_2_6';root=w.parent
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory')
def read(p,index=None):return pd.read_csv(p,sep='\t',index_col=index)
a=read(b/'metabolism_v1/source_data/four_pathway_AUCell_cell_scores_v1.tsv.gz',0)
u=read(b/'trajectory_v1/input/scvi_umap_v1.tsv.gz',0).loc[a.index]
f=read(b/'metabolism_v1/source_data/scFEA_four_axis_cell_scores_v1.tsv.gz',0).loc[a.index]
ds=read(w/'data/PATHWAY_DOTPLOT_SUMMARY.tsv');valid=json.loads((q/'SECTION26_VALIDATION.json').read_text());assert valid['status']=='PASS'
factors=['F'+str(i) for i in range(1,8)];npats=[valid['paired_patient_counts_per_program'][s] for s in factors]
pathnames=['Pyrimidine metabolism','Cysteine and methionine metabolism','Arginine and proline metabolism','Glycerophospholipid metabolism']
auc=['AUC_Pyrimidine','AUC_CysMet','AUC_ArgPro','AUC_Glycerophospholipid'];flux=['Pyrimidine_flux','CysMet_flux','ArgPro_flux','Choline_proxy']
short=['Pyrimidine','Cys / Met','Arg / Pro','Glycerophospholipid'];fshort=['Pyrimidine proxy','Cys / Met proxy','Arg / Pro proxy','Choline proxy']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.7,'pdf.fonttype':42,'savefig.facecolor':'white'})
registry=[]
def save(fig,name,title,foot,reference):
 fig.text(.035,.943 if fig.get_figheight()<8 else .963,title,ha='left',fontsize=16,fontweight='bold')
 fig.text(.035,.025,foot,ha='left',va='bottom',fontsize=8,linespacing=1.45)
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=360)
 plt.close(fig)
 subprocess.run(['pdftoppm','-r','115','-singlefile','-png',str(out/(name+'.pdf')),str(q/(name+'_pdf'))],check=True)
 registry.append({'name':name,'reference':reference,'title':title});print('SAVED',name,flush=True)
def umap_row(frame,cols,names,name,title,isflux=False):
 fig,axs=plt.subplots(1,4,figsize=(17,4.9));fig.subplots_adjust(left=.05,right=.985,top=.79,bottom=.20,wspace=.46)
 for ax,col,label in zip(axs,cols,names):
  v=frame[col].values;order=np.argsort(v,kind='stable')
  ax.set_axisbelow(True);ax.grid(color='#E5E5E5',lw=.45)
  im=ax.scatter(u.iloc[order,0],u.iloc[order,1],c=v[order],s=2.8,cmap='magma',vmin=v.min(),vmax=v.max(),linewidths=0,rasterized=True,zorder=3)
  ax.set_title(label,loc='left',fontsize=11);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.set_aspect('equal');ax.set_xlim(u.UMAP1.min()-.35,u.UMAP1.max()+.35);ax.set_ylim(u.UMAP2.min()-.35,u.UMAP2.max()+.35)
  ax.xaxis.set_major_locator(MaxNLocator(4));ax.yaxis.set_major_locator(MaxNLocator(4));ax.tick_params(labelsize=8,length=2)
  cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.025,shrink=.35);cb.ax.set_title('Proxy' if isflux else 'AUC',fontsize=8,pad=5);cb.ax.tick_params(labelsize=7,length=2);cb.locator=MaxNLocator(3);cb.update_ticks()
 save(fig,name,title,'11,950 malignant cells | 26 patients | original UMAP and scores | full per-pathway color ranges\n'+('scFEA: mean module-z proxies; choline is partial and MTA is not an explicit model output.' if isflux else 'Color indicates original AUCell pathway score. Higher values are shown in yellow.'),'image16: four horizontal feature panels, visible axes/grid, magma palette, small side color bars')
umap_row(a,auc,short,'01_four_pathway_UMAP_reference_v2','Four metabolic pathways in malignant epithelial cells')
def dot_facets(paths,name,title):
 full=len(paths)>10;fig=plt.figure(figsize=(16,22) if full else (13,5.7))
 left=.405 if full else .31;pw=.22 if full else .245;gap=.012;bottom=.075 if full else .19;top=.93 if full else .83;height=top-bottom
 vmax=float(ds[ds.pathway.isin(paths)].mean_AUC.max());maxsize=118 if full else 290;sc=None
 for side,group in enumerate(['TYMP_negative','TYMP_positive']):
  ax=fig.add_axes([left+side*(pw+gap),bottom,pw,height]);sub=ds[ds.group==group].set_index(['pathway','factor'])
  for i,path in enumerate(paths):
   for j,fac in enumerate(factors):
    if (path,fac) not in sub.index:
     ax.plot(j,i,marker='x',ms=3,color='#C2C7CC',mew=.7);continue
    r=sub.loc[(path,fac)];sc=ax.scatter(j,i,s=r.fraction*maxsize,c=[r.mean_AUC],cmap='viridis',vmin=0,vmax=vmax,edgecolors='#34363B',linewidths=.35,zorder=3)
  ax.set_xlim(-.55,6.55);ax.set_ylim(len(paths)-.5,-.5)
  ax.set_xticks(range(7),[f'{fac}\nn={nn}' for fac,nn in zip(factors,npats)],fontsize=8 if full else 9)
  ax.set_yticks(range(len(paths)),paths if side==0 else ['']*len(paths),fontsize=7 if full else 9)
  ax.tick_params(length=2,pad=5);ax.set_xlabel('cNMF program | paired patients',fontsize=9)
  ax.set_axisbelow(True);ax.grid(color='#D7D9DC',lw=.55)
  ax.add_patch(Rectangle((0,1.002),1,.02 if full else .085,transform=ax.transAxes,clip_on=False,facecolor='#E5E5E5',edgecolor='#333333',lw=.6))
  ax.text(.5,1.012 if full else 1.044,'TYMP-' if side==0 else 'TYMP+',transform=ax.transAxes,ha='center',va='center',fontsize=11)
 lx=left+2*pw+gap+.035
 cax=fig.add_axes([lx,.52 if full else .39,.013,.15 if full else .24]);cb=fig.colorbar(sc,cax=cax);cb.set_label('Patient-equal mean AUC',fontsize=9);cb.ax.tick_params(labelsize=8)
 lax=fig.add_axes([lx-.005,.29 if full else .19,.115,.17 if full else .17]);lax.axis('off')
 handles=[plt.scatter([],[],s=v*maxsize,facecolors='#A2B7C7',edgecolors='#34363B',linewidths=.4,label=f'{int(v*100)}%') for v in [.1,.3,.6,1]]
 lax.legend(handles=handles,title='Cells with\nAUC > 0',frameon=False,loc='upper left',fontsize=8,title_fontsize=9,labelspacing=1.3)
 save(fig,name,title,'Same patients in the TYMP- and TYMP+ facets within each program; >=10 cells per arm per patient. Each patient contributes equally.\nDot area = fraction of cells with AUC >0; color = mean AUC; x = no eligible pair (F3). Descriptive display; dot size is not a P value.','image15: two group facets, pathways as rows, cNMF groups as columns, filled circles, grid and separate legends')
dot_facets(pathnames,'02_four_pathway_TYMP_dotplot_reference_v2','Metabolic pathway profiles in TYMP- and TYMP+ cells')
allpaths=list(read(b/'metabolism_v1/source_data/KEGG_AUCell_all_pathways_v1.tsv.gz',0).columns)
dot_facets(allpaths,'03_all85_pathways_TYMP_dotplot_reference_v2','All 85 metabolic pathways in TYMP- and TYMP+ cells')
# A true circle corrplot, with patient-level effects, not rectangular heatmap cells.
cs=read(w/'data/CORRPLOT_SUMMARY.tsv');cols=['F1']+auc;names=['F1','Pyrimidine','Cys / Met','Arg / Pro','Glycerophospholipid']
fig=plt.figure(figsize=(12,7.2));ax=fig.add_axes([.20,.24,.44,.61]);mat=cs.pivot(index='x',columns='y',values='rho').loc[cols,cols].values
for i in range(5):
 for j in range(5):
  r=mat[i,j];sc=ax.scatter(j,i,s=abs(r)*1250,c=[r],cmap='RdBu_r',vmin=-1,vmax=1,edgecolors='none',zorder=3)
  ax.text(j,i,f'{r:.2f}',ha='center',va='center',fontsize=9,color='white' if abs(r)>.65 else '#30353A',zorder=4)
ax.set_xlim(-.5,4.5);ax.set_ylim(4.5,-.5);ax.set_aspect('equal');ax.set_xticks(range(5),names,rotation=38,ha='right',fontsize=9);ax.set_yticks(range(5),names,fontsize=9);ax.tick_params(length=0)
ax.set_xticks(np.arange(-.5,5,.1)[::10],minor=True);ax.set_yticks(np.arange(-.5,5,.1)[::10],minor=True);ax.grid(which='minor',color='#E1E4E7',lw=.6);ax.tick_params(which='minor',length=0)
cb=fig.colorbar(sc,ax=ax,fraction=.04,pad=.025,shrink=.8);cb.set_label('Median patient Spearman rho',fontsize=9)
summ=read(b/'metabolism_v1/tables/four_pathway_AUCell_correlation_summary_v1.tsv');txt=fig.add_axes([.73,.23,.25,.62]);txt.axis('off');txt.text(0,1,'F1 associations',fontsize=11,fontweight='bold',va='top')
for i,(col,label) in enumerate(zip(auc,short)):
 r=summ[(summ.pathway==col)&(summ.predictor=='F1')].iloc[0]
 txt.text(0,.87-i*.2,f'{label}\nrho={r.median_patient_rho:+.3f}; positive {int(round(r.fraction_positive*r.n_patients))}/{int(r.n_patients)}\nP={r.signed_rank_PValue:.3g}; FDR={r.FDR:.3g}',fontsize=9,va='top',linespacing=1.6)
save(fig,'04_pathway_F1_circle_corrplot_v2','Correlations among metabolic pathways and the F1 program','Circle area = absolute rho; color = direction and magnitude | 20 patients; correlations calculated within each patient\nF1 P/FDR retain the original four-pathway family. Pathway-pathway entries are descriptive summaries.','section2.6 text: Spearman Corrplot; circle encoding of correlation, complete side statistics')
# Transposed gene DotPlot follows the example's row-gene/column-group orientation.
dot=read(root/'phase129_sections2_6_2_7_redraw_remote_v1/qc/ENZYME_DOTPLOT_PATIENT_PROFILES.tsv');genes=['TYMP','TK1','SLC29A1','SLC29A2','ODC1','SRM','SMS','MTAP','CHKA','PCYT1A','PCYT2','CEPT1']
gs=dot.groupby(['gene','quartile']).agg(mu=('mean_log1p','mean'),fraction=('fraction_detected','mean')).reset_index();vmax=gs.mu.max()
fig=plt.figure(figsize=(8.5,9));ax=fig.add_axes([.24,.17,.39,.69])
for i,gene in enumerate(genes):
 for j in range(4):
  r=gs[(gs.gene==gene)&(gs.quartile==j+1)].iloc[0]
  sc=ax.scatter(j,i,s=r.fraction*460,c=[r.mu],cmap='viridis',vmin=0,vmax=vmax,edgecolors='#32343A',linewidths=.5,zorder=3)
ax.set_xlim(-.5,3.5);ax.set_ylim(11.5,-.5);ax.set_xticks(range(4),['Q1\nlow F1','Q2','Q3','Q4\nhigh F1']);ax.set_yticks(range(12),genes,fontstyle='italic');ax.tick_params(length=2)
ax.set_axisbelow(True);ax.grid(color='#D9DCE0',lw=.55)
for y in [3.5,7.5]:ax.axhline(y,color='#66727D',lw=1)
ax.add_patch(Rectangle((0,1.004),1,.06,transform=ax.transAxes,clip_on=False,facecolor='#E5E5E5',edgecolor='#333333',lw=.6));ax.text(.5,1.034,'Within-patient F1 quartile',transform=ax.transAxes,ha='center',va='center',fontsize=10)
cax=fig.add_axes([.74,.50,.025,.23]);cb=fig.colorbar(sc,cax=cax);cb.set_label('Mean log1p CP10K',fontsize=9)
leg=fig.add_axes([.73,.2,.2,.23]);leg.axis('off');handles=[plt.scatter([],[],s=v*460,c='#A0B5C5',edgecolors='#333333',linewidths=.4,label=f'{int(v*100)}%') for v in [.1,.3,.6,1]]
leg.legend(handles=handles,title='Detected cells',frameon=False,loc='upper left',fontsize=9,labelspacing=1.4)
save(fig,'05_enzyme_F1_dotplot_reference_v2','Enzyme and transporter expression along F1','17 patients | >=5 cells in each F1 quartile; ties retained | patient-equal summaries\nDot area = RNA detection fraction; color = unscaled mean log1p expression.','image15 DotPlot grammar and section2.6 fixed 12-gene requirement')
umap_row(f,flux,fshort,'06_scFEA_proxy_UMAP_reference_v2','scFEA metabolic proxies in malignant epithelial cells',True)
# Reference-style statistical bubbles with valid biological-sample P values.
s=read(b/'metabolism_v1/tables/four_pathway_AUCell_paired_statistics_v1.tsv');fig=plt.figure(figsize=(10,6.8));maxabs=float(abs(s.median_patient_delta).max());norm=TwoSlopeNorm(vmin=-maxabs,vcenter=0,vmax=maxabs)
for j,contrast in enumerate(['TYMP-positive_vs_negative','F1-high_clone_vs_other']):
 ax=fig.add_axes([.32+j*.19,.22,.18,.58]);sub=s[s.contrast==contrast].set_index('pathway')
 for i,path in enumerate(auc):
  r=sub.loc[path];pval=-np.log10(r.paired_wilcox_PValue);sc=ax.scatter(0,i,s=pval*180,c=[r.median_patient_delta],cmap='RdBu_r',norm=norm,edgecolors='black' if r.FDR<.05 else '#AEB6BE',linewidths=1.4 if r.FDR<.05 else .55,zorder=3)
 ax.set_xlim(-.5,.5);ax.set_ylim(3.5,-.5);ax.set_xticks([0],['High minus low'],rotation=30,ha='right',fontsize=8);ax.set_yticks(range(4),pathnames if j==0 else ['']*4,fontsize=9);ax.grid(color='#D9DCE0',lw=.6);ax.set_axisbelow(True)
 ax.add_patch(Rectangle((0,1.004),1,.13,transform=ax.transAxes,clip_on=False,facecolor='#E5E5E5',edgecolor='#333333',lw=.6));ax.text(.5,1.069,'TYMP+ vs TYMP-' if j==0 else 'F1-high vs other\nprimary clones',transform=ax.transAxes,ha='center',va='center',fontsize=9)
cax=fig.add_axes([.79,.48,.02,.23]);cb=fig.colorbar(sc,cax=cax);cb.set_label('Median paired AUC delta',fontsize=9)
leg=fig.add_axes([.77,.17,.21,.23]);leg.axis('off');handles=[plt.scatter([],[],s=v*180,c='#B3C3CE',edgecolors='black',linewidths=.5,label=str(v)) for v in [.5,1,2,3]]
leg.legend(handles=handles,title='-log10(P)',frameon=False,loc='upper left',fontsize=8,labelspacing=1.2)
save(fig,'07_four_pathway_statistical_dotplot_v2','Patient-paired differences in four metabolic pathways','14 paired patients per contrast | >=20 cells per arm | size = patient-level -log10(P); color = median paired effect\nBlack outline = FDR <0.05; BH across four pathways separately within each contrast.','image15: two facets and significance-sized bubbles, using valid paired-patient tests')
(q/'SECTION26_FIGURE_REGISTRY.json').write_text(json.dumps(registry,indent=2))
print('SECTION26_REFERENCE_COMPLETE',len(registry),flush=True)
