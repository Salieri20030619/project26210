from pathlib import Path
import json,hashlib,sys,subprocess
import numpy as np,pandas as pd,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap,Normalize
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch,Rectangle
from matplotlib.ticker import MaxNLocator
from scipy.cluster.hierarchy import dendrogram

w=Path(__file__).resolve().parents[1];r=w.parent
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
sources=[]
def read(p,index=False):
 p=Path(p);sources.append(p);return pd.read_csv(p,sep='\t',index_col=0 if index else None)
def npz(p):
 p=Path(p);sources.append(p);return np.load(p)
gate=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(w/'qc/MEMORY_GATE_DRAW.txt').write_text(gate.stdout+gate.stderr);gate.check_returncode()
m=read(b/'05_cnmf/input_v1/cnmf_cell_metadata_v1.tsv.gz',True)
u=read(b/'05_cnmf/full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz',True)
z=npz(r/'phase123_section2_2_reference_redraw_remote_v1/data/plot_values_server_only.npz')
cells=z['cells'].astype(str);xy=z['xy'];expr=z['expression'];genes=z['genes'].tolist()
m=m.loc[cells];u=u.loc[cells]
assert len(cells)==11950 and len(set(cells))==11950
calls=read(b/'03_malignancy_cnv/malignancy_cnv_summary_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz')
hpath=b/'02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad';sources.append(hpath)
with h5py.File(hpath,'r')as h:
 idxkey=h['obs'].attrs.get('_index','_index')
 if isinstance(idxkey,bytes):idxkey=idxkey.decode()
 ids=h['obs'][idxkey].asstr()[:];fullxy=h['obsm/X_umap'][:]
index=pd.Index(ids).get_indexer(calls.exact_cell_id);assert (index>=0).all()
cxy=fullxy[index]
assert np.allclose(fullxy[pd.Index(ids).get_indexer(cells)],xy,atol=1e-7)
statuses=np.where(calls.malignant_conservative_v1,'Malignant',np.where(calls.malignancy_status_v1=='non_malignant_consensus','Nonmalignant','Unresolved'))
assert len(calls)==38670 and sum(statuses=='Malignant')==11950 and sum(statuses=='Nonmalignant')==498
pt=read(r/'phase123_section2_2_reference_redraw_remote_v1/qc/F1_TYMP_patients.tsv')
pc=read(b/'05_cnmf/tymp_rank_evidence_v1/factor_patient_stratified_correlations_v1.tsv.gz')
cnv=pc[(pc['rank']==7)&(pc.factor=='F1')&(pc.outcome=='infercnv_burden_over_threshold_v1')].dropna(subset=['rho'])
reg=read(r/'phase127_section2_4_pyscenic_redraw_remote_v1/qc/CANDIDATE_CONTINUOUS_DISPLAY.tsv')
act=read(r/'phase128_section2_5_metaviper_redraw_remote_v1/qc/m12_f1_activity_patient_summary_v1.tsv')
tya=read(r/'phase133_sections2_5_2_7_targeted_support_remote_v1/qc/TYMP_ACTIVITY_SUMMARY.tsv')
auc=read(r/'phase129_sections2_6_2_7_redraw_remote_v1/qc/original_four_pathway_AUCell_paired_statistics_v1.tsv')
tr=read(b/'09_metabolism_trajectory/trajectory_v1/source_data/trajectory_all_methods_cell_scores_v1.tsv.gz',True).loc[cells]
seg=read(b/'09_metabolism_trajectory/trajectory_v1/source_data/monocle3_principal_graph_segments_v1.tsv')
heat=npz(r/'phase132_section2_7_gradient_heatmaps_remote_v1/data/GRADIENT_HEATMAP_VALUES.npz')
states=read(r/'phase132_section2_7_gradient_heatmaps_remote_v1/data/PATIENT_STATE_COMPOSITION.tsv.gz')
term=read(r/'phase133_sections2_5_2_7_targeted_support_remote_v1/qc/TERMINAL_ENZYME_SUMMARY.tsv')
assert len(pt)==20 and len(cnv)==19 and len(term)==39 and len(reg)==6
W,H=190,270
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,'axes.titlesize':8,
 'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5,
 'xtick.major.size':2,'ytick.major.size':2,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white',
 'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
blue='#387CA0';red='#C45569';gray='#BDC4CA';dark='#24333D'
bluecm=LinearSegmentedColormap.from_list('F1',['#E9EDF2','#B5C6DD','#396B9C'])
redcm=LinearSegmentedColormap.from_list('TYMP',['#E7EAED','#E5B2B2','#A62C50'])
purple=LinearSegmentedColormap.from_list('density',['#FFFFFF','#EEE9F5','#D4C8E8','#A991CB','#684D9E'])
div=LinearSegmentedColormap.from_list('effects',['#327CA2','#F7F7F7','#C65564'])
gradient=LinearSegmentedColormap.from_list('dynamics',[(0,'#130078'),(.23,'#5400A8'),(.48,'#8900AC'),(.70,'#EF7949'),(1,'#F5F520')])
panel_records=[];text_records=[]
def txt(x,y,text,**kw):
 defaults=dict(fontsize=7,ha='left',va='top',color=dark)
 defaults.update(kw)
 artist=fig.text(x/W,1-y/H,text,**defaults);text_records.append(artist);return artist
def axmm(x,y,ww,hh):return fig.add_axes([x/W,1-(y+hh)/H,ww/W,hh/H])
def panel(letter,x,y,title):
 txt(x,y,letter,fontsize=11,fontweight='bold',color='black')
 txt(x+5,y+.3,title,fontsize=8.1,fontweight='bold',color='black')
 panel_records.append({'panel':letter,'title':title,'x_mm':x,'top_mm':y})
def clean(ax):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(pad=1.5)
def mapax(ax):
 ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
 for sp in ax.spines.values():sp.set_visible(False)
def bar(im,x,y,ww,label,ticks=None):
 ca=axmm(x,y,ww,1.7);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=ticks)
 cb.ax.tick_params(labelsize=6.5,length=1.8,pad=1)
 cb.set_label(label,fontsize=6.5,labelpad=1.1);return cb

txt(4,2.5,'TYMP-associated malignant epithelial states and metabolic programs',fontsize=10,fontweight='bold',color='black')
txt(4,8,'75-library cohort  |  Frozen analysis and original patient-level statistics',fontsize=7,color='#56636D')

# A: the entire tumor-epithelial CNV scope, no omission of unresolved cells.
panel('A',4,15,'CNV-defined epithelial cells')
ax=axmm(6,23,55,34)
for st,col,size in [('Unresolved','#D2D6D9',.25),('Nonmalignant',blue,.65),('Malignant',red,.5)]:
 ix=statuses==st;ax.scatter(cxy[ix,0],cxy[ix,1],s=size,c=col,linewidths=0,rasterized=True,label=st)
mapax(ax)
for i,(st,col,n)in enumerate([('Malignant',red,11950),('Nonmalignant',blue,498),('Unresolved','#B6BEC5',26222)]):
 fig.add_artist(Line2D([7/W],[1-(59+i*3.2)/H],marker='o',markersize=3,color=col,linestyle='none',transform=fig.transFigure))
 txt(10,57.9+i*3.2,f'{st}: {n:,}',fontsize=6.6)
txt(43,55.5,'UMAP1 / UMAP2',fontsize=6.5)

# B/C: identical original malignant-cell coordinates.
for letter,x,val,title,cm,label in [('B',67,u.F1.to_numpy(),'F1 program usage',bluecm,'F1 usage'),('C',130,expr[:,genes.index('TYMP')],'TYMP RNA',redcm,'Log1p RNA')]:
 panel(letter,x,15,title);ax=axmm(x+2,23,54,34)
 vmax=1.0 if letter=='B' else float(np.quantile(val[val>0],.99))
 order=np.argsort(val,kind='stable')
 im=ax.scatter(xy[order,0],xy[order,1],c=val[order],cmap=cm,s=.75,linewidths=0,vmin=0,vmax=vmax,rasterized=True)
 mapax(ax);bar(im,x+15,59,28,label,[0,vmax]);txt(x+2,68,'11,950 cells; 26 patients',fontsize=6.5)

# D: patient effects, not cell-level inferential tests.
panel('D',4,76,'F1 patient associations')
ax=axmm(17,87,42,27)
for i,(vals,col)in enumerate([(pt.effect.to_numpy(),red),(cnv.rho.to_numpy(),blue)]):
 jitter=np.linspace(-.13,.13,len(vals))
 ax.scatter(vals,i+jitter,s=9,color=col,alpha=.58,linewidths=0)
 med=np.median(vals);ax.scatter([med],[i],s=22,marker='D',color=dark,edgecolors='white',linewidths=.4,zorder=4)
ax.axvline(0,c='#A2ABB3',ls='--',lw=.55);ax.set_xlim(-1,1);ax.set_ylim(1.55,-.55)
ax.set_yticks([0,1],['TYMP','CNV']);ax.set_xticks([-1,-.5,0,.5,1]);ax.set_xlabel('Within-patient rho',labelpad=2);clean(ax)
txt(5,121,'TYMP: n=20, median=0.115, q=0.0197',fontsize=6.55)
txt(5,125,'CNV: n=19, median=-0.033, q=0.923',fontsize=6.55)

# E: all three candidate TFs, both methods.
panel('E',67,76,'F1-associated regulons')
ax=axmm(83,87,39,27)
for j,method in enumerate(['GENIE3','GRNBoost2']):
 for i,tf in enumerate(['ILF2','HOXD3','HNF4A']):
  v=reg[(reg.method==method)&(reg.tf==tf)].iloc[0]
  ax.scatter(v.median_patient_effect,i+(-.12 if j==0 else .12),s=24,c=[blue if j==0 else red],marker='o'if j==0 else'^',edgecolors='white',lw=.3)
ax.axvline(0,color='#A2ABB3',ls='--',lw=.55);ax.set_xlim(-.015,.36);ax.set_ylim(2.5,-.5)
ax.set_yticks(range(3),['ILF2','HOXD3','HNF4A']);ax.set_xticks([0,.1,.2,.3]);ax.set_xlabel('Median patient rho',labelpad=2);clean(ax)
fig.legend(handles=[Line2D([],[],marker='o',c=blue,ls='',label='GENIE3',markersize=3),Line2D([],[],marker='^',c=red,ls='',label='GRNBoost2',markersize=3)],loc='upper left',bbox_to_anchor=(68/W,1-120/H),borderaxespad=0,frameon=False,ncol=2,fontsize=6.5,handletextpad=.3,columnspacing=.8)
txt(68,125,'n=20; all six associations q<0.05',fontsize=6.55)

# F: exact saved KDE grids, four prespecified proteins.
panel('F',130,76,'F1-linked protein activity')
densstats=[]
for k,gene in enumerate(['CASP9','CHMP1A','PRSS8','PXN']):
 xx=135+(k%2)*27;yy=88+(k//2)*21
 dz=npz(r/'phase131_sections2_5_2_6_reference_correction_remote_v1/data'/('F1_'+gene+'_DENSITY_GRID.npz'))
 a=axmm(xx,yy,20,13)
 a.contourf(dz['x'],dz['y'],dz['relative_density'],levels=np.linspace(0,1,9),cmap=purple)
 a.contour(dz['x'],dz['y'],dz['relative_density'],levels=[.10,.25,.45,.65,.85],colors='#6E549F',linewidths=.35)
 a.axvline(0,c='#888888',ls='--',lw=.35);a.axvline(1,c='#888888',ls='--',lw=.35)
 a.set_xlim(-.2,1.2);a.set_xticks([0,1]);a.yaxis.set_major_locator(MaxNLocator(2));a.tick_params(labelsize=6.5,length=1.5,pad=1);clean(a)
 v=act[(act.protein==gene)&(act.endpoint=='continuous_F1_spearman')].iloc[0]
 a.set_title(gene,fontsize=7,pad=1.5,loc='left')
 a.text(.98,.96,f'rho {v.median_patient_effect:+.2f}',transform=a.transAxes,ha='right',va='top',fontsize=6.5,bbox=dict(facecolor='white',alpha=.80,edgecolor='none',pad=.3))
 densstats.append(v.to_dict())
txt(131,104,'NES',fontsize=6.5,rotation=90,ha='center',va='center')
txt(160,127,'F1 usage',fontsize=6.5,ha='center')

# G: all four activity associations and their mandatory technical sensitivity.
panel('G',4,138,'TYMP-protein activity')
ax=axmm(19,150,40,28)
for j,model in enumerate(['unadjusted','technical_adjusted']):
 for i,gene in enumerate(['CASP9','CHMP1A','PRSS8','PXN']):
  v=tya[(tya.target==gene)&(tya.model==model)].iloc[0];y=i+(-.13 if j==0 else .13);col=blue if j==0 else red
  ax.plot([v.median_CI_low,v.median_CI_high],[y,y],c=col,lw=.65)
  ax.scatter(v.median_effect,y,s=19,marker='o'if j==0 else's',facecolors=col if v.FDR<.05 else'white',edgecolors=col,lw=.65,zorder=4)
ax.axvline(0,color='#A2ABB3',ls='--',lw=.55);ax.set_xlim(-.10,.215);ax.set_ylim(3.55,-.55)
ax.set_xticks([-.1,0,.1,.2]);ax.set_yticks(range(4),['CASP9','CHMP1A','PRSS8','PXN']);ax.set_xlabel('Median patient rho',labelpad=2);clean(ax)
fig.legend(handles=[Line2D([],[],marker='o',c=blue,label='Primary',markersize=3,lw=.5),Line2D([],[],marker='s',c=red,label='Adjusted',markersize=3,lw=.5)],loc='upper left',bbox_to_anchor=(6/W,1-185/H),borderaxespad=0,frameon=False,ncol=2,fontsize=6.5,handlelength=1.2,handletextpad=.4,columnspacing=1)
txt(6,190,'n=20; filled symbols: q<0.05',fontsize=6.5)

# H: both prespecified comparisons and all four pathways retained.
panel('H',67,138,'Paired metabolic activity')
ax=axmm(88,150,33,27)
paths=['AUC_Pyrimidine','AUC_CysMet','AUC_ArgPro','AUC_Glycerophospholipid']
contrasts=['TYMP-positive_vs_negative','F1-high_clone_vs_other']
norm=Normalize(-.03,.03)
for i,path in enumerate(paths):
 for j,ct in enumerate(contrasts):
  v=auc[(auc.pathway==path)&(auc.contrast==ct)].iloc[0]
  size=10+24*(-np.log10(v.paired_wilcox_PValue))
  ax.scatter(j,i,s=size,c=[v.median_patient_delta],cmap=div,norm=norm,edgecolors=dark if v.FDR<.05 else'#C2C7CB',linewidths=.75 if v.FDR<.05 else.35)
ax.set_xlim(-.5,1.5);ax.set_ylim(3.5,-.5);ax.set_xticks([0,1],['TYMP+','F1 clone']);ax.set_yticks(range(4),['Pyrimidine','Cys/Met','Arg/Pro','GPL']);ax.tick_params(length=0,pad=2)
for sp in ax.spines.values():sp.set_visible(False)
ax.grid(axis='y',lw=.4,color='#E8EBEE');ax.set_axisbelow(True)
txt(70,182,'n=14; area = -log10(P)',fontsize=6.5)
bar(plt.cm.ScalarMappable(norm=norm,cmap=div),81,187,34,'Median delta AUCell',[-.03,0,.03])

# I: actual learned graph, same frozen UMAP as B/C.
panel('I',130,138,'Monocle3 trajectory')
ax=axmm(132,148,54,34)
vv=tr.Monocle3_pseudotime_patient_v1.to_numpy();order=np.argsort(vv,kind='stable')
im=ax.scatter(xy[order,0],xy[order,1],c=vv[order],s=.7,cmap='viridis',vmin=0,vmax=1,linewidths=0,rasterized=True)
ax.add_collection(LineCollection(np.stack([seg[['x','y']].values,seg[['xend','yend']].values],axis=1),colors='#333D47',linewidths=.36,alpha=.85,zorder=4))
mapax(ax);bar(im,146,186,27,'Patient-scaled pseudotime',[0,1])

# J: all 2511 rows, original Monocle3 gene ordering and clusters.
panel('J',4,202,'Gene dynamics along Monocle3')
grid=heat['grid'];leaves=heat['leaves'];n=len(heat['genes']);clusters=heat['clusters']
assert n==2511
composition=states[states.method=='Monocle3'].groupby(['rank','state']).fraction.mean().unstack().reindex(columns=['Base','Hypoxia','Cycle','TYMP_state']).sort_index()
assert composition.shape==(200,4) and np.allclose(composition.sum(axis=1),1)
statecols=['#8CA8B6','#CB8B4B','#8766AB','#C94C64'];clcols=['#9B7890','#47879A','#87A86F','#C3A987']
for j,(name,col)in enumerate(zip(['Base','Hypoxia','Cycle','TYMP'],statecols)):
 fig.add_artist(Rectangle(((15+j*23)/W,1-211/H),2/W,1.7/H,transform=fig.transFigure,color=col,ec='none'))
 txt(18+j*23,209.1,name,fontsize=6.5)
ca=axmm(14,213,60,4);ca.stackplot(grid,composition.values.T,colors=statecols,lw=0);ca.set_xlim(grid[0],grid[-1]);ca.set_ylim(0,1);ca.axis('off')
a=axmm(14,219,60,35)
ordered=heat['Monocle3'][leaves];fine=np.linspace(grid[0],grid[-1],1000);show=np.stack([np.interp(fine,grid,row)for row in ordered])
im=a.imshow(show,aspect='auto',interpolation='nearest',cmap=gradient,vmin=-2.5,vmax=2.5,extent=[grid[0],grid[-1],n-.5,-.5],rasterized=True)
a.set_yticks([]);a.set_xticks([.025,.5,.975],['0.025','0.50','0.975']);a.tick_params(pad=1,length=1.5);a.set_xlabel('Within-patient pseudotime rank',fontsize=6.5,labelpad=1.5)
for sp in a.spines.values():sp.set_visible(False)
for boundary in np.flatnonzero(np.diff(clusters)!=0)+.5:a.axhline(boundary,c='white',lw=.45)
strip=axmm(11,219,1.8,35);strip.imshow((clusters-1)[:,None],aspect='auto',cmap=ListedColormap(clcols),vmin=0,vmax=3,interpolation='nearest');strip.axis('off')
ta=axmm(4.5,219,5.4,35);dendrogram(heat['tree'],orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#8A9299',ax=ta)
for coll in ta.collections:coll.set_linewidth(.15)
ta.invert_yaxis();ta.axis('off')
ordergenes=heat['genes'][leaves].tolist();selected=['TYMP','TK1','SLC29A1','SLC29A2','ODC1','SRM','SMS','MTAP','CHKA','PCYT1A','PCYT2','CEPT1']
labels=[]
for gene in selected:
 row=ordergenes.index(gene);y=219+35*(row+.5)/n;labels.append([gene,row,y,y])
labels.sort(key=lambda a:a[2]);gap=2.85
for k in range(1,len(labels)):labels[k][3]=max(labels[k][3],labels[k-1][3]+gap)
labels[-1][3]=min(labels[-1][3],254)
for k in range(len(labels)-2,-1,-1):labels[k][3]=min(labels[k][3],labels[k+1][3]-gap)
assert labels[0][3]>=218
for gene,row,y,yp in labels:
 fig.add_artist(Line2D([74.5/W,79/W,84/W],[1-y/H,1-y/H,1-yp/H],transform=fig.transFigure,c='#66727D',lw=.4))
 txt(85,yp,gene,fontsize=6.5,va='center',fontweight='bold'if gene=='TYMP'else'normal')
bar(im,16,260,25,'RNA row z-score',[-2.5,0,2.5])
txt(50,260,'2,511 genes; 20 patients',fontsize=6.5)
for jj,col in enumerate(clcols):
 fig.add_artist(Rectangle(((50+jj*14)/W,1-265.8/H),2/W,1.7/H,transform=fig.transFigure,color=col,ec='none'))
 txt(53+jj*14,264.1,f'C{jj+1}',fontsize=6.5)

# K: all39 prespecified endpoints; star denotes nominal P only.
panel('K',119,202,'Terminal enzyme RNA')
genesK=['TYMP','RRM2B','TK1','SLC29A1','SLC29A2','ODC1','SRM','SMS','MTAP','CHKA','PCYT1A','PCYT2','CEPT1'];methods=['Monocle3','Slingshot','scTour']
matrix=term.pivot(index='gene',columns='method',values='median_effect').loc[genesK,methods]
ax=axmm(140,216,44,37)
im=ax.imshow(matrix,cmap=div,vmin=-.7,vmax=.7,aspect='auto',interpolation='nearest')
ax.set_yticks(range(13),genesK,fontsize=6.5)
ax.set_xticks(range(3),methods,rotation=35,ha='right',fontsize=6.5);ax.tick_params(length=0,pad=1.5)
for sp in ax.spines.values():sp.set_visible(False)
stars=[]
for ii,gene in enumerate(genesK):
 for jj,method in enumerate(methods):
  v=term[(term.gene==gene)&(term.method==method)].iloc[0]
  if v.P<.05:
   ax.text(jj,ii,'*',fontsize=10,fontweight='bold',ha='center',va='center',color='white'if abs(v.median_effect)>.35 else'#343941');stars.append((method,gene))
assert len(stars)==4 and (term.FDR>=.05).all()
bar(im,150,208.5,28,'Delta log2(CPM+1)',[-.7,0,.7])
txt(121,264,'n=20; * P<0.05; all q>=0.05',fontsize=6.5)

fig.canvas.draw()
from matplotlib.text import Text
renderer=fig.canvas.get_renderer()
outside=[]
font_sizes=[]
for text in fig.findobj(match=Text):
 if not text.get_visible() or not text.get_text().strip():continue
 # axis tick artists hidden by an invisible axis are not rendered
 if text.axes is not None and not text.axes.axison:continue
 box=text.get_window_extent(renderer)
 font_sizes.append(float(text.get_fontsize()))
 if box.x0<0 or box.y0<0 or box.x1>fig.bbox.width or box.y1>fig.bbox.height:
  outside.append({'text':text.get_text(),'bbox':[float(v)for v in box.bounds]})
(w/'qc/RENDER_TEXT_BOUNDS.json').write_text(json.dumps({'outside_page':outside,'minimum_font_size':min(font_sizes),'canvas_pixels':[fig.bbox.width,fig.bbox.height]},indent=2))
# Save independently checkable mappings and display statistics.
pd.DataFrame(densstats).to_csv(w/'source_data/F1_ACTIVITY_PANEL_STATISTICS.tsv',sep='\t',index=False)
reg.to_csv(w/'source_data/REGULON_PANEL_STATISTICS.tsv',sep='\t',index=False)
tya[tya.type=='metaVIPER_activity'].to_csv(w/'source_data/TYMP_ACTIVITY_PANEL_STATISTICS.tsv',sep='\t',index=False)
auc.to_csv(w/'source_data/METABOLIC_PANEL_STATISTICS.tsv',sep='\t',index=False)
term.to_csv(w/'source_data/TERMINAL_PANEL_STATISTICS.tsv',sep='\t',index=False)
pd.DataFrame(labels,columns=['gene','matrix_row','original_y_mm','label_y_mm']).to_csv(w/'source_data/HEATMAP_LABEL_POSITIONS.tsv',sep='\t',index=False)
pd.DataFrame({'cell_id':calls.exact_cell_id,'class':statuses,'UMAP1':cxy[:,0],'UMAP2':cxy[:,1]}).to_csv(w/'source_data/CNV_UMAP_POINTS.tsv.gz',sep='\t',index=False)
metadata={'width_mm':W,'height_mm':H,'panels':panel_records,'no_new_tests':True,'libraries':75,'cnv_scope_cells':len(calls),'malignant_cells':len(cells),'malignant_patients':int(m.curated_patient_id.nunique()),'all39_terminal_endpoints':True,'nominal_star_endpoints':stars,'heatmap_genes':n,'font_min_pt':6.5,'full_source_statistics_preserved':True,'tymp_display_vmax_positive99percentile':float(np.quantile(expr[:,genes.index('TYMP')][expr[:,genes.index('TYMP')]>0],.99))}
(w/'qc/FIGURE_METADATA.json').write_text(json.dumps(metadata,indent=2))
for ext in ['pdf','png','svg']:
 path=w/'figures'/('Figure1_TYMP_integrated_A4_v1.'+ext)
 fig.savefig(path,dpi=600)
 print('SAVED',path.name,path.stat().st_size,flush=True)
plt.close(fig)
manifest=[]
for path in sorted(set(sources)):
 manifest.append({'path':str(path),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
(w/'qc/SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
(w/'qc/SOFTWARE_VERSIONS.json').write_text(json.dumps({'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__,'h5py':h5py.__version__},indent=2))
print('DRAW_COMPLETE',flush=True)
