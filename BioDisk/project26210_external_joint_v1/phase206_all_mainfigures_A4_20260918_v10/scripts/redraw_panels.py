from pathlib import Path
import json,ast
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Wedge,FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap,to_rgb
from matplotlib.lines import Line2D
from matplotlib.path import Path as MPath
r=Path.cwd();dest=r/'phase206_all_mainfigures_A4_20260918_v10'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':7,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
cm=LinearSegmentedColormap.from_list('div',['#568AA7','#FAF8F0','#A5417B'])
def save(fig,key):
 fig.savefig(dest/'panels'/f'{key}.pdf',dpi=400);plt.close(fig);print(key,flush=True)
# Fig2A: frozen displayed values, designed at 96 mm width.
s=pd.read_csv(r/'phase198_fig2_contrast_corr_revision_20260918_v5/tables/TOP20_DISPLAY_VALUES.tsv',sep='\t');paths=json.loads((r/'phase198_fig2_contrast_corr_revision_20260918_v5/tables/TOP20_SELECTION.json').read_text())['paths'];lim=s.row_z.abs().max()
fig=plt.figure(figsize=(96/25.4,141/25.4))
short={'Pentose and glucuronate interconversions':'Pentose / glucuronate','Valine, leucine and isoleucine degradation':'Branched-chain AA degradation','Glycine, serine and threonine metabolism':'Glycine / serine / threonine','Alanine, aspartate and glutamate metabolism':'Alanine / aspartate / glutamate','Drug metabolism - other enzymes':'Drug metabolism: other enzymes','Drug metabolism - cytochrome P450':'Drug metabolism: CYP450','Citrate cycle (TCA cycle)':'TCA cycle','Glycolysis / Gluconeogenesis':'Glycolysis / gluconeogenesis'}
for side,group in enumerate(['TYMP_negative','TYMP_positive']):
 ax=fig.add_axes([.47+side*.25,.14,.23,.80]);sub=s[s.group==group].set_index(['pathway','factor'])
 for i,path in enumerate(paths):
  for j in range(7):
   key=(path,'F'+str(j+1))
   if key not in sub.index:ax.plot(j,i,'x',color='#BBBBBB',ms=2);continue
   v=sub.loc[key];sc=ax.scatter(j,i,s=v.fraction*17,c=[v.row_z],cmap='RdBu_r',vmin=-lim,vmax=lim,lw=.25,edgecolors='#555555')
 ax.set(xlim=(-.55,6.55),ylim=(19.5,-.5));ax.set_xticks(range(7),['F'+str(i+1) for i in range(7)],fontsize=6);ax.set_yticks(range(20),[short.get(p,p) for p in paths] if side==0 else ['']*20,fontsize=6);ax.tick_params(length=0,pad=2);ax.grid(color='#E9E9E9',lw=.4);ax.set_axisbelow(True);ax.set_title(['TYMP-','TYMP+'][side],fontsize=7,pad=5)
ca=fig.add_axes([.50,.075,.22,.015]);cb=fig.colorbar(sc,cax=ca,orientation='horizontal',ticks=[-lim,0,lim]);cb.ax.set_xticklabels([f'{-lim:.1f}','0',f'{lim:.1f}']);cb.ax.tick_params(labelsize=6,length=2);ca.set_title('Row z',fontsize=6,pad=2)
fig.legend(handles=[Line2D([],[],ls='',marker='o',ms=np.sqrt(17*v),mfc='#92AABB',mec='#555555',mew=.3,label=f'{v:.0%}') for v in [.1,.5,1]],loc='lower right',bbox_to_anchor=(1,.014),frameon=False,ncol=3,fontsize=6,handletextpad=.2,columnspacing=.5)
save(fig,'2A')
# Fig2C: same patient median rho and direction counts; detailed P/FDR stay in caption.
q=r/'phase131_sections2_5_2_6_reference_correction_remote_v1/data';cs=pd.read_csv(q/'CORRPLOT_SUMMARY.tsv',sep='\t');eff=pd.read_csv(q/'CORRPLOT_PATIENT_EFFECTS.tsv',sep='\t');cols=['F1','AUC_Pyrimidine','AUC_CysMet','AUC_ArgPro','AUC_Glycerophospholipid'];names=['F1','Pyrimidine','Cys/Met','Arg/Pro','GPL'];mat=cs.pivot(index='x',columns='y',values='rho').loc[cols,cols].values
fig=plt.figure(figsize=(88/25.4,70/25.4));ax=fig.add_axes([.20,.20,.62,.72])
for i in range(5):
 for j in range(5):
  if i>j:
   v=mat[i,j];ax.scatter(j,i,s=abs(v)*400,c=[v],cmap='RdBu_r',vmin=-.5,vmax=.5);ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=6,color='white' if abs(v)>.34 else '#333333')
  elif i<j:
   vv=eff[(eff.x==cols[i])&(eff.y==cols[j])].rho;assert len(vv)==20;k=int((vv>0).sum());ax.add_patch(Rectangle((j-.34,i+.12),.68,.10,fc='#E6E9EA'));ax.add_patch(Rectangle((j-.34,i+.12),.68*k/20,.10,fc='#408E91'));ax.text(j,i-.10,f'{k}/20',ha='center',va='center',fontsize=6)
  else:ax.text(j,i,'-',ha='center',va='center',fontsize=6,color='#AAAAAA')
ax.set(xlim=(-.5,4.5),ylim=(4.5,-.5),aspect='equal');ax.set_xticks(range(5),names,rotation=35,ha='right',fontsize=6);ax.set_yticks(range(5),names,fontsize=6);ax.tick_params(length=0);ax.set_xticks(np.arange(-.5,5),minor=True);ax.set_yticks(np.arange(-.5,5),minor=True);ax.grid(which='minor',color='#E8E8E8',lw=.4)
for sp in ax.spines.values():sp.set_visible(False)
ca=fig.add_axes([.87,.40,.022,.38]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-.5,.5),cmap='RdBu_r'),cax=ca,ticks=[-.5,0,.5]);cb.ax.tick_params(labelsize=6,length=2);ca.set_title('rho',fontsize=6)
save(fig,'2C')
# Fig2E: same concentric glyphs, no paragraph-style annotations.
v=pd.read_csv(r/'phase202_candidate_labels_comparison_20260918_v6/tables/CANDIDATE_VS_WHOLE_SUMMARY.tsv',sep='\t')
fig=plt.figure(figsize=(190/25.4,45/25.4));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,190),ylim=(0,45),aspect='equal');ax.axis('off');colors=['#993F80','#146A98','#397348','#B9852D','#7D8993']
for x,z,col,name in zip(np.linspace(17,173,5),v.itertuples(),colors,['C4','C10','C15','Union','All macrophages']):
 ax.text(x,42,name,ha='center',fontsize=7,fontweight='bold',color=col);ax.text(x,37.5,f'n={int(z.n_expression)}',ha='center',fontsize=6,color='#606F7B')
 for radius,width,val,cc in [(11.3,2,z.contribution,col),(8.4,1.8,z.detection,tuple(.52*np.array(to_rgb(col))+.48))]:
  ax.add_patch(Wedge((x,23),radius,0,360,width=width,fc='#E7ECEF'))
  if val>0:ax.add_patch(Wedge((x,23),radius,90-3.6*val,90,width=width,fc=cc))
 ax.text(x,23,f'{z.expression:.2f}',ha='center',va='center',fontsize=10,fontweight='bold',color='#293C49');ax.text(x,7.5,f'{z.detection:.1f}% / {z.contribution:.1f}%',ha='center',fontsize=6.5,color=col)
ax.text(95,1,'Center: SPP1   |   Inner: positive cells   |   Outer: RNA share (n=16)',ha='center',fontsize=6.5,color='#53616D');save(fig,'2E')
# Fig2F: unchanged tissue order and row standardization.
bulk=pd.read_csv(r/'phase192_sections28_212_figures_20260917_v1/tables/TISSUE_MTA_ORDERED_PATIENTS.tsv',sep='\t',index_col=0);features=['MET:5-methylthioadenosine (MTA)','SPP1','LST1','AIF1','CSF1R','CD68','FCER1G','SMS','MTAP'];fig=plt.figure(figsize=(190/25.4,51/25.4));start=0
for i,(co,dat) in enumerate(bulk.groupby('cohort',sort=False)):
 ax=fig.add_axes([.075+.86*start/255,.07,.86*len(dat)/255-.008,.78]);vv=dat[features].values.T.copy();vv[0]=np.log2(vv[0]+1);zz=(vv-vv.mean(1,keepdims=True))/vv.std(1,keepdims=True);im=ax.imshow(zz,aspect='auto',cmap=cm,vmin=-2.5,vmax=2.5);ax.set_title(co.replace('_',' ')+f' (n={len(dat)})',fontsize=7);ax.set_yticks(range(9),['MTA']+features[1:] if i==0 else ['']*9,fontsize=6.5);ax.set_xticks([]);ax.tick_params(length=0);ax.add_patch(Rectangle((-.5,.5),len(dat),1,fill=False,ec='#993F80',lw=.9));start+=len(dat)
ca=fig.add_axes([.953,.19,.009,.56]);cb=fig.colorbar(im,cax=ca,ticks=[-2,0,2]);cb.ax.tick_params(labelsize=6,length=2);ca.set_title('Row z',fontsize=6);save(fig,'2F')
# Fig3A marker heatmap: same row z across the full 12-state matrix.
v=pd.read_csv(r/'phase192_sections28_212_figures_20260917_v1/tables/MARKER_PATIENT_EQUAL_LOG2CPM.tsv',sep='\t',index_col=0);zz=(v.values-v.values.mean(1,keepdims=True))/v.values.std(1,keepdims=True)
fig=plt.figure(figsize=(90/25.4,115/25.4));ax=fig.add_axes([.22,.065,.66,.91]);im=ax.imshow(zz,aspect='auto',cmap=cm,vmin=-2.5,vmax=2.5);ax.set_yticks(range(len(v)),v.index,fontsize=6,fontstyle='italic');ax.set_xticks(range(12),v.columns,fontsize=6,rotation=45,ha='right');ax.tick_params(length=0,pad=2);ax.add_patch(Rectangle((list(v.columns).index('C4')-.5,-.5),1,len(v),fill=False,ec='#993F80',lw=1))
ca=fig.add_axes([.92,.43,.020,.30]);cb=fig.colorbar(im,cax=ca,ticks=[-2,0,2]);cb.ax.tick_params(labelsize=6,length=2);ca.set_title('Row z',fontsize=6);save(fig,'3A')
# Fig1E frozen 24-regulon matrix, same order, all genes retained.
v=pd.read_csv(r/'phase130_section2_4_BE_correction_remote_v1/data/HEATMAP_VALUES_UNCHANGED.tsv',sep='\t',index_col=0);fig=plt.figure(figsize=(90/25.4,101/25.4));ax=fig.add_axes([.13,.12,.62,.79]);im=ax.imshow(v.values,aspect='auto',cmap=LinearSegmentedColormap.from_list('reg',['#5E936F','#F7F6F2','#927093']),vmin=-2,vmax=2)
ax.set_yticks(range(24),[z.replace('(+)','') for z in v.index],fontsize=6);ax.yaxis.tick_right();ax.set_xticks(np.arange(7)*2+.5,['F'+str(i+1) for i in range(7)],fontsize=6);ax.tick_params(length=0,pad=2)
for t in ax.get_yticklabels():
 if t.get_text() in ['ILF2','HOXD3','HNF4A']:t.set_fontweight('bold')
for k in np.arange(1.5,13,2):ax.axvline(k,color='white',lw=.8)
raw=(r/'phase130_section2_4_BE_correction_remote_v1/scripts/redraw_B_annotations.py').read_text();prefix=raw[:raw.index('plt.rcParams.update')];scope={'__file__':str(r/'phase130_section2_4_BE_correction_remote_v1/scripts/redraw_B_annotations.py')};exec(compile(prefix,'frozen_regulon_tree','exec'),scope);assert np.allclose(v.values,scope['Zall'][scope['leaves']]);ta=fig.add_axes([.005,.12,.105,.79]);scope['dendrogram'](scope['tree'],orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#888888',ax=ta);ta.invert_yaxis();ta.axis('off')
ca=fig.add_axes([.20,.955,.37,.018]);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-2,0,2]);cb.ax.tick_params(labelsize=6,length=2);fig.text(.025,.97,'Row z',fontsize=6);fig.text(.025,.025,'G3 / GB2',fontsize=6);save(fig,'1E')
# Fig1G: original KDE grids.
fig,axs=plt.subplots(2,2,figsize=(90/25.4,82/25.4));fig.subplots_adjust(left=.16,right=.86,bottom=.14,top=.92,hspace=.7,wspace=.55);pc=LinearSegmentedColormap.from_list('protein',['#FFFFFF','#D4C8E8','#7D63B1','#2C155C'])
for ax,g in zip(axs.ravel(),['CASP9','CHMP1A','PRSS8','PXN']):
 z=np.load(r/'phase131_sections2_5_2_6_reference_correction_remote_v1/data'/f'F1_{g}_DENSITY_GRID.npz');lev=[.025,.06,.12,.22,.36,.52,.70,.86,1.000001];im=ax.contourf(z['x'],z['y'],z['relative_density'],levels=lev,cmap=pc,vmin=0,vmax=1);ax.contour(z['x'],z['y'],z['relative_density'],levels=lev[:-1],colors='#695279',linewidths=.25);ax.set_xticks([0,.5,1]);ax.set_title(g,fontsize=7,fontweight='bold',pad=3);ax.tick_params(labelsize=6,pad=1,length=2)
fig.text(.02,.5,'metaVIPER NES',rotation=90,va='center',fontsize=7);fig.text(.49,.02,'F1 usage',ha='center',fontsize=7);ca=fig.add_axes([.91,.2,.016,.35]);cb=fig.colorbar(im,cax=ca,ticks=[.025,1]);cb.ax.tick_params(labelsize=6,length=2);save(fig,'1G')
# Fig1F: preserve all 71 nodes and edges; final-size font set explicitly.
o=r/'phase130_section2_4_BE_correction_remote_v1';raw=(o/'scripts/redraw_E_gene_network.py').read_text();fn=next(n for n in ast.parse(raw).body if isinstance(n,ast.FunctionDef) and n.name=='network');tfs=['ILF2','HOXD3','HNF4A'];colors={'Both methods':'#806080','GENIE3 only':'#B57756','GRNBoost2 only':'#397C9C'};nodecolor={'TF':'#6A263E','Target':'#8898A3','Shared target':'#C5913D'};shared={'PTER','PTN','TRPV4'};edges=pd.read_csv(o/'qc/GENE_NETWORK_DISPLAY_EDGES.tsv',sep='\t');nodes=pd.read_csv(o/'qc/GENE_NETWORK_NODE_POSITIONS.tsv',sep='\t');exec(compile(ast.Module(body=[fn],type_ignores=[]),'network','exec'))
fig,ax=plt.subplots(figsize=(94/25.4,78/25.4));fig.subplots_adjust(left=.02,right=.98,top=.98,bottom=.04);pos=network(ax,edges,nodes.gene.tolist(),'',size=6);ax.set_xlim(-1.65,1.65);ax.set_ylim(-1.52,1.52)
for coll in ax.collections:coll.set_sizes(coll.get_sizes()*.22)
for patch in ax.patches:patch.set_linewidth(.45)
save(fig,'1F')
print('CUSTOM_PANELS_READY')
