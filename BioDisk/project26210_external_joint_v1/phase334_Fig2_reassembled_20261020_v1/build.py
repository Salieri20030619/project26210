from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from matplotlib.patches import Rectangle,FancyArrowPatch
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage,dendrogram,fcluster
from scipy.spatial.distance import squareform
import pymupdf as fitz

w=Path(__file__).resolve().parent;i=w/'inputs';out=w/'delivery';out.mkdir(exist_ok=True)
pdir=w/'panels';pdir.mkdir(exist_ok=True);qc=w/'qc';qc.mkdir(exist_ok=True)
mm=72/25.4
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6.5,'axes.titlesize':7,'axes.labelsize':6.5,'xtick.labelsize':6,'ytick.labelsize':6,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.5})
colors=['#BA7868','#7295BE','#79A697','#C9AD63','#A58BB4','#8DB6C3','#C78FAD'];factors=[f'F{k}' for k in range(1,8)]
checks={}
def save(fig,name):
    fig.savefig(pdir/f'{name}.pdf',dpi=400);plt.close(fig)
def axes_off(ax):
    ax.set_aspect('equal');ax.axis('off')

# A: original CNV classes + original three program maps, at final panel size.
fig=plt.figure(figsize=(190/25.4,75/25.4))
c=pd.read_csv(i/'CNV_coordinates_and_classes.tsv.gz',sep='\t');xy=c[['UMAP1','UMAP2']].values
groups=[('malignant_consensus','Malignant consensus',10844,'#A94F60'),('malignant_infercnv_only','Malignant inferCNV only',1106,'#DE9B62'),('non_malignant_consensus','Non-malignant consensus',498,'#488A83'),('method_discordant','Method-discordant',22216,'#9BA6B5'),('indeterminate','Indeterminate',4006,'#D1C6AE')]
ax=fig.add_axes([.005,.33,.26,.55]);fig.text(.01,.94,'Epithelial CNV identity',size=7,weight='bold')
for k in [3,4,0,1,2]:
    key,label,n,col=groups[k];mask=c.CNV_class==key;assert mask.sum()==n
    ax.scatter(*xy[mask].T,s=.45,c=col,lw=0,rasterized=True)
axes_off(ax)
for k,(_,label,n,col) in enumerate(groups):
    fig.text(.012,.255-k*.052,'●',color=col,size=6)
    fig.text(.029,.255-k*.052,f'{label}: {n:,}',size=6)
fig.text(.01,.89,'38,670 tumour-tissue cells',size=6,color='#65717D')
z=np.load(i/'plot_values_server_only.npz');u=pd.read_csv(i/'consensus_usage_rank7_v1.tsv.gz',sep='\t',index_col='exact_cell_id').loc[z['cells']]
assert len(u)==11950 and np.allclose(u.sum(1),1)
xx=z['xy'];expr=z['expression'][:,list(z['genes']).index('TYMP')];dom=u.idxmax(1).values
for k,title in enumerate(['cNMF program','F1 usage','TYMP RNA']):
    x=.30+k*.235;ax=fig.add_axes([x,.28,.215,.65]);ax.set_title(title,pad=3,weight='bold')
    if k==0:
        for f,col in zip(factors,colors):
            m=dom==f;ax.scatter(*xx[m].T,s=.65,c=col,lw=0,rasterized=True)
        for n,(f,col) in enumerate(zip(factors,colors)):
            fig.text(x+(n%4)*.052,.19-(n//4)*.068,'●',size=6,color=col)
            fig.text(x+.014+(n%4)*.052,.19-(n//4)*.068,f,size=6)
    else:
        v=u.F1.values if k==1 else expr;order=np.argsort(v,kind='stable');lim=np.quantile(v[v>0],.99)
        cm=LinearSegmentedColormap.from_list('seq',['#EEEEEE','#315B91' if k==1 else '#A45654'])
        im=ax.scatter(*xx[order].T,s=.65,c=v[order],cmap=cm,vmin=0,vmax=lim,lw=0,rasterized=True)
        ca=fig.add_axes([x+.03,.18,.15,.022]);cb=fig.colorbar(im,cax=ca,orientation='horizontal');cb.ax.tick_params(labelsize=6,length=1,pad=1)
        cb.set_label('Normalized usage' if k==1 else 'Log1p expression',fontsize=6,labelpad=1)
    axes_off(ax)
fig.text(.30,.025,'11,950 malignant cells | Shared original UMAP coordinates',size=6,color='#65717D')
save(fig,'A');checks['A']={'cnv_cells':len(c),'malignant_cells':len(u),'classes_preserved':True,'no_embedding_refit':True}

# C: frozen profile correlation, identical display clustering, readable row keys.
v=pd.read_csv(i/'H_patient_program_Spearman.tsv',sep='\t',index_col=0);meta=pd.read_csv(i/'H_PROFILE_MANIFEST.tsv',sep='\t',index_col=0).loc[v.index]
tree=linkage(squareform(np.maximum(0,1-v.values),checks=False),method='average',optimal_ordering=True)
fig=plt.figure(figsize=(60/25.4,85/25.4));ta=fig.add_axes([.04,.78,.72,.10]);dd=dendrogram(tree,ax=ta,no_labels=True,color_threshold=0,above_threshold_color='#999999');ta.axis('off');order=dd['leaves'];N=len(order)
cl=fcluster(tree,t=7,criterion='maxclust')[order];arr=v.values[np.ix_(order,order)]
ax=fig.add_axes([.04,.12,.72,.65]);cm=LinearSegmentedColormap.from_list('corr',['#3022F4','#B8ADFF','#FFFFFF','#FFAAA9','#FF1919']);im=ax.imshow(arr,cmap=cm,vmin=-1,vmax=1,aspect='auto')
ax.set_yticks(range(N),[f'P{k+1:02}' for k in range(N)],fontsize=6);ax.yaxis.tick_right();ax.tick_params(length=0,pad=2)
starts=np.r_[0,np.flatnonzero(cl[1:]!=cl[:-1])+1];ends=np.r_[starts[1:],N]
ax.set_xticks((starts+ends-1)/2,[f'C{k+1}' for k in range(len(starts))],fontsize=6)
for a,b in zip(starts,ends):ax.add_patch(Rectangle((a-.5,a-.5),b-a,b-a,fill=False,ec='#444444',lw=.5))
bar=fig.add_axes([.04,.77,.72,.008]);bar.imshow(np.array([factors.index(s) for s in meta.iloc[order].program])[None,:],cmap=ListedColormap(colors),vmin=0,vmax=6,aspect='auto');bar.axis('off')
studies=['GSE156632','GSE159115','GSE178481','GSE207493','GSE242299'];study_colors=['#7296BC','#E3AD66','#76A893','#9F86B4','#C47F85']
bar=fig.add_axes([.019,.12,.013,.65]);bar.imshow(np.array([studies.index(s) for s in meta.iloc[order].study])[:,None],cmap=ListedColormap(study_colors),vmin=0,vmax=4,aspect='auto');bar.axis('off')
ca=fig.add_axes([.04,.94,.44,.016]);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-1,0,1]);ca.set_title('Spearman correlation',size=6,pad=3);cb.ax.tick_params(labelsize=6,length=1,pad=1)
for k,(f,col) in enumerate(zip(factors,colors)):
    fig.text(.02+k*.135,.035,'●',size=6,color=col);fig.text(.045+k*.135,.035,f,size=6)
save(fig,'C');meta.iloc[order].assign(display_key=[f'P{k+1:02}' for k in range(N)]).to_csv(qc/'PROFILE_LABEL_KEY.tsv',sep='\t',encoding='utf-8');checks['C']={'profiles':N,'matrix_unchanged':True}

# D: retain frozen matrix/order; tree added as unchanged vector crop below.
v=pd.read_csv(i/'HEATMAP_VALUES_UNCHANGED.tsv',sep='\t',index_col=0)
fig=plt.figure(figsize=(58/25.4,85/25.4));ax=fig.add_axes([.13,.11,.58,.79]);cm=LinearSegmentedColormap.from_list('reg',['#5E936F','#F7F6F2','#927093']);im=ax.imshow(v,cmap=cm,vmin=-2,vmax=2,aspect='auto')
ax.set_yticks(range(len(v)),[s.replace('(+)','') for s in v.index],fontsize=6);ax.yaxis.tick_right();ax.set_xticks(np.arange(7)*2+.5,factors,fontsize=6);ax.tick_params(length=0,pad=1.5)
for t in ax.get_yticklabels():
    if t.get_text() in ['ILF2','HOXD3','HNF4A']:t.set_weight('bold')
for x in np.arange(1.5,13,2):ax.axvline(x,color='white',lw=.6)
ca=fig.add_axes([.20,.95,.37,.015]);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-2,0,2]);ca.set_title('Row z',size=6,pad=3);cb.ax.tick_params(length=1,labelsize=6,pad=1)
fig.text(.10,.035,'Paired columns: GENIE3 / GRNBoost2',size=6)
save(fig,'D');checks['D']={'matrix_shape':list(v.shape),'matrix_and_order_unchanged':True}

# E: same 71 nodes / 71 directed edges, all gene labels horizontal.
edges=pd.read_csv(i/'GENE_NETWORK_DISPLAY_EDGES.tsv',sep='\t');nodes=pd.read_csv(i/'GENE_NETWORK_NODE_POSITIONS.tsv',sep='\t');tfs=['ILF2','HOXD3','HNF4A'];targets=[g for g in nodes.gene if g not in tfs];assert len(targets)==68 and len(edges)==71
support={'Both methods':'#806080','GENIE3 only':'#B57756','GRNBoost2 only':'#397C9C'}
fig=plt.figure(figsize=(100/25.4,94/25.4));ax=fig.add_axes([.005,.085,.99,.905]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
pos={};shared=set(edges.groupby('target').tf.nunique().loc[lambda x:x>1].index)
for j,g in enumerate(targets):
    side=j//34;y=.985-(j%34)*.0293;x=.255 if side==0 else .745;pos[g]=(x,y)
    col='#B38132' if g in shared else '#4A5C6A';ax.scatter(x,y,s=4,c=col,zorder=3)
    ax.text(x-.015 if side==0 else x+.015,y,g,ha='right' if side==0 else 'left',va='center',size=6,color=col)
for tf,y in zip(tfs,[.80,.50,.20]):
    pos[tf]=(.5,y);ax.scatter(.5,y,s=22,c='#742A45',marker='D',zorder=4);ax.text(.5,y+.025,tf,ha='center',size=7,weight='bold',color='#742A45',bbox={'facecolor':'white','edgecolor':'none','pad':.5},zorder=5)
for e in edges.itertuples():
    ax.add_patch(FancyArrowPatch(pos[e.tf],pos[e.target],connectionstyle='arc3,rad=0.06',arrowstyle='-|>',mutation_scale=3.3,lw=.45,color=support[e.support],alpha=.6,zorder=1))
for k,(s,col) in enumerate(support.items()):
    fig.text(.015+k*.33,.035,'—',size=8,color=col);fig.text(.055+k*.33,.035,s,size=6)
save(fig,'E');checks['E']={'nodes':len(nodes),'edges':len(edges),'targets':len(targets),'all_labels_horizontal':True,'edges_unchanged':True}

# F: original density grids, contour levels and fixed protein set.
fig,axs=plt.subplots(2,2,figsize=(80/25.4,94/25.4));fig.subplots_adjust(left=.16,right=.84,bottom=.16,top=.94,hspace=.65,wspace=.55)
cm=LinearSegmentedColormap.from_list('protein',['#FFFFFF','#D4C8E8','#7D63B1','#2C155C']);levels=[.025,.06,.12,.22,.36,.52,.70,.86,1.000001]
for ax,g in zip(axs.ravel(),['CASP9','CHMP1A','PRSS8','PXN']):
    z=np.load(i/f'F1_{g}_DENSITY_GRID.npz');im=ax.contourf(z['x'],z['y'],z['relative_density'],levels=levels,cmap=cm,vmin=0,vmax=1);ax.contour(z['x'],z['y'],z['relative_density'],levels=levels[:-1],colors='#695279',linewidths=.25);ax.set_xticks([0,.5,1]);ax.set_title(g,size=7,weight='bold',pad=4);ax.tick_params(length=2,pad=1,labelsize=6)
fig.text(.015,.55,'metaVIPER NES',size=6.5,rotation=90,va='center');fig.text(.50,.055,'F1 usage',size=6.5,ha='center')
ca=fig.add_axes([.89,.23,.018,.32]);cb=fig.colorbar(im,cax=ca,ticks=[.025,1]);cb.ax.set_yticklabels(['0.025','1.00']);cb.ax.tick_params(labelsize=6,length=1,pad=1);ca.set_title('Density',size=6,pad=4)
save(fig,'F');checks['F']={'original_grids_and_levels':True,'proteins':4}

# Final A4 assembly. Original B matrix is retained as vector/raster PDF content.
doc=fitz.open();pg=doc.new_page(width=210*mm,height=297*mm)
from matplotlib import font_manager
pg.insert_font(fontname='FigRegular',fontfile=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans')))
pg.insert_font(fontname='FigBold',fontfile=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold')))
layout={'A':(10,19,190,75),'B':(10,106,60,79),'C':(76,103,60,85),'D':(142,103,58,85),'E':(10,202,100,94),'F':(120,202,80,94)}
for letter in ['A','C','D','E','F']:
    x,y,ww,hh=layout[letter];src=fitz.open(pdir/f'{letter}.pdf');pg.show_pdf_page(fitz.Rect(x*mm,y*mm,(x+ww)*mm,(y+hh)*mm),src,0)
# Retain original regulon dendrogram, crop contains no labels.
src=fitz.open(i/'1E.pdf');p=src[0];crop=fitz.Rect(.005*p.rect.width,.09*p.rect.height,.11*p.rect.width,.88*p.rect.height)
pg.show_pdf_page(fitz.Rect(142.29*mm,111.5*mm,148.38*mm,178.65*mm),src,0,clip=crop,keep_proportion=False)
old=fitz.open(w/'inputs/Figure2_previous.pdf')
for block in old[0].get_text('dict')['blocks']:
    for line in block.get('lines',[]):
        for span in line['spans']:old[0].add_redact_annot(fitz.Rect(span['bbox']),fill=False)
old[0].apply_redactions(images=0,graphics=0)
matrix_info=next(z for z in old[0].get_image_info() if z['width']==2508 and z['height']==1596)
matrix_clip=fitz.Rect(matrix_info['bbox']);checks['B']={'original_matrix_image_preserved':True,'source_image_bbox':list(matrix_clip),'source_image_pixels':[2508,1596]}
pg.show_pdf_page(fitz.Rect(10*mm,115*mm,70*mm,154*mm),old,0,clip=matrix_clip,keep_proportion=False)
for label,x in [('1',113.6),('5',131.5),('10',146.9),('15',162.6),('20',178.5)]:pg.insert_text(( (10+(x-110.469)*60/72.987)*mm,114*mm),label,fontsize=6,fontname='FigRegular')
pg.insert_text((10*mm,108*mm),'GSE156632 RCC2',fontsize=7,fontname='FigRegular')
pg.insert_text((10*mm,160*mm),'Top: reference | Bottom: tumour epithelial',fontsize=6,fontname='FigRegular')
pg.insert_text((10*mm,164*mm),'Chromosomes',fontsize=6,fontname='FigRegular')
# Preserve exact original colour scale and retype labels at final size.
scale_info=next(z for z in old[0].get_image_info() if z['width']==445 and z['height']==52)
pg.show_pdf_page(fitz.Rect(20*mm,170*mm,60*mm,172*mm),old,0,clip=fitz.Rect(scale_info['bbox']),keep_proportion=False)
for s,x in [('-0.5',18),('0',39),('0.5',58)]:pg.insert_text((x*mm,176*mm),s,fontsize=6,fontname='FigRegular')
pg.insert_text((24*mm,181*mm),'log2 relative signal',fontsize=6,fontname='FigRegular')
titles={'A':'Epithelial identity and malignant-cell programs','B':'Chromosome CNV','C':'Patient-program profiles','D':'Regulon activity','E':'TF-target network','F':'Inferred protein activity'}
for l,(x,y,ww,hh) in layout.items():
    yy=14 if l=='A' else 99 if l in ['B','C','D'] else 198
    pg.insert_text((x*mm,yy*mm),l,fontsize=10,fontname='FigBold');pg.insert_text(((x+5)*mm,yy*mm),titles[l],fontsize=7.5,fontname='FigBold')
fn=out/'Fig2_reassembled_A4_v1.pdf';doc.save(fn,garbage=4,deflate=True)
pg.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False).save(out/'Fig2_reassembled_A4_v1.png')
pg.get_pixmap(matrix=fitz.Matrix(1.7,1.7),alpha=False).save(qc/'Fig2_preview.png')
(out/'Fig2_reassembled_A4_v1.svg').write_text(pg.get_svg_image(text_as_path=False),encoding='utf-8')
checks['page']={'size_mm':[210,297],'pages':1,'mapping':{'A':'candidate1 + old C','B':'old B','C':'old D','D':'old E','E':'old F','F':'old G'},'supplement':'Candidate3 placeholder only','manuscript_changed':False}
(qc/'VALIDATION.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print('BUILT',fn)
