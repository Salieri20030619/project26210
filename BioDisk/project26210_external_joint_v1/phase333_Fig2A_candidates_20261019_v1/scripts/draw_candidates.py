from pathlib import Path
import json, hashlib, sys
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

root=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
out=root/'phase333_Fig2A_candidates_20261019_v1'
for s in ['delivery','source_data','qc','scripts']:(out/s).mkdir(parents=True,exist_ok=True)
base=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
hpath=base/'02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad'
calls_path=base/'03_malignancy_cnv/malignancy_cnv_summary_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz'
loading_path=base/'05_cnmf/full_ranks4_8_v1/consensus_loading_rank7_v1.tsv.gz'
fixed_path=root/'phase123_section2_2_reference_redraw_remote_v1/qc/I_gene_contribution_aggregate.tsv'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
checks={}
def save(fig,name):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    bad=[]
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(renderer)
        if b.x0 < -1 or b.y0 < -1 or b.x1>fig.bbox.width+1 or b.y1>fig.bbox.height+1:bad.append(t.get_text())
    assert not bad,bad
    for ext in ['pdf','png','svg']:fig.savefig(out/'delivery'/f'{name}.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
    checks[name]={'text_inside_canvas':True}

calls=pd.read_csv(calls_path,sep='\t',index_col='exact_cell_id')
with h5py.File(hpath,'r') as h:
    idxkey=h['obs'].attrs.get('_index','_index')
    ids=pd.Index(h['obs'][idxkey].asstr()[:])
    pos=ids.get_indexer(calls.index)
    assert (pos>=0).all() and ids.is_unique and calls.index.is_unique
    xy=h['obsm/X_umap'][:][pos]
    c=h['obs/cell_type_broad_v1']
    labels=c['categories'].asstr()[:][c['codes'][:][pos]]
    assert (labels=='Epithelial cell').all()
groups=[('malignant_consensus','Malignant: consensus',10844,'#A94F60'),('malignant_infercnv_only','Malignant: inferCNV only',1106,'#DE9B62'),('non_malignant_consensus','Non-malignant: consensus',498,'#488A83'),('method_discordant','Method-discordant',22216,'#9BA6B5'),('indeterminate','Indeterminate',4006,'#D1C6AE')]
status=calls.malignancy_status_v1.to_numpy()
assert set(status)==set(x[0] for x in groups) and len(calls)==38670 and np.isfinite(xy).all()
fig=plt.figure(figsize=(190/25.4,116/25.4))
fig.text(.04,.945,'Epithelial malignancy classification',size=12,weight='bold',color='#263A49')
fig.text(.04,.895,'38,670 tumour-tissue epithelial cells | Frozen UMAP coordinates',size=8,color='#52616E')
ax=fig.add_axes([.035,.12,.50,.71])
# Background classes first; confirmed classes last. All cells retained.
for i in [3,4,0,1,2]:
    key,title,n,color=groups[i];mask=status==key;assert int(mask.sum())==n
    ax.scatter(xy[mask,0],xy[mask,1],s=1.4 if i in [3,4] else 2,c=color,lw=0,rasterized=True,alpha=.8)
ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
for sp in ax.spines.values():sp.set_visible(False)
ax.annotate('',xy=(.20,.025),xytext=(.045,.025),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.7,'color':'#46525D'})
ax.annotate('',xy=(.045,.18),xytext=(.045,.025),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.7,'color':'#46525D'})
ax.text(.12,-.035,'UMAP1',transform=ax.transAxes,ha='center',size=7)
ax.text(-.005,.18,'UMAP2',transform=ax.transAxes,ha='left',size=7)
leg=fig.add_axes([.565,.17,.40,.65]);leg.axis('off')
for i,(key,title,n,col) in enumerate(groups):
    y=.92-i*.18
    leg.scatter(.015,y,s=35,color=col,edgecolor='none',clip_on=False)
    leg.text(.065,y,title,va='center',size=8)
    leg.text(.065,y-.065,f'{n:,} cells  |  {n/len(calls)*100:.2f}%',size=7.5,color='#61717E')
leg.set(xlim=(0,1),ylim=(0,1))
fig.text(.04,.04,'Malignant analysis set: 11,950 cells = 10,844 consensus + 1,106 inferCNV-only',size=7.5,color='#52616E')
save(fig,'Candidate_1_epithelial_malignancy')
pd.DataFrame({'exact_cell_id':calls.index,'UMAP1':xy[:,0],'UMAP2':xy[:,1],'CNV_class':status}).to_csv(out/'source_data/CNV_coordinates_and_classes.tsv.gz',sep='\t',index=False,encoding='utf-8')

fixed=pd.read_csv(fixed_path,sep='\t',index_col=0)
load=pd.read_csv(loading_path,sep='\t',index_col='gene')
factors=[f'F{i}' for i in range(1,8)]
assert fixed.shape==(35,7) and list(fixed.columns)==factors and fixed.index.is_unique
expected=load.loc[fixed.index,factors].div(load.loc[fixed.index,factors].max(axis=1),axis=0)
assert np.allclose(fixed,expected,atol=1e-12)
colors=['#BA7868','#7295BE','#79A697','#C9AD63','#A58BB4','#8DB6C3','#C78FAD']
cm=LinearSegmentedColormap.from_list('loading',['#F6F8FB','#CEDBEA','#6F96BD','#245383'])
fig=plt.figure(figsize=(190/25.4,116/25.4))
fig.text(.04,.945,'Representative genes of cNMF programs',size=12,weight='bold',color='#263A49')
fig.text(.04,.895,'F1-F7 | Five frozen representative genes per program',size=8,color='#52616E')
# Two aligned blocks keep all gene names horizontal at readable size.
for block,(lo,hi,x,w) in enumerate([(0,20,.18,.27),(20,35,.68,.27)]):
    n=hi-lo;height=n*.029
    ax=fig.add_axes([x,.80-height,w,height])
    im=ax.imshow(fixed.iloc[lo:hi].values,cmap=cm,vmin=0,vmax=1,aspect='auto',interpolation='nearest')
    ax.set_yticks(np.arange(n),fixed.index[lo:hi],fontsize=7)
    ax.set_xticks(range(7),factors,fontsize=8);ax.xaxis.tick_top();ax.tick_params(length=0,pad=3)
    for j in range(n):
        ax.add_patch(plt.Rectangle((6.62,j-.5),.18,1,facecolor=colors[(j+lo)//5],edgecolor='none',clip_on=False))
    for j in range(5,n,5):ax.axhline(j-.5,color='white',lw=1.7)
    for sp in ax.spines.values():sp.set_visible(False)
    ax.set_xticks(np.arange(-.5,7,1),minor=True);ax.set_yticks(np.arange(-.5,n,1),minor=True)
    ax.grid(which='minor',color='white',lw=.35);ax.tick_params(which='minor',length=0)
ca=fig.add_axes([.68,.21,.27,.025]);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[0,.5,1]);cb.ax.tick_params(length=2,labelsize=7)
ca.set_title('Relative loading (within gene)',size=8,pad=6)
for i,c in enumerate(colors):
    x=.055+i*.073
    fig.patches.append(plt.Rectangle((x,.115),.014,.018,transform=fig.transFigure,facecolor=c,edgecolor='none'))
    fig.text(x+.019,.115,f'F{i+1}',size=7,va='bottom')
fig.text(.04,.06,'Colour strip: source program of each representative gene.',size=7.2,color='#52616E')
fig.text(.04,.025,'Each gene is scaled to its maximum loading across F1-F7; values are not differential-expression statistics.',size=7.2,color='#52616E')
save(fig,'Candidate_3_cNMF_representative_genes')
fixed.to_csv(out/'source_data/frozen_gene_relative_loadings.tsv',sep='\t',encoding='utf-8')
sources=[]
for p in [calls_path,loading_path,fixed_path]:sources.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
sources.append({'path':str(hpath),'bytes':hpath.stat().st_size,'mtime_ns':hpath.stat().st_mtime_ns,'read_scope':'obs IDs/cell type and obsm/X_umap only'})
meta={'status':'DATA_CHECKS_PASS_VISUAL_REVIEW_PENDING','class_counts':{g[0]:g[2] for g in groups},'cells':len(calls),'cnmf_genes':35,'cnmf_programs':7,'max_loading_difference':float(np.abs(fixed.values-expected.values).max()),'no_new_models':True,'no_new_tests':True,'manuscript_modified':False,'fig2_modified':False,'sources':sources,'checks':checks,'versions':{'python':sys.version,'matplotlib':matplotlib.__version__,'numpy':np.__version__,'pandas':pd.__version__}}
(out/'qc/VALIDATION.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(meta,ensure_ascii=False))
