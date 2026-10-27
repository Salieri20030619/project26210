from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad,matplotlib,hashlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap
from matplotlib.patches import Patch,Rectangle
from scipy.stats import spearmanr,rankdata
from scipy.cluster.hierarchy import linkage,dendrogram,fcluster
from scipy.spatial.distance import squareform
w=Path(__file__).resolve().parents[1];base=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1');r=base/'05_cnmf';factors=[f'F{i}' for i in range(1,8)]
u=pd.read_csv(r/'full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz',sep='\t').set_index('exact_cell_id');m=pd.read_csv(r/'input_v1/cnmf_cell_metadata_v1.tsv.gz',sep='\t').set_index('exact_cell_id').loc[u.index];genes=pd.read_csv(r/'full_ranks4_8_v1/consensus_loading_rank7_v1.tsv.gz',sep='\t').gene.tolist();assert len(u)==11950 and len(genes)==2500
m['program']=u.idxmax(axis=1);m['profile_id']=m.curated_patient_id+' / '+m.program;counts=m.groupby('profile_id',sort=True).size();eligible=counts[counts>=20].index.tolist();profile_index={s:i for i,s in enumerate(eligible)}
records=[]
for pid in eligible:
 mm=m[m.profile_id==pid];records.append({'profile_id':pid,'patient':mm.curated_patient_id.iloc[0],'study':mm.study.iloc[0],'program':mm.program.iloc[0],'cells':len(mm)})
meta=pd.DataFrame(records).set_index('profile_id');sums=np.zeros((len(eligible),len(genes)),dtype=float);ns=np.zeros(len(eligible),dtype=int)
a=ad.read_h5ad(base/'02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad',backed='r');idx=a.obs_names.get_indexer(m.index);assert (idx>=0).all();order=np.argsort(idx);seen=[]
for offset in range(0,len(order),1200):
 ii=order[offset:offset+1200];x=a[idx[ii],genes].X.toarray();assert np.isfinite(x).all()
 pg=m.iloc[ii].profile_id.to_numpy()
 for pid in np.unique(pg):
  if pid not in profile_index:continue
  jj=profile_index[pid];mask=pg==pid;sums[jj]+=x[mask].sum(axis=0,dtype=float);ns[jj]+=int(mask.sum())
 seen.extend(ii.tolist())
a.file.close();assert len(set(seen))==11950 and np.array_equal(ns,meta.cells.to_numpy())
means=sums/ns[:,None];scale=means.std(axis=0);keep=scale>1e-12;z=(means[:,keep]-means[:,keep].mean(axis=0))/scale[keep]
# Profiles carry equal weight in the gene-wise scaling. No cell-level P values.
corr=spearmanr(z,axis=1).statistic;independent=np.corrcoef(rankdata(z,axis=1));assert np.max(np.abs(corr-independent))<1e-12;assert np.isfinite(corr).all() and np.allclose(corr,corr.T) and np.allclose(np.diag(corr),1)
tree=linkage(squareform(np.maximum(0,1-corr),checks=False),method='average',optimal_ordering=True);display_clusters=fcluster(tree,t=7,criterion='maxclust')
pd.DataFrame(means,index=eligible,columns=genes).to_csv(w/'data/H_patient_program_mean_expression.tsv.gz',sep='\t');pd.DataFrame(z,index=eligible,columns=np.array(genes)[keep]).to_csv(w/'data/H_gene_standardized_profiles.tsv.gz',sep='\t');pd.DataFrame(corr,index=eligible,columns=eligible).to_csv(w/'data/H_patient_program_Spearman.tsv',sep='\t');meta.to_csv(w/'qc/H_PROFILE_MANIFEST.tsv',sep='\t')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'axes.linewidth':.75});program_colors=['#BA7868','#7295BE','#79A697','#C9AD63','#A58BB4','#8DB6C3','#C78FAD'];studies=['GSE156632','GSE159115','GSE178481','GSE207493','GSE242299'];study_colors=['#7296BC','#E3AD66','#76A893','#9F86B4','#C47F85'];cmap=LinearSegmentedColormap.from_list('reference_H',['#3022F4','#B8ADFF','#FFFFFF','#FFAAA9','#FF1919'])
fig=plt.figure(figsize=(10,9.2));ta=fig.add_axes([.10,.785,.57,.12]);dd=dendrogram(tree,ax=ta,no_labels=True,color_threshold=0,above_threshold_color='#999999',link_color_func=lambda k:'#999999');ta.axis('off');leaves=np.array(dd['leaves']);ordered=meta.iloc[leaves];cl=display_clusters[leaves];N=len(meta)
ax=fig.add_axes([.10,.155,.57,.62]);im=ax.imshow(corr[np.ix_(leaves,leaves)],cmap=cmap,vmin=-1,vmax=1,aspect='auto',interpolation='nearest');ax.set_xticks(np.arange(-.5,N,1),minor=True);ax.set_yticks(np.arange(-.5,N,1),minor=True);ax.grid(which='minor',color='#FFFFFF',lw=.16);ax.tick_params(which='minor',length=0)
ax.set_yticks(range(N),[p.replace('__','  ') for p in ordered.index],fontsize=7.2);ax.yaxis.tick_right();ax.tick_params(axis='y',length=0,pad=4);starts=np.r_[0,np.flatnonzero(cl[1:]!=cl[:-1])+1];ends=np.r_[starts[1:],N];centers=(starts+ends-1)/2;ax.set_xticks(centers,['C'+str(i+1) for i in range(len(starts))]);ax.tick_params(axis='x',length=0,pad=5);ax.set_xlabel('Display clusters of patient × program profiles',fontsize=9,labelpad=9)
for lo,hi in zip(starts,ends):ax.add_patch(Rectangle((lo-.5,lo-.5),hi-lo,hi-lo,fill=False,edgecolor='#111111',lw=1.05))
program_bar=np.array([factors.index(s) for s in ordered.program]);study_bar=np.array([studies.index(s) for s in ordered.study]);at=fig.add_axes([.10,.775,.57,.009]);at.imshow(program_bar[None,:],cmap=ListedColormap(program_colors),vmin=0,vmax=6,aspect='auto',interpolation='nearest');at.axis('off');al=fig.add_axes([.087,.155,.011,.62]);al.imshow(study_bar[:,None],cmap=ListedColormap(study_colors),vmin=0,vmax=4,aspect='auto',interpolation='nearest');al.axis('off')
ca=fig.add_axes([.10,.946,.22,.015]);fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-1,0,1]);ca.set_title('Spearman correlation',fontsize=9,pad=5);fig.text(.018,.98,'H',fontsize=18,fontweight='bold',va='top')
fig.legend(handles=[Patch(facecolor=c,label=f) for f,c in zip(factors,program_colors)],loc='upper right',bbox_to_anchor=(.965,.98),frameon=False,ncol=7,title='Original cNMF program',fontsize=7,title_fontsize=8,handlelength=.85,columnspacing=.7)
fig.legend(handles=[Patch(facecolor=c,label=s) for s,c in zip(studies,study_colors) if s in ordered.study.values],loc='lower left',bbox_to_anchor=(.10,.050),frameon=False,ncol=5,fontsize=7,handlelength=.9,columnspacing=.9)
fig.text(.10,.025,f'{N} profiles from {meta.patient.nunique()} patients; ≥20 cells/profile; {keep.sum():,} variable genes. Gene-wise standardized mean RNA.',fontsize=8,color='#444444')
for ext in ['png','pdf']:fig.savefig(w/'figures'/('H_patient_program_expression_blocks_v2.'+ext),dpi=350,bbox_inches='tight',pad_inches=.1)
plt.close(fig)
result={'status':'PASS','profiles':N,'patients':int(meta.patient.nunique()),'studies':int(meta.study.nunique()),'cells_represented':int(ns.sum()),'source_cells':11950,'minimum_cells_per_profile':20,'genes_input':2500,'genes_variable':int(keep.sum()),'fixed_program_labels':True,'profile_definition':'within-patient arithmetic mean log1p RNA of cells with each frozen dominant cNMF factor','scaling':'each gene z-scored across eligible patient-program mean profiles; constant genes omitted','correlation':'Spearman across variable genes; independent rank-Pearson check passed','clustering':'average linkage on 1-r; seven display clusters C1-C7, not refitted cNMF or new biological types','no_correlation_interpolation_or_smoothing':True,'new_inferential_P_or_FDR':False,'reference_relation':'sample-associated expression profiles, not KINOMO independently inferred sample programs; old H was random-start loadings'}
(w/'qc/H_VALIDATION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
