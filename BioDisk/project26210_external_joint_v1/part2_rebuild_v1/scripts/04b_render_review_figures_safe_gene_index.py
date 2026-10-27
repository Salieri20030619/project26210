from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');C=R/(R/'provenance/CURRENT_CLUSTERING_RUN.txt').read_text().strip();M=R/(R/'provenance/CURRENT_MARKER_RUN.txt').read_text().strip();V=R/(R/'provenance/CURRENT_VALIDATION_RUN.txt').read_text().strip();O=R/'10_figures/run_v2';O.mkdir()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
z=ad.read_h5ad(C/'CLUSTERING_GRAPH.h5ad');ann=pd.read_csv(M/'CLUSTER_ANNOTATIONS.tsv',sep='\t',dtype={'cluster':str}).set_index('cluster');labels=z.obs.cluster.astype(str);clusters=sorted(labels.unique(),key=int);colors={c:plt.get_cmap('tab20')(i/20) for i,c in enumerate(clusters)};xy=z.obsm['X_umap'];rng=np.random.default_rng(262147);order=rng.permutation(len(z));files=[]
def save(fig,name):
 for ext in ['png','svg']:
  p=O/(name+'.'+ext);fig.savefig(p,dpi=300,bbox_inches='tight',facecolor='white');files.append(p.name)
 plt.close(fig)
def plain(ax):ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.set_xticks([]);ax.set_yticks([])
fig,axs=plt.subplots(1,3,figsize=(17,5),layout='constrained')
for c in clusters:
 ii=order[labels.to_numpy()[order]==c];axs[0].scatter(xy[ii,0],xy[ii,1],s=1.2,color=colors[c],rasterized=True,linewidths=0);med=np.median(xy[ii],axis=0);axs[0].text(*med,'C'+c,fontsize=10,weight='bold',ha='center',bbox=dict(facecolor='white',alpha=.75,edgecolor='none',pad=1))
axs[0].set_title('A  Exploratory fine clusters')
studies=sorted(z.obs.study.astype(str).unique());pal=plt.get_cmap('Set2')
for i,s in enumerate(studies):
 ii=order[z.obs.study.astype(str).to_numpy()[order]==s];axs[1].scatter(xy[ii,0],xy[ii,1],s=1.2,color=pal(i),label=s,rasterized=True,linewidths=0)
axs[1].legend(markerscale=5,frameon=False,fontsize=8,loc='upper right');axs[1].set_title('B  Study identity')
flag=z.obs.lineage_review_flag.to_numpy(bool);axs[2].scatter(xy[~flag,0],xy[~flag,1],s=1,color='#cbd5df',rasterized=True,linewidths=0);axs[2].scatter(xy[flag,0],xy[flag,1],s=4,color='#c45143',rasterized=True,linewidths=0);axs[2].set_title('C  Lineage review flags (retained)')
for ax in axs:plain(ax)
fig.suptitle(f'Part 2 rebuild | {len(z):,} candidate cells | 40 patients | 75 source libraries',fontsize=14);save(fig,'01_unsupervised_atlas')
comp=pd.read_csv(V/'CLUSTER_COMPOSITION_QC.tsv',sep='\t',dtype={'cluster':str}).set_index('cluster').loc[clusters];study=pd.crosstab(labels,z.obs.study).reindex(clusters);frac=study.div(study.sum(axis=1),axis=0)
fig,axs=plt.subplots(1,2,figsize=(13,max(4,len(clusters)*.45)),gridspec_kw={'width_ratios':[1,1.35]},layout='constrained');axs[0].barh(np.arange(len(clusters)),comp.n_cells,color=[colors[c] for c in clusters]);axs[0].set_yticks(np.arange(len(clusters)),['C'+c+'  '+ann.loc[c,'short_label'] for c in clusters]);axs[0].invert_yaxis();axs[0].set_xlabel('Cells (descriptive count)');axs[0].set_title('A  Candidate cluster sizes')
im=axs[1].imshow(frac.to_numpy(),vmin=0,vmax=1,cmap='Blues',aspect='auto');axs[1].set_xticks(range(len(study.columns)),study.columns,rotation=35,ha='right');axs[1].set_yticks(range(len(clusters)),['C'+c for c in clusters]);axs[1].set_title('B  Study composition within each cluster')
for i in range(len(clusters)):
 for j in range(len(study.columns)):axs[1].text(j,i,f'{100*frac.iloc[i,j]:.0f}%',ha='center',va='center',fontsize=8,color='white' if frac.iloc[i,j]>.6 else '#152a3a')
fig.colorbar(im,ax=axs[1],label='Fraction of cluster cells',shrink=.75);save(fig,'02_cluster_study_composition')
can=pd.read_csv(M/'CANONICAL_PANEL_ALL_CLUSTERS.tsv',sep='\t',dtype={'cluster':str});genes=['LST1','TYROBP','CSF1R','FCN1','S100A8','S100A9','C1QA','C1QC','APOE','APOC1','TREM2','SPP1','GPNMB','FABP5','CD9','LPL','FOLR2','LYVE1','MRC1','IL1B','CXCL8','CCL3','HLA-DRA','CD74','ISG15','IFIT1','MKI67','TOP2A','EPCAM','KRT18','CD3D','GNLY','MS4A1','VWF'];genes=list(dict.fromkeys(genes+ann.main_marker.astype(str).tolist()+['CCL18','CYP27A1']))
extra=[]
for c in clusters:
 t=pd.read_csv(M/'full_marker_tables'/('cluster_'+c+'_all_genes.tsv.gz'),sep='\t',dtype={'cluster':str});extra.append(t[t.gene.isin(genes)])
can=pd.concat(extra,ignore_index=True);genes=[g for g in genes if g in set(can.gene)]
mu=can.pivot(index='cluster',columns='gene',values='mean_log1pCP10K_in').reindex(index=clusters,columns=genes);pct=can.pivot(index='cluster',columns='gene',values='pct_in').reindex(index=clusters,columns=genes);scaled=(mu-mu.mean(axis=0))/mu.std(axis=0).replace(0,1);xx,yy=np.meshgrid(np.arange(len(genes)),np.arange(len(clusters)))
fig,ax=plt.subplots(figsize=(max(14,len(genes)*.38),max(4.5,len(clusters)*.5)),layout='constrained');scat=ax.scatter(xx.ravel(),yy.ravel(),s=(pct.to_numpy().ravel()*110),c=scaled.to_numpy().ravel(),cmap='RdBu_r',vmin=-2,vmax=2,edgecolors='none');ax.set_xticks(range(len(genes)),genes,rotation=65,ha='right');ax.set_yticks(range(len(clusters)),['C'+c+' '+ann.loc[c,'short_label'] for c in clusters]);ax.invert_yaxis();ax.set_xlim(-.7,len(genes)-.3);ax.set_ylim(len(clusters)-.3,-.7);ax.set_title('Post-cluster marker panel | dot size: fraction expressing; color: gene-scaled mean log1p CP10K');fig.colorbar(scat,ax=ax,label='Mean expression z score across clusters',shrink=.7);save(fig,'03_marker_dotplot')
# SPP1 examined only after clustering; no threshold determines labels.
a=ad.read_h5ad(R/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad');a=a[z.obs_names].copy();idx=a.var_names.get_loc('SPP1');raw=np.asarray(a.X[:,idx].toarray()).ravel();total=np.asarray(a.X.sum(axis=1)).ravel();expr=np.log1p(raw/total*10000);vmax=float(np.quantile(expr[expr>0],.99)) if (expr>0).any() else 1
fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained');oo=np.argsort(expr);s=axs[0].scatter(xy[oo,0],xy[oo,1],c=expr[oo],s=1.4,cmap='magma_r',vmin=0,vmax=vmax,rasterized=True,linewidths=0);plain(axs[0]);axs[0].set_title('A  SPP1 on candidate-cell UMAP');fig.colorbar(s,ax=axs[0],label='log1p CP10K',shrink=.7)
MM=R/'04_markers_annotation/run_v4_macrophage_only';macro_clusters=[c for c in clusters if str(ann.loc[c,'macrophage_comparator_eligible']).lower()=='true'];pair_ns=[]
for i,c in enumerate(macro_clusters):
 d=np.load(MM/'paired_patient_effects'/('cluster_'+c+'.npz'));gi=list(pd.read_csv(MM/'PSEUDOBULK_GENE_INDEX.tsv',sep='\t').gene).index('SPP1');v=d['effects'][:,gi];pair_ns.append(len(v));axs[1].scatter(i+rng.uniform(-.17,.17,len(v)),v,s=18,color=colors[c],alpha=.7);axs[1].plot([i-.25,i+.25],[np.median(v),np.median(v)],color='#172b40',lw=2)
axs[1].axhline(0,color='#777',lw=.7);axs[1].set_xticks(range(len(macro_clusters)),['C'+c+'\nn='+str(n) for c,n in zip(macro_clusters,pair_ns)],rotation=0);axs[1].set_ylabel('Patient-paired log2(CPM+1) difference\ncluster minus other macrophages');axs[1].set_title('B  SPP1 within macrophage catalogue (39 patients)');save(fig,'04_spp1_postcluster_check')
m=pd.read_csv(C/'RESOLUTION_STABILITY.tsv',sep='\t');fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained');ax.plot(m.resolution,m.minimum_between_seed_ARI,'o-',label='Minimum seed-pair ARI',color='#355f8d');ax.plot(m.resolution,m.subsample_regraph_ARI,'s-',label='80% regraph ARI',color='#3a9d8f');ax.axhline(.85,color='#355f8d',ls=':',lw=1);ax.axhline(.8,color='#3a9d8f',ls=':',lw=1);ax.axvline(z.uns['selected_resolution'],color='#a6514c',ls='--',label='Selected resolution');ax.set_xlabel('Leiden resolution');ax.set_ylabel('Adjusted Rand index');ax.set_ylim(0,1.02);ax.legend(frameon=False);ax.set_title('Purified candidate grid: global stability gates unmet');save(fig,'05_clustering_stability')
core=pd.read_csv(C/'CLUSTER_CORE_STABILITY.tsv',sep='\t',dtype={'cluster':str}).set_index('cluster').loc[clusters]
fig,ax=plt.subplots(figsize=(11,4.5),layout='constrained');ii=np.arange(len(clusters));ax.bar(ii-.18,core.minimum_seed_Jaccard,width=.36,label='Minimum seed Jaccard',color='#355f8d');ax.bar(ii+.18,core.subsample_Jaccard,width=.36,label='80% regraph Jaccard',color='#3a9d8f');ax.axhline(.75,color='#355f8d',ls=':',lw=1);ax.axhline(.70,color='#3a9d8f',ls=':',lw=1);ax.set_xticks(ii,['C'+c for c in clusters],rotation=45);ax.set_ylim(0,1.05);ax.set_ylabel('Best-matched cluster Jaccard');ax.legend(frameon=False);ax.set_title('Algorithmic stability of individual candidate cores');save(fig,'06_candidate_core_stability')
d=np.load(MM/'paired_patient_effects/cluster_4.npz');keygenes=['SPP1','GPNMB','CCL18','CYP27A1','LPL','APOE'];ev=pd.read_csv(MM/'C4_KEY_MARKER_EVIDENCE.tsv',sep='\t').set_index('gene');fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
for i,g in enumerate(keygenes):
 v=d['effects'][:,list(pd.read_csv(MM/'PSEUDOBULK_GENE_INDEX.tsv',sep='\t').gene).index(g)];ax.scatter(i+rng.uniform(-.17,.17,len(v)),v,s=25,color='#377e89',alpha=.65);ax.plot([i-.25,i+.25],[np.median(v)]*2,color='#233444',lw=2)
ax.axhline(0,color='#888',lw=.8);ax.set_xticks(range(len(keygenes)),[g+'\nq='+format(ev.loc[g,'BH_FDR_cluster_gene_family'],'.2g') for g in keygenes]);ax.set_ylabel('Patient-paired log2(CPM+1) difference\nC4 minus other macrophage candidates');ax.set_title('C4: GPNMB/CCL18-high, SPP1-enriched candidate | 26 paired patients');save(fig,'07_c4_multigene_patient_evidence')
(O/'FIGURE_FILES.json').write_text(json.dumps(files,indent=2));print('FIGURES',files)
