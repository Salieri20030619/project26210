from pathlib import Path
import gzip,json,hashlib
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata
from scipy.cluster.hierarchy import linkage,fcluster,dendrogram
from scipy.sparse import diags
w=Path(__file__).resolve().parents[1];q=w/'qc';root=w.parent
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory/trajectory_v1')
old=root/'phase129_sections2_6_2_7_redraw_remote_v1'
genes=(b/'input/genes_v1.tsv').read_text().splitlines();cells=(b/'input/cells_v1.tsv').read_text().splitlines()
m=pd.read_csv(b/'input/cell_metadata_v1.tsv.gz',sep='\t',index_col=0).loc[cells]
tr=pd.read_csv(b/'source_data/trajectory_all_methods_cell_scores_v1.tsv.gz',sep='\t',index_col=0).loc[cells]
lib=pd.read_csv(old/'qc/RNA_EXTRACTION_ROWS.tsv',sep='\t',index_col=0).loc[cells]
enz=pd.read_csv(old/'data/ENZYME_LOG1P_CP10K.tsv.gz',sep='\t',index_col=0).loc[cells]
assert len(genes)==2500 and len(cells)==11950 and len(set(genes))==2500
with gzip.open(b/'input/raw_counts_genes_by_cells_v1.mtx.gz','rb') as f:raw=mmread(f).tocsr().T.tocsr()
assert raw.shape==(11950,2500) and (lib.total_retained_UMI>0).all()
X=(diags(1e4/lib.total_retained_UMI.values)@raw).tocsr();X.data=np.log1p(X.data)
X=X.toarray();extra=[g for g in enz.columns if g not in genes]
for gene in enz.columns:
 if gene in genes:assert np.max(abs(X[:,genes.index(gene)]-enz[gene].values))<2e-6
if extra:X=np.column_stack([X,enz[extra].values])
genes=genes+extra;assert len(set(genes))==len(genes)
patient_ids=[pid for pid,sub in m.groupby('curated_patient_id') if len(sub)>=30];assert len(patient_ids)==20
keep=m.curated_patient_id.isin(patient_ids);assert keep.sum()==11921
methods=['Monocle3','Slingshot','scTour'];grid=np.linspace(.025,.975,200);bin_rows=[];states=[];curves={};patient_curves={};max_recheck=0
for method in methods:
 fits=[]
 for pid in patient_ids:
  ix=np.flatnonzero(m.curated_patient_id.values==pid);t=tr.iloc[ix][method+'_pseudotime_patient_v1'].values
  ranks=rankdata(t)/len(t);ib=np.clip(np.ceil(ranks*20),1,20).astype(int);ub=np.unique(ib)
  nn=np.array([(ib==k).sum() for k in ub]);bx=(ub-.5)/20
  means=np.array([X[ix[ib==k]].mean(axis=0) for k in ub])
  weights=np.exp(-.5*((grid[:,None]-bx[None,:])/.075)**2)*nn[None,:];weights/=weights.sum(axis=1,keepdims=True)
  fit=weights@means;fits.append(fit)
  # Direct sum formula, independently checked at fixed grid positions across every gene.
  for j in [0,50,100,150,199]:
   ww=np.exp(-((grid[j]-bx)**2)/(2*.075**2))*nn;direct=np.sum(ww[:,None]*means,axis=0)/ww.sum()
   max_recheck=max(max_recheck,float(np.max(abs(direct-fit[j]))))
  for k,n in zip(ub,nn):bin_rows.append({'method':method,'patient':pid,'bin':int(k),'n_cells':int(n)})
  for st in ['Base','Hypoxia','Cycle','TYMP_state']:
   v=(tr.iloc[ix].malignant_state_v1==st).to_numpy(float);vals=np.array([v[ib==k].mean() for k in ub]);sf=weights@vals
   for x,y in zip(grid,sf):states.append({'method':method,'patient':pid,'rank':x,'state':st,'fraction':y})
 fits=np.stack(fits,axis=0)
 patient_curves[method]=fits.astype('float32');curves[method]=fits.mean(axis=0).T
 print(method,'patients',len(fits),'genes',len(genes),flush=True)
assert max_recheck<1e-11
# Fixed gene registry; only numerical-flat Monocle3 display rows are omitted, and explicitly registered.
span=np.ptp(curves['Monocle3'],axis=1);included=span>1e-10
registry=pd.DataFrame({'gene':genes,'source':['frozen_trajectory_HVG' if g not in extra else 'prespecified_enzyme_addition' for g in genes],'Monocle3_fitted_range':span,'included':included})
registry.to_csv(q/'GENE_INCLUSION_REGISTRY.tsv',sep='\t',index=False)
selected=np.array(genes)[included];assert set(enz.columns)<=set(selected)
zs={};row_sd={}
for method in methods:
 c=curves[method][included];sd=c.std(axis=1,ddof=0);row_sd[method]=sd
 z=np.divide(c-c.mean(axis=1,keepdims=True),sd[:,None],out=np.zeros_like(c),where=sd[:,None]>1e-10)
 assert np.isfinite(z).all();zs[method]=z
# Ward clustering of the original Monocle3 fitted profiles; same row order retained for all methods.
tree=linkage(zs['Monocle3'],method='ward',metric='euclidean')
# Rotate existing branches using weighted peak positions; this changes only leaf orientation.
peaks=np.argmax(zs['Monocle3'],axis=1).astype(float);scores={i:peaks[i] for i in range(len(selected))};sizes={i:1 for i in range(len(selected))}
for i,row in enumerate(tree):
 a,bn=map(int,row[:2])
 if scores[a]>scores[bn]:tree[i,0],tree[i,1]=bn,a
 scores[len(selected)+i]=(scores[a]*sizes[a]+scores[bn]*sizes[bn])/(sizes[a]+sizes[bn]);sizes[len(selected)+i]=sizes[a]+sizes[bn]
leaves=np.array(dendrogram(tree,no_plot=True)['leaves']);rawcl=fcluster(tree,4,criterion='maxclust')[leaves];ids=list(dict.fromkeys(rawcl));clusters=np.array([ids.index(c)+1 for c in rawcl])
for c in np.unique(clusters):assert np.all(np.diff(np.flatnonzero(clusters==c))==1)
np.savez_compressed(w/'data/GRADIENT_HEATMAP_VALUES.npz',genes=selected,all_genes=np.array(genes),grid=grid,tree=tree,leaves=leaves,clusters=clusters,**{method:zs[method] for method in methods},**{method+'_mean_log1p':curves[method][included] for method in methods})
np.savez_compressed(w/'data/PATIENT_GENE_CURVES.npz',genes=np.array(genes),patients=np.array(patient_ids),grid=grid,**patient_curves)
sf=pd.DataFrame(states);assert np.allclose(sf.groupby(['method','patient','rank']).fraction.sum(),1)
sf.to_csv(w/'data/PATIENT_STATE_COMPOSITION.tsv.gz',sep='\t',index=False)
pd.DataFrame(bin_rows).to_csv(q/'PATIENT_BIN_COUNTS.tsv',sep='\t',index=False)
pd.DataFrame({'gene':selected[leaves],'cluster':['C'+str(x) for x in clusters]}).to_csv(q/'DISPLAY_ROW_ORDER.tsv',sep='\t',index=False)
validation={'status':'PASS','source_cells':11950,'display_cells':int(keep.sum()),'independent_patients':20,'excluded_small_patient_cells':int((~keep).sum()),'input_HVGs':2500,'added_prespecified_enzymes':extra,'input_genes':len(genes),'display_genes':len(selected),'flat_rows_excluded':registry.loc[~included,'gene'].tolist(),'cluster_counts':{str(c):int(np.sum(clusters==c)) for c in np.unique(clusters)},'fixed_row_order_across_methods':True,'all_12_enzymes_included':True,'smoothing':'20 within-patient rank bins; Gaussian h=.075; mean across 20 equally weighted patients; 200 display rank positions','independent_fit_max_error':max_recheck,'state_fractions_sum_to_one':True,'new_inferential_tests':False,'models_rerun':False}
(q/'VALIDATION.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation),flush=True)
