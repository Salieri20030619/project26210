from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad,scipy.io as sio
P=Path(__file__).resolve().parents[1];S=P.parent/'phase305_larger_immune_sequential_20261002_v1'
a=ad.read_h5ad(P/'inputs/OBRADOVIC_macrophage_raw.h5ad')
symbol_n=np.asarray((a.X>0).sum(axis=1)).ravel();delta=a.obs.n_genes.to_numpy()-symbol_n
rows=[]
samples=json.loads((S/'metadata/OBRADOVIC_AUTHOR_SAMPLE_MAP.json').read_text())
lookup={Path(s['path']).parts[2]:s for s in samples if s['tissue']=='Tumor'}
for lib,g in a.obs.groupby('library',observed=True):
 ix=a.obs.library.eq(lib).to_numpy();diff=delta[ix]
 rec={'library':lib,'cells':int(ix.sum()),'cells_with_difference':int((diff!=0).sum()),'min_difference':int(diff.min()),'max_difference':int(diff.max())}
 if (diff!=0).any():
  d=S/'inputs/OBRADOVIC_raw'/lookup[str(lib)]['path']
  ft=pd.read_csv(next(iter(list(d.glob('genes.tsv*'))+list(d.glob('features.tsv*')))),sep='\t',header=None)
  bc=pd.read_csv(next(d.glob('barcodes.tsv*')),sep='\t',header=None).iloc[:,0].astype(str)
  rr=sio.mmread(next(d.glob('matrix.mtx*'))).T.tocsr()
  if ft.shape[1]>=3:
   keep=ft.iloc[:,2].eq('Gene Expression').to_numpy();rr=rr[:,keep].tocsr();ft=ft.loc[keep]
  jj=pd.Index(str(lib)+':'+bc).get_indexer(g.index);assert (jj>=0).all()
  rr=rr[jj];ng=np.diff(rr.indptr)
  assert np.array_equal(ng,g.n_genes.to_numpy())
  gene=ft.iloc[:,1].astype(str).to_numpy();un,inv=np.unique(gene,return_inverse=True)
  import scipy.sparse as sp
  collapsed=rr@sp.csr_matrix((np.ones(len(inv)),(np.arange(len(inv)),inv)),shape=(len(inv),len(un)))
  assert np.array_equal(np.asarray((collapsed>0).sum(axis=1)).ravel(),symbol_n[ix])
  rec.update(original_features=len(gene),unique_symbols=len(un),duplicate_symbol_features=len(gene)-len(un),explicit_zero_entries=int((rr.data==0).sum()),original_premerge_complexity_exact=True,postmerge_complexity_exact=True)
 rows.append(rec);print(json.dumps(rec),flush=True)
pd.DataFrame(rows).to_csv(P/'qc/COMPLEXITY_DEFINITION_BY_LIBRARY.tsv',sep='\t',index=False)
q={'status':'PASS','affected_cells':int((delta!=0).sum()),'max_difference':int(delta.max()),'definition':'Inherited n_genes is original pre-symbol-collapse feature complexity, while matrix columns are merged gene symbols. Original feature counts independently reproduced for every affected library. Existing model intentionally retains original QC covariate; no source or model edits.'}
(P/'qc/COMPLEXITY_DEFINITION_VALIDATED.json').write_text(json.dumps(q,indent=2))
print(json.dumps(q),flush=True)
