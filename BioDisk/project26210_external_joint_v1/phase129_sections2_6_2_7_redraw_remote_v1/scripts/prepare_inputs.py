from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd,h5py
from scipy.sparse import csr_matrix
w=Path(__file__).resolve().parents[1]
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
t=b/'09_metabolism_trajectory'
genes=['TYMP','TK1','SLC29A1','SLC29A2','ODC1','SRM','SMS','MTAP','CHKA','PCYT1A','PCYT2','CEPT1']
m=pd.read_csv(t/'trajectory_v1/input/cell_metadata_v1.tsv.gz',sep='\t').set_index('exact_cell_id')
p=b/'03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad'
with h5py.File(p) as f:
 def decode(a):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in a])
 obs=decode(f['obs'][f['obs'].attrs['_index']][:]);var=decode(f['var'][f['var'].attrs['_index']][:])
 oi=pd.Index(obs).get_indexer(m.index);gi=pd.Index(var).get_indexer(genes)
 assert (oi>=0).all() and (gi>=0).all() and len(m)==11950
 g=f['layers/counts'];assert g.attrs['encoding-type']=='csr_matrix'
 ptr=g['indptr'][:];raw=np.zeros((len(m),len(genes)));tot=np.zeros(len(m))
 # Read exact source rows only; library sizes include all retained genes.
 for j,i in enumerate(oi):
  a,z=ptr[i:i+2];cols=g['indices'][a:z];vals=g['data'][a:z]
  tot[j]=vals.sum();lookup=dict(zip(cols,vals));raw[j]=[lookup.get(k,0) for k in gi]
  if j%3000==0:print('source_rows',j,flush=True)
 expr=np.log1p(raw/tot[:,None]*1e4)
 assert np.array_equal(raw[:,0],m.TYMP_raw_count_v1.values)
 err=float(np.max(np.abs(expr[:,0]-m.TYMP_log1p_CP10K_v1.values)))
 assert err<2e-6,err
 pd.DataFrame(expr,index=m.index,columns=genes).to_csv(w/'data/ENZYME_LOG1P_CP10K.tsv.gz',sep='\t')
 pd.DataFrame(raw,index=m.index,columns=genes).to_csv(w/'data/ENZYME_RAW_COUNTS.tsv.gz',sep='\t')
 pd.DataFrame({'exact_cell_id':m.index,'source_row':oi,'total_retained_UMI':tot}).to_csv(w/'qc/RNA_EXTRACTION_ROWS.tsv',sep='\t',index=False)
 (w/'qc/RNA_EXTRACTION_VALIDATION.json').write_text(json.dumps({'status':'PASS','cells':len(m),'genes':genes,'TYMP_max_error':err,'source':str(p),'source_bytes':p.stat().st_size,'source_mtime_ns':p.stat().st_mtime_ns},indent=2))
print('INPUT_EXTRACTION_PASS',flush=True)
