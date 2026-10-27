from common import *
import h5py
from anndata.io import read_elem
m=read(P/'inputs/STROMAL_INPUT_METADATA.tsv.gz')
with h5py.File(H,'r') as h:
 ids=h['obs/cell_id'].asstr()[:];ix=pd.Index(ids).get_indexer(m.cell_id);assert (ix>=0).all()
 d=pd.DataFrame({'cell_id':m.cell_id})
 for k in ['malignancy_status_v1','malignant_conservative_v1','cnv_assessed_v1','cnv_analysis_scope_v1']:d[k]=np.asarray(read_elem(h['obs'][k]))[ix]
d.to_csv(P/'inputs/CANDIDATE_CNV_FLAGS.tsv.gz',sep='\t',index=False)
