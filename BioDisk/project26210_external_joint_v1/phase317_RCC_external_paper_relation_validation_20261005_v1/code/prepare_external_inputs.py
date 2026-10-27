from pathlib import Path
import json,hashlib,subprocess,gc,sys
import numpy as np,pandas as pd,h5py
from scipy import sparse
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent
reg=pd.read_csv(R/'phase234_part4_replan_execution_20260921_v1/tables/INPUT_SECTION_REGISTER.tsv',sep='\t');reg.to_csv(P/'resources/EXTERNAL_SECTION_REGISTER.tsv',sep='\t',index=False)
records=[]
def sha(f):
 h=hashlib.sha256()
 with open(f,'rb') as x:
  for block in iter(lambda:x.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
for a in reg.itertuples():
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True)
 d=P/'inputs'/a.sample;d.mkdir(exist_ok=True)
 old=R/'phase289_state_source_boundary_process_20260929_v1/inputs'/f'{a.sample}_SPOT_PROGRAMS.tsv.gz';meta=pd.read_csv(old,sep='\t',float_precision='round_trip')
 with h5py.File(a.input_path,'r') as h:
  genes=np.asarray(read_elem(h['var/gene_symbol_original'])) if 'gene_symbol_original' in h['var'] else np.asarray(read_elem(h['var'][h['var'].attrs['_index']]))
  genes=genes.astype(str);ids=np.asarray(read_elem(h['obs'][h['obs'].attrs['_index']])).astype(str);assert pd.Index(ids).is_unique and meta.spot_id.is_unique
  ind=pd.Index(ids).get_indexer(meta.spot_id);assert (ind>=0).all()
  layer='layers/counts' if 'layers/counts' in h else 'X';x=read_elem(h[layer]).tocsr()[ind].astype(np.float64)
  source_total=np.asarray(read_elem(h['obs/total_counts']))[ind] if 'total_counts' in h['obs'] else None
 assert np.isfinite(x.data).all() and (x.data>=0).all() and np.all(x.data==np.floor(x.data))
 before=x.shape[1]
 if pd.Index(genes).has_duplicates:
  ug,inv=np.unique(genes,return_inverse=True);S=sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),inv)),shape=(len(genes),len(ug)));x=(x@S).tocsr();genes=ug
 x.eliminate_zeros();x.sort_indices();total=np.asarray(x.sum(1)).ravel();assert (total>0).all()
 x.data.astype('<f8').tofile(d/'values.bin');x.indices.astype('<i4').tofile(d/'indices.bin');x.indptr.astype('<i4').tofile(d/'indptr.bin')
 (d/'genes.txt').write_text('\n'.join(genes)+'\n');(d/'spots.txt').write_text('\n'.join(meta.spot_id)+'\n')
 cols=['spot_id','SPP1_share','total_macro_fraction','malignant_fraction','boundary_input_valid','inside_morphology_supported','inside_expression_sensitivity']
 z=meta[cols].copy();z=z.rename(columns={'SPP1_share':'mapped_SPP1_share'});z['library_total']=total;z['log1p_UMI']=np.log1p(total);z['detected_genes']=np.diff(x.indptr)
 for g in ['FAP','SPP1','CD3D','CD8A']:
  if g in genes:
   vv=x[:,np.flatnonzero(genes==g)[0]].toarray().ravel();z['g_'+g]=np.log1p(vv/total*1e4);z['det_'+g]=vv>0
  else:z['g_'+g]=np.nan;z['det_'+g]=np.nan
 z.to_csv(d/'METADATA.tsv',sep='\t',index=False)
 info={'sample':a.sample,'dataset':a.dataset,'patient':a.patient,'preparation':a.preparation,'input_path':a.input_path,'input_bytes':Path(a.input_path).stat().st_size,'input_sha256':sha(a.input_path),'source_gene_columns':before,'n_genes':len(genes),'n_spots':len(total),'nnz':len(x.data),'layer':layer,'count_dtype':'float64 little endian exact nonnegative integers','normalization':'sum of actual measured genes in count matrix','obs_total_mismatch_spots':int(np.sum(total!=source_total)) if source_total is not None else None,'old_metadata_sha256':sha(old),'genes_sha256':sha(d/'genes.txt'),'spots_sha256':sha(d/'spots.txt'),'metadata_sha256':sha(d/'METADATA.tsv')}
 (d/'INPUT_QC.json').write_text(json.dumps(info,indent=2));records.append(info);print('PREPARED',a.sample,len(total),len(genes),len(x.data),flush=True)
 del x,z,meta;gc.collect()
(P/'qc/EXTERNAL_INPUT_MANIFEST.json').write_text(json.dumps(records,indent=2));print('ALL_PREPARED',len(records),flush=True)
