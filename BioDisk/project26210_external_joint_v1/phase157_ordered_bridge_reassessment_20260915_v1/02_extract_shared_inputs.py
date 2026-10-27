from pathlib import Path
import json,subprocess,datetime
import numpy as np,pandas as pd,h5py
from scipy import sparse
from anndata.io import read_elem
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase157_ordered_bridge_reassessment_20260915_v1';D=O/'02_inputs';D.mkdir(exist_ok=False)
assert json.loads((O/'01_public_TYMP_MTA/STATUS.json').read_text())['status']=='COMPLETE_PUBLIC_MEASURED_MTA_REASSESSMENT'
with (D/'memory.log').open('w') as log:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=log,stderr=subprocess.STDOUT,check=True)
plan={'sequence':['public_measured_TYMP_MTA','MTA_proxies_to_new_C4','other_public_X'],'scope':'Frozen 16 matched tumour-library patients for cross-compartment links; all current mature macrophages for C4 versus others response compatibility. Client results not used.','MTA_metrics':['signed 89-gene public fibroblast MTA response transcript proxy','ODC1/SRM/SMS production transcript proxy','MTAP clearance RNA; separate indicator, not measured MTA'],'other_X':['single-gene Secreted Signaling ligands from frozen public DB','50 fixed Hallmark programmes excluding TYMP/SPP1','13 predefined metabolism/transport genes from Phase155'],'QC':'RNA genes need >=10 pseudobulk counts in >=12/16 malignant libraries; programme >=10 covered genes; fixed filters independent of C4 association. Multi-subunit ligands remain in catalogue but not primary single-gene screen.','stats':'Patient only; raw and within-study associations; one-X conditional C4 association and standardized rank indirect effect; separate BH by family and endpoint; no FDR action gate. Preserve all negative and opposing results.','figures':False}
(O/'ORDERED_PLAN.json').write_text(json.dumps(plan,indent=2))
H=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
m=pd.read_csv(R/'08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str},low_memory=False);m=m[~m.candidate.isin(['1','Monocyte'])]
with h5py.File(H,'r') as f:
 ids=f['obs']['cell_id'].asstr()[:];genes=f['var']['gene'].asstr()[:]
 cols=['study','curated_patient_id','sample_id','source_tissue','total_counts_raw','eligible_for_subclone_selection_v1']
 obs=pd.DataFrame({k:np.asarray(read_elem(f['obs'][k])) for k in cols},index=ids)
mal=obs[obs.eligible_for_subclone_selection_v1==True].copy();mal['role']='Malignant'
mac=obs.loc[m.index].copy();mac['role']=np.where(m.candidate=='4','C4','OtherMacrophage')
meta=pd.concat([mal,mac]).loc[lambda z:~z.index.duplicated()].copy();meta=meta.loc[obs.index[obs.index.isin(meta.index)]]
keys=['study','curated_patient_id','sample_id','source_tissue','role']
meta['group_key']=meta[keys].astype(str).agg('|'.join,axis=1)
labels=sorted(meta.group_key.unique());lookup={x:i for i,x in enumerate(labels)}
meta['group_index']=meta.group_key.map(lookup)
group=meta.groupby('group_index').agg(study=('study','first'),patient=('curated_patient_id','first'),sample_id=('sample_id','first'),source_tissue=('source_tissue','first'),role=('role','first'),cells=('study','size')).reindex(range(len(labels))).reset_index()
counts=np.zeros((len(labels),len(genes)),dtype=np.float64);detected=np.zeros_like(counts,dtype=np.int32)
mask=obs.index.isin(meta.index);n=len(ids)
sig=pd.read_csv('/mnt/BioDisk/wangmingyuan/project26210_phase77_mta_perturbation_signature_v1/training_v1/phase77_weighted_signatures_v1.tsv',sep='\t');sig=sig[sig.signature_id=='DIRECT_MTA_GSE249082'].copy()
sig=sig[~sig.human_symbol.isin(['TYMP','SPP1'])]
sig=sig[sig.human_symbol.isin(genes)]
si=pd.Index(genes).get_indexer(sig.human_symbol)
sig_expr=np.zeros((len(meta),len(si)),dtype=np.float64)
with h5py.File(H,'r',rdcc_nbytes=128*1024*1024) as f:
 X=f['layers/counts'];ptr=X['indptr'][:]
 for st in range(0,n,8192):
  en=min(st+8192,n);sel=mask[st:en]
  if not sel.any():continue
  low,high=int(ptr[st]),int(ptr[en]);mat=sparse.csr_matrix((X['data'][low:high],X['indices'][low:high],ptr[st:en+1]-low),shape=(en-st,len(genes)))[sel]
  cell_ids=ids[st:en][sel];lib=np.asarray(mat.sum(axis=1,dtype=np.float64)).ravel()
  assert np.allclose(lib,obs.loc[cell_ids,'total_counts_raw'],rtol=0,atol=1e-5)
  codes=meta.loc[cell_ids,'group_index'].to_numpy();G=sparse.csr_matrix((np.ones(len(codes)),(codes,np.arange(len(codes)))),shape=(len(labels),len(codes)))
  counts+=(G@mat).toarray();detected+=(G@(mat>0)).toarray().astype(np.int32)
  loc=meta.index.get_indexer(cell_ids);sig_expr[loc]=np.log1p(mat[:,si].toarray()*10000/lib[:,None])
assert (counts>=0).all() and (detected<=group.cells.to_numpy()[:,None]).all()
np.savez_compressed(D/'GROUP_PSEUDOBULK.npz',counts=counts,detected=detected,genes=genes)
group.to_csv(D/'GROUP_INDEX.tsv',sep='\t',index=False)
meta[keys+['group_index']].to_csv(D/'CELL_IDENTITY.tsv.gz',sep='\t',index_label='cell_id')
sig.to_csv(D/'MTA_SIGNATURE_GENE_REGISTRY.tsv',sep='\t',index=False)
# Public response score: within study and compartment gene z-scores, original signed weights.
weights=sig.weight.to_numpy()*np.where(sig.direction=='UP',1,-1);den=float(sig.weight.sum())
scores=np.full(len(meta),np.nan);comp=np.where(meta.role=='Malignant','Malignant','Macrophage')
for study in meta.study.unique():
 for compartment in ['Malignant','Macrophage']:
  ix=np.flatnonzero((meta.study.to_numpy()==study)&(comp==compartment))
  if len(ix)<2:continue
  x=sig_expr[ix];sd=x.std(axis=0,ddof=1);z=np.divide(x-x.mean(axis=0),sd,out=np.zeros_like(x),where=sd>0)
  scores[ix]=z@weights/den
meta['MTA_response_transcript_score']=scores
meta[keys+['MTA_response_transcript_score']].to_csv(D/'CELL_MTA_RESPONSE_SCORES.tsv.gz',sep='\t',index_label='cell_id')
status={'status':'INPUTS_COMPLETE','cells':len(meta),'groups':len(group),'genes':len(genes),'MTA_signature_genes':len(sig),'MTA_up':int((sig.direction=='UP').sum()),'MTA_down':int((sig.direction=='DOWN').sum()),'raw_library_sums_checked':True,'figures_created':False}
(D/'STATUS.json').write_text(json.dumps(status,indent=2));print(json.dumps(status))
