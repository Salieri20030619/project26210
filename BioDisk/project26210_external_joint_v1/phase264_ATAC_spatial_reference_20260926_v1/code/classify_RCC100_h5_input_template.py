from pathlib import Path
import json,re,hashlib,traceback
import numpy as np,pandas as pd,h5py
from scipy import sparse,io
P=Path(__file__).resolve().parents[1];D=P/'downloads/RCC100'
assert (P/'DOWNLOAD_EXIT_STATUS').read_text().strip()=='0'
plan=json.load(open(P/'metadata/LINEAGE_TYPING_PLAN.json'))
h5s=sorted(D.rglob('*filtered*peak*matrix*.h5'))
if not h5s:raise RuntimeError('No complete filtered ATAC peak-by-cell matrix found; do not infer labels from fragment filenames.')
assert len(h5s)==1,h5s
path=h5s[0]
with h5py.File(path,'r') as f:
 g=f['matrix'];shape=tuple(g['shape'][:])
 peaks=g['features/name'].asstr()[:] if 'features/name' in g else g['features/id'].asstr()[:]
 barcodes=g['barcodes'].asstr()[:]
 assert shape==(len(peaks),len(barcodes)) and shape[1]<100000
 counts=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=shape)
assert pd.Index(barcodes).is_unique and pd.Index(peaks).is_unique
assert np.isfinite(counts.data).all() and counts.data.min()>=0 and np.equal(counts.data,np.floor(counts.data)).all()
parsed=[re.fullmatch(r'(chr[^:]+):(\d+)-(\d+)',x) or re.fullmatch(r'(chr[^-]+)-(\d+)-(\d+)',x) for x in peaks]
assert all(parsed),'Unrecognized peak coordinate format'
peakdf=pd.DataFrame([(m[1],int(m[2]),int(m[3]),i) for i,m in enumerate(parsed)],columns=['chrom','start','end','peak_index'])
regions=pd.read_csv(P/'reference/protein_coding_gene_activity_regions.tsv',sep='\t')
centroid=pd.read_csv(P/'reference/broad_lineage_RNA_centroids.tsv',sep='\t',index_col=0)
genes=sorted(set(regions.gene)&set(centroid.index));gindex={g:i for i,g in enumerate(genes)}
ri=[];ci=[]
for chrom,pr in peakdf.groupby('chrom'):
 pr=pr.sort_values('start');starts=pr.start.to_numpy();ends=pr.end.to_numpy();indices=pr.peak_index.to_numpy()
 maxlen=int((ends-starts).max())
 for row in regions[(regions.chrom==chrom)&regions.gene.isin(gindex)].itertuples():
  lo=np.searchsorted(starts,row.start-maxlen,side='left');hi=np.searchsorted(starts,row.end,side='left')
  take=indices[lo:hi][ends[lo:hi]>row.start]
  ri.extend([gindex[row.gene]]*len(take));ci.extend(take.tolist())
overlap=sparse.coo_matrix((np.ones(len(ri)),(ri,ci)),shape=(len(genes),counts.shape[0])).tocsr()
overlap.data[:]=1
activity=(overlap@counts).tocsr()
total=np.asarray(activity.sum(axis=0)).ravel()
n_genes=np.asarray((activity>0).sum(axis=0)).ravel()
peakfrags=np.asarray(counts.sum(axis=0)).ravel()
qc=(peakfrags>=1000)&(n_genes>=500)&(total>0)
assert qc.sum()>=100,'Fewer than 100 called nuclei pass basic gene-activity eligibility'
features=pd.read_csv(P/'reference/broad_lineage_training_features.tsv',sep='\t')
feats=sorted(set(features.gene)&set(genes));assert len(feats)>=240
ix=[gindex[g] for g in feats]
q=np.log1p(activity[ix].toarray()/np.maximum(total,1)[None,:]*1e4)
q-=q[:,qc].mean(axis=1,keepdims=True)
q/=np.maximum(np.linalg.norm(q,axis=0),1e-12)
ref=np.log1p(centroid.loc[feats].to_numpy())
ref-=ref.mean(axis=1,keepdims=True);ref/=np.maximum(np.linalg.norm(ref,axis=0),1e-12)
score=q.T@ref;names=centroid.columns.to_numpy()
best=score.argmax(axis=1);margin=np.sort(score,axis=1)[:,-1]-np.sort(score,axis=1)[:,-2]
allhold=sorted(set(sum(plan['canonical_holdout'].values(),[]))&set(genes))
hraw=activity[[gindex[g] for g in allhold]].toarray()
h=np.log1p(hraw/np.maximum(total,1)[None,:]*1e4)
h=(h-h[:,qc].mean(axis=1,keepdims=True))/np.maximum(h[:,qc].std(axis=1,keepdims=True),1e-6)
hs=[];nz=[]
for name in names:
 inds=[allhold.index(g) for g in plan['canonical_holdout'][name] if g in allhold]
 assert len(inds)>=3,name
 hs.append(h[inds].mean(axis=0));nz.append((hraw[inds]>0).sum(axis=0))
hs=np.array(hs).T;nz=np.array(nz).T
hold=hs.argmax(axis=1);support=nz[np.arange(len(best)),best]
qualified=qc&(margin>=.05)&(best==hold)&(support>=2)
labels=np.where(qualified,names[best],'Unassigned')
out=pd.DataFrame(dict(barcode=barcodes,peak_fragments=peakfrags,gene_activity_features=n_genes,basic_qc=qc,reference_label=names[best],cosine=score.max(axis=1),margin=margin,heldout_label=names[hold],heldout_nonzero=support,qualified=qualified,lineage=labels))
out.to_csv(P/'tables/RCC100_ATAC_LINEAGE_ASSIGNMENTS.tsv',sep='\t',index=False)
for i,n in enumerate(names):out[n+'_score']=score[:,i]
out.to_csv(P/'tables/RCC100_ATAC_CLASSIFIER_DETAIL.tsv.gz',sep='\t',index=False)
summary=out.groupby('lineage').agg(n_cells=('barcode','size'),median_peak_fragments=('peak_fragments','median'),median_cosine=('cosine','median'),median_margin=('margin','median')).reset_index()
summary['n_eligible_for_pseudobulk']=summary.n_cells.ge(100)&summary.lineage.ne('Unassigned')
summary.to_csv(P/'tables/RCC100_LINEAGE_QC.tsv',sep='\t',index=False)
peakdf[peakdf.chrom=='chr1'][['chrom','start','end']].to_csv(P/'tables/RCC100_chr1_all_peaks.bed',sep='\t',header=False,index=False)
sparse.save_npz(P/'tables/RCC100_gene_activity.npz',activity)
pd.Series(genes).to_csv(P/'tables/RCC100_gene_activity_genes.txt',index=False,header=False)
pd.Series(barcodes).to_csv(P/'tables/RCC100_gene_activity_barcodes.txt',index=False,header=False)
status=dict(status='LINEAGE_TYPING_COMPLETE',matrix_source=str(path),matrix_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),peaks=counts.shape[0],called_nuclei=counts.shape[1],basic_eligible=int(qc.sum()),confident_assignments=int(qualified.sum()),typing_features=len(feats),heldout_genes=len(allhold),qualified_lineages=summary.loc[summary.n_eligible_for_pseudobulk,'lineage'].tolist())
(P/'qc/LINEAGE_TYPING_STATUS.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status),flush=True);print(summary.to_string(index=False),flush=True)

