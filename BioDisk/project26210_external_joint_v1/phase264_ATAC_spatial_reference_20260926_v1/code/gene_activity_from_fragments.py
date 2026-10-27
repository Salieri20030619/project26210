from pathlib import Path
import json,time
import numpy as np,pandas as pd
from scipy import sparse
from ncls import NCLS
P=Path(__file__).resolve().parents[1];D=P/'downloads/RCC100';starttime=time.time()
assert (P/'DOWNLOAD_EXIT_STATUS').read_text().strip()=='0'
meta=pd.read_csv(D/'singlecell.csv');meta=meta[meta.is__cell_barcode.eq(1)].copy().reset_index(drop=True)
assert meta.barcode.is_unique
barcodes=meta.barcode.to_numpy();lookup=dict(zip(barcodes,range(len(meta))))
regions=pd.read_csv(P/'reference/protein_coding_gene_activity_regions.tsv',sep='\t')
cent=pd.read_csv(P/'reference/broad_lineage_RNA_centroids.tsv',sep='\t',index_col=0)
regions=regions[regions.gene.isin(cent.index)]
ambiguous=regions.groupby('gene').chrom.nunique();bad=ambiguous[ambiguous>1].index.tolist()
regions=regions[~regions.gene.isin(bad)].groupby(['gene','chrom'],as_index=False).agg(start=('start','min'),end=('end','max'))
genes=sorted(regions.gene.unique());gi={g:i for i,g in enumerate(genes)};trees={}
for chrom,t in regions.groupby('chrom'):
 trees[chrom]=NCLS(t.start.to_numpy(np.int64),t.end.to_numpy(np.int64),np.array([gi[g] for g in t.gene],dtype=np.int64))
rows=[];cols=[];frgcount=np.zeros(len(meta),dtype=np.int64);nread=0;nused=0
# Independent brute-force overlap check on the first 1000 called-cell fragments.
checked=0;maxcheck=1000
for chunk in pd.read_csv(D/'fragments.tsv.gz',sep='\t',header=None,names=['chrom','start','end','barcode','support'],usecols=[0,1,2,3],dtype={'chrom':str,'start':np.int64,'end':np.int64,'barcode':str},compression='gzip',chunksize=500000,comment='#'):
 nread+=len(chunk);idx=chunk.barcode.map(lookup)
 chunk=chunk[idx.notna()].copy();chunk['cell_idx']=idx[idx.notna()].to_numpy(dtype=np.int64)
 nused+=len(chunk);frgcount+=np.bincount(chunk.cell_idx,minlength=len(meta))
 for chrom,z in chunk.groupby('chrom',sort=False):
  if chrom not in trees:continue
  z=z.reset_index(drop=True);rid,gid=trees[chrom].all_overlaps_both(z.start.to_numpy(np.int64),z.end.to_numpy(np.int64),np.arange(len(z),dtype=np.int64))
  if checked<maxcheck:
   sub=regions[regions.chrom==chrom]
   n=min(maxcheck-checked,len(z))
   for j in range(n):
    expected=set(sub.loc[(sub.start<z.end.iloc[j])&(sub.end>z.start.iloc[j]),'gene'].map(gi))
    observed=set(gid[rid==j].tolist());assert expected==observed,(chrom,j)
   checked+=n
  if len(gid):
   rows.append(gid.astype(np.int32));cols.append(z.cell_idx.to_numpy()[rid].astype(np.int32))
 if nread%5000000==0:print('FRAGMENTS_READ',nread,'CALLED_CELL_FRAGMENTS',nused,'SECONDS',round(time.time()-starttime,1),flush=True)
r=np.concatenate(rows);c=np.concatenate(cols)
act=sparse.coo_matrix((np.ones(len(r),dtype=np.float32),(r,c)),shape=(len(genes),len(meta))).tocsr()
assert np.isfinite(act.data).all() and np.equal(act.data,np.floor(act.data)).all()
sparse.save_npz(P/'tables/RCC100_gene_activity.npz',act)
pd.Series(genes).to_csv(P/'tables/RCC100_gene_activity_genes.txt',index=False,header=False)
pd.Series(barcodes).to_csv(P/'tables/RCC100_gene_activity_barcodes.txt',index=False,header=False)
meta['observed_unique_fragments']=frgcount;meta.to_csv(P/'tables/RCC100_CALLED_CELL_METADATA.tsv',sep='\t',index=False)
peakdf=pd.read_csv(D/'peaks.bed',sep='\t',header=None,names=['chrom','start','end'],usecols=[0,1,2])
peakdf[peakdf.chrom.eq('chr1')].to_csv(P/'tables/RCC100_chr1_all_peaks.bed',sep='\t',header=False,index=False)
result=dict(status='PASS',source=str(D/'fragments.tsv.gz'),read_all_gzip_to_end=True,gzip_crc_verified=True,called_cells=len(meta),genes=len(genes),all_fragment_rows=nread,called_cell_fragments=nused,activity_nonzero=int(act.nnz),direct_overlap_checks=checked,ambiguous_gene_symbols_excluded=bad,elapsed_seconds=time.time()-starttime,peak_count=len(peakdf),chr1_peaks=int(peakdf.chrom.eq('chr1').sum()))
(P/'qc/GENE_ACTIVITY_BUILD.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)

