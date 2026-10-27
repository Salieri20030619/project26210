from pathlib import Path
import json,collections,datetime,hashlib
import numpy as np,pandas as pd
from scipy import sparse
P=Path(__file__).resolve().parents[1];R=P/'resources/region1';R.mkdir(exist_ok=True)
assert (P/'logs/acquire_region1.exit').read_text().strip()=='0'
src=P/'downloads/tumor1_region1_transcripts.csv.gz';cm={};A=collections.Counter();N=collections.Counter();xyz=collections.defaultdict(lambda:np.zeros(4));features=collections.Counter();audit=collections.Counter()
for d in pd.read_csv(src,chunksize=300000,dtype={'cell_id':str,'feature_name':str}):
 audit['raw_transcripts']+=len(d);features.update(d.feature_name.value_counts().to_dict())
 assigned=~d.cell_id.isin(['-1','UNASSIGNED','nan','None','']);quality=d.qv.ge(20);control=d.feature_name.str.match(r'(?i)^(NegControl|NegativeControl|Blank|Unassigned|antisense)')
 audit['unassigned']+=int((~assigned).sum());audit['qv_below20']+=int((~quality).sum());audit['control']+=int(control.sum())
 k=d[assigned&quality&~control].copy();audit['assigned_qv20']+=len(k)
 assert np.isfinite(k[['x_location','y_location','z_location']].to_numpy()).all()
 for c in k.cell_id.unique():
  if c not in cm:cm[c]=len(cm)
 A.update(k.groupby(['cell_id','feature_name'],sort=False).size().to_dict())
 n=k[k.overlaps_nucleus.eq(1)];N.update(n.groupby(['cell_id','feature_name'],sort=False).size().to_dict());audit['nuclear_qv20']+=len(n)
 s=n.groupby('cell_id')[['x_location','y_location','z_location']].sum();ct=n.groupby('cell_id').size()
 for c,row in s.iterrows():xyz[c]+=np.array([row.x_location,row.y_location,row.z_location,ct[c]])
genes=sorted(features);gm={g:i for i,g in enumerate(genes)}
def build(c):
 return sparse.coo_matrix((np.fromiter(c.values(),np.int64),(np.fromiter((cm[x] for x,g in c),np.int64),np.fromiter((gm[g] for x,g in c),np.int64))),shape=(len(cm),len(genes))).tocsr()
a=build(A);n=build(N);assert int(a.sum())==audit['assigned_qv20'];assert int(n.sum())==audit['nuclear_qv20'];assert (a-n).min()>=0
sparse.save_npz(R/'assigned_counts.npz',a);sparse.save_npz(R/'nuclear_counts.npz',n)
rows=[]
for c,i in cm.items():
 v=xyz[c];coord=v[:3]/v[3] if v[3] else [np.nan]*3
 rows.append({'cell_id':c,'matrix_row':i,'assigned_RNA':int(a[i].sum()),'assigned_genes':int(a[i].getnnz()),'nuclear_RNA':int(n[i].sum()),'nuclear_genes':int(n[i].getnnz()),'nucleus_x_proxy_um':coord[0],'nucleus_y_proxy_um':coord[1],'nucleus_z_proxy_um':coord[2]})
cells=pd.DataFrame(rows);cells['technical_qc_pass']=(cells.assigned_RNA>=20)&(cells.assigned_genes>=5)&(cells.nuclear_RNA>=5);cells.to_csv(R/'CELL_METADATA.tsv',sep='\t',index=False)
pd.DataFrame({'gene':genes,'matrix_col':np.arange(len(genes)),'raw_transcripts':[features[g] for g in genes],'assigned_transcripts':np.asarray(a.sum(axis=0)).ravel(),'nuclear_transcripts':np.asarray(n.sum(axis=0)).ravel()}).to_csv(R/'FEATURE_AUDIT.tsv',sep='\t',index=False)
rec={'status':'RECONSTRUCTED','time':datetime.datetime.now().astimezone().isoformat(),'region':'tumor1_region1','raw_file':str(src),'cells_with_assigned_counts':len(cells),'QC_cells':int(cells.technical_qc_pass.sum()),'features_observed':len(genes),'counts':dict(audit),'gzip_full_read':'PASS','raw_download_validation':'qc/REGION1_DOWNLOAD_VALIDATION.json'}
(P/'qc/REGION1_RECONSTRUCTION.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec,indent=2),flush=True)
