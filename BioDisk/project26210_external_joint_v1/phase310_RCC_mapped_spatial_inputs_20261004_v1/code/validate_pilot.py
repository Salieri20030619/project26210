from pathlib import Path
import csv,gzip,json,re,collections,hashlib,sys
import numpy as np,pandas as pd
from scipy import sparse
P=Path(__file__).resolve().parents[1]
cells=pd.read_csv(P/'tables/PILOT_CELL_METADATA.tsv',sep='\t',dtype={'cell_id':str});genes=pd.read_csv(P/'tables/PILOT_FEATURE_AUDIT.tsv',sep='\t')
A=sparse.load_npz(P/'resources/pilot_assigned_counts_qv20.npz');N=sparse.load_npz(P/'resources/pilot_nuclear_counts_qv20.npz')
cm={str(c):int(i) for c,i in zip(cells.cell_id,cells.matrix_row)};gm=dict(zip(genes.gene,genes.matrix_col))
# Independent csv reader and scalar accumulation; selected complete sparse rows for count verification.
selected=set(sorted(cm,key=lambda x:hashlib.sha256(x.encode()).hexdigest())[:64])
ca=collections.Counter();cn=collections.Counter();ga=collections.Counter();gn=collections.Counter();coords=collections.defaultdict(lambda:[0.,0.,0.,0])
selected_a=collections.Counter();selected_n=collections.Counter();raw_gene=collections.Counter();rawrows=unassigned=qualified=nuclear=0
ctrl=re.compile(r'^(NegControl|NegativeControl|Blank|Unassigned|antisense)',re.I)
with gzip.open(P/'downloads/pilot_tumor1_region6.csv.gz','rt') as f:
 for r in csv.DictReader(f):
  rawrows+=1;raw_gene[r['feature_name']]+=1;c=r['cell_id'];g=r['feature_name']
  if c in ['-1','UNASSIGNED','nan','None','']:
   unassigned+=1;continue
  if float(r['qv'])<20 or ctrl.match(g):continue
  qualified+=1;ca[c]+=1;ga[g]+=1
  if c in selected:selected_a[(c,g)]+=1
  if r['overlaps_nucleus']=='1':
   nuclear+=1;cn[c]+=1;gn[g]+=1;s=coords[c]
   for k,v in enumerate(['x_location','y_location','z_location']):s[k]+=float(r[v])
   s[3]+=1
   if c in selected:selected_n[(c,g)]+=1
checks=[]
def check(name,ok):
 checks.append({'check':name,'pass':bool(ok)})
 if not ok:raise AssertionError(name)
check('raw_rows',rawrows==3194464)
check('assigned_conservation',qualified==int(A.sum())==sum(ca.values()))
check('nuclear_conservation',nuclear==int(N.sum())==sum(cn.values()))
check('cell_ids',set(ca)==set(cm))
check('all_cell_assigned_totals',all(ca[c]==int(cells.iloc[i].assigned_RNA) for c,i in cm.items()))
check('all_cell_nuclear_totals',all(cn[c]==int(cells.iloc[i].nuclear_RNA) for c,i in cm.items()))
check('all_gene_raw_totals',all(raw_gene[g]==int(genes.iloc[i].all_transcripts_raw) for g,i in gm.items()))
check('all_gene_assigned_totals',all(ga[g]==int(A[:,i].sum()) for g,i in gm.items()))
check('all_gene_nuclear_totals',all(gn[g]==int(N[:,i].sum()) for g,i in gm.items()))
centroid_error=0.
for c,s in coords.items():
 i=cm[c];ref=np.array(s[:3])/s[3];actual=cells.iloc[i][['nucleus_x_proxy_um','nucleus_y_proxy_um','nucleus_z_proxy_um']].to_numpy(float);centroid_error=max(centroid_error,float(np.max(np.abs(ref-actual))))
check('all_nuclear_centroid_proxies',centroid_error<1e-8)
for c in selected:
 i=cm[c]
 check('assigned_sparse_row_'+c,np.array_equal(A[i].toarray().ravel(),[selected_a[c,g] for g in genes.gene]))
 check('nuclear_sparse_row_'+c,np.array_equal(N[i].toarray().ravel(),[selected_n[c,g] for g in genes.gene]))
check('control_columns_zero',all(int(A[:,i].sum())==0 and int(N[:,i].sum())==0 for g,i in gm.items() if ctrl.match(g)))
check('nested_counts',(A-N).min()>=0)
expected=((cells.assigned_RNA>=20)&(cells.assigned_genes>=5)&(cells.nuclear_RNA>=5))
check('technical_qc',np.array_equal(expected.to_numpy(),cells.technical_qc_pass.to_numpy()))
meta=json.loads((P/'metadata/figshare_current.json').read_text());filemap=pd.read_csv(P/'tables/FILE_REGION_MAP.tsv',sep='\t')
check('all_10_file_paths',len(filemap)==10 and all(meta['folder_structure'][str(r.file_id)].endswith(r.region) for r in filemap.itertuples()))
check('three_labels',set(filemap.tumor_label)=={'tumor1','tumor3','tumor4'})
check('original_download_md5',hashlib.md5((P/'downloads/pilot_tumor1_region6.csv.gz').read_bytes()).hexdigest()=='cfeea875a325f82b340c2a192e1a9016')
legacy=json.loads((P/'qc/PILOT_INPUT_VALIDATION.json').read_text());legacy['legacy_unassigned_flag_counts_not_for_use']=legacy.pop('unassigned_flag_counts');legacy['unassigned_transcripts_corrected']=unassigned;legacy['correction']='Raw numeric -1 marks unassigned in this older Xenium schema; initial reader only tested UNASSIGNED. Independent reconstruction already handled both; no counts/coordinates/effects changed.'
(P/'qc/PILOT_INPUT_VALIDATION_v3.json').write_text(json.dumps(legacy,indent=2))
summary={'status':'PASS','named_checks':len(checks),'checks':checks,'independent_csv_full_read':rawrows,'selected_complete_sparse_rows':len(selected),'all_cells_checked':len(cells),'all_features_checked':len(genes),'maximum_centroid_error_um':centroid_error,'noncontrol_observed_genes':sum(not bool(ctrl.match(g)) for g in gm),'unassigned_transcripts':unassigned,'no_P_or_FDR':True,'scope':'technical input correctness, not biological validation'}
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='checks'},indent=2),flush=True)
