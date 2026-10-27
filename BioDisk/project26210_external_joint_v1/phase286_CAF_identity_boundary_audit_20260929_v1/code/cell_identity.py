from common import *
import h5py
from scipy import sparse
from anndata.io import read_elem
PANELS={
'fibro_core':['DCN','LUM','COL1A1','COL1A2','COL3A1','PDGFRA'],
'fibro_support':['DPT','COL6A1','COL6A2','FBLN1','PI16','CCDC80'],
'pericyte':['RGS5','CSPG4','MCAM','ABCC9','KCNJ8','NOTCH3','PDGFRB'],
'smc':['MYH11','CNN1','ACTA2','TAGLN','MYL9'],
'endothelial':['PECAM1','VWF','EMCN','KDR','CDH5'],
'immune':['PTPRC','LST1','TYROBP','FCER1G','CD3D','CD3E'],
'epithelial':['EPCAM','KRT8','KRT18','KRT19','PAX8'],
'activation':['FAP','POSTN','COL11A1','ACTA2','TAGLN','CTGF']}
def tobool(s):
 return s.astype(str).str.lower().map({'true':True,'false':False,'1':True,'0':False,'1.0':True,'0.0':False})
def main():
 gate('cell_audit_start');assert sha(P/'PROTOCOL_BEFORE_AUDIT.md')==json.loads((P/'FREEZE.json').read_text())['protocol_sha256']
 stat=H.stat();finger=dict(path=str(H),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns)
 with h5py.File(H,'r') as f:
  ids=f['obs']['cell_id'].asstr()[:]
  fields=['study','curated_patient_id','sample_id','source_tissue','cell_type_broad_v1','annotation_cluster_v1','annotation_confidence_v1','annotation_evidence_v1','total_counts_raw','n_genes_raw','low_quality','analysis_keep','scrublet_predicted_doublet','pct_counts_mt_raw']
  m=pd.DataFrame({k:np.asarray(read_elem(f['obs'][k])) for k in fields});m.insert(0,'cell_id',ids)
  genes=f['var']['gene'].asstr()[:].tolist();assert len(genes)==len(set(genes)) and m.cell_id.is_unique and len(m)==258120
 if 'CTGF' not in genes and 'CCN2' in genes:PANELS['activation'][-1]='CCN2'
 panelmeta=[dict(panel=k,requested=len(v),present=sum(g in genes for g in v),missing=';'.join(g for g in v if g not in genes),genes=';'.join(v)) for k,v in PANELS.items()]
 assert all(x in genes for x in PANELS['fibro_core'])
 union=sorted(set(sum(list(PANELS.values()),[]))&set(genes));gi={g:genes.index(g) for g in union};index={g:i for i,g in enumerate(union)}
 raw=np.zeros((len(m),len(union)),dtype=np.float64);tot=np.zeros(len(m));ng=np.zeros(len(m),int)
 with h5py.File(H,'r',rdcc_nbytes=64*1024*1024) as f:
  x=f['layers/counts'];ptr=x['indptr'][:]
  for start in range(0,len(m),4096):
   end=min(start+4096,len(m));lo,hi=int(ptr[start]),int(ptr[end])
   block=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(genes)))
   assert np.isfinite(block.data).all() and (block.data>=0).all()
   raw[start:end]=block[:,[gi[g] for g in union]].toarray()
   tot[start:end]=np.asarray(block.sum(axis=1,dtype=np.float64)).ravel();ng[start:end]=np.asarray((block>0).sum(axis=1)).ravel()
   if start//4096%10==0:print('EXTRACTED',end,'/',len(m),flush=True)
 assert np.all(tot>0)
 np.savez_compressed(P/'inputs/IDENTITY_PANEL_COUNTS.npz',cell_id=ids.astype(str),genes=np.array(union),counts=raw,total_counts=tot,n_genes_current=ng)
 m.to_csv(P/'inputs/ORIGINAL_CELL_METADATA.tsv.gz',sep='\t',index=False)
 log=np.log1p(raw/tot[:,None]*1e4)
 for name,gs in PANELS.items():
  jj=[index[g] for g in gs if g in index];m['det_'+name]=(raw[:,jj]>0).sum(1);m['score_'+name]=log[:,jj].mean(1)
 for g in union:m['expr_'+g]=log[:,index[g]]
 anchors=(raw[:,[index[g] for g in ['DCN','LUM','PDGFRA']]]>0).sum(1)
 coll=(raw[:,[index[g] for g in ['COL1A1','COL1A2','COL3A1']]]>0).sum(1)
 m['fibro_candidate']=(m.det_fibro_core>=3)&(anchors>=1)&(coll>=1)
 m['fibro_strict']=(m.det_fibro_core>=4)&(anchors>=2)&(coll>=1)
 keep=tobool(m.analysis_keep);low=tobool(m.low_quality);doublet=tobool(m.scrublet_predicted_doublet)
 m['qc_flags_known']=keep.notna()&low.notna()&doublet.notna()
 m['qc_pass']=m.qc_flags_known&keep.fillna(False)&~low.fillna(True)&~doublet.fillna(True)&m.n_genes_raw.ge(500)&m.total_counts_raw.ge(1000)
 m['competing_immune']=m.det_immune>=2;m['competing_endo']=m.det_endothelial>=2;m['competing_epi']=m.det_epithelial>=2
 m['strong_mural']=((raw[:,[index[g] for g in ['RGS5','CSPG4','ABCC9','KCNJ8']]]>0).sum(1)>=2)|((raw[:,index['MYH11']]>0)&(raw[:,index['CNN1']]>0))
 m['high_priority_review']=m.fibro_strict&m.qc_pass&~m.competing_immune&~m.competing_endo&~m.competing_epi&~m.strong_mural
 m['tumor']=m.source_tissue.astype(str).str.lower().eq('tumor')
 m.to_csv(P/'inputs/CELL_IDENTITY_AUDIT.tsv.gz',sep='\t',index=False)
 cols=['cell_id','study','curated_patient_id','sample_id','source_tissue','cell_type_broad_v1','annotation_cluster_v1','annotation_confidence_v1','fibro_candidate','fibro_strict','qc_pass','competing_immune','competing_endo','competing_epi','strong_mural','high_priority_review']+[c for c in m if c.startswith(('det_','score_','expr_'))]
 m.loc[m.fibro_candidate|m.cell_type_broad_v1.eq('Fibroblast'),cols].to_csv(P/'tables/CANDIDATE_CELL_REVIEW.tsv.gz',sep='\t',index=False)
 groups=[]
 for k,g in m.groupby(['study','source_tissue','cell_type_broad_v1'],observed=True):
  groups.append(dict(zip(['study','tissue','original_label'],k),cells=len(g),patients=g.curated_patient_id.nunique(),fibro_candidate=int(g.fibro_candidate.sum()),fibro_strict=int(g.fibro_strict.sum()),strict_QC=int((g.fibro_strict&g.qc_pass).sum()),high_priority_review=int(g.high_priority_review.sum()),candidate_mural=int((g.fibro_candidate&g.strong_mural).sum()),candidate_other_lineage=int((g.fibro_candidate&(g.competing_epi|g.competing_endo|g.competing_immune)).sum())))
 tab(groups,'ORIGINAL_LABEL_SCREEN_COUNTS.tsv')
 patients=[]
 for k,g in m.groupby(['study','curated_patient_id','source_tissue'],observed=True):
  patients.append(dict(zip(['study','patient','tissue'],k),cells=len(g),original_fibro=int(g.cell_type_broad_v1.eq('Fibroblast').sum()),candidate=int(g.fibro_candidate.sum()),strict=int(g.fibro_strict.sum()),strict_QC=int((g.fibro_strict&g.qc_pass).sum()),high_priority_review=int(g.high_priority_review.sum()),high_priority_nonFibro=int((g.high_priority_review&g.cell_type_broad_v1.ne('Fibroblast')).sum()),libraries=';'.join(sorted(set(g.sample_id.astype(str))))))
 pt=pd.DataFrame(patients);tab(pt,'PATIENT_CAF_REVIEW_CAPACITY.tsv')
 capacity=[]
 for tissue,z in pt.groupby('tissue'):
  for threshold in [10,20,50]:
   for field in ['original_fibro','candidate','strict_QC','high_priority_review']:
    ss=z[z[field]>=threshold];capacity.append(dict(tissue=tissue,threshold=threshold,screen=field,patients=len(ss),studies=ss.study.nunique(),cells=int(ss[field].sum())))
 tab(capacity,'REFERENCE_CAPACITY_DESCRIPTIVE.tsv')
 marker=[]
 for name,mask in [('original_fibro',m.cell_type_broad_v1.eq('Fibroblast')),('high_priority_review',m.high_priority_review),('candidate_with_mural',m.fibro_candidate&m.strong_mural)]:
  for (study,tissue),gg in m[mask].groupby(['study','source_tissue'],observed=True):
   jj=gg.index.to_numpy(int)
   for gene in union:marker.append(dict(group=name,study=study,tissue=tissue,gene=gene,cells=len(jj),positive=int((raw[jj,index[gene]]>0).sum()),detection_fraction=float((raw[jj,index[gene]]>0).mean()),mean_log1pCP10K=float(log[jj,index[gene]].mean())))
 tab(marker,'IDENTITY_MARKER_DESCRIPTIVE.tsv');tab(panelmeta,'IDENTITY_PANEL_COVERAGE.tsv')
 ann=m.groupby(['cell_type_broad_v1','annotation_cluster_v1','annotation_confidence_v1','annotation_evidence_v1'],observed=True).size().reset_index(name='cells');tab(ann,'ORIGINAL_ANNOTATION_RULES.tsv')
 now=H.stat();assert (now.st_size,now.st_mtime_ns)==(stat.st_size,stat.st_mtime_ns)
 dump(P/'resources/IDENTITY_PANELS_FROZEN.json',PANELS)
 dump(P/'qc/CELL_AUDIT_COMPLETE.json',dict(status='COMPUTED',source=finger,cells=len(m),patients=m.curated_patient_id.nunique(),libraries=m.sample_id.nunique(),count_total_max_difference=float(np.max(np.abs(tot-m.total_counts_raw.to_numpy()))),n_genes_max_difference=int(np.max(np.abs(ng-m.n_genes_raw.to_numpy()))),tumor_original_fibro=int((m.tumor&m.cell_type_broad_v1.eq('Fibroblast')).sum()),tumor_candidates=int((m.tumor&m.fibro_candidate).sum()),tumor_high_priority=int((m.tumor&m.high_priority_review).sum()),labels_changed=False,new_P_values=False,figures=0))
 print('CELL_AUDIT_COMPLETE',flush=True)
if __name__=='__main__':main()
