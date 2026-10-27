
from pathlib import Path
import json,gzip,hashlib,datetime,sys
import numpy as np,pandas as pd,h5py,scipy
from anndata.io import read_elem
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix,save_npz
from scipy.sparse.csgraph import connected_components
P=Path('/mnt/BioDisk/wangmingyuan'); O=P/'project26210_external_joint_v1/phase166_A_technical_pilot_20260915_v1';O.mkdir(exist_ok=True)
S=P/'project26210_external_joint_v1_phase54_gse175540_spatial_v1';F=P/'project26210_phase78_fine_state_spatial_v1/run_v1/sections_v1'
cfg=dict(sections=['GSM5924034','GSM5924045'],k=6,max_edge_over_median_nearest=1.3,min_valid_neighbors=3,region_quantile=.5,min_region_spots=50,condition_number_review_threshold=30,candidates=['SRM','SMS','SLC16A3'],association_tests=False)
(O/'PILOT_CONFIG.json').write_text(json.dumps(cfg,indent=2))
genes='TYMP SRM SMS SLC16A3 ODC1 MTAP AMD1 MAT2A SPP1 COL1A1 COL1A2 DCN LUM CA9 EGLN3 VEGFA BNIP3 PDK1'.split()
rows=[];cover=[];design=[]
def z(v):
 v=np.asarray(v);sd=v.std(axis=0);return (v-v.mean(axis=0))/np.where(sd>0,sd,1)
for gsm in cfg['sections']:
 a=S/'run_v1/sections_v1'/gsm/'processed_section_v1.h5ad'
 with h5py.File(a) as f:
  obs=read_elem(f['obs']);var=read_elem(f['var']);xy=read_elem(f['obsm/spatial']);c=read_elem(f['layers/counts'])
 names=list(var.index); inds=[names.index(g) for g in genes];raw=c[:,inds].toarray();tot=np.asarray(c.sum(axis=1)).ravel()
 expr=np.log1p(raw/tot[:,None]*1e4);e=pd.DataFrame(expr,index=obs.index,columns=genes)
 ab=pd.read_csv(F/gsm/'fine_state_abundance_v1.tsv.gz',sep='\t',index_col=0).loc[obs.index]
 q=ab.filter(regex='^q05__');fra=q.div(q.sum(axis=1),axis=0)
 mal=fra['q05__Malignant'].to_numpy();mac=fra[['q05__C1QC-state Macrophage','q05__Other Macrophage','q05__SPP1-state Macrophage']].sum(axis=1).to_numpy()
 dist,idx=cKDTree(xy).query(xy,k=7);spacing=float(np.median(dist[:,1]));rr=np.repeat(np.arange(len(obs)),6);cc=idx[:,1:].ravel();dd=dist[:,1:].ravel();keep=dd<=spacing*1.3
 ar=obs.array_row.to_numpy();ac=obs.array_col.to_numpy();dr=np.abs(ar[rr]-ar[cc]);dc=np.abs(ac[rr]-ac[cc]);lattice=((dr==0)&(dc==2))|((dr==1)&(dc==1))
 # require both physical first-ring proximity and Visium array adjacency
 keep=keep&lattice
 A=csr_matrix((np.ones(keep.sum()),(rr[keep],cc[keep])),shape=(len(obs),len(obs)));A=A.maximum(A.T);degree=np.asarray(A.sum(axis=1)).ravel()
 ncomp,lab=connected_components(A);W=csr_matrix(A.multiply(1/np.maximum(degree,1)[:,None]));save_npz(O/(gsm+'_adjacency.npz'),A)
 mh=mac>=np.median(mac);sh=mal>=np.median(mal);recA=A.multiply(mh[None,:]).tocsr();recn=np.asarray(recA.sum(axis=1)).ravel();R=csr_matrix(recA.multiply(1/np.maximum(recn,1)[:,None]))
 masks={'source':sh,'receiver_link':sh&(recn>=3)}
 stroma=z(e[['COL1A1','COL1A2','DCN','LUM']]).mean(axis=1);hyp=z(e[['CA9','EGLN3','VEGFA','BNIP3','PDK1']]).mean(axis=1)
 basic=np.column_stack([np.log1p(tot),mal,mac]);extended=np.column_stack([basic,stroma,hyp])
 for role,mask in masks.items():
  for tier,C in [('composition',basic),('extended',extended)]:
   if role=='receiver_link':
    C=np.column_stack([C,R@np.log1p(tot),R@mal,R@mac])
    if tier=='extended':C=np.column_stack([C,R@stroma,R@hyp])
   for candidate in cfg['candidates']:
    cols=np.column_stack([e[candidate],C])
    if role=='receiver_link':cols=np.column_stack([cols,e['TYMP']])
    X=np.column_stack([np.ones(mask.sum()),z(cols[mask])]);rank=int(np.linalg.matrix_rank(X));cond=float(np.linalg.cond(X))
    design.append(dict(section=gsm,role=role,tier=tier,candidate=candidate,spots=int(mask.sum()),columns=X.shape[1],rank=rank,condition_number=cond,full_rank=rank==X.shape[1],passes_pilot_rule=bool(mask.sum()>=50 and rank==X.shape[1] and cond<30)))
 for j,g in enumerate(genes):
  for role,mask in {'all':np.ones(len(obs),bool),**masks}.items():
   cover.append(dict(section=gsm,gene=g,region=role,spots=int(mask.sum()),detected=int((raw[mask,j]>0).sum()),fraction=float((raw[mask,j]>0).mean()) if mask.any() else None))
 scale_files=list((S/'input_extract_v1').glob(gsm+'*scalefactors*'))
 scale=json.loads(gzip.open(scale_files[0],'rt').read()) if scale_files else {}
 B=np.column_stack([np.ones(len(obs)),ar,ac]);pred=B@np.linalg.lstsq(B,xy,rcond=None)[0];res=np.sqrt(((pred-xy)**2).sum(axis=1))
 rows.append(dict(section=gsm,spots=len(obs),spacing_pixels=spacing,edge_cutoff_pixels=spacing*1.3,knn_directed_edges=len(keep),retained_directed_edges=int(keep.sum()),removed_fraction=float(1-keep.mean()),lattice_and_distance_disagreement=int(((dd<=spacing*1.3)!=lattice).sum()),undirected_edges=A.nnz//2,components=ncomp,largest_component_fraction=float(np.bincount(lab).max()/len(obs)),isolated=int((degree==0).sum()),spots_ge3_neighbors=int((degree>=3).sum()),source_spots=int(sh.sum()),receiver_link_spots=int(masks['receiver_link'].sum()),affine_coordinate_max_residual_pixels=float(res.max()),scalefactors=json.dumps(scale),malignant_threshold=float(np.median(mal)),macrophage_threshold=float(np.median(mac))))
 pd.DataFrame(dict(spot_id=obs.index,degree=degree,component=lab,source_region=sh,macrophage_region=mh,eligible_receiver_link=masks['receiver_link'],eligible_macrophage_neighbors=recn)).to_csv(O/(gsm+'_SPOT_TECHNICAL_QC.tsv.gz'),sep='\t',index=False)
 print('completed',gsm,flush=True)
pd.DataFrame(rows).to_csv(O/'GRAPH_QC.tsv',sep='\t',index=False);pd.DataFrame(cover).to_csv(O/'REGION_GENE_COVERAGE.tsv',sep='\t',index=False);pd.DataFrame(design).to_csv(O/'DESIGN_IDENTIFIABILITY.tsv',sep='\t',index=False)
ref=P/'project26210_part4_tymp_reset_phase17_20260812_v1/reference_tymp_continuous_v1/reference_tymp_continuous_v1_signatures_study_balanced_mean_counts.tsv.gz'
r=pd.read_csv(ref,sep='\t',index_col=0); r.reindex(genes).to_csv(O/'REFERENCE_TARGET_GENE_AUDIT.tsv',sep='\t')
status=dict(stage='Phase166',technical_pilot_completed=True,all_designs_pass=bool(pd.DataFrame(design).passes_pilot_rule.all()),association_tests=0,new_figures=0,patient_mapping='UNRESOLVED',python=sys.version,numpy=np.__version__,scipy=scipy.__version__,created_at=datetime.datetime.now().astimezone().isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print(json.dumps(rows));print(json.dumps(status));print(pd.DataFrame(design).to_string(index=False))

