from pathlib import Path
import json,sys,hashlib,datetime,platform
import numpy as np,pandas as pd,h5py,scipy,anndata
from anndata.io import read_elem
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix,save_npz
from scipy.sparse.csgraph import connected_components
from scipy.stats import rankdata
P=Path('/mnt/BioDisk/wangmingyuan')
ROOT=P/'project26210_external_joint_v1'
O=Path(__file__).resolve().parent
IN=ROOT/'phase165_A_input_confirmation_20260915_v1'
old=ROOT/'phase166_A_technical_pilot_20260915_v1/strict_region_v2'
sections=pd.read_csv(IN/'SECTION_INPUT_QC.tsv',sep='\t')
assert len(sections)==24 and sections.section.is_unique
cfg=json.loads((old/'PILOT_CONFIG.json').read_text())
cfg.update(sections=sections.section.tolist(),stage='Phase167',rank_transform='average ranks within each endpoint mask, separately for every non-intercept column and outcome',outcome_gate='finite, nonzero variance; no association fitting',overlap_gate='descriptive only; no invented independence threshold')
(O/'FROZEN_CONFIG.json').write_text(json.dumps(cfg,indent=2))
genes='TYMP SRM SMS SLC16A3 ODC1 MTAP AMD1 MAT2A SPP1 COL1A1 COL1A2 DCN LUM CA9 EGLN3 VEGFA BNIP3 PDK1'.split()
checks=[]; inputs=[]; grows=[]; designs=[]; coverage=[]; overlaps=[]
def check(section,name,value):
 checks.append(dict(section=section,check=name,passed=bool(value)))
 if not value: raise RuntimeError(f'{section}: {name}')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def record(p,role):
 p=Path(p); st=p.stat()
 inputs.append(dict(path=str(p),role=role,bytes=st.st_size,mtime_ns=st.st_mtime_ns,sha256=sha(p)))
def z(v):
 v=np.asarray(v,dtype=float)
 if len(v)==0:return v
 sd=v.std(axis=0)
 return (v-v.mean(axis=0))/np.where(sd>0,sd,1)
for p in [IN/'SECTION_INPUT_QC.tsv',IN/'metadata/GSE175540_family.soft',old/'PILOT_CONFIG.json',old.parent/'run_strict_pilot.py']:
 record(p,'provenance')
for s in sections.itertuples(index=False):
 gsm=s.section; prep=s.preparation
 record(s.input_path,'counts_coordinates');record(s.abundance_path,'frozen_abundance')
 with h5py.File(s.input_path,'r') as f:
  obs=read_elem(f['obs']);var=read_elem(f['var']);xy=read_elem(f['obsm/spatial']);c=read_elem(f['layers/counts'])
 names=list(var.index);n=len(obs)
 check(gsm,'unique_gene_and_spot_ids',len(set(names))==len(names) and obs.index.is_unique)
 check(gsm,'all_18_genes_present',all(g in names for g in genes))
 raw=c[:,[names.index(g) for g in genes]].toarray().astype(float);tot=np.asarray(c.sum(axis=1)).ravel()
 check(gsm,'finite_nonnegative_integer_counts',(np.isfinite(c.data)&(c.data>=0)&(c.data==np.floor(c.data))).all())
 check(gsm,'positive_depth_matches_obs',(tot>0).all() and np.allclose(tot,obs.total_counts))
 check(gsm,'finite_unique_coordinates',np.isfinite(xy).all() and len(np.unique(xy,axis=0))==n)
 e=pd.DataFrame(np.log1p(raw/tot[:,None]*1e4),index=obs.index,columns=genes)
 ab=pd.read_csv(s.abundance_path,sep='\t',index_col=0)
 check(gsm,'abundance_exact_barcode_order',list(obs.index)==list(ab.index))
 q=ab.filter(regex='^q05__')
 check(gsm,'abundance_valid',q.shape[1]==13 and np.isfinite(q.to_numpy()).all() and (q.to_numpy()>=0).all() and (q.sum(axis=1)>0).all())
 fra=q.div(q.sum(axis=1),axis=0)
 mal=fra['q05__Malignant'].to_numpy()
 maccols=['q05__C1QC-state Macrophage','q05__Other Macrophage','q05__SPP1-state Macrophage']
 mac=fra[maccols].sum(axis=1).to_numpy()
 other=fra.drop(columns=['q05__Malignant']+maccols).max(axis=1).to_numpy()
 sh=(mal>=.2)&(mal>=np.maximum(mac,other));mh=(mac>=.1)&(mac>=np.median(mac))
 dist,idx=cKDTree(xy).query(xy,k=7);spacing=float(np.median(dist[:,1]))
 rr=np.repeat(np.arange(n),6);cc=idx[:,1:].ravel();dd=dist[:,1:].ravel()
 ar=obs.array_row.to_numpy();ac=obs.array_col.to_numpy()
 dr=np.abs(ar[rr]-ar[cc]);dc=np.abs(ac[rr]-ac[cc])
 lattice=((dr==0)&(dc==2))|((dr==1)&(dc==1))
 keep=(dd<=spacing*1.3)&lattice
 A=csr_matrix((np.ones(keep.sum()),(rr[keep],cc[keep])),shape=(n,n));A=A.maximum(A.T)
 degree=np.asarray(A.sum(axis=1)).ravel()
 ni,nj=A.nonzero()
 check(gsm,'symmetric_no_self_degree_le6',(A!=A.T).nnz==0 and np.count_nonzero(A.diagonal())==0 and (degree<=6).all())
 check(gsm,'all_edges_distance_lattice',(((np.abs(ar[ni]-ar[nj])==0)&(np.abs(ac[ni]-ac[nj])==2))|((np.abs(ar[ni]-ar[nj])==1)&(np.abs(ac[ni]-ac[nj])==1))).all() and (np.linalg.norm(xy[ni]-xy[nj],axis=1)<=spacing*1.3+1e-8).all())
 nc,lab=connected_components(A)
 recA=A.multiply(mh[None,:]).tocsr();recn=np.asarray(recA.sum(axis=1)).ravel()
 R=csr_matrix(recA.multiply(1/np.maximum(recn,1)[:,None]))
 masks={'source':sh,'receiver_link':sh&(recn>=3)}
 check(gsm,'receiver_normalization',np.allclose(np.asarray(R.sum(axis=1)).ravel()[recn>0],1))
 save_npz(O/(gsm+'_adjacency.npz'),A)
 pd.DataFrame(dict(spot_id=obs.index,degree=degree,component=lab,source_region=sh,macrophage_region=mh,eligible_receiver_link=masks['receiver_link'],eligible_macrophage_neighbors=recn)).to_csv(O/(gsm+'_SPOT_TECHNICAL_QC.tsv.gz'),sep='\t',index=False)
 stroma=z(e[['COL1A1','COL1A2','DCN','LUM']]).mean(axis=1)
 hyp=z(e[['CA9','EGLN3','VEGFA','BNIP3','PDK1']]).mean(axis=1)
 basic=np.column_stack([np.log1p(tot),mal,mac]);extended=np.column_stack([basic,stroma,hyp])
 for role,mask in masks.items():
  m=int(mask.sum())
  for tier,C0 in [('composition',basic),('extended',extended)]:
   C=C0.copy()
   if role=='receiver_link':
    C=np.column_stack([C,R@np.log1p(tot),R@mal,R@mac])
    if tier=='extended':C=np.column_stack([C,R@stroma,R@hyp])
   for candidate in cfg['candidates']:
    pred=e['TYMP'].to_numpy() if role=='source' else e[candidate].to_numpy()
    outcome=e[candidate].to_numpy() if role=='source' else R@e['SPP1'].to_numpy()
    cols=np.column_stack([pred,C])
    if role=='receiver_link':cols=np.column_stack([cols,e['TYMP']])
    for transform in ['linear','rank']:
     v=cols[mask];y=outcome[mask]
     if transform=='rank' and m:
      v=rankdata(v,axis=0,method='average');y=rankdata(y,method='average')
     X=np.column_stack([np.ones(m),z(v)])
     finite=bool(np.isfinite(X).all() and np.isfinite(y).all())
     r=int(np.linalg.matrix_rank(X)) if m and finite else 0
     cond=float(np.linalg.cond(X)) if m and finite else np.inf
     yvar=float(np.var(y)) if m else 0.0
     reasons=[]
     if m<50:reasons.append('fewer_than_50_source_spots')
     if not finite:reasons.append('nonfinite')
     if r!=X.shape[1]:reasons.append('rank_deficient')
     if not np.isfinite(cond) or cond>=30:reasons.append('condition_number_ge30')
     if yvar<=0:reasons.append('constant_outcome')
     designs.append(dict(section=gsm,preparation=prep,role=role,tier=tier,candidate=candidate,transform=transform,spots=m,columns=X.shape[1],rank=r,condition_number=cond,outcome_variance=yvar,outcome_unique=len(np.unique(y)),predictor_unique=len(np.unique(pred[mask])),passes_frozen_design_rule=bool(m>=50 and finite and r==X.shape[1] and cond<30),eligible=not reasons,exclusion_reasons=';'.join(reasons)))
 for j,g in enumerate(genes):
  for role,mask in {'all':np.ones(n,bool),**masks,'macrophage_region':mh}.items():
   coverage.append(dict(section=gsm,preparation=prep,gene=g,region=role,spots=int(mask.sum()),detected=int((raw[mask,j]>0).sum()),fraction=float((raw[mask,j]>0).mean()) if mask.any() else None))
 M=recA[masks['receiver_link']].tocsr();M.eliminate_zeros()
 uses=np.asarray(M.sum(axis=0)).ravel();active=uses[uses>0]
 shared=(M@M.T).tocsr();shared.setdiag(0);shared.eliminate_zeros()
 m=M.shape[0];paircount=shared.nnz//2
 overlaps.append(dict(section=gsm,preparation=prep,source_endpoints=m,unique_receiver_spots=int(len(active)),total_receiver_uses=int(active.sum()),receiver_reuse_max=float(active.max()) if len(active) else 0,receiver_reuse_median=float(np.median(active)) if len(active) else 0,receivers_reused_fraction=float((active>1).mean()) if len(active) else None,source_pairs_sharing_receiver=int(paircount),fraction_source_pairs_sharing_receiver=float(paircount/(m*(m-1)/2)) if m>1 else None,source_endpoints_with_shared_receiver=int((np.diff(shared.indptr)>0).sum())))
 B=np.column_stack([np.ones(n),ar,ac]);res=np.linalg.norm(B@np.linalg.lstsq(B,xy,rcond=None)[0]-xy,axis=1)
 grows.append(dict(section=gsm,preparation=prep,spots=n,source_spots=int(sh.sum()),receiver_region_spots=int(mh.sum()),receiver_link_spots=int(masks['receiver_link'].sum()),spacing_pixels=spacing,edge_cutoff_pixels=spacing*1.3,removed_knn_fraction=float(1-keep.mean()),lattice_distance_disagreement=int(((dd<=spacing*1.3)!=lattice).sum()),undirected_edges=A.nnz//2,components=nc,largest_component_fraction=float(np.bincount(lab).max()/n),isolated=int((degree==0).sum()),coordinate_affine_max_residual=float(res.max()),source_min_fraction=.2,source_dominant_required=True,receiver_fraction_threshold=max(.1,float(np.median(mac)))))
 print('completed',gsm,prep,'source',int(sh.sum()),'receiver_link',int(masks['receiver_link'].sum()),flush=True)
for name,rows in [('INPUT_MANIFEST',inputs),('GRAPH_QC',grows),('DESIGN_ENDPOINT_ELIGIBILITY',designs),('REGION_GENE_COVERAGE',coverage),('NEIGHBOR_OVERLAP',overlaps),('VALIDATION',checks)]:
 pd.DataFrame(rows).to_csv(O/(name+'.tsv'),sep='\t',index=False)
status=dict(stage='Phase167',status='COMPUTATION_COMPLETE_PENDING_REVIEW',sections=24,association_tests=0,new_figures=0,patient_mapping='UNDER_REVIEW',validation_checks=len(checks),validation_passed=sum(x['passed'] for x in checks),python=sys.version,executable=sys.executable,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,h5py=h5py.__version__,anndata=anndata.__version__,created_at=datetime.datetime.now().astimezone().isoformat())
(O/'COMPUTE_STATUS.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status),flush=True)
