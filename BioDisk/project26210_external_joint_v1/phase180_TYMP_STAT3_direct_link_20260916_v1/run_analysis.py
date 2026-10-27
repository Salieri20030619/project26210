from pathlib import Path
import pandas as pd,numpy as np,json,gzip,hashlib,shutil
from scipy.io import mmread
from scipy.stats import rankdata,spearmanr,wilcoxon,norm,t
from threadpoolctl import threadpool_limits
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase180_TYMP_STAT3_direct_link_20260916_v1";prev=r/"phase179_STAT3_external_activity_validation_20260916_v1/run_v2"
K=pd.read_csv(prev/"FROZEN_STAT3_TARGETS.tsv",sep="\t");shutil.copy2(prev/"FROZEN_STAT3_TARGETS.tsv",d/"FROZEN_STAT3_TARGETS.tsv")
files=[prev/"FROZEN_STAT3_TARGETS.tsv"];checks=[]
def rz(v):
 v=rankdata(v);return (v-v.mean())/v.std(ddof=1)
def score(E,sign,sc):
 A=rankdata(E,axis=0) if sc=="rankmean" else E
 sd=A.std(0,ddof=1);ok=sd>1e-12
 return np.round((((A[:,ok]-A[:,ok].mean(0))/sd[ok])*sign[ok]).mean(1),12),int(ok.sum())
def corr(x,y,covs):
 C=np.column_stack([np.ones(len(x))]+[rankdata(v) for v in covs]);a=rankdata(x);b=rankdata(y)
 aa=a-C@np.linalg.lstsq(C,a,rcond=None)[0];bb=b-C@np.linalg.lstsq(C,b,rcond=None)[0];rh=float(np.corrcoef(aa,bb)[0,1]);k=np.linalg.matrix_rank(C)-1;df=len(a)-k-2
 q,_=np.linalg.qr(C);ar=a-q@(q.T@a);br=b-q@(q.T@b)
 # Only compare QR when full rank.
 if np.linalg.matrix_rank(C)==C.shape[1]:checks.append({"check":"partial_QR","pass":abs(rh-np.corrcoef(ar,br)[0,1])<1e-10})
 return rh,float(2*t.sf(abs(rh)*np.sqrt(df/max(1e-15,1-rh*rh)),df)),k
def summary(v):
 v=np.asarray(v);rng=np.random.default_rng(180);b=np.median(v[rng.integers(0,len(v),(5000,len(v)))],axis=1)
 return dict(n=len(v),median_rho=np.median(v),positive=int((v>0).sum()),negative=int((v<0).sum()),P=wilcoxon(v,method="approx",correction=True).pvalue if len(v)>=10 else np.nan,CI_low=np.quantile(b,.025),CI_high=np.quantile(b,.975))
# A malignant cells.
q=Path("/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/07_pyscenic_v1/input")
gf=q/"genes_v1.tsv";cf=q/"cells_v1.tsv";xf=q/"raw_counts_genes_by_cells_v1.mtx.gz";mf=r/"phase162_MTA_polyamine_part1_reassessment_20260915_v1/PART1_TARGET_CELL_METADATA.tsv.gz";files +=[gf,cf,xf,mf]
genes=pd.read_csv(gf,sep="\t").gene.tolist();ids=pd.read_csv(cf,sep="\t").exact_cell_id.tolist();m=pd.read_csv(mf,sep="\t").set_index("exact_cell_id").loc[ids]
with gzip.open(xf,"rb") as f:
 with threadpool_limits(limits=1):X=mmread(f).tocsr()
assert X.shape==(len(genes),11950)
tot=np.asarray(X.sum(0)).ravel();ng=np.asarray((X>0).sum(0)).ravel();yy=X[genes.index("TYMP")].toarray().ravel()
assert np.array_equal(yy,m["count:TYMP"]) and np.array_equal(tot,m.total_counts)
assert set(K.gene)<=set(genes)
target=sorted(set(K.gene));E=np.log1p(X[[genes.index(g) for g in target]].toarray().T/tot[:,None]*1e4)
gmt=r/"part2_rebuild_v1/06_functional_programs/run_v1/resources/MSigDB_Hallmark_2020.gmt";files.append(gmt);prog={}
for line in gmt.read_text().splitlines():
 fields=line.split("\t")
 if fields[0] in ["G2-M Checkpoint","Hypoxia"]:
  gg=[g for g in fields[2:] if g in genes and g not in {"TYMP","SPP1","STAT3"}]
  prog[fields[0]]=np.log1p(X[[genes.index(g) for g in gg]].toarray().T/tot[:,None]*1e4)
del X
rows=[];elig=[];cache={}
for patient,mm in m.groupby("curated_patient_id"):
 ix=m.index.get_indexer(mm.index);dep=np.log1p(tot[ix]);gene=ng[ix];cy=prog["G2-M Checkpoint"][ix];hy=prog["Hypoxia"][ix]
 cyc=score(cy,np.ones(cy.shape[1]),"zmean")[0];hyp=score(hy,np.ones(hy.shape[1]),"zmean")[0]
 for res,kk in K.groupby("resource"):
  ee=E[ix][:,[target.index(g) for g in kk.gene]]
  for sc in ["zmean","rankmean"]:
   s,nv=score(ee,kk.sign.to_numpy(),sc);ok=len(ix)>=50 and (yy[ix]>0).sum()>=10 and nv>=(20 if res=="collectri" else 10)
   elig.append(dict(patient=patient,resource=res,scoring=sc,cells=len(ix),TYMP_detected=int((yy[ix]>0).sum()),variable_targets=nv,eligible=ok))
   if not ok:continue
   for adj,cv in [("raw",[]),("depth",[dep]),("depth_genes",[dep,gene]),("depth_cycle_hypoxia",[dep,cyc,hyp]),("depth_genes_cycle_hypoxia",[dep,gene,cyc,hyp])]:
    rr,_,_=corr(yy[ix],s,cv);rows.append(dict(patient=patient,study=mm.study.iloc[0],resource=res,scoring=sc,adjustment=adj,cells=len(ix),rho=rr))
   if sc=="zmean":cache[res+"__"+patient]=np.column_stack([yy[ix],s,dep,gene,cyc,hyp])
PD=pd.DataFrame(rows);PD.to_csv(d/"MALIGNANT_PATIENT_EFFECTS.tsv",sep="\t",index=False);pd.DataFrame(elig).to_csv(d/"MALIGNANT_ELIGIBILITY.tsv",sep="\t",index=False);np.savez_compressed(d/"MALIGNANT_VALIDATION_CACHE.npz",**cache)
sm=[];los=[]
for (res,sc,adj),p in PD.groupby(["resource","scoring","adjustment"]):
 sm.append(dict(resource=res,scoring=sc,adjustment=adj,**summary(p.rho.to_numpy())))
 for st in p.study.unique():
  v=p[p.study!=st].rho;los.append(dict(resource=res,scoring=sc,adjustment=adj,omitted=st,n=len(v),median_rho=v.median(),positive=int((v>0).sum())))
SM=pd.DataFrame(sm);SM.to_csv(d/"MALIGNANT_SUMMARY.tsv",sep="\t",index=False);pd.DataFrame(los).to_csv(d/"MALIGNANT_LEAVE_STUDY.tsv",sep="\t",index=False)
print("MALIGNANT\n"+SM[(SM.resource=="collectri")&(SM.scoring=="zmean")].to_csv(sep="\t",index=False),flush=True)
# B Same-library sender to receiver: aggregate first, score across patients.
bf=r/"phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv";gf=r/"phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_INDEX.tsv";zf=r/"phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_PSEUDOBULK.npz";files +=[bf,gf,zf]
b=pd.read_csv(bf,sep="\t");G=pd.read_csv(gf,sep="\t");Z=np.load(zf,allow_pickle=True);gn=list(Z["genes"]);ex=[];audit=[]
for _,row in b.iterrows():
 gg=G[(G.patient==row.patient)&(G.sample_id==row.sample_id)&G.role.isin(["C4","OtherMacrophage"])]
 assert gg.cells.sum()==row.macrophage_cells
 c=Z["counts"][gg.group_index.to_numpy(int)].sum(0);assert c[gn.index("SPP1")]==row.macrophage_SPP1_counts
 ex.append(np.log2(c[[gn.index(g) for g in target]]/c.sum()*1e6+1))
 audit.append(dict(patient=row.patient,sample_id=row.sample_id,cells=int(gg.cells.sum()),counts=float(c.sum()),SPP1_counts=float(c[gn.index("SPP1")])))
ex=np.array(ex);pd.DataFrame(audit).to_csv(d/"SAME_LIBRARY_AGGREGATION_AUDIT.tsv",sep="\t",index=False);br=[]
groups=[np.asarray(i) for i in b.groupby("study").groups.values()]
D=pd.get_dummies(b.study,drop_first=True,dtype=float).to_numpy()
for res,kk in K.groupby("resource"):
 for sc in ["zmean","rankmean"]:
  s,nv=score(ex[:,[target.index(g) for g in kk.gene]],kk.sign.to_numpy(),sc);b[res+"_"+sc]=s
  x=rankdata(b.TYMP_pseudobulk);y=rankdata(s)
  for adj in ["raw","study","study_C4"]:
   cv=[] if adj=="raw" else [D[:,j] for j in range(D.shape[1])]
   if adj=="study_C4":cv.append(rankdata(b.C4_fraction))
   C=np.column_stack([np.ones(len(b))]+cv);Q,_=np.linalg.qr(C);a=x-Q@(Q.T@x);y0=y-Q@(Q.T@y);obs=np.corrcoef(a,y0)[0,1]
   if adj=="raw":pv=spearmanr(x,y).pvalue
   else:
    rng=np.random.default_rng(262180);hits=0
    for _ in range(9999):
     xp=a.copy()
     for idx in groups:xp[idx]=rng.permutation(xp[idx])
     xp-=Q@(Q.T@xp);rr=np.dot(xp,y0)/np.sqrt(np.dot(xp,xp)*np.dot(y0,y0))
     hits+=abs(rr)>=abs(obs)-1e-14
    pv=(hits+1)/10000
   br.append(dict(resource=res,scoring=sc,adjustment=adj,n=len(b),rho=obs,P=pv,method="Spearman" if adj=="raw" else "within-study residual permutation",variable_targets=nv))
b.to_csv(d/"SAME_LIBRARY_PATIENT_VALUES.tsv",sep="\t",index=False);BR=pd.DataFrame(br);BR.to_csv(d/"SAME_LIBRARY_TYMP_MACRO_STAT3.tsv",sep="\t",index=False)
print("BRIDGE\n"+BR.to_csv(sep="\t",index=False),flush=True)
# C paired bulk, frozen previous final scores retained unchanged.
bo=[]
model={"raw":[],"MTAP_myeloid":["MTAP","MYELOID_PROXY"],"MTAP_myeloid_MTA":["MTAP","MYELOID_PROXY","MET:5-methylthioadenosine (MTA)"],"MTAP_myeloid_hypoxia":["MTAP","MYELOID_PROXY","HYPOXIA_PROXY"]}
for co in ["RC18","RC20","CPTAC","CPTAC_val"]:
 for res in ["collectri","dorothea"]:
  for sc in ["zmean","rankmean"]:
   fpath=prev/f"BULK_PATIENTS_{co}_{res}_{sc}.tsv";files.append(fpath);f=pd.read_csv(fpath,sep="\t",index_col=0)
   for adj,cc in model.items():
    ff=f[["TYMP","activity"]+cc].dropna();rr,pv,k=corr(ff.TYMP,ff.activity,[ff[c] for c in cc]);bo.append(dict(cohort=co,resource=res,scoring=sc,adjustment=adj,n=len(ff),rho=rr,P=pv,k=k))
BO=pd.DataFrame(bo);BO.to_csv(d/"BULK_TYMP_STAT3_COHORT.tsv",sep="\t",index=False)
def meta(p):
 zz=np.arctanh(p.rho.to_numpy());var=1/(p.n-p.k-3).to_numpy();w=1/var;mu=np.average(zz,weights=w);Q=np.sum(w*(zz-mu)**2);df=len(p)-1;tau=max(0,(Q-df)/(w.sum()-(w*w).sum()/w.sum()));w=1/(var+tau);mu=np.average(zz,weights=w);se=1/np.sqrt(w.sum())
 return dict(n=int(p.n.sum()),positive_cohorts=int((p.rho>0).sum()),rho=np.tanh(mu),P=2*norm.sf(abs(mu/se)),CI_low=np.tanh(mu-1.96*se),CI_high=np.tanh(mu+1.96*se),I2=max(0,(Q-df)/Q)*100 if Q else 0)
bm=[];bl=[]
for (res,sc,adj),p in BO.groupby(["resource","scoring","adjustment"]):
 bm.append(dict(resource=res,scoring=sc,adjustment=adj,**meta(p)))
 for co in p.cohort:bl.append(dict(resource=res,scoring=sc,adjustment=adj,omitted=co,**meta(p[p.cohort!=co])))
BM=pd.DataFrame(bm);BM.to_csv(d/"BULK_TYMP_STAT3_META.tsv",sep="\t",index=False);pd.DataFrame(bl).to_csv(d/"BULK_TYMP_STAT3_LEAVE_COHORT.tsv",sep="\t",index=False)
# Explicitly retain old endpoints with byte-identical copies and provenance.
keep=d/"RETAINED_PHASE179";keep.mkdir(exist_ok=True)
for name in ["SC_ACTIVITY_SPP1_SUMMARY.tsv","BULK_EDGES_META.tsv","BULK_SHORTPATH_META.tsv"]:
 shutil.copy2(prev/name,keep/name);checks.append({"check":"retained_"+name,"pass":(prev/name).read_bytes()==(keep/name).read_bytes()});files.append(prev/name)
manifest=[]
for f in files:
 h=hashlib.sha256()
 with f.open("rb") as fh:
  for block in iter(lambda:fh.read(8*1024*1024),b""):h.update(block)
 manifest.append({"path":str(f),"bytes":f.stat().st_size,"sha256":h.hexdigest()})
(d/"INPUT_MANIFEST.json").write_text(json.dumps(manifest,indent=2));pd.DataFrame(checks).to_csv(d/"VALIDATION.tsv",sep="\t",index=False);assert all(v["pass"] for v in checks)
print("BULK\n"+BM.to_csv(sep="\t",index=False),flush=True);print("CHECKS",len(checks),flush=True)
