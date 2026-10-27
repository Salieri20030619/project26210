from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,sys
from scipy.sparse import load_npz
from scipy.stats import rankdata,wilcoxon,norm,t
from statsmodels.stats.multitest import multipletests
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase179_STAT3_external_activity_validation_20260916_v1"
K=pd.read_csv(d/"FROZEN_STAT3_TARGETS.tsv",sep="\t")
audit=pd.read_csv(d/"TARGET_AUDIT.tsv",sep="\t");assert audit[audit.decision=="keep"].groupby(["resource","gene"]).sign.nunique().max()==1
assert not set(K.gene)&{"TYMP","SPP1","STAT3"}
resources=["collectri","dorothea"];scorings=["zmean","rankmean"]
covsets={"raw":[],"MTAP_myeloid":["MTAP","MYELOID_PROXY"],"MTAP_myeloid_hypoxia":["MTAP","MYELOID_PROXY","HYPOXIA_PROXY"]}
inputs=[d/"FROZEN_STAT3_TARGETS.tsv"];checks=[]
def z(v):
 v=np.asarray(v,float);sd=v.std(axis=0,ddof=1);ok=sd>1e-12
 out=np.zeros_like(v);out[:,ok]=(v[:,ok]-v[:,ok].mean(0))/sd[ok] if v.ndim==2 else 0
 return out,ok
def rz(v):
 a=rankdata(v);return (a-a.mean())/a.std(ddof=1)
def score(E,sign,scoring):
 A=rankdata(E,axis=0) if scoring=="rankmean" else E
 zz,ok=z(A);return (zz[:,ok]*sign[ok]).mean(1),int(ok.sum())
def partial(x,y,covs):
 X=np.column_stack([np.ones(len(x))]+[rankdata(v) for v in covs]);a=rankdata(x);b=rankdata(y)
 a-=X@np.linalg.lstsq(X,a,rcond=None)[0];b-=X@np.linalg.lstsq(X,b,rcond=None)[0]
 if np.std(a)<1e-10 or np.std(b)<1e-10:return np.nan,np.nan
 rr=float(np.corrcoef(a,b)[0,1]);df=len(a)-np.linalg.matrix_rank(X)-1
 return rr,float(2*t.sf(abs(rr)*np.sqrt(df/max(1e-15,1-rr*rr)),df))
def wil(v):
 v=np.asarray(v);return float(wilcoxon(v,method="approx",correction=True).pvalue) if len(v)>=10 and np.any(v) else (1. if len(v)>=10 else np.nan)
mfile=r/"part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz"
gfile=r/"part2_rebuild_v1/08_regulons/input_v1/GENES.tsv"
xfile=r/"part2_rebuild_v1/08_regulons/input_v1/GENES_BY_CELLS_RAW.npz"
ptfile=r/"part2_rebuild_v1/07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz"
inputs +=[mfile,gfile,xfile,ptfile]
m=pd.read_csv(mfile,sep="\t",index_col=0,dtype={"candidate":str},low_memory=False)
genes=pd.read_csv(gfile,sep="\t").gene.tolist();X=load_npz(xfile).tocsr();assert X.shape==(len(genes),len(m))
allg=sorted(set(K.gene)|{"SPP1"});gi={g:i for i,g in enumerate(allg)}
raw=X[[genes.index(g) for g in allg]].toarray().T.astype(float);E=np.log1p(raw/m.total_counts_raw.to_numpy()[:,None]*1e4)
sel=(~m.candidate.isin(["1","Monocyte"]))&(m.source_tissue.str.lower()=="tumor")
tm=pd.read_csv(ptfile,sep="\t",index_col=0,dtype={"partition":str},low_memory=False);tm=tm[tm["mode"]=="primary"];assert not tm.index.duplicated().any()
rows=[];coverage=[];c4=[];cellscore=[]
for patient,mm in m[sel].groupby("curated_patient_id"):
 ids=m.index.get_indexer(mm.index);yy=E[ids,gi["SPP1"]];dep=np.log1p(mm.total_counts_raw.to_numpy());ng=mm.n_genes_raw.to_numpy();cc=(mm.candidate=="4").to_numpy(float)
 for res in resources:
  kk=K[K.resource==res];ix=[gi[g] for g in kk.gene];sgn=kk.sign.to_numpy(float)
  for sc in scorings:
   ss,nvar=score(E[ids][:,ix],sgn,sc);minvar=20 if res=="collectri" else 10
   ok=len(mm)>=50 and (yy>0).sum()>=10 and nvar>=minvar and np.std(ss)>1e-12
   coverage.append(dict(patient=patient,resource=res,scoring=sc,cells=len(mm),SPP1_detected=int((yy>0).sum()),variable_targets=nvar,eligible=ok))
   if not ok:continue
   common=dict(patient=patient,study=mm.study.iloc[0],resource=res,scoring=sc)
   for adj,cv in [("raw",[]),("depth",[dep]),("depth_genes",[dep,ng]),("depth_C4",[dep,cc])]:
    rho,pv=partial(ss,yy,cv);rows.append({**common,"adjustment":adj,"cells":len(mm),"rho":rho})
   q=tm.reindex(mm.index);pr=[]
   for part in q.partition.dropna().unique():
    mask=(q.partition==part).to_numpy()&np.isfinite(q.pt_primary.to_numpy())
    if mask.sum()<50 or (yy[mask]>0).sum()<10 or len(np.unique(q.pt_primary.to_numpy()[mask]))<10:continue
    rr,_=partial(ss[mask],yy[mask],[dep[mask],q.pt_primary.to_numpy()[mask]]);pr.append(rr)
   if pr:rows.append({**common,"adjustment":"depth_pseudotime","cells":len(mm),"rho":float(np.nanmedian(pr))})
   if (cc==1).sum()>=20 and (cc==0).sum()>=20:c4.append({**common,"C4_minus_other_score":float(ss[cc==1].mean()-ss[cc==0].mean())})
   if sc=="zmean":
    cellscore.append(pd.DataFrame({"cell_id":mm.index,"patient":patient,"resource":res,"score":ss,"SPP1":yy,"depth":dep,"candidate":mm.candidate.to_numpy()}))
P=pd.DataFrame(rows);P.to_csv(d/"SC_PATIENT_EFFECTS.tsv",sep="\t",index=False)
pd.DataFrame(coverage).to_csv(d/"SC_ELIGIBILITY.tsv",sep="\t",index=False);pd.DataFrame(c4).to_csv(d/"C4_ACTIVITY_PATIENTS.tsv",sep="\t",index=False)
pd.concat(cellscore).to_csv(d/"SC_PRIMARY_SCORE_CACHE.tsv.gz",sep="\t",index=False,compression="gzip")
summ=[];loo=[];rng=np.random.default_rng(179)
for (res,sc,adj),p in P.groupby(["resource","scoring","adjustment"]):
 vals=p.rho.dropna().to_numpy();bs=np.median(vals[rng.integers(0,len(vals),(5000,len(vals)))],axis=1)
 summ.append(dict(resource=res,scoring=sc,adjustment=adj,n=len(vals),median_rho=np.median(vals),positive=int((vals>0).sum()),negative=int((vals<0).sum()),P=wil(vals),median_CI_low=np.quantile(bs,.025),median_CI_high=np.quantile(bs,.975)))
 for st in p.study.unique():
  v=p.loc[p.study!=st,"rho"].dropna();loo.append(dict(resource=res,scoring=sc,adjustment=adj,omitted_study=st,n=len(v),median_rho=v.median(),positive=int((v>0).sum())))
S=pd.DataFrame(summ);S.to_csv(d/"SC_ACTIVITY_SPP1_SUMMARY.tsv",sep="\t",index=False);pd.DataFrame(loo).to_csv(d/"SC_LEAVE_STUDY.tsv",sep="\t",index=False)
cs=[]
for (res,sc),p in pd.DataFrame(c4).groupby(["resource","scoring"]):
 v=p.C4_minus_other_score.to_numpy();cs.append(dict(resource=res,scoring=sc,n=len(v),median_delta=np.median(v),positive=int((v>0).sum()),negative=int((v<0).sum()),P=wil(v)))
pd.DataFrame(cs).to_csv(d/"C4_ACTIVITY_SUMMARY.tsv",sep="\t",index=False)
print("SC_COMPLETE\n"+S[(S.resource=="collectri")&(S.scoring=="zmean")].to_csv(sep="\t",index=False),flush=True)
# Bulk: frozen score calibration; rerank each patient bootstrap, not causal identification.
pub=Path("/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1")
mapfile=pub/"MasterMapping_ccRCC.csv";inputs.append(mapfile);mapping=pd.read_csv(mapfile,dtype=str).fillna("")
edges=[];paths=[];boots={};B=2000
def regfit(a,mm,y,cv):
 C=np.column_stack([np.ones(len(a))]+[rz(v) for v in cv]);a=rz(a);mm=rz(mm);y=rz(y)
 ba=np.linalg.lstsq(np.column_stack([C,a]),mm,rcond=None)[0][-1]
 beta=np.linalg.lstsq(np.column_stack([C,a,mm]),y,rcond=None)[0];total=np.linalg.lstsq(np.column_stack([C,a]),y,rcond=None)[0][-1]
 return ba,beta[-1],beta[-2],total
for ci,cohort in enumerate(["RC18","RC20","CPTAC","CPTAC_val"]):
 fp=r/"phase162_MTA_polyamine_part1_reassessment_20260915_v1"/("PUBLIC_PATIENTS_"+cohort+".tsv");rp=pub/"downloads_v1"/("matched_tpm_"+cohort+".csv");inputs +=[fp,rp]
 f=pd.read_csv(fp,sep="\t",index_col=0);first=pd.read_csv(rp,nrows=0).columns[0]
 rr=pd.read_csv(rp,index_col=0,usecols=[first]+sorted(set(K.gene)))
 mp=mapping[(mapping.Dataset==cohort)&(mapping.TN=="Tumor")].copy()
 mp["patient"]=mp.CommonID.where(mp.CommonID!="",mp.RNAID.str.replace(r"-T$","",regex=True)).str.replace(r"^RC(?:18|20):","",regex=True)
 expr=np.log2(rr.loc[mp.RNAID]+1);expr.index=mp.patient;expr=expr.loc[f.index];assert expr.index.is_unique
 for ri,res in enumerate(resources):
  kk=K[K.resource==res];ee=expr[kk.gene].to_numpy();sign=kk.sign.to_numpy(float)
  for si,sc in enumerate(scorings):
   f["activity"],nvar=score(ee,sign,sc);assert nvar>= (20 if res=="collectri" else 10)
   f.to_csv(d/f"BULK_PATIENTS_{cohort}_{res}_{sc}.tsv",sep="\t")
   for ai,(adj,cn) in enumerate(covsets.items()):
    ff=f[["MET:5-methylthioadenosine (MTA)","activity","SPP1"]+cn].dropna();a,mm,y=[ff.iloc[:,j].to_numpy() for j in range(3)];cv=[ff[c].to_numpy() for c in cn]
    aa,bb,cp,tt=regfit(a,mm,y,cv);assert abs(tt-cp-aa*bb)<1e-10
    base=dict(cohort=cohort,resource=res,scoring=sc,adjustment=adj,n=len(ff),variable_targets=nvar)
    for edge,xv,yv,covs in [("MTA_activity",a,mm,cv),("activity_SPP1_given_MTA",mm,y,cv+[a]),("MTA_SPP1_total",a,y,cv)]:
     rh,pv=partial(xv,yv,covs);edges.append({**base,"edge":edge,"rho":rh,"P":pv,"k":len(covs)})
    rng=np.random.default_rng(262179+ci*1000+ri*100+si*10+ai);ii=rng.integers(0,len(ff),(B,len(ff)))
    ar=[]
    for v in [a,mm,y]+cv:
     q=rankdata(v[ii],axis=1);q=(q-q.mean(1)[:,None])/q.std(1,ddof=1)[:,None];ar.append(q)
    av,mv,yv=ar[:3];cc=[np.ones_like(av)]+ar[3:]
    def batreg(Y,cols):
     des=np.stack(cols,axis=2);return np.einsum("bij,bj->bi",np.linalg.pinv(np.einsum("bni,bnj->bij",des,des)),np.einsum("bni,bn->bi",des,Y))
    ba=batreg(mv,cc+[av])[:,-1];by=batreg(yv,cc+[av,mv])[:,-1];draw=ba*by
    assert np.isfinite(draw).all()
    for j in range(10):
     ii0=ii[j];aa0,bb0,_,_=regfit(a[ii0],mm[ii0],y[ii0],[v[ii0] for v in cv]);checks.append({"check":f"boot_{cohort}_{res}_{sc}_{adj}_{j}","pass":bool(np.isclose(aa0*bb0,draw[j],atol=1e-9))})
    key="__".join([cohort,res,sc,adj]);boots[key]=draw
    paths.append({**base,"a_beta":aa,"b_beta":bb,"ab":aa*bb,"cprime":cp,"total":tt,"identity_error":tt-cp-aa*bb,"bootstrap_SE":draw.std(ddof=1),"bootstrap_low":np.quantile(draw,.025),"bootstrap_high":np.quantile(draw,.975),"bootstrap_positive_fraction":np.mean(draw>0)})
 print("BULK_COMPLETE",cohort,flush=True)
ED=pd.DataFrame(edges);PA=pd.DataFrame(paths)
for _,ix in ED.groupby(["cohort","resource","scoring","adjustment"]).groups.items():ED.loc[ix,"FDR_three_edges"]=multipletests(ED.loc[ix,"P"],method="fdr_bh")[1]
ED.to_csv(d/"BULK_EDGES_COHORT.tsv",sep="\t",index=False);PA.to_csv(d/"BULK_SHORTPATH_COHORT.tsv",sep="\t",index=False);np.savez_compressed(d/"BULK_SHORTPATH_BOOTSTRAP.npz",**boots)
def meta(v,var):
 w=1/var;mu=np.average(v,weights=w);Q=np.sum(w*(v-mu)**2);df=len(v)-1;den=w.sum()-(w*w).sum()/w.sum();tau=max(0,(Q-df)/den);ww=1/(var+tau);mu=np.average(v,weights=ww);se=1/np.sqrt(ww.sum())
 return mu,se,max(0,(Q-df)/Q)*100 if Q else 0,tau,ww/ww.sum()
em=[];pm=[];lc=[]
for (res,sc,adj,edge),p in ED.groupby(["resource","scoring","adjustment","edge"]):
 mu,se,i2,tau,w=meta(np.arctanh(p.rho.to_numpy()),1/(p.n-p.k-3).to_numpy())
 em.append(dict(resource=res,scoring=sc,adjustment=adj,edge=edge,n=int(p.n.sum()),positive_cohorts=int((p.rho>0).sum()),rho=np.tanh(mu),P=2*norm.sf(abs(mu/se)),CI_low=np.tanh(mu-1.96*se),CI_high=np.tanh(mu+1.96*se),I2=i2))
for (res,sc,adj),p in PA.groupby(["resource","scoring","adjustment"]):
 mu,se,i2,tau,w=meta(p.ab.to_numpy(),p.bootstrap_SE.to_numpy()**2)
 bs=np.sum(np.stack([boots["__".join([c,res,sc,adj])] for c in p.cohort])*w[:,None],axis=0)
 pm.append(dict(resource=res,scoring=sc,adjustment=adj,n=int(p.n.sum()),positive_cohorts=int((p.ab>0).sum()),ab=mu,SE=se,P=2*norm.sf(abs(mu/se)),CI_low=mu-1.96*se,CI_high=mu+1.96*se,I2=i2,bootstrap_low=np.quantile(bs,.025),bootstrap_high=np.quantile(bs,.975)))
 for co in p.cohort:
  pp=p[p.cohort!=co];vv,ss,_,_,_=meta(pp.ab.to_numpy(),pp.bootstrap_SE.to_numpy()**2);lc.append(dict(resource=res,scoring=sc,adjustment=adj,omitted=co,ab=vv,CI_low=vv-1.96*ss,CI_high=vv+1.96*ss))
EM=pd.DataFrame(em)
for _,ix in EM.groupby(["resource","scoring","adjustment"]).groups.items():EM.loc[ix,"FDR_three_edges"]=multipletests(EM.loc[ix,"P"],method="fdr_bh")[1]
EM.to_csv(d/"BULK_EDGES_META.tsv",sep="\t",index=False);pd.DataFrame(pm).to_csv(d/"BULK_SHORTPATH_META.tsv",sep="\t",index=False);pd.DataFrame(lc).to_csv(d/"BULK_SHORTPATH_LEAVE_COHORT.tsv",sep="\t",index=False)
pd.DataFrame(checks).to_csv(d/"BOOTSTRAP_VALIDATION.tsv",sep="\t",index=False);assert all(c["pass"] for c in checks)
manifest=[]
for f in inputs:
 h=hashlib.sha256()
 with f.open("rb") as fh:
  for block in iter(lambda:fh.read(8*1024*1024),b""):h.update(block)
 manifest.append({"path":str(f),"bytes":f.stat().st_size,"sha256":h.hexdigest()})
(d/"INPUT_MANIFEST.json").write_text(json.dumps(manifest,indent=2))
(d/"COMPUTE_STATUS.json").write_text(json.dumps({"status":"computed","bootstrap_checks":len(checks),"bootstrap_fits":len(PA)*B,"sc_patients":int(P.patient.nunique()),"input_files":len(manifest)},indent=2))
print("MAIN_EDGES\n"+EM[(EM.resource=="collectri")&(EM.scoring=="zmean")].to_csv(sep="\t",index=False),flush=True)
print("SHORT\n"+pd.DataFrame(pm).to_csv(sep="\t",index=False),flush=True)
