from pathlib import Path
import pandas as pd,numpy as np,json,sys,subprocess,hashlib,re,datetime
from scipy.stats import rankdata,t,norm,chi2,spearmanr
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase160_external_bridge_investigation_20260915_v1';Q=Path('/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1')
g=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(O/'ANALYSIS_MEMORY.log').write_text(g.stdout+g.stderr);assert g.returncode==0
inv=json.loads((O/'INPUT_INVENTORY.json').read_text());mapping=pd.read_csv(Q/'MasterMapping_ccRCC.csv',dtype=str).fillna('')
gene_primary=['MSMO1','LDHA','SLC16A3','CYP27A1'];gene_context=['HK2','LDHB']
myeloid=['LST1','AIF1','CSF1R','CD68','FCER1G']
met_primary=['lactate','cholesterol'];met_context=['5-methylthioadenosine (MTA)']
met_aux=['pyruvate','desmosterol','7-hydroxycholesterol (alpha or beta)','cholesterol sulfate']
adjustments={'raw':[],'myeloid':['MYELOID_PROXY'],'myeloid_hypoxia':['MYELOID_PROXY','HYPOXIA_PROXY']}
links=[('LDHA','lactate'),('LDHB','lactate'),('SLC16A3','lactate'),('HK2','lactate'),('MSMO1','cholesterol'),('MSMO1','desmosterol'),('CYP27A1','cholesterol')]
plan=dict(gene_primary=gene_primary,gene_context=gene_context,met_primary=met_primary,met_context=met_context,met_aux=met_aux,myeloid_proxy=myeloid,adjustments=adjustments,gene_metabolite_links=links,unit='one patient, Tumor only; exact mapped patient IDs; any repeated tumor specimen per patient aggregated as median log2TPM+1 / median provided metabolite abundance',missing='no imputation; per-test complete-case; >=10 independent patients for inferential meta component',scope='External bulk RNA and observed reference metabolomics; not inferred metabolites, not isolated macrophages or C4; no client data; no figures',meta='Fisher z, approximate variance 1/(n-k-3); fixed plus DerSimonian-Laird random effects normal intervals; random effects primary; Q and I2 reported; no pooled metabolite concentrations',fdr='Within gene_primary4/gene_context2/met_primary2/met_context1/met_aux4 families, separately for edge, cohort/meta, and adjustment; gene-metabolite links separate. TYMP-SPP1 one planned relation unadjusted for multiplicity.')
(O/'PLAN_BEFORE_ASSOCIATIONS.json').write_text(json.dumps(plan,indent=2))
gmt=P/'part2_rebuild_v1/06_functional_programs/run_v1/resources/MSigDB_Hallmark_2020.gmt'
hyp=[]
for line in gmt.read_text().splitlines():
 parts=line.split('\t')
 if parts[0].lower()=='hypoxia':hyp=parts[2:]
assert len(hyp)>10
excluded=set(['TYMP','SPP1']+gene_primary+gene_context+myeloid)
hyp=[x for x in hyp if x not in excluded]
def rz(x):
 r=rankdata(np.asarray(x,float));sd=r.std(ddof=1)
 return (r-r.mean())/sd if sd>0 else r*0
def corr(x,y,C):
 a=x-C@np.linalg.lstsq(C,x,rcond=None)[0];b=y-C@np.linalg.lstsq(C,y,rcond=None)[0]
 k=int(np.linalg.matrix_rank(C)-1);n=len(x)
 if a.std()<1e-12 or b.std()<1e-12:return np.nan,np.nan,k
 r=float(np.corrcoef(a,b)[0,1]);df=n-k-2
 return r,float(2*t.sf(abs(r)*np.sqrt(df/max(1-r*r,1e-15)),df)),k
cohort_frames={};meta_frames={};sample=[];missing=[];inputmanifest=[];allpatient=[]
for it in inv:
 name=it['cohort'];rawrna=pd.read_csv(it['rna_file'],index_col=0);met=pd.read_csv(it['met_file'],index_col=0);mm=mapping[mapping.Dataset==name].copy()
 assert not mm.RNAID.duplicated().any() and set(mm.RNAID)==set(rawrna.index.astype(str))
 assert set(mm.MetabID)==set(met.index.astype(str))
 mm=mm[mm.TN=='Tumor'].copy()
 mm['patient']=mm.CommonID.where(mm.CommonID!='',mm.RNAID.str.replace(r'-T$','',regex=True)).str.replace(r'^RC(?:18|20):','',regex=True)
 assert mm.patient.ne('').all()
 genes=['TYMP','SPP1']+gene_primary+gene_context+myeloid
 h=[x for x in hyp if x in rawrna.columns]
 assert len(h)>=10 and all(x in rawrna for x in genes)
 numeric=rawrna.loc[mm.RNAID,sorted(set(genes+h))].apply(pd.to_numeric,errors='coerce')
 assert (numeric.fillna(0)>=0).all().all()
 exp=np.log2(numeric+1);exp.index=mm.patient.to_numpy();exp=exp.groupby(level=0,sort=False).median()
 mt=met.loc[mm.MetabID].apply(pd.to_numeric,errors='coerce');mt.index=mm.patient.to_numpy();mt=mt.groupby(level=0,sort=False).median()
 assert exp.index.equals(mt.index)
 f=exp[genes].copy()
 for score,gg in [('MYELOID_PROXY',myeloid),('HYPOXIA_PROXY',h)]:
  sub=exp[gg];sd=sub.std(ddof=1);valid=sd[sd>0].index;assert len(valid)>=3
  f[score]=((sub[valid]-sub[valid].mean())/sd[valid]).mean(1)
 for key in met_primary+met_context+met_aux:
  if key not in mt:missing.append(dict(cohort=name,metabolite=key,reason='not_measured_in_reference_file'))
  else:f['MET:'+key]=mt[key]
 f.replace([np.inf,-np.inf],np.nan,inplace=True);cohort_frames[name]=f
 f.to_csv(O/('PATIENT_VALUES_'+name+'.tsv'),sep='\t')
 sample.append(dict(cohort=name,tumor_specimens=len(mm),independent_mapped_patients=len(f),hypoxia_genes=len(h),gene_expression_rows=len(exp)))
 for pat in f.index:allpatient.append(dict(cohort=name,patient=pat))
 for fn in [it['rna_file'],it['met_file']]:
  p=Path(fn);inputmanifest.append(dict(path=fn,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
pd.DataFrame(sample).to_csv(O/'SAMPLE_AUDIT.tsv',sep='\t',index=False)
ids=pd.DataFrame(allpatient);ids.to_csv(O/'PATIENT_ID_AUDIT.tsv',sep='\t',index=False)
assert not ids.patient.duplicated().any(),'Cross-cohort mapped patient duplicates require resolution'
pd.DataFrame(missing).to_csv(O/'UNMEASURED_METABOLITES.tsv',sep='\t',index=False)
rows=[];direct=[];linkrows=[];checks=[]
for name,f in cohort_frames.items():
 candidates=[(x,'gene_primary') for x in gene_primary]+[(x,'gene_context') for x in gene_context]+[('MET:'+x,'met_primary') for x in met_primary if 'MET:'+x in f]+[('MET:'+x,'met_context') for x in met_context if 'MET:'+x in f]+[('MET:'+x,'met_aux') for x in met_aux if 'MET:'+x in f]
 for candidate,family in candidates:
  for adj,cn in adjustments.items():
   d=f[['TYMP','SPP1',candidate]+cn].dropna();n=len(d)
   if n<10:continue
   A=rz(d.TYMP);M=rz(d[candidate]);Y=rz(d.SPP1);C=np.column_stack([np.ones(n)]+[rz(d[x]) for x in cn])
   a,pa,ka=corr(A,M,C);xy,pxy,kxy=corr(M,Y,C);b,pb,kb=corr(M,Y,np.column_stack([C,A]))
   if not np.isfinite([a,xy,b]).all():continue
   ca=np.linalg.lstsq(np.column_stack([C,A]),M,rcond=None)[0][-1];cb=np.linalg.lstsq(np.column_stack([C,A,M]),Y,rcond=None)[0];ct=np.linalg.lstsq(np.column_stack([C,A]),Y,rcond=None)[0][-1]
   row=dict(cohort=name,candidate=candidate,family=family,adjustment=adj,n=n,a_rho=a,a_P=pa,a_k=ka,XY_rho=xy,XY_P=pxy,XY_k=kxy,b_rho=b,b_P=pb,b_k=kb,indirect_ab=ca*cb[-1],direct_cprime=cb[-2],total_c=ct,identity_error=ct-cb[-2]-ca*cb[-1])
   rows.append(row);checks.append(dict(check=name+candidate+adj+'_identity',passed=abs(row['identity_error'])<1e-10))
   if adj=='raw':
    s,p=spearmanr(d[candidate],d.SPP1);checks.append(dict(check=name+candidate+'_scipy',passed=abs(s-xy)<1e-10 and abs(p-pxy)<1e-10))
 for adj,cn in adjustments.items():
  d=f[['TYMP','SPP1']+cn].dropna();C=np.column_stack([np.ones(len(d))]+[rz(d[x]) for x in cn]);rho,p,k=corr(rz(d.TYMP),rz(d.SPP1),C)
  direct.append(dict(cohort=name,candidate='TYMP--SPP1',family='direct',adjustment=adj,n=len(d),rho=rho,P=p,k=k))
 for gene,met in links:
  if 'MET:'+met not in f:continue
  for adj,cn in adjustments.items():
   d=f[[gene,'MET:'+met]+cn].dropna();C=np.column_stack([np.ones(len(d))]+[rz(d[x]) for x in cn]);rho,p,k=corr(rz(d[gene]),rz(d['MET:'+met]),C)
   linkrows.append(dict(cohort=name,candidate=gene+'--'+met,family='gene_metabolite',adjustment=adj,n=len(d),rho=rho,P=p,k=k))
R=pd.DataFrame(rows)
for _,ix in R.groupby(['cohort','family','adjustment']).groups.items():
 for e in ['a','XY','b']:R.loc[ix,e+'_FDR']=multipletests(R.loc[ix,e+'_P'],method='fdr_bh')[1]
R.to_csv(O/'COHORT_BRIDGE_RESULTS.tsv',sep='\t',index=False)
DR=pd.DataFrame(direct);DR.to_csv(O/'COHORT_DIRECT_TYMP_SPP1.tsv',sep='\t',index=False)
LR=pd.DataFrame(linkrows)
for _,ix in LR.groupby(['cohort','adjustment']).groups.items():LR.loc[ix,'FDR']=multipletests(LR.loc[ix,'P'],method='fdr_bh')[1]
LR.to_csv(O/'COHORT_GENE_METABOLITE_LINKS.tsv',sep='\t',index=False)
def meta(g,rcol,pcol,kcol):
 n=g.n.to_numpy(float);k=g[kcol].to_numpy(float);r=np.clip(g[rcol].to_numpy(float),-.999999,.999999);z=np.arctanh(r);w=n-k-3
 assert (w>0).all()
 mu=(w*z).sum()/w.sum();se=1/np.sqrt(w.sum());q=(w*(z-mu)**2).sum();df=len(g)-1;c=w.sum()-(w*w).sum()/w.sum();tau=max(0,(q-df)/c) if c>0 else 0
 wr=1/(1/w+tau);mr=(wr*z).sum()/wr.sum();sr=1/np.sqrt(wr.sum())
 return dict(cohorts=len(g),total_n=int(n.sum()),positive_cohorts=int((r>0).sum()),negative_cohorts=int((r<0).sum()),fixed_rho=np.tanh(mu),fixed_P=2*norm.sf(abs(mu/se)),random_rho=np.tanh(mr),random_P=2*norm.sf(abs(mr/sr)),random_low=np.tanh(mr-1.96*sr),random_high=np.tanh(mr+1.96*sr),Q=q,Q_P=chi2.sf(q,df) if df else np.nan,I2=max(0,(q-df)/q)*100 if q else 0,tau2=tau)
meta_rows=[]
for (candidate,family,adj),g in R.groupby(['candidate','family','adjustment']):
 for e in ['a','XY','b']:meta_rows.append(dict(candidate=candidate,family=family,adjustment=adj,edge=e,**meta(g,e+'_rho',e+'_P',e+'_k')))
MR=pd.DataFrame(meta_rows)
for _,ix in MR.groupby(['family','adjustment','edge']).groups.items():
 for m in ['fixed','random']:MR.loc[ix,m+'_FDR']=multipletests(MR.loc[ix,m+'_P'],method='fdr_bh')[1]
MR.to_csv(O/'META_BRIDGE_RESULTS.tsv',sep='\t',index=False)
othermeta=[]
for frame in [DR,LR]:
 for (candidate,family,adj),g in frame.groupby(['candidate','family','adjustment']):
  othermeta.append(dict(candidate=candidate,family=family,adjustment=adj,**meta(g,'rho','P','k')))
OM=pd.DataFrame(othermeta);OM['random_FDR']=np.nan;OM['fixed_FDR']=np.nan
for (fam,adj),ix in OM.groupby(['family','adjustment']).groups.items():
 for m in ['fixed','random']:OM.loc[ix,m+'_FDR']=multipletests(OM.loc[ix,m+'_P'],method='fdr_bh')[1] if fam!='direct' else OM.loc[ix,m+'_P']
OM.to_csv(O/'META_DIRECT_AND_GENE_METABOLITE.tsv',sep='\t',index=False)
# Leave-one-cohort meta (direction/stability, no cohort removal selected by outcome).
lo=[]
for (candidate,family,adj),g in R.groupby(['candidate','family','adjustment']):
 if len(g)<3:continue
 for omitted in g.cohort:
  h=g[g.cohort!=omitted]
  for e in ['a','XY','b']:lo.append(dict(candidate=candidate,family=family,adjustment=adj,omitted=omitted,edge=e,**meta(h,e+'_rho',e+'_P',e+'_k')))
pd.DataFrame(lo).to_csv(O/'LEAVE_ONE_COHORT_META.tsv',sep='\t',index=False)
V=pd.DataFrame(checks);assert V.passed.all();V.to_csv(O/'VALIDATION.tsv',sep='\t',index=False)
inputmanifest += [dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in [Q/'MasterMapping_ccRCC.csv',gmt]]
(O/'SOURCE_MANIFEST.json').write_text(json.dumps(inputmanifest,indent=2))
status=dict(status='COMPUTATION_COMPLETE_REPORT_PENDING',cohorts=4,patients=len(ids),tests=len(R),checks=len(V),passed=int(V.passed.sum()),figures=False,client_data=False,metabolites='observed reference data, no imputation',time=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print('STATUS',json.dumps(status))
print('META PRIMARY\n'+MR[MR.family.isin(['gene_primary','met_primary','met_context'])][['candidate','adjustment','edge','cohorts','total_n','positive_cohorts','random_rho','random_P','random_FDR','I2']].to_string(index=False))
print('DIRECT AND GENE-METABOLITE\n'+OM[OM.adjustment=='raw'][['candidate','cohorts','total_n','positive_cohorts','random_rho','random_P','random_FDR','I2']].to_string(index=False))

