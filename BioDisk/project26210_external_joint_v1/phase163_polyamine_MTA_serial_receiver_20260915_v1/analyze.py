from pathlib import Path
import pandas as pd,numpy as np,json,subprocess,sys,datetime,hashlib
from scipy.stats import rankdata,norm,chi2,wilcoxon,spearmanr
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');S=P/'phase162_MTA_polyamine_part1_reassessment_20260915_v1';O=P/'phase163_polyamine_MTA_serial_receiver_20260915_v1';O.mkdir(exist_ok=True)
gate=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(O/'MEMORY_GATE.log').write_text(gate.stdout+gate.stderr);assert gate.returncode==0
sys.path.insert(0,str(P/'phase157_ordered_bridge_reassessment_20260915_v1'));from bridge_stats import rankz,corrp
Xs=['SRM','SMS','production_core'];adjs={'raw':[],'MTAP':['MTAP'],'MTAP_myeloid':['MTAP','MYELOID_PROXY']};B=2000
plan=dict(retained='Phase162 retained unchanged; production SRM/SMS + MTAP clearance background; C4 localization, whole macrophage functional layer',candidates=Xs,serial='TYMP(A) -> candidate X -> measured MTA(Z) -> bulk SPP1(Y)',adjustments=adjs,bootstrap=B,bootstrap_method='Patient resampling within each cohort; rerank each draw; fixed previously defined production RNA score; percentile CI of a*d*b; descriptive association decomposition, no causal mediation proportion',primary='MTAP_myeloid',fdr='3 X candidates per edge/cohort/adjustment; separate pooled family; no significance requirement on total effect',receptor_panel=['ADORA2B','ADORA2A','STAT3','SOCS3','IL10','VEGFA','ARG1'],receptor_scope='Frozen tumor C4 vs other macrophage matched patients >=20 cells each; RNA availability and expression only, not phosphorylation or validated MTA response',figures=False,client_data_in_models=False)
(O/'PLAN_BEFORE_FIT.json').write_text(json.dumps(plan,indent=2))
agreement='\n\n## 2026-09-15 用户确认保留MTA多胺轴并继续\n用户确认保留Phase162部分并继续分析。保留第一部分SRM/多胺状态与公共SMS—MTA/SPP1分段证据，MTAP单列清除及遗传背景，SLC16A3并行。不得因后续串联模型不支持而删除既有分段事实；不得把分段相关拼成已确认因果中介。整体巨噬与C4双层分工、无图及客户实验仅作背景继续适用。\n'
p=P/'AGENTS.md'
if '用户确认保留MTA多胺轴并继续' not in p.read_text():p.write_text(p.read_text()+agreement,encoding='utf-8')
def fit(a,x,z,y,c):
 a=rankz(a);x=rankz(x);z=rankz(z);y=rankz(y);C=np.column_stack([np.ones(len(a))]+[rankz(v) for v in c])
 aa=np.linalg.lstsq(np.column_stack([C,a]),x,rcond=None)[0][-1]
 mz=np.linalg.lstsq(np.column_stack([C,a,x]),z,rcond=None)[0]
 yy=np.linalg.lstsq(np.column_stack([C,a,x,z]),y,rcond=None)[0]
 total=np.linalg.lstsq(np.column_stack([C,a]),y,rcond=None)[0][-1]
 paths=dict(a=aa,d=mz[-1],a_to_MTA_direct=mz[-2],b=yy[-1],X_to_SPP1_direct=yy[-2],A_to_SPP1_direct=yy[-3],serial_adb=aa*mz[-1]*yy[-1],X_direct_path=aa*yy[-2],MTA_direct_path=mz[-2]*yy[-1],total=total)
 paths['identity_error']=total-paths['A_to_SPP1_direct']-paths['X_direct_path']-paths['MTA_direct_path']-paths['serial_adb']
 edges=[]
 for e,v,w,cc in [('a',a,x,C),('d',x,z,np.column_stack([C,a])),('b',z,y,np.column_stack([C,a,x])),('X_direct',x,y,np.column_stack([C,a,z])),('MTA_without_X',z,y,np.column_stack([C,a]))]:
  rho,p=corrp(v,w,cc);edges.append(dict(edge=e,rho=rho,P=p,k=np.linalg.matrix_rank(cc)-1))
 return paths,edges
def boot(a,x,z,y,c,seed):
 rng=np.random.default_rng(seed);n=len(a);ii=rng.integers(0,n,(B,n))
 vals=[np.asarray(a),np.asarray(x),np.asarray(z),np.asarray(y)]+[np.asarray(v) for v in c]
 ar=[rankdata(v[ii],axis=1) for v in vals]
 ar=[(v-v.mean(1)[:,None])/v.std(1,ddof=1)[:,None] for v in ar]
 aa,xx,zz,yy=ar[:4];cov=[np.ones((B,n))]+ar[4:]
 def reg(response,preds):
  design=np.stack(preds,axis=2);gram=np.einsum('bni,bnj->bij',design,design);rhs=np.einsum('bni,bn->bi',design,response)
  return np.einsum('bij,bj->bi',np.linalg.pinv(gram),rhs)
 ac=reg(xx,cov+[aa])[:,-1];zc=reg(zz,cov+[aa,xx]);yc=reg(yy,cov+[aa,xx,zz])
 prod=ac*zc[:,-1]*yc[:,-1]
 return prod
rows=[];edge_rows=[];draws={};files=[]
for ci,cohort in enumerate(['RC18','RC20','CPTAC','CPTAC_val']):
 file=S/('PUBLIC_PATIENTS_'+cohort+'.tsv');files.append(file);f=pd.read_csv(file,sep='\t',index_col=0)
 for xi,x in enumerate(Xs):
  for ai,(adj,cn) in enumerate(adjs.items()):
   d=f[['TYMP',x,'MET:5-methylthioadenosine (MTA)','SPP1']+cn].dropna()
   a=d.TYMP;xv=d[x];z=d['MET:5-methylthioadenosine (MTA)'];y=d.SPP1;c=[d[v] for v in cn]
   result,edges=fit(a,xv,z,y,c);bb=boot(a,xv,z,y,c,262163+ci*100+xi*10+ai);valid=bb[np.isfinite(bb)];lo,hi=np.quantile(valid,[.025,.975])
   rows.append(dict(cohort=cohort,candidate=x,adjustment=adj,n=len(d),**result,bootstrap_valid=len(valid),serial_low=lo,serial_high=hi,bootstrap_positive_fraction=(valid>0).mean()))
   for e in edges:edge_rows.append(dict(cohort=cohort,candidate=x,adjustment=adj,n=len(d),**e))
   draws[cohort+'__'+x+'__'+adj]=bb
R=pd.DataFrame(rows);E=pd.DataFrame(edge_rows)
for _,ix in E.groupby(['cohort','adjustment','edge']).groups.items():E.loc[ix,'FDR']=multipletests(E.loc[ix,'P'],method='fdr_bh')[1]
R.to_csv(O/'COHORT_SERIAL_DECOMPOSITION.tsv',sep='\t',index=False);E.to_csv(O/'COHORT_CONDITIONAL_EDGES.tsv',sep='\t',index=False);np.savez_compressed(O/'PATIENT_BOOTSTRAP_SERIAL_DRAWS.npz',**draws)
def mc(d):
 z=np.arctanh(d.rho.to_numpy());w=(d.n-d.k-3).to_numpy();m=np.average(z,weights=w);q=np.sum(w*(z-m)**2);tau=max(0,(q-len(d)+1)/(w.sum()-(w*w).sum()/w.sum()));wr=1/(1/w+tau);mu=np.average(z,weights=wr);se=1/np.sqrt(wr.sum())
 return dict(cohorts=len(d),total_n=int(d.n.sum()),positive_cohorts=int((d.rho>0).sum()),random_rho=np.tanh(mu),random_P=2*norm.sf(abs(mu/se)),random_low=np.tanh(mu-1.96*se),random_high=np.tanh(mu+1.96*se),I2=max(0,(q-len(d)+1)/q)*100 if q else 0)
mm=[];loo=[]
for (x,adj,e),d in E.groupby(['candidate','adjustment','edge']):
 mm.append(dict(candidate=x,adjustment=adj,edge=e,**mc(d)))
 for omit in d.cohort:loo.append(dict(candidate=x,adjustment=adj,edge=e,omitted=omit,**mc(d[d.cohort!=omit])))
M=pd.DataFrame(mm)
for _,ix in M.groupby(['adjustment','edge']).groups.items():M.loc[ix,'random_FDR']=multipletests(M.loc[ix,'random_P'],method='fdr_bh')[1]
M.to_csv(O/'META_CONDITIONAL_EDGES.tsv',sep='\t',index=False);pd.DataFrame(loo).to_csv(O/'LEAVE_ONE_COHORT_EDGES.tsv',sep='\t',index=False)
# Public mechanistic receptor panel: no assignment from receptor expression alone.
G=pd.read_csv(P/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_INDEX.tsv',sep='\t');Z=np.load(P/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_PSEUDOBULK.npz',allow_pickle=True);genes=list(Z['genes'])
rr=[]
for pat,gg in G[G.source_tissue.str.lower()=='tumor'].groupby('patient'):
 c=gg[gg.role=='C4'];o=gg[gg.role=='OtherMacrophage']
 if c.cells.sum()<20 or o.cells.sum()<20:continue
 ci=c.group_index.to_numpy(int);oi=o.group_index.to_numpy(int);cc=Z['counts'][ci].sum(0);oc=Z['counts'][oi].sum(0)
 for gene in plan['receptor_panel']:
  j=genes.index(gene);ce=np.log2(cc[j]/cc.sum()*1e6+1);oe=np.log2(oc[j]/oc.sum()*1e6+1)
  rr.append(dict(patient=pat,gene=gene,C4_cells=int(c.cells.sum()),other_cells=int(o.cells.sum()),C4_detection=Z['detected'][ci,j].sum()/c.cells.sum(),other_detection=Z['detected'][oi,j].sum()/o.cells.sum(),C4_log2CPM=ce,other_log2CPM=oe,delta=ce-oe))
RP=pd.DataFrame(rr);RP.to_csv(O/'MTA_RECEIVER_PANEL_PATIENTS.tsv',sep='\t',index=False)
rs=[]
for gene,d in RP.groupby('gene'):
 v=d.delta;rs.append(dict(gene=gene,n_patients=len(d),C4_detection_median=d.C4_detection.median(),C4_ge10pct=int((d.C4_detection>=.1).sum()),other_detection_median=d.other_detection.median(),other_ge10pct=int((d.other_detection>=.1).sum()),C4_higher=int((v>0).sum()),C4_lower=int((v<0).sum()),median_delta=v.median(),P=wilcoxon(v).pvalue if (v!=0).any() else 1.))
RS=pd.DataFrame(rs);RS['FDR']=multipletests(RS.P,method='fdr_bh')[1];RS.to_csv(O/'MTA_RECEIVER_PANEL_SUMMARY.tsv',sep='\t',index=False)
checks=[]
for i,r in R.iterrows():checks.append(dict(check='path_identity_'+str(i),passed=abs(r.identity_error)<1e-10 and r.bootstrap_valid>=1900))
for cohort,g in E.groupby('cohort'):
 f=pd.read_csv(S/('PUBLIC_PATIENTS_'+cohort+'.tsv'),sep='\t',index_col=0)
 for i,r in g.iterrows():
  cov=adjs[r.adjustment];a='TYMP';x=r.candidate;z='MET:5-methylthioadenosine (MTA)';y='SPP1'
  cols={'a':[a,x]+cov,'d':[x,z,a]+cov,'b':[z,y,a,x]+cov,'X_direct':[x,y,a,z]+cov,'MTA_without_X':[z,y,a]+cov}[r.edge]
  f2=f[cols].dropna();pr=np.linalg.inv(np.corrcoef([rankdata(f2[v]) for v in cols]));rho=-pr[0,1]/np.sqrt(pr[0,0]*pr[1,1])
  checks.append(dict(check='independent_precision_'+str(i),passed=abs(rho-r.rho)<1e-9))
V=pd.DataFrame(checks);V.to_csv(O/'VALIDATION.tsv',sep='\t',index=False);assert V.passed.all()
status=dict(status='COMPUTATION_COMPLETE_REPORT_PENDING',public_patients=255,cohorts=4,candidates=3,serial_models=len(R),bootstrap_draws_per_model=B,checks=len(V),passed=int(V.passed.sum()),figures=False,client_data_in_models=False,time=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
(O/'SOURCE_MANIFEST.json').write_text(json.dumps([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files+[P/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_PSEUDOBULK.npz',P/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_INDEX.tsv']],indent=2))
print('META PRIMARY\n'+M[M.adjustment=='MTAP_myeloid'].to_string(index=False))
print('SERIAL PRIMARY\n'+R[R.adjustment=='MTAP_myeloid'][['cohort','candidate','a','d','b','serial_adb','serial_low','serial_high','X_to_SPP1_direct']].to_string(index=False))
print('RECEIVER\n'+RS.to_string(index=False));print('STATUS',json.dumps(status))

