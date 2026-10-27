from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import numpy as np,pandas as pd,anndata as ad,json,subprocess,hashlib,gc
from scipy.stats import rankdata,binomtest,false_discovery_control
P=Path(__file__).resolve().parents[1];R=P.parent;O=R/'phase314_TYMP_participation_immune_exploration_20261004_v1';pr=json.loads((P/'resources/ANALYSIS_PROTOCOL.json').read_text());panels=pr['panels'];allgenes=sorted(set(sum(panels.values(),[])));q=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(P/'logs/memory_gate.log').write_text(q.stdout+q.stderr);assert q.returncode==0
rawpath=O/'inputs/OBRADOVIC_macrophage_raw.h5ad';a=ad.read_h5ad(rawpath);x=a.X.tocsr().astype(float);umi=np.asarray(x.sum(axis=1)).ravel();norm=x.multiply((1e4/umi)[:,None]).tocsr();norm.data=np.log1p(norm.data);indices=a.var_names.get_indexer(allgenes);allY=np.zeros((a.n_obs,len(allgenes)));okidx=np.flatnonzero(indices>=0);allY[:,okidx]=norm[:,indices[okidx]].toarray();allDet=(allY>0)
base_effects=pd.read_csv(O/'tables/ALL_PATIENT_EFFECTS.tsv',sep='\t');oldsummary=pd.read_csv(O/'tables/ALL_EFFECT_SUMMARIES.tsv',sep='\t');full=[];los=[];coverage=[];checks=[];caches=[]
def matrix_rho(source,Y,D):
 if Y.ndim==1:Y=Y[:,None]
 Q,sv,_=np.linalg.svd(D,full_matrices=False);rk=int(np.linalg.matrix_rank(D));Q=Q[:,:rk];s=rankdata(np.round(source,12));s=s-Q@(Q.T@s);yn=rankdata(np.round(Y,12),axis=0);yn=yn-Q@(Q.T@yn);ss=np.linalg.norm(s);ys=np.linalg.norm(yn,axis=0);rr=s@yn/(ss*ys);rr[(ys<1e-9)|(ss<1e-9)]=np.nan;return rr,rk
for condition,definition in pr['memberships'].items():
 for f in sorted((O/'inputs').glob(condition+'__*.npz')):
  c=np.load(f,allow_pickle=True);rows=c['raw_rows'];n=len(rows);patient=str(a.obs.iloc[rows[0]].patient);names=c['value_names'].tolist();Y=allY[rows];det=allDet[rows].sum(axis=0);measured=np.array([g in set(c['available']) for g in allgenes]);sources=np.column_stack([c['values'][:,names.index(s)] for s in ['TYMP','SPP1']]);baseline=np.column_stack([c['values'][:,names.index(k)] for k in panels]);D=c['D_core'];outname=f.name;np.savez_compressed(P/'inputs'/outname,genes=np.array(allgenes,dtype=str),values=Y,detected=det,measured=measured,sources=sources,baseline=baseline,D_core=D,raw_rows=rows);caches.append(dict(condition=condition,patient=patient,file=outname,cells=n))
  for source,target in pr['edges']:
   sg=['TYMP','SPP1'].index(source);other=1-sg;v=sources[:,sg];Dc=np.column_stack([D,rankdata(np.round(sources[:,other],12))]);gs=panels[target];expected=len(gs);js=[allgenes.index(g) for g in gs if measured[allgenes.index(g)]];total=Y[:,js].sum(axis=1);rawscore=total/len(js);bc=baseline[:,list(panels).index(target)];err=float(np.max(np.abs(rawscore-bc)));assert err<1e-12
   baseok=n>=definition['minimum_cells'] and (v>0).sum()>=10 and (sources[:,other]>0).sum()>=10
   def valid(geneset,rank):
    jj=[allgenes.index(g) for g in geneset if measured[allgenes.index(g)]];seen=int((det[jj]>=max(3,int(np.ceil(.01*n)))).sum());good=baseok and len(geneset)>=5 and len(jj)>=.8*len(geneset) and seen>=.2*len(geneset) and n-rank>=8
    return jj,good,seen
   for label,dd in [('original_core',D),('conditional',Dc)]:
    rr,rk=matrix_rho(v,bc,dd);_,good,seen=valid(gs,rk);value=float(rr[0]) if good else np.nan
    full.append(dict(condition=condition,patient=patient,source=source,target=target,variant=label,rho=value,cells=n,design_rank=rk,genes_original=expected,genes_remaining=expected,genes_measured=len(js),genes_detected=seen,removed_genes='',eligible=bool(good)))
    if label=='original_core':
     old=base_effects[(base_effects.condition==condition)&(base_effects.patient==patient)&(base_effects.source==source)&(base_effects.target==target)&(base_effects.branch=='immune')].iloc[0];checks.append(dict(condition=condition,patient=patient,source=source,target=target,raw_score_max_error=err,old_rho=float(old.rho),rounded_rho=value,difference=value-old.rho,sign_unchanged=bool(np.sign(value)==np.sign(old.rho))))
   for label,remain in pr['variants'][target].items():
    jj,good,seen=valid(remain,np.linalg.matrix_rank(Dc));score=Y[:,jj].mean(axis=1) if jj else np.zeros(n);rr,rk=matrix_rho(v,score,Dc);full.append(dict(condition=condition,patient=patient,source=source,target=target,variant=label,rho=float(rr[0]) if good else np.nan,cells=n,design_rank=rk,genes_original=expected,genes_remaining=len(remain),genes_measured=len(jj),genes_detected=seen,removed_genes=';'.join(sorted(set(gs)-set(remain))),eligible=bool(good)))
   # Every original gene is deleted once, including measured-zero genes. Missing genes do not enter means.
   for start in range(0,len(gs),40):
    block=gs[start:start+40];vv=[];goods=[]
    for gene in block:
     j=allgenes.index(gene);remain=[g for g in gs if g!=gene];jj,good,seen=valid(remain,np.linalg.matrix_rank(Dc));vv.append((total-Y[:,j])/(len(js)-1) if measured[j] else total/len(js));goods.append(good)
    rr,rk=matrix_rho(v,np.column_stack(vv),Dc)
    for gene,value,good in zip(block,rr,goods):los.append(dict(condition=condition,patient=patient,source=source,target=target,removed_gene=gene,rho=float(value) if good else np.nan,eligible=bool(good)))
  print('PATIENT_DONE',condition,patient,n,flush=True)
 del c;gc.collect()
f=pd.DataFrame(full);lo=pd.DataFrame(los);f.to_csv(P/'tables/FORMAL_PATIENT_EFFECTS.tsv',sep='\t',index=False);lo.to_csv(P/'tables/LEAVE_ONE_GENE_PATIENT_EFFECTS.tsv.gz',sep='\t',index=False);pd.DataFrame(checks).to_csv(P/'qc/ORIGINAL_BASELINE_REPRODUCTION.tsv',sep='\t',index=False);(P/'inputs/CACHE_MANIFEST.json').write_text(json.dumps(caches,indent=2))
def summary(df,groupcols,formal):
 out=[]
 for key,z in df.groupby(groupcols,sort=False):
  z=z[z.rho.notna()];n=len(z);positive=int((z.rho>0).sum());negative=int((z.rho<0).sum());nn=positive+negative
  out.append(dict(zip(groupcols,key))|dict(patients=n,positive=positive,negative=negative,median_rho=float(z.rho.median()) if n else np.nan,min_rho=float(z.rho.min()) if n else np.nan,max_rho=float(z.rho.max()) if n else np.nan,P=float(binomtest(positive,nn,.5).pvalue) if formal and nn>=10 else np.nan))
 return pd.DataFrame(out)
su=summary(f,['condition','source','target','variant'],True);su['FDR']=np.nan;su['family']=su.variant.map(lambda v:'original224_inherited' if v=='original_core' else 'conditional20' if v=='conditional' else 'nonoverlap36');su['family_planned_n']=su.family.map({'original224_inherited':224,'conditional20':20,'nonoverlap36':36})
for family,dd in su[su.variant!='original_core'].groupby('family'):
 ids=dd[dd.P.notna()].index;su.loc[dd.index,'family_tested_n']=len(ids)
 if len(ids):su.loc[ids,'FDR']=false_discovery_control(su.loc[ids,'P'].to_numpy())
for idx,z in su[su.variant=='original_core'].iterrows():
 orig=oldsummary[(oldsummary.condition==z.condition)&(oldsummary.source==z.source)&(oldsummary.target==z.target)&(oldsummary.branch=='immune')].iloc[0];assert z.P==orig.P;su.loc[idx,'FDR']=orig.FDR;su.loc[idx,'family_tested_n']=orig.family_tested_n
su.to_csv(P/'tables/FORMAL_SUMMARIES.tsv',sep='\t',index=False);ls=summary(lo,['condition','source','target','removed_gene'],False);ls.to_csv(P/'tables/LEAVE_ONE_GENE_SUMMARIES.tsv',sep='\t',index=False)
# Select a single broad-derived five-gene stress list; carry unchanged into the strict membership.
five={};selection=[]
for source,target in pr['edges']:
 b=su[(su.condition=='L2_broad')&(su.source==source)&(su.target==target)&(su.variant=='conditional')].iloc[0];z=ls[(ls.condition=='L2_broad')&(ls.source==source)&(ls.target==target)].copy();z['signed_attenuation']=np.sign(b.median_rho)*(b.median_rho-z.median_rho);z=z.sort_values(['signed_attenuation','removed_gene'],ascending=[False,True]);chosen=z.head(5).removed_gene.tolist();five[source+'|'+target]=chosen
 for row in z.to_dict('records'):selection.append(row|dict(selected_for_joint5=row['removed_gene'] in chosen,baseline_rho=b.median_rho))
(P/'resources/FIVE_GENE_STRESS_LIST.json').write_text(json.dumps(five,indent=2));pd.DataFrame(selection).to_csv(P/'tables/FIVE_GENE_SELECTION_AUDIT.tsv',sep='\t',index=False)
rows=[]
for cm in caches:
 c=np.load(P/'inputs'/cm['file']);Y=c['values'];det=c['detected'];measured=c['measured'];genes=c['genes'].tolist();n=cm['cells'];D=c['D_core'];sources=c['sources']
 for source,target in pr['edges']:
  sg=['TYMP','SPP1'].index(source);v=sources[:,sg];Dc=np.column_stack([D,rankdata(np.round(sources[:,1-sg],12))]);drop=five[source+'|'+target];remain=[g for g in panels[target] if g not in drop];jj=[genes.index(g) for g in remain if measured[genes.index(g)]];score=Y[:,jj].mean(axis=1);rr,rk=matrix_rho(v,score,Dc);good=n>=pr['memberships'][cm['condition']]['minimum_cells'] and (v>0).sum()>=10 and (sources[:,1-sg]>0).sum()>=10 and len(remain)>=5 and len(jj)>=.8*len(remain) and (det[jj]>=max(3,int(np.ceil(.01*n)))).sum()>=.2*len(remain) and n-rk>=8
  rows.append(dict(condition=cm['condition'],patient=cm['patient'],source=source,target=target,variant='joint_five_gene_deletion',rho=float(rr[0]) if good else np.nan,removed_genes=';'.join(drop),patients_are_independent=True))
fv=pd.DataFrame(rows);fv.to_csv(P/'tables/FIVE_GENE_PATIENT_EFFECTS.tsv',sep='\t',index=False);fs=summary(fv,['condition','source','target','variant'],False);fs.to_csv(P/'tables/FIVE_GENE_SUMMARIES.tsv',sep='\t',index=False)
(P/'qc/ANALYSIS_COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','formal_patient_rows':len(f),'formal_summary_rows':len(su),'single_gene_patient_rows':len(lo),'single_gene_summary_rows':len(ls),'five_gene_rows':len(fv)},indent=2));print('CONDITIONAL_SUMMARY');print(su[su.variant=='conditional'][['condition','source','target','patients','positive','negative','median_rho','P','FDR']].to_string(index=False));print('ANALYSIS_COMPLETE',flush=True)
