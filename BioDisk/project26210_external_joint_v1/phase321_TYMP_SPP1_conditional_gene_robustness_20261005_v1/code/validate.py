from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import json,hashlib,subprocess,math
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata
from scipy.linalg import qr
P=Path(__file__).resolve().parents[1];R=P.parent;O=R/'phase314_TYMP_participation_immune_exploration_20261004_v1';pr=json.loads((P/'resources/ANALYSIS_PROTOCOL.json').read_text());checks=[]
def ck(name,ok):
 checks.append(dict(check=name,pass_=bool(ok)))
 if not ok:raise AssertionError(name)
q=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(P/'logs/validate_memory_gate.log').write_text(q.stdout+q.stderr);ck('memory_gate',q.returncode==0)
ck('protocol_hash',hashlib.sha256((P/'resources/ANALYSIS_PROTOCOL.json').read_bytes()).hexdigest()==(P/'resources/ANALYSIS_PROTOCOL.sha256').read_text())
a=ad.read_h5ad(O/'inputs/OBRADOVIC_macrophage_raw.h5ad');raw=a.X.tocsr().astype(float);umi=np.asarray(raw.sum(axis=1)).ravel();yn=raw.multiply((10000/umi)[:,None]).tocsr();yn.data=np.log1p(yn.data)
formal=pd.read_csv(P/'tables/FORMAL_PATIENT_EFFECTS.tsv',sep='\t');loo=pd.read_csv(P/'tables/LEAVE_ONE_GENE_PATIENT_EFFECTS.tsv.gz',sep='\t');fv=pd.read_csv(P/'tables/FIVE_GENE_PATIENT_EFFECTS.tsv',sep='\t');su=pd.read_csv(P/'tables/FORMAL_SUMMARIES.tsv',sep='\t');ls=pd.read_csv(P/'tables/LEAVE_ONE_GENE_SUMMARIES.tsv',sep='\t');fs=pd.read_csv(P/'tables/FIVE_GENE_SUMMARIES.tsv',sep='\t');five=json.loads((P/'resources/FIVE_GENE_STRESS_LIST.json').read_text());manifest=json.loads((P/'inputs/CACHE_MANIFEST.json').read_text());flookup=formal.set_index(['condition','patient','source','target','variant']);llookup=loo.set_index(['condition','patient','source','target','removed_gene']);vlookup=fv.set_index(['condition','patient','source','target'])
for cm in manifest:
 c=np.load(P/'inputs'/cm['file']);old=np.load(O/'inputs'/cm['file'],allow_pickle=True);genes=c['genes'].tolist();idx=a.var_names.get_indexer(genes);expected=np.zeros(c['values'].shape);valididx=np.flatnonzero(idx>=0);expected[:,valididx]=yn[c['raw_rows']][:,idx[valididx]].toarray();ck(cm['file']+'_raw_gene_values',np.allclose(expected,c['values'],atol=1e-13));ck(cm['file']+'_rows',np.array_equal(c['raw_rows'],old['raw_rows']));ck(cm['file']+'_detection',np.array_equal((expected>0).sum(axis=0),c['detected']));ck(cm['file']+'_core_design',np.array_equal(c['D_core'],old['D_core']));av=set(old['available']);ck(cm['file']+'_measurement',np.array_equal(c['measured'],[g in av for g in genes]));Y=c['values'];n=Y.shape[0];sources=c['sources'];base=c['baseline'];D=c['D_core'];det=c['detected'];measured=c['measured'];condition=cm['condition'];patient=cm['patient']
 for k,g in enumerate(['TYMP','SPP1']):ck(cm['file']+'_'+g+'_raw_source',np.allclose(sources[:,k],yn[c['raw_rows'],a.var_names.get_loc(g)].toarray().ravel(),atol=1e-13))
 bases={}
 for source in ['TYMP','SPP1']:
  sg=['TYMP','SPP1'].index(source)
  for variant,dd in [('original_core',D),('conditional',np.column_stack([D,rankdata(np.round(sources[:,1-sg],12))]))]:
   Q,RR,piv=qr(dd,mode='economic',pivoting=True);rank=int(np.linalg.matrix_rank(dd));Q=Q[:,:rank];vr=rankdata(np.round(sources[:,sg],12));vr=vr-Q@(Q.T@vr);bases[(source,variant)]=(Q,vr,rank)
 def corr(source,variant,ys):
  if ys.ndim==1:ys=ys[:,None]
  Q,vr,rank=bases[(source,variant)];ry=rankdata(np.round(ys,12),axis=0);res=ry-Q@(Q.T@ry);nr=np.linalg.norm(res,axis=0);nv=np.linalg.norm(vr);out=np.einsum('i,ij->j',vr,res)/(nv*nr);out[(nr<1e-9)|(nv<1e-9)]=np.nan;return out,rank
 def eligibility(source,gs,rank):
  js=[genes.index(g) for g in gs if measured[genes.index(g)]];seen=int((det[js]>=max(3,math.ceil(.01*n))).sum());ok=n>=pr['memberships'][condition]['minimum_cells'] and all((sources[:,k]>0).sum()>=10 for k in [0,1]) and len(gs)>=5 and len(js)>=.8*len(gs) and seen>=.2*len(gs) and n-rank>=8;return js,ok
 for source,target in pr['edges']:
  gs=pr['panels'][target];jj=[genes.index(g) for g in gs if measured[genes.index(g)]];total=Y[:,jj].sum(axis=1);bs=base[:,list(pr['panels']).index(target)];ck(cm['file']+'_'+target+'_baseline',np.allclose(bs,total/len(jj),atol=1e-12))
  for variant in ['original_core','conditional']+list(pr['variants'][target]):
   model=variant if variant in ['original_core','conditional'] else 'conditional';remain=gs if variant in ['original_core','conditional'] else pr['variants'][target][variant];rank=bases[(source,model)][2];js,ok=eligibility(source,remain,rank);ys=bs if variant in ['original_core','conditional'] else Y[:,js].mean(axis=1);v,_=corr(source,model,ys);v=float(v[0]) if ok else np.nan;record=flookup.loc[(condition,patient,source,target,variant)];ck(str((condition,patient,source,target,variant)),np.isclose(v,record.rho,atol=1e-9,equal_nan=True) and bool(ok)==bool(record.eligible))
  for start in range(0,len(gs),60):
   block=gs[start:start+60];ys=[];good=[]
   for gene in block:
    j=genes.index(gene);remain=[g for g in gs if g!=gene];js,ok=eligibility(source,remain,bases[(source,'conditional')][2]);ys.append((total-Y[:,j])/(len(jj)-1) if measured[j] else total/len(jj));good.append(ok)
   vals,_=corr(source,'conditional',np.column_stack(ys))
   for gene,v,ok in zip(block,vals,good):
    v=float(v) if ok else np.nan;record=llookup.loc[(condition,patient,source,target,gene)];ck(str((condition,patient,source,target,gene)),np.isclose(v,record.rho,atol=1e-9,equal_nan=True) and bool(ok)==bool(record.eligible))
  remain=[g for g in gs if g not in five[source+'|'+target]];js,ok=eligibility(source,remain,bases[(source,'conditional')][2]);val,_=corr(source,'conditional',Y[:,js].mean(axis=1));val=float(val[0]) if ok else np.nan;ck(str((condition,patient,source,target,'joint5')),np.isclose(val,vlookup.loc[(condition,patient,source,target)].rho,atol=1e-9,equal_nan=True))
 print('VALIDATED_PATIENT',condition,patient,flush=True)
def sums(pat,summary,keys,test):
 lookup=summary.set_index(keys)
 for key,g in pat.groupby(keys,sort=False):
  z=g.rho.dropna().to_numpy();n=len(z);pos=int((z>0).sum());neg=int((z<0).sum());m=pos+neg;med=float(np.median(z)) if n else np.nan;row=lookup.loc[key];ck(str(key)+'_summary',n==row.patients and pos==row.positive and neg==row.negative and np.isclose(med,row.median_rho,atol=1e-12,equal_nan=True));pv=min(1,2*sum(math.comb(m,k) for k in range(min(pos,neg)+1))/2**m) if test and m>=10 else np.nan;ck(str(key)+'_P',np.isclose(pv,row.P,atol=1e-12,equal_nan=True))
sums(formal,su,['condition','source','target','variant'],True);sums(loo,ls,['condition','source','target','removed_gene'],False);sums(fv,fs,['condition','source','target','variant'],False)
for fam,g in su[su.variant!='original_core'].groupby('family'):
 x=g[g.P.notna()].sort_values('P');m=len(x);bh=np.minimum.accumulate((x.P.to_numpy()*m/np.arange(1,m+1))[::-1])[::-1].clip(0,1);ck(fam+'_BH',np.allclose(bh,x.FDR,atol=1e-12))
for source,target in pr['edges']:
 baseline=su[(su.condition=='L2_broad')&(su.source==source)&(su.target==target)&(su.variant=='conditional')].iloc[0].median_rho;z=ls[(ls.condition=='L2_broad')&(ls.source==source)&(ls.target==target)].copy();z['a']=np.sign(baseline)*(baseline-z.median_rho);top=z.sort_values(['a','removed_gene'],ascending=[False,True]).head(5).removed_gene.tolist();ck(source+'|'+target+'_top5',top==five[source+'|'+target])
ck('formal_count',len(su)==76);ck('baseline_directions_unchanged',pd.read_csv(P/'qc/ORIGINAL_BASELINE_REPRODUCTION.tsv',sep='\t').sign_unchanged.all());ck('no_small_n_P',su[su.patients<10].P.isna().all())
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'raw_matrix_values_checked':True,'all_formal_and_deletion_models_independent_QR':True,'details':checks},indent=2));print('FINAL_VALIDATION_PASS',len(checks),flush=True)
