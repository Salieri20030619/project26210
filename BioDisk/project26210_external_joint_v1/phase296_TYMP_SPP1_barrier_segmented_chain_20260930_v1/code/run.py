from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,subprocess,importlib.util
from scipy import stats
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=lambda z:z.item() if hasattr(z,'item') else str(z)))
def gate(stage):
 with (P/'logs/memory_gate.log').open('a') as f:
  f.write(stage+'\n');f.flush();subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,stderr=f,check=True)
gate('freeze_load')
paths=[R/'phase292_barrier_execution_and_external_feasibility_20260930_v1/tables/PATIENT_INPUTS.tsv',R/'phase280_TYMP_source_niche_function_20260928_v1/tables/PAIRED_SOURCE_PATIENT_INPUTS.tsv',R/'phase280_TYMP_source_niche_function_20260928_v1/code/source.py']
dump(P/'qc/INPUT_FREEZE.json',dict(protocol_sha256=sha(P/'PROTOCOL_FROZEN.md'),inputs=[dict(path=str(p),sha256=sha(p)) for p in paths]))
a=pd.read_csv(paths[0],sep='\t');b=pd.read_csv(paths[1],sep='\t')
d=a[a.eligible].merge(b[['patient','study','library','C4_fraction']],on=['patient','study','library'],validate='one_to_one')
assert len(d)==13 and d.patient.is_unique
d.to_csv(P/'tables/FROZEN_PATIENT_INPUTS.tsv',sep='\t',index=False)
d.groupby('study').size().rename('patients').to_csv(P/'tables/STUDY_COMPOSITION.tsv',sep='\t')
spec=importlib.util.spec_from_file_location('prior_source',paths[2]);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
Y='collagen_maturation';M='macro_SPP1';F='C4_fraction';sources=['malignant_TYMP','macro_TYMP']
edges=[]
for x in sources:edges.append(('core_edges',x,M,[]))
edges.append(('core_edges',M,Y,[]))
for x in sources:edges.append(('core_edges',x,Y,[]))
for x in sources:
 edges.append(('conditional_edges',x,Y,[M]));edges.append(('conditional_edges',M,Y,[x]))
for x in sources:edges.append(('state_proxy',x,F,[]))
edges.append(('state_proxy',F,Y,[]))
main=list(edges)
for fam,x,y,cov in main:
 if y==Y:edges.append(('depth_sensitivity',x,y,cov+['log_receiver_depth']))
assert len(edges)==20
def fit(z,x,y,cov,inference):
 ext=z[cov].to_numpy() if cov else None
 n=len(z);D=np.column_stack([np.ones(n)]+[(z.study.to_numpy()==st).astype(float) for st in sorted(z.study.unique())[1:]]+[stats.rankdata(z[c]) for c in cov])
 rk=np.linalg.matrix_rank(D);df=n-rk-1
 do=inference and n>=10 and df>=5
 o=mod.cross(z[x],z[y],z.study,ext,seed=296,B=99999 if do else 0)
 if not do:o['P']=np.nan;o['method']='DESCRIPTIVE_OR_DF_INSUFFICIENT'
 return o
out=[];leave=[]
for i,(family,x,y,cov) in enumerate(edges):
 gate('edge_'+str(i))
 z=d.dropna(subset=[x,y]+cov);res=fit(z,x,y,cov,True);eid='E'+str(i+1).zfill(2)
 out.append(dict(edge_id=eid,family=family,source=x,endpoint=y,controls='study'+(';'+';'.join(cov) if cov else ''),**res))
 for typ in ['patient','study']:
  for value in sorted(z[typ].unique()):
   w=z[z[typ]!=value];q=fit(w,x,y,cov,False)
   leave.append(dict(edge_id=eid,omission_type=typ,omitted=value,n=q['n'],rho=q['rho']))
o=pd.DataFrame(out);o['FDR']=np.nan
for family,ix in o.groupby('family').groups.items():
 vals=o.loc[ix,'P'];ok=vals.notna()
 # Prespecified full family retained; unavailable tests assigned1 for adjustment then restoredNA.
 q=stats.false_discovery_control(vals.fillna(1).to_numpy())
 o.loc[ix,'FDR']=np.where(ok,q,np.nan);o.loc[ix,'family_size']=len(ix)
o['evidence']=np.where(o.P.isna(),'NOT_ESTIMABLE',np.where(o.FDR.lt(.05),'FDR_SIGNIFICANT',np.where(o.P.lt(.05),'NOMINAL','CURRENT_DATA_NOT_SUPPORTED')))
o.to_csv(P/'tables/SEGMENTED_ASSOCIATIONS.tsv',sep='\t',index=False)
pd.DataFrame(leave).to_csv(P/'tables/LEAVE_OUT_EFFECTS.tsv',sep='\t',index=False)
checks=[]
for row,(family,x,y,cov) in zip(o.itertuples(),edges):
 z=d.dropna(subset=[x,y]+cov);D=np.column_stack([np.ones(len(z))]+[(z.study.to_numpy()==st).astype(float) for st in sorted(z.study.unique())[1:]]+[stats.rankdata(z[c]) for c in cov])
 u=stats.rankdata(z[x]);v=stats.rankdata(z[y]);u-=D@np.linalg.lstsq(D,u,rcond=None)[0];v-=D@np.linalg.lstsq(D,v,rcond=None)[0]
 manual=np.corrcoef(u,v)[0,1];assert np.isclose(manual,row.rho,atol=1e-12)
 checks.append(dict(check=row.edge_id+'_independent_lstsq_rho',passed=True,error=abs(manual-row.rho)))
# Audit endpoints against prior same-definition main correlations: effects must match, MC P may vary by seed.
prior=pd.read_csv(R/'phase292_barrier_execution_and_external_feasibility_20260930_v1/tables/ASSOCIATIONS_v2.tsv',sep='\t')
for x in sources+[M]:
 old=prior[(prior.source==x)&(prior.endpoint==Y)&(prior.model=='study_only')].iloc[0]
 new=o[(o.family=='core_edges')&(o.source==x)&(o.endpoint==Y)].iloc[0]
 assert np.isclose(old.rho,new.rho,atol=1e-12)
 checks.append(dict(check=x+'_prior_same_edge_rho',passed=True,error=abs(old.rho-new.rho)))
assert o.groupby('family').size().to_dict()==dict(core_edges=5,conditional_edges=4,state_proxy=3,depth_sensitivity=8)
checks.append(dict(check='prespecified_family_sizes',passed=True,error=0))
pd.DataFrame(checks).to_csv(P/'qc/INDEPENDENT_CHECKS.tsv',sep='\t',index=False)
dump(P/'qc/COMPLETION.json',dict(status='COMPUTED_VALIDATED',patients=len(d),edges=len(o),checks_passed=len(checks),protocol_sha256=sha(P/'PROTOCOL_FROZEN.md'),all_inputs_unchanged=all(sha(p)==f['sha256'] for p,f in zip(paths,json.loads((P/'qc/INPUT_FREEZE.json').read_text())['inputs']))))
print(o[['edge_id','family','source','endpoint','controls','n','rho','P','FDR','evidence']].to_string(index=False));print('PHASE296_ANALYSIS_COMPLETE',flush=True)
