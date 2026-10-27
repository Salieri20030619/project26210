from pathlib import Path
import numpy as np,pandas as pd,json,hashlib,subprocess
from scipy import stats
P=Path(__file__).resolve().parents[1];R=P.parent;A=R/'phase291_spatial_null_TYMP_matched_regions_20260930_v1';B=R/'phase298_matrix_to_immune_entry_20260930_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gate(s):
 with (P/'logs/memory_gate.log').open('a') as f:
  f.write(s+'\n');f.flush();subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,check=True)
s=pd.read_csv(B/'tables/SECTION_RESULTS.tsv',sep='\t')
q=s[(s.dataset=='TJRCC')&(s.radius==200)&(s.model=='vascular')&(s.program=='ECM_organization')&(s.endpoint=='T')&s.reason.eq('OK')]
assert len(q)==9 and q.patient.is_unique
freeze={'protocol_sha256':sha(P/'PROTOCOL_FROZEN.md'),'patients':q.patient.tolist(),'sample_list':q['sample'].tolist()}
(P/'qc/FREEZE.json').write_text(json.dumps(freeze,indent=2))
rows=[];checks=[];man={}
for r in q.itertuples():
 gate(r.sample)
 ep=A/'inputs'/f'{r.sample}_200_REGIONS.tsv.gz';e=pd.read_csv(ep,sep='\t').set_index('edge');man[str(ep)]=sha(ep)
 for model in ['core','vascular']:
  f=A/'inputs'/f'{r.sample}_200_{model}_joint_to_T_distribution_SPECTRAL.npz';g=A/'inputs'/f'{r.sample}_200_{model}_joint_to_T_distribution_EXACT_NULL.npz'
  man[str(f)]=sha(f);man[str(g)]=sha(g)
  c=np.load(f);v=np.load(g);z=e.loc[c['edge']]
  D=c['D'];W=c['W'];V=v['V']
  names={'M':'M_state','E':'out_p_ECM_organization','T':'T_delta'}
  assert np.isfinite(z[list(names.values())]).all().all()
  rank={k:stats.rankdata(z[col]) for k,col in names.items()}
  residual={k:rr-D@np.linalg.lstsq(D,rr,rcond=None)[0] for k,rr in rank.items()}
  for edge,x,y in [('M_E','M','E'),('E_T','E','T'),('M_T','M','T')]:
   a=residual[x];b=residual[y];den=np.linalg.norm(a)*np.linalg.norm(b);assert den>1e-10
   rho=float(a@b/den);aa=V.T@a;assert np.max(abs(V@aa-a))<1e-7
   seed=int(hashlib.sha256((r.sample+model+edge+'phase300').encode()).hexdigest()[:8],16)
   signs=np.random.default_rng(seed).choice([-1.,1.],size=(len(aa),999));null=(V.T@b)@(aa[:,None]*signs)/den
   p=float((1+(abs(null)>=abs(rho)-1e-12).sum())/1000)
   assert abs(np.corrcoef(a,b)[0,1]-rho)<1e-10
   mor=(a@W@a)/(a@a)
   for k in [0,499,998]:
    sim=V@(aa*signs[:,k]);assert abs(D.T@sim).max()<1e-6
    assert abs(sim@sim/(a@a)-1)<1e-8 and abs((sim@W@sim)/(sim@sim)-mor)<1e-7
    assert abs(b@sim/den-null[k])<1e-10
   old=s[(s['sample']==r.sample)&(s.radius==200)&(s.model==model)&(s.program=='ECM_organization')&(s.endpoint=='T')].iloc[0]
   if edge=='E_T':assert abs(rho-old.effect)<1e-10
   rows.append(dict(patient=r.patient,sample=r.sample,model=model,edge=edge,regions=len(z),rho=rho,spatial_P=p,old_E_T_spatial_P=old.spatial_P if edge=='E_T' else np.nan,old_E_T_spatial_FDR=old.spatial_FDR if edge=='E_T' else np.nan))
   checks.append(dict(patient=r.patient,model=model,edge=edge,passed=True))
o=pd.DataFrame(rows);o['spatial_FDR']=np.nan
for _,ix in o.groupby(['sample','model']).groups.items():o.loc[ix,'spatial_FDR']=stats.false_discovery_control(o.loc[ix,'spatial_P'])
o.to_csv(P/'tables/SAME_REGION_EDGES.tsv',sep='\t',index=False)
d=o.pivot(index=['patient','sample','model'],columns='edge',values='rho').reset_index();d['two_links_compatible']=d.M_E.gt(0)&d.E_T.lt(0);d['all_three_compatible']=d.two_links_compatible&d.M_T.lt(0)
d.to_csv(P/'tables/PATIENT_LINK_ALIGNMENT.tsv',sep='\t',index=False)
su=[]
for model,z in d.groupby('model'):
 su.append(dict(model=model,n=len(z),M_E_positive=int(z.M_E.gt(0).sum()),M_E_median=float(z.M_E.median()),E_T_negative=int(z.E_T.lt(0).sum()),E_T_median=float(z.E_T.median()),M_T_negative=int(z.M_T.lt(0).sum()),M_T_median=float(z.M_T.median()),two_links_compatible=int(z.two_links_compatible.sum()),all_three_compatible=int(z.all_three_compatible.sum())))
pd.DataFrame(su).to_csv(P/'tables/ALIGNMENT_SUMMARY.tsv',sep='\t',index=False)
assert all(sha(Path(p))==h for p,h in man.items())
(P/'qc/COMPLETION.json').write_text(json.dumps(dict(status='COMPUTED_VALIDATED',n=9,edge_models=len(o),checked=len(checks),new_patient_P=False,source_hashes=man,summary=su),indent=2))
print(pd.DataFrame(su).to_string(index=False));print(d.to_string(index=False));print('PHASE300_ANALYSIS_COMPLETE')
