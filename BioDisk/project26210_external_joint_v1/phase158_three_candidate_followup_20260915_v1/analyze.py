from pathlib import Path
import sys,json,hashlib,subprocess,datetime
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
S=P/'phase157_ordered_bridge_reassessment_20260915_v1'
O=P/'phase158_three_candidate_followup_20260915_v1'
O.mkdir(exist_ok=True)
gate=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(O/'MEMORY_GATE.log').write_text(gate.stdout+gate.stderr);assert gate.returncode==0
sys.path.insert(0,str(S))
from bridge_stats import bridge,rankz,corrp
B=pd.read_csv(P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t').set_index('patient')
X=pd.read_csv(S/'04_other_X_screen/PATIENT_CANDIDATE_VALUES.tsv',sep='\t',index_col=0).loc[B.index]
selected=['ligand:IGFBP3','program:Cholesterol Homeostasis','program:Glycolysis']
assert len(B)==16 and B.index.is_unique and np.isfinite(X.to_numpy()).all()
plan={'scope':selected,'primary_endpoint':'C4_fraction','additional_endpoints':['SPP1_pseudobulk','SPP1_detection_fraction'],'multiplicity':'New SPP1 tests retain original 56 ligand/50 program/10 metabolism families separately by endpoint, edge and adjustment; selected three do not reset discovery FDR. C4 FDR inherited Phase157. Sensitivity families descriptive and explicitly separate.','post_selection':True,'C4_frozen':True,'figures':False,'client_data':False}
(O/'PLAN.json').write_text(json.dumps(plan,indent=2))
results=[]
for endpoint in ['C4_fraction','SPP1_pseudobulk','SPP1_detection_fraction']:
 for candidate in X:
  for adjustment in ['raw','study']:
   r=bridge(B.TYMP_pseudobulk,X[candidate],B[endpoint],B.study,adjustment)
   results.append(dict(candidate=candidate,family=candidate.split(':')[0],endpoint=endpoint,adjustment=adjustment,**r))
R=pd.DataFrame(results)
for _,idx in R.groupby(['family','endpoint','adjustment']).groups.items():
 for edge in ['a','XY','b']:R.loc[idx,edge+'_FDR']=multipletests(R.loc[idx,edge+'_P'],method='fdr_bh')[1]
R.to_csv(O/'ALL_FAMILY_ENDPOINT_RESULTS.tsv',sep='\t',index=False)
F=R[R.candidate.isin(selected)].copy()
F.to_csv(O/'THREE_CANDIDATE_ENDPOINT_RESULTS.tsv',sep='\t',index=False)
sens=[]
for candidate in selected:
 for covname in ['program:Hypoxia','program:G2-M Checkpoint']:
  for adj in ['raw','study']:
   C=np.ones((16,1)) if adj=='raw' else pd.get_dummies(B.study,dtype=float).to_numpy()
   C=np.column_stack([C,rankz(X[covname])])
   A=rankz(B.TYMP_pseudobulk);M=rankz(X[candidate]);Y=rankz(B.C4_fraction)
   ar,ap=corrp(A,M,C);xr,xp=corrp(M,Y,C);br,bp=corrp(M,Y,np.column_stack([C,A]))
   sens.append(dict(candidate=candidate,covariate=covname,adjustment=adj,n=16,a_rho=ar,a_P=ap,XY_rho=xr,XY_P=xp,b_rho=br,b_P=bp))
T=pd.DataFrame(sens)
for _,idx in T.groupby(['covariate','adjustment']).groups.items():
 for edge in ['a','XY','b']:T.loc[idx,edge+'_FDR_selected_sensitivity']=multipletests(T.loc[idx,edge+'_P'],method='fdr_bh')[1]
T.to_csv(O/'HYPOXIA_PROLIFERATION_SENSITIVITY.tsv',sep='\t',index=False)
loo=[]
for candidate in selected:
 for endpoint in ['C4_fraction','SPP1_pseudobulk','SPP1_detection_fraction']:
  for kind in ['patient','study']:
   groups=B.index if kind=='patient' else sorted(B.study.unique())
   for omit in groups:
    mask=B.index!=omit if kind=='patient' else B.study!=omit
    r=bridge(B.loc[mask,'TYMP_pseudobulk'],X.loc[mask,candidate],B.loc[mask,endpoint],B.loc[mask,'study'],'raw')
    loo.append(dict(candidate=candidate,endpoint=endpoint,kind=kind,omitted=omit,**r))
pd.DataFrame(loo).to_csv(O/'LEAVE_ONE_RESULTS.tsv',sep='\t',index=False)
# Exact same-library C4/other expression, retaining C4 zeros for proportion; >=20 cells only for compartment expression.
G=pd.read_csv(S/'02_inputs/GROUP_INDEX.tsv',sep='\t')
Z=np.load(S/'02_inputs/GROUP_PSEUDOBULK.npz',allow_pickle=True)
counts=Z['counts'];genes=list(Z['genes']);spp=genes.index('SPP1')
comp=[]
for pat,b in B.iterrows():
 for role in ['C4','OtherMacrophage']:
  g=G[(G.patient==pat)&(G.sample_id==b.sample_id)&(G.role==role)]
  if len(g)!=1:continue
  row=g.iloc[0];i=int(row.group_index)
  comp.append(dict(patient=pat,role=role,cells=int(row.cells),SPP1_log2CPM=float(np.log2(counts[i,spp]/counts[i].sum()*1e6+1)),SPP1_detection=float(Z['detected'][i,spp]/row.cells)))
Cdf=pd.DataFrame(comp);Cdf.to_csv(O/'MATCHED_COMPARTMENT_SPP1.tsv',sep='\t',index=False)
extra=[]
for role in ['C4','OtherMacrophage']:
 c=Cdf[(Cdf.role==role)&(Cdf.cells>=20)].set_index('patient')
 for candidate in selected:
  for adj in ['raw','study']:
   r=bridge(B.loc[c.index,'TYMP_pseudobulk'],X.loc[c.index,candidate],c.SPP1_log2CPM,B.loc[c.index,'study'],adj)
   extra.append(dict(candidate=candidate,compartment=role,adjustment=adj,**r))
E=pd.DataFrame(extra)
for _,idx in E.groupby(['compartment','adjustment']).groups.items():
 for edge in ['a','XY','b']:E.loc[idx,edge+'_FDR_selected_secondary']=multipletests(E.loc[idx,edge+'_P'],method='fdr_bh')[1]
E.to_csv(O/'COMPARTMENT_SPP1_SECONDARY.tsv',sep='\t',index=False)
# Independent direct Spearman audit across all raw endpoint fits plus path decomposition.
checks=[]
for i,r in R.iterrows():
 checks.append(dict(check='path_identity_'+str(i),passed=abs(r.identity_error)<1e-10))
 if r.adjustment=='raw':
  q,p=spearmanr(X[r.candidate],B[r.endpoint])
  checks.append(dict(check='scipy_XY_'+str(i),passed=bool(abs(q-r.XY_rho)<1e-10 and abs(p-r.XY_P)<1e-10)))
checks.append(dict(check='patient_count',passed=len(B)==16))
checks.append(dict(check='no_figures',passed=not any(O.glob('*.png')) and not any(O.glob('*.pdf'))))
V=pd.DataFrame(checks);V.to_csv(O/'VALIDATION.tsv',sep='\t',index=False);assert V.passed.all()
inputs=[P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',S/'04_other_X_screen/PATIENT_CANDIDATE_VALUES.tsv',S/'02_inputs/GROUP_PSEUDOBULK.npz',S/'02_inputs/GROUP_INDEX.tsv',S/'bridge_stats.py']
(O/'SOURCE_MANIFEST.json').write_text(json.dumps([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputs],indent=2))
(O/'STATUS.json').write_text(json.dumps(dict(status='COMPUTATION_COMPLETE_REPORT_PENDING',checks=len(V),passed=int(V.passed.sum()),time=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2))
print(F[['candidate','endpoint','adjustment','a_rho','a_P','XY_rho','XY_P','XY_FDR','b_partial_rho','b_P','b_FDR']].to_string(index=False))
print('SENSITIVITY\n'+T.to_string(index=False))
print('COMPARTMENT\n'+E[['candidate','compartment','adjustment','n_patients','XY_rho','XY_P','b_partial_rho','b_P']].to_string(index=False))
print('VALIDATION',len(V),int(V.passed.sum()))

