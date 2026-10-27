from pathlib import Path
import pandas as pd,numpy as np,sys,json
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');S=P/'phase157_ordered_bridge_reassessment_20260915_v1';O=P/'phase159_metabolic_parts_C4_value_20260915_v1'
sys.path.insert(0,str(S));from bridge_stats import rankz,corrp
B=pd.read_csv(P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t').set_index('patient')
C=pd.read_csv(O/'MACROPHAGE_SIDE_PATIENT_VALUES.tsv',sep='\t')
F=pd.read_csv(O/'MACROPHAGE_SIDE_TYMP_X_SPP1.tsv',sep='\t')
D=pd.read_csv(O/'SAME_LIBRARY_C4_CONTRIBUTION.tsv',sep='\t').set_index('patient')
out=[];checks=[]
for scope,g in C.groupby('scope'):
 X=g.pivot(index='patient',columns='gene',values='log2CPM')
 Y=B.loc[X.index,'SPP1_pseudobulk'] if scope=='all_macrophage' else D.loc[X.index, 'C4_SPP1_log2CPM' if scope=='C4' else 'Other_SPP1_log2CPM']
 for i,r in F[(F.scope==scope)&(F.adjustment=='raw')].iterrows():
  rho,p=spearmanr(X[r.gene],Y)
  checks.append(dict(check='receiver_scipy_'+str(i),passed=abs(rho-r.XY_rho)<1e-10 and abs(p-r.XY_P)<1e-10))
 if scope!='all_macrophage':continue
 for gene in ['CYP27A1','LDHA','SLC16A3']:
  for adj in ['raw','study']:
   A=rankz(B.loc[X.index,'TYMP_pseudobulk']);M=rankz(X[gene]);YY=rankz(Y)
   cov=np.ones((len(X),1)) if adj=='raw' else pd.get_dummies(B.loc[X.index,'study'],dtype=float).to_numpy()
   cov=np.column_stack([cov,rankz(B.loc[X.index,'C4_fraction'])])
   a,pa=corrp(A,M,cov);xy,pxy=corrp(M,YY,cov);b,pb=corrp(M,YY,np.column_stack([cov,A]))
   out.append(dict(gene=gene,adjustment=adj,additional_covariate='C4_fraction',a_rho=a,a_P=pa,XY_rho=xy,XY_P=pxy,b_rho=b,b_P=pb))
T=pd.DataFrame(out)
for _,ix in T.groupby('adjustment').groups.items():
 for e in ['a','XY','b']:T.loc[ix,e+'_FDR_selected_sensitivity']=multipletests(T.loc[ix,e+'_P'],method='fdr_bh')[1]
T.to_csv(O/'MACROPHAGE_C4_COMPOSITION_SENSITIVITY.tsv',sep='\t',index=False)
v=pd.DataFrame(checks);v.to_csv(O/'RECEIVER_VALIDATION.tsv',sep='\t',index=False);assert v.passed.all()
status=json.loads((O/'STATUS.json').read_text());status.update(receiver_checks=len(v),receiver_checks_passed=int(v.passed.sum()),total_checks=status['checks']+status['raw_checks']+len(v))
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print(T.to_string(index=False))
R=pd.read_csv(O/'ALL_COMPONENT_BRIDGE_RESULTS.tsv',sep='\t')
print('DETAIL\n'+R[(R.candidate.isin(['gene:MSMO1','gene:HK2','gene:LDHA','gene:SLC16A3','gene:CYP27A1','module:glycolysis_lactate_conversion']))].to_string(index=False))
print('RECEIVER A Q\n'+F[(F.gene=='CYP27A1')&(F.scope=='all_macrophage')].to_string(index=False))
print('STATUS',json.dumps(status))

