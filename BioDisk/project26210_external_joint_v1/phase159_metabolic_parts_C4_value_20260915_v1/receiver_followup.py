from pathlib import Path
import numpy as np,pandas as pd,sys,json
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');S=P/'phase157_ordered_bridge_reassessment_20260915_v1';O=P/'phase159_metabolic_parts_C4_value_20260915_v1'
sys.path.insert(0,str(S));from bridge_stats import bridge,rankz,corrp
B=pd.read_csv(P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t').set_index('patient')
G=pd.read_csv(S/'02_inputs/GROUP_INDEX.tsv',sep='\t');Z=np.load(S/'02_inputs/GROUP_PSEUDOBULK.npz',allow_pickle=True);genes=list(Z['genes']);cnt=Z['counts'];dt=Z['detected']
plan=json.loads((O/'PLAN_BEFORE_FIT.json').read_text())
rows=[];coverage=[];vv=[]
for scope,roles in [('all_macrophage',['C4','OtherMacrophage']),('C4',['C4']),('OtherMacrophage',['OtherMacrophage'])]:
 c=[];pp=[];st=[]
 for pat,b in B.iterrows():
  gg=G[(G.patient==pat)&(G.sample_id==b.sample_id)&G.role.isin(roles)]
  if gg.cells.sum()<20:continue
  pp.append(pat);c.append(cnt[gg.group_index.to_numpy(int)].sum(0));st.append(b.study)
 c=np.array(c);ex=pd.DataFrame(np.log2(c/c.sum(1)[:,None]*1e6+1),index=pp,columns=genes)
 for family,gg in plan['receiver_genes'].items():
  for gene in gg:
   j=genes.index(gene);k=int((c[:,j]>=10).sum());eligible=k>=np.ceil(.75*len(pp))
   coverage.append(dict(scope=scope,family=family,gene=gene,n=len(pp),libraries_ge10=k,eligible=eligible))
   for pat in pp:vv.append(dict(scope=scope,patient=pat,gene=gene,log2CPM=ex.loc[pat,gene]))
   if not eligible:continue
   for adj in ['raw','study']:
    rows.append(dict(scope=scope,family=family,gene=gene,adjustment=adj,**bridge(B.loc[pp,'TYMP_pseudobulk'],ex[gene],ex.SPP1,st,adj)))
R=pd.DataFrame(rows)
for _,idx in R.groupby(['scope','family','adjustment']).groups.items():
 for edge in ['a','XY','b']:R.loc[idx,edge+'_FDR']=multipletests(R.loc[idx,edge+'_P'],method='fdr_bh')[1]
R.to_csv(O/'MACROPHAGE_SIDE_TYMP_X_SPP1.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(O/'MACROPHAGE_SIDE_COVERAGE.tsv',sep='\t',index=False);pd.DataFrame(vv).to_csv(O/'MACROPHAGE_SIDE_PATIENT_VALUES.tsv',sep='\t',index=False)
# All sender gene/module fits additional hypoxia sensitivity, a distinct predefined followup family.
X=pd.read_csv(O/'PATIENT_SENDER_VALUES.tsv',sep='\t',index_col=0).loc[B.index]
REG=pd.read_csv(O/'CANDIDATE_REGISTRY.tsv',sep='\t');old=pd.read_csv(S/'04_other_X_screen/PATIENT_CANDIDATE_VALUES.tsv',sep='\t',index_col=0).loc[B.index]
sr=[]
for _,r in REG.iterrows():
 for endpoint in ['C4_fraction','SPP1_pseudobulk','SPP1_detection_fraction']:
  for adj in ['raw','study']:
   C=np.ones((16,1)) if adj=='raw' else pd.get_dummies(B.study,dtype=float).to_numpy()
   C=np.column_stack([C,rankz(old['program:Hypoxia'])])
   A=rankz(B.TYMP_pseudobulk);M=rankz(X[r.candidate]);Y=rankz(B[endpoint])
   a,pa=corrp(A,M,C);xy,pxy=corrp(M,Y,C);b,pb=corrp(M,Y,np.column_stack([C,A]))
   sr.append(dict(candidate=r.candidate,family=r.family,level=r.level,endpoint=endpoint,adjustment=adj,a_rho=a,a_P=pa,XY_rho=xy,XY_P=pxy,b_rho=b,b_P=pb))
SR=pd.DataFrame(sr)
for _,idx in SR.groupby(['family','level','endpoint','adjustment']).groups.items():
 for edge in ['a','XY','b']:SR.loc[idx,edge+'_FDR']=multipletests(SR.loc[idx,edge+'_P'],method='fdr_bh')[1]
SR.to_csv(O/'ALL_SENDER_HYPOXIA_SENSITIVITY.tsv',sep='\t',index=False)
print('MACRO BOTH PATHS POSITIVE (all values retained)')
print(R[(R.a_rho>0)&(R.b_partial_rho>0)&(R.a_P<.2)&(R.b_P<.1)][['scope','gene','adjustment','n_patients','a_rho','a_P','a_FDR','XY_rho','XY_P','b_partial_rho','b_P','b_FDR']].to_string(index=False))
print('CYP27A1, NR1H3, LDHA, SLC16A3\n'+R[R.gene.isin(['CYP27A1','NR1H3','LDHA','SLC16A3'])][['scope','gene','adjustment','n_patients','a_rho','a_P','XY_rho','XY_P','b_partial_rho','b_P','b_FDR']].to_string(index=False))
print('HYP SENDER\n'+SR[(SR.candidate.isin(['gene:MSMO1','gene:HK2','gene:PFKP','gene:FDPS','gene:LDHA','gene:SLC16A3']))&(SR.adjustment=='raw')].to_string(index=False))

