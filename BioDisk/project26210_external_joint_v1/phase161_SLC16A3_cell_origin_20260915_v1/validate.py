from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import rankdata,wilcoxon
from statsmodels.stats.multitest import multipletests
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase161_SLC16A3_cell_origin_20260915_v1';S=P/'phase157_ordered_bridge_reassessment_20260915_v1'
C=pd.read_csv(O/'CELL_TARGET_EXPRESSION.tsv.gz',sep='\t');H=np.load(O/'HYPoxia_CELL_INPUT.npz',allow_pickle=True)['expression'];R=pd.read_csv(O/'PATIENT_WITHIN_CELL_CORRELATIONS.tsv',sep='\t')
checks=[]
for i,r in R[R.adjustment!='raw'].iterrows():
 mask=(C.curated_patient_id==r.patient)&(C.role.isin(['C4','OtherMacrophage']) if r.compartment=='Macrophage' else C.role==r.compartment)
 d=C[mask];hs=H[d.index];sd=hs.std(0,ddof=1);ok=sd>0;score=((hs[:,ok]-hs[:,ok].mean(0))/sd[ok]).mean(1)
 vals=[d['expr:'+r.x],d['expr:'+r.y],np.log1p(d.total_counts)]
 if r.adjustment=='depth_hypoxia':vals.append(score)
 mat=np.array([rankdata(v) for v in vals]);prec=np.linalg.inv(np.corrcoef(mat));rho=-prec[0,1]/np.sqrt(prec[0,0]*prec[1,1])
 checks.append(dict(check='precision_'+str(i),passed=abs(rho-r.rho)<1e-7))
Z=np.load(S/'02_inputs/GROUP_PSEUDOBULK.npz',allow_pickle=True);G=pd.read_csv(S/'02_inputs/GROUP_INDEX.tsv',sep='\t');genes=list(Z['genes'])
SRC=pd.read_csv(O/'PAIRED_COMPARTMENT_SOURCE.tsv',sep='\t')
for i,r in SRC.iterrows():
 samples=C.loc[C.curated_patient_id==r.patient,'sample_id'].unique();assert len(samples)==1
 roles=['Malignant'] if r.compartment=='Malignant' else ['C4','OtherMacrophage']
 ids=G[(G.patient==r.patient)&(G.sample_id==samples[0])&G.role.isin(roles)].group_index.to_numpy(int);cnt=Z['counts'][ids].sum(0)
 expected=np.log2(cnt[genes.index(r.gene)]/cnt.sum()*1e6+1)
 checks.append(dict(check='cache_count_'+str(i),passed=abs(expected-r.log2CPM)<1e-6))
out=[]
for (x,y,adj),gg in R.groupby(['x','y','adjustment']):
 c=gg[gg.compartment=='C4'].set_index('patient');o=gg[gg.compartment=='OtherMacrophage'].set_index('patient')
 common=c.index.intersection(o.index)
 if not len(common):continue
 delta=c.loc[common,'rho']-o.loc[common,'rho'];p=wilcoxon(delta).pvalue if len(delta)>=3 else np.nan
 out.append(dict(x=x,y=y,adjustment=adj,n_patients=len(common),C4_greater=int((delta>0).sum()),other_greater=int((delta<0).sum()),median_C4_rho=c.loc[common,'rho'].median(),median_other_rho=o.loc[common,'rho'].median(),median_paired_difference=delta.median(),P=p,exploratory=True))
D=pd.DataFrame(out)
for _,ix in D.groupby('adjustment').groups.items():D.loc[ix,'FDR_exploratory_pair_family']=multipletests(D.loc[ix,'P'],method='fdr_bh')[1]
D.to_csv(O/'C4_OTHER_COMMON_PATIENT_CORRELATION_DIFFERENCE.tsv',sep='\t',index=False)
V=pd.DataFrame(checks);V.to_csv(O/'INDEPENDENT_VALIDATION.tsv',sep='\t',index=False)
print('FAILED',V[~V.passed].to_dict('records'));assert V.passed.all()
status=json.loads((O/'STATUS.json').read_text());status.update(independent_checks=len(V),independent_passed=int(V.passed.sum()),total_checks=status['checks']+len(V));(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print(D.to_string(index=False));print(json.dumps(status))

