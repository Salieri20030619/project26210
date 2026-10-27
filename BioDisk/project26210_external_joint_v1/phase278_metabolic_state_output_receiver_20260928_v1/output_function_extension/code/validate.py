from pathlib import Path
import json,datetime,ast,math
import numpy as np,pandas as pd
D=Path(__file__).resolve().parents[1];P=D.parent;R=P.parent
read=lambda f:pd.read_csv(f,sep='\t',float_precision='round_trip')
tree=ast.parse((R/'phase277_metabolic_macrophage_barrier_cooperation_20260928_v1/code/validate.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['eq','bhcheck']],type_ignores=[]),'helpers','exec'),globals())
s=read(D/'tables/OUTPUT_CD8_FUNCTION_SUMMARY.tsv');p=read(D/'tables/PATIENT_INPUTS.tsv');l=read(D/'tables/LIBRARY_INPUTS.tsv')
assert len(s)==60 and s.n.eq(11).all()
for row in p.itertuples():
 g=l[(l.patient==row.patient)&(l.gene==row.gene)&(l.endpoint==row.endpoint)]
 eq(row.x,g.x.mean());eq(row.y,g.y.mean());assert row.n_libraries==len(g)==g.library.nunique()
old=read(R/'phase213_immune_restart_execution_20260919_v1/I1_patient_inputs/LIBRARY_METRICS_ALL_VARIANTS.tsv')
for row in l[l.gene=='SPP1'].itertuples():
 g=old[(old.patient==row.patient)&(old.library==row.library)&(old.endpoint==row.endpoint)&(old.variant=='primary')]
 assert len(g)==1;eq(row.x,g.iloc[0].C4_SPP1_pseudobulk,1e-10);eq(row.y,g.iloc[0].outcome,1e-10)
errors=[]
for row in s.itertuples():
 g=p[(p.gene==row.gene)&(p.endpoint==row.endpoint)];vals=g[['x','y']].rank().to_numpy()
 st=g.study.to_numpy() if row.model=='study_adjusted' else np.repeat('all',len(g))
 Dm=np.column_stack([np.ones(len(g)),pd.get_dummies(st,drop_first=True,dtype=float).to_numpy()]);Q=np.linalg.qr(Dm,mode='reduced')[0];res=vals-Q@(Q.T@vals)
 rho=float(np.corrcoef(res.T)[0,1]);eq(rho,row.rho,1e-9);errors.append(abs(rho-row.rho))
 null=np.load(D/'inputs'/f'NULL_{row.gene}_{row.endpoint}_{row.model}.npz');v=null['permutation']
 pv=(1+int((np.abs(v)>=abs(rho)-1e-12).sum()))/(len(v)+1) if row.method.startswith('MC') else float((np.abs(v)>=abs(rho)-1e-12).mean());eq(pv,row.P)
b=s.copy();b['P_wilcoxon']=b.P;b['FDR_wilcoxon']=b.FDR;bhcheck(b)
assert (D/'logs/exit_status.txt').read_text()=='0'
out=dict(status='PASS',at=datetime.datetime.now().astimezone().isoformat(),endpoints=60,patients=11,independent_rank_error=max(errors),SPP1_source_reproduces_prior_frozen_library_values=True,library_and_patient_aggregation=True,permutation_P_recomputed=True,BH60_verified=True,figures=0)
(D/'qc/FINAL_VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
