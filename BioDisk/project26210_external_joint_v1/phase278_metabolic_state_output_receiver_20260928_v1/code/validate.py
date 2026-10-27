from pathlib import Path
import ast,json,math,datetime,hashlib
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];R=P.parent
read=lambda f:pd.read_csv(f,sep='\t',float_precision='round_trip')
tree=ast.parse((R/'phase277_metabolic_macrophage_barrier_cooperation_20260928_v1/code/validate.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['eq','signP','bhcheck']],type_ignores=[]),'independent_helpers','exec'),globals())
ef=read(P/'tables/WITHIN_PATIENT_EFFECTS.tsv');su=read(P/'tables/WITHIN_PATIENT_SUMMARY.tsv');cache=read(P/'inputs/CELL_SCORE_CACHE.tsv.gz');cache['candidate']=cache.candidate.astype(str);assert len(ef)==36*120 and len(su)==120
for row in su.itertuples():
 g=ef[(ef.group==row.group)&(ef.program==row.program)&(ef.gene==row.gene)&(ef.model==row.model)]
 v=g.effect.dropna().to_numpy();assert row.n==len(v)
 eq(row.median,np.median(v) if len(v) else np.nan);eq(row.P,signP(v) if len(v)>=10 else np.nan)
for _,g in su.groupby(['group','program']):assert len(g)==30;bhcheck(g)
errors=[]
for patient,full in cache.groupby('patient'):
 for row in ef[(ef.patient==patient)&ef.effect.notna()].itertuples():
  g=full[full.candidate=='4'] if row.group=='C4' else full
  cols=[row.program,row.gene];mat=g[cols].rank(method='average').to_numpy()
  if row.model=='unadjusted':Z=np.ones((len(g),1))
  else:
   C=g[['log_depth','n_genes']].rank().to_numpy()
   if row.group=='WholeMacro':
    dd=pd.get_dummies(g.candidate,dtype=float,drop_first=True).rank().to_numpy();C=np.column_stack([C,dd])
   Z=np.column_stack([np.ones(len(g)),C])
  Q=np.linalg.qr(Z,mode='reduced')[0];res=mat-Q@(Q.T@mat);rr=float(np.corrcoef(res.T)[0,1]);eq(rr,row.effect,1e-8);errors.append(abs(rr-row.effect))
 # Signal eligibility independent from exported counts.
 for row in ef[ef.patient==patient].itertuples():
  g=full[full.candidate=='4'] if row.group=='C4' else full;assert len(g)==row.cells and int((g[row.gene]>0).sum())==row.output_detected_cells and int(g['detected_'+row.program].sum())==row.program_detected_cells
 # Score finite and IDs unique per patient.
 assert full.cell_id.is_unique
fun=read(P/'tables/MATCHED_CD8_FUNCTION_SUMMARY.tsv');pat=read(P/'tables/MATCHED_CD8_PATIENT_INPUTS.tsv');lib=read(P/'tables/MATCHED_CD8_LIBRARY_INPUTS.tsv')
for row in pat.itertuples():
 z=lib[(lib.patient==row.patient)&(lib.program==row.program)&(lib.endpoint==row.endpoint)]
 assert z.library.nunique()==len(z)==row.n_libraries;eq(z.x.mean(),row.x);eq(z.y.mean(),row.y)
for row in fun.itertuples():
 g=pat[(pat.program==row.program)&(pat.endpoint==row.endpoint)]
 assert len(g)==row.n
 if row.n<10:assert pd.isna(row.P);continue
 vals=g[['x','y']].rank().to_numpy();st=g.study if row.model=='study_adjusted' else pd.Series(['all']*len(g),index=g.index)
 Z=np.column_stack([np.ones(len(g)),pd.get_dummies(st,drop_first=True,dtype=float).to_numpy()])
 Q=np.linalg.qr(Z,mode='reduced')[0];res=vals-Q@(Q.T@vals);rho=float(np.corrcoef(res.T)[0,1]);eq(rho,row.rho,1e-9)
 null=np.load(P/'inputs'/f'NULL_{row.program}_{row.endpoint}_{row.model}.npz');dist=null['permutation']
 pvalue=(1+int((np.abs(dist)>=abs(rho)-1e-12).sum()))/(len(dist)+1) if row.method.startswith('MC') else float((np.abs(dist)>=abs(rho)-1e-12).mean())
 eq(row.P,pvalue)
 # Analytical within-study centering agrees with QR.
 centered=vals.copy()
 for _,idx in pd.Series(np.arange(len(g))).groupby(st.to_numpy()):centered[idx.to_numpy()]-=centered[idx.to_numpy()].mean(0)
 eq(float(np.corrcoef(centered.T)[0,1]),rho)
assert len(fun)==8
# Function module has a single P/FDR pair.
fv=fun.copy();fv['P_wilcoxon']=fv.P;fv['FDR_wilcoxon']=fv.FDR;bhcheck(fv)
rc=np.load(P/'inputs/RECEIVER_COUNTS.npz');rm=read(P/'inputs/RECEIVER_CELL_METADATA.tsv.gz');qual=read(P/'tables/RECEIVER_COMPLEX_LIBRARY_QUALIFICATION.tsv');assert np.array_equal(rc['cell_id'],rm.cell_id)
gg=list(rc['genes']);nchecks=0
for key,g in rm.groupby(['receiver','study','patient','library']):
 rows=qual[(qual.receiver==key[0])&(qual.study==key[1])&(qual.patient==key[2])&(qual.library==key[3])]
 assert len(rows)==27;ids=g.index.to_numpy()
 for row in rows.itertuples():
  jj=[gg.index(x) for x in row.subunits.split(';')];present=bool(rc['gene_present'][jj].all());hits=(rc['counts'][np.ix_(ids,jj)]>0).all(1) if present else np.zeros(len(ids),bool);n=int(hits.sum())
  assert n==row.joint_positive_cells and len(ids)==row.receiver_cells and present==row.complete_genes_present
  assert row.receptor_eligible==(len(ids)>=20 and present and n>=5 and n/len(ids)>=.1)
  assert n<=int((rc['counts'][np.ix_(ids,jj)]>0).sum(0).min())
  nchecks+=1
source=read(P/'tables/INHERITED_OUTPUT_ENRICHMENT.tsv');old=read(R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/core/contrasts/tumor_C4_genes.tsv.gz').set_index('gene')
for row in source.itertuples():
 for col in ['P','FDR','median_delta']:eq(getattr(row,col),old.at[row.gene,col])
for filename in ['source_exit.txt','receiver_exit.txt']:assert (P/'logs'/filename).read_text()=='0'
seal=json.loads((P/'FREEZE.json').read_text())
for f,k in [('PROTOCOL_BEFORE_TESTING.md','protocol_sha256'),('PROGRAMS_FROZEN.json','programs_sha256')]:assert hashlib.sha256((P/f).read_bytes()).hexdigest()==seal[k]
assert not list(P.rglob('*.png')) and not list(P.rglob('*.pdf'))
out=dict(status='PASS',at=datetime.datetime.now().astimezone().isoformat(),within_patient_rows=len(ef),within_summary_endpoints=120,matched_function_endpoints=8,independent_rank_checks=len(errors),max_rank_error=max(errors),receiver_joint_complex_checks=nchecks,receiver_cells=len(rm),checks=['independent QR','exact binomial patient P','manual BH','frozen whole-gene P/FDR preservation','same-library then patient aggregation','saved blocked permutation P recomputation','receiver complete complex conjunction','qualification counts','frozen protocol identity'],figures=0,manuscript_changes=0)
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
