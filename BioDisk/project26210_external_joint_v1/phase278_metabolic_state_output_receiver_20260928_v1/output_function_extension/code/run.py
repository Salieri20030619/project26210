from pathlib import Path
import json,importlib.util,traceback
import numpy as np,pandas as pd
from scipy import sparse,stats
D=Path(__file__).resolve().parents[1];P=D.parent;R=P.parent
spec=importlib.util.spec_from_file_location('base',P/'code/source_analysis.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b);b.P=D
def main():
 b.gate('output_function');genes_out=json.loads((P/'PROGRAMS_FROZEN.json').read_text())['outputs']
 src=R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs'
 m=pd.read_csv(src/'MACRO_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str});X=sparse.load_npz(src/'MACRO_COUNTS.npz').tocsr();genes=pd.read_csv(src/'GENES.tsv',sep='\t').gene.tolist()
 total=np.asarray(X.sum(1)).ravel();raw=X[:,[genes.index(g) for g in genes_out]].toarray();assert np.allclose(total,m.library_counts,atol=1e-6,rtol=0)
 f=R/'phase213_immune_restart_execution_20260919_v1/I1_patient_inputs/LIBRARY_METRICS_ALL_VARIANTS.tsv'
 libs=pd.read_csv(f,sep='\t');libs=libs[(libs.variant=='primary')&libs.endpoint.isin(['CD8_cytotoxicity','CD8_dysfunction'])&libs.C4_intrinsic_eligible];rows=[]
 for r in libs.itertuples():
  mask=m.source_tissue.str.lower().eq('tumor')&m.curated_patient_id.eq(r.patient)&m.sample_id.eq(r.library)&m.candidate.eq('4')
  assert mask.sum()==r.n_C4 and mask.sum()>=20
  xx=np.log2(raw[mask].sum(0)/total[mask].sum()*1e6+1)
  for j,g in enumerate(genes_out):rows.append(dict(patient=r.patient,study=r.study,library=r.library,gene=g,endpoint=r.endpoint,x=float(xx[j]),y=r.outcome,C4_cells=int(mask.sum()),receiver_cells=r.receiver_cells))
 lib=pd.DataFrame(rows);b.tab(lib,'LIBRARY_INPUTS.tsv')
 pat=lib.groupby(['patient','study','gene','endpoint']).agg(x=('x','mean'),y=('y','mean'),n_libraries=('library','size'),libraries=('library',lambda x:';'.join(sorted(x)))).reset_index();b.tab(pat,'PATIENT_INPUTS.tsv')
 out=[]
 for i,((gene,ep),g) in enumerate(pat.groupby(['gene','endpoint'],sort=True)):
  for model in ['unadjusted','study_adjusted']:
   st=g.study.to_numpy() if model=='study_adjusted' else np.repeat('one_block',len(g))
   row=dict(gene=gene,endpoint=ep,model=model,n=len(g),rho=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan,method='DESCRIPTIVE',permutations=0)
   if len(g)>=10:
    res,null=b.main_test(g.x.to_numpy(),g.y.to_numpy(),st,seed=262378+i,B=99999);lo,hi,valid,boot=b.bootstrap_effect(g.x,g.y,st,B=5000,seed=262478+i)
    row.update(rho=res['rho'],P=res['p'],CI_low=lo,CI_high=hi,method=res['method'],permutations=res['permutations'],possible_permutations=res['possible_permutations'],residual_df=res['residual_df'],valid_bootstraps=valid)
    np.savez_compressed(D/'inputs'/f'NULL_{gene}_{ep}_{model}.npz',permutation=null,bootstrap=boot)
   out.append(row)
  print('FUNCTION_GENE_ENDPOINT',i+1,gene,ep,flush=True)
 su=pd.DataFrame(out);assert len(su)==60
 su['FDR']=np.where(su.P.notna(),stats.false_discovery_control(su.P.fillna(1)),np.nan);su['family_size']=60;b.tab(su,'OUTPUT_CD8_FUNCTION_SUMMARY.tsv')
 b.dump(D/'qc/COMPUTATION_COMPLETED.json',dict(at=b.now(),endpoints=60,patients=pat.patient.nunique(),studies=pat.study.nunique(),figures=0))
 print('OUTPUT_FUNCTION_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception:b.status('FAILED',traceback=traceback.format_exc());raise
