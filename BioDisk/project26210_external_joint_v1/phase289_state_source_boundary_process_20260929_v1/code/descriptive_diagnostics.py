from common import *
from scipy import stats
def main():
 d=read(P/'tables/MATCHED_SOURCE_RECEIVER_INPUTS.tsv');s=read(P/'tables/SOURCE_RECEIVER_ASSOCIATIONS.tsv')
 rows=[]
 for r in s.itertuples():
  if r.model!='study_adjusted':continue
  for st in sorted(d.loc[d.eligible,'study'].unique()):
   for scope,mask in [('within_study',d.eligible&d.study.eq(st)),('leave_study',d.eligible&d.study.ne(st))]:
    x=d.loc[mask,r.predictor].to_numpy();y=d.loc[mask,r.endpoint].to_numpy();n=len(x)
    effect=np.nan
    if n>=3 and np.std(x)>0 and np.std(y)>0:
     X=np.column_stack([np.ones(n)]+[(d.loc[mask,'study'].to_numpy()==t).astype(float) for t in sorted(d.loc[mask,'study'].unique())[1:]])
     a=stats.rankdata(x);b=stats.rankdata(y);a-=X@np.linalg.lstsq(X,a,rcond=None)[0];b-=X@np.linalg.lstsq(X,b,rcond=None)[0]
     if n-np.linalg.matrix_rank(X)>=3 and np.std(a)>1e-10 and np.std(b)>1e-10:effect=float(np.corrcoef(a,b)[0,1])
    rows.append(dict(family=r.family,predictor=r.predictor,endpoint=r.endpoint,study=st,scope=scope,patients=n,rho=effect,P=np.nan,interpretation='descriptive_cohort_concentration_check'))
 tab(rows,'PROCESS_STUDY_CONCENTRATION_DESCRIPTIVE.tsv')
 src=read(P/'tables/SOURCE_GENE_PATIENT_LIBRARY.tsv')
 # Report patient tissue aggregated source expression with at least20 cells, no cell-level inference.
 agg=[]
 for key,g in src.groupby(['study','patient','tissue','source_group','gene']):
  cells=int(g.cells.sum());total=g.total_UMI.sum();count=g.raw_counts.sum()
  agg.append(dict(zip(['study','patient','tissue','source_group','gene'],key),cells=cells,positive_cells=int(g.positive_cells.sum()),detection_fraction=float(g.positive_cells.sum()/cells),logCP10K=float(np.log1p(count/total*1e4)),eligible_descriptive=cells>=20))
 tab(agg,'SOURCE_GENE_PATIENT_DESCRIPTIVE.tsv')
 cap=[]
 for pop in ['activated','transition','conservative']:
  for scope,base in [('all_16_matched',np.ones(len(d),bool)),('eligible_13_source_receiver',d.eligible)]:
   g=d[base&d[pop+'_cells'].ge(20)];cap.append(dict(population=pop,scope=scope,patients=len(g),studies=g.study.nunique(),cells=int(g[pop+'_cells'].sum())))
 tab(cap,'SPECIFIC_STATE_CAPACITY_BY_SCOPE.tsv')
 source=[R/'phase287_stromal_refinement_image_morphology_20260929_v1/tables/STROMAL_ANNOTATION_V1.tsv.gz',R/'phase287_stromal_refinement_image_morphology_20260929_v1/inputs/STROMAL_COUNTS.npz',R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/MACRO_METADATA.tsv.gz']
 dump(P/'qc/STATE_INPUT_MANIFEST.json',[dict(path=str(f),sha256=sha(f)) for f in source])
 print('DESCRIPTIVE_DIAGNOSTICS_COMPLETE',flush=True)
if __name__=='__main__':main()
