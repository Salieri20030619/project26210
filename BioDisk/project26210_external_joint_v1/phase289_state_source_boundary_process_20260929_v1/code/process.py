from common import *
import importlib.util
spec=importlib.util.spec_from_file_location('source280',R/'phase280_TYMP_source_niche_function_20260928_v1/code/source.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
def main():
 gate('process_linkage')
 panels=json.loads((P/'PANELS_FROZEN.json').read_text())
 pair=read(R/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv')
 group=read(R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_INDEX.tsv')
 z=np.load(R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_PSEUDOBULK.npz',allow_pickle=True)
 genes=z['genes'].tolist();gi={g:i for i,g in enumerate(genes)};counts=z['counts']
 stromal=read(P/'tables/STROMAL_STATE_PSEUDOBULK.tsv')
 rows=[];qs=[];allsource=[]
 for r in pair.itertuples():
  rec=dict(patient=r.patient,study=r.study,library=r.sample_id)
  for label,roles in [('malignant',['Malignant']),('macro',['C4','OtherMacrophage']),('C4',['C4'])]:
   use=group.patient.eq(r.patient)&group.sample_id.eq(r.sample_id)&group.source_tissue.eq('Tumor')&group.role.isin(roles)
   nc=int(group.loc[use,'cells'].sum());cc=counts[use.to_numpy()].sum(0);total=float(cc.sum());log=np.log1p(cc/max(total,1)*1e4)
   rec[label+'_cells']=nc
   for gene in ['TYMP','SPP1','TGFB1','CTSB']:rec[label+'_'+gene]=float(log[gi[gene]]) if nc else np.nan
   for prog,gs in panels.items():
    inds=[gi[g] for g in gs if g in gi];value=float(log[inds].mean()) if len(inds)>=.8*len(gs) and nc else np.nan
    allsource.append(dict(patient=r.patient,study=r.study,library=r.sample_id,source=label,cells=nc,program=prog,score=value))
  for pop in ['stromal_pool','activated','transition','conservative']:
   hit=stromal.patient.eq(r.patient)&stromal.library.eq(r.sample_id)&stromal.tissue.eq('Tumor')&stromal.population.eq(pop)
   assert hit.sum()<=1
   rec[pop+'_cells']=int(stromal.loc[hit,'cells'].iloc[0]) if hit.any() else 0
   for prog in panels:
    rec[pop+'_'+prog]=float(stromal.loc[hit,prog].iloc[0]) if hit.any() else np.nan
   if pop=='stromal_pool':
    for gene in ['CD44','ITGAV','ITGB1','TGFBR1','TGFBR2']:
     rec[gene+'_positive_cells']=int(stromal.loc[hit,gene+'_positive_cells'].iloc[0]) if hit.any() else 0
     rec[gene+'_logCP10K']=float(stromal.loc[hit,gene+'_logCP10K'].iloc[0]) if hit.any() else np.nan
  rec['eligible']=rec['malignant_cells']>=50 and rec['macro_cells']>=100 and rec['stromal_pool_cells']>=50
  rec['eligibility_reason']='OK' if rec['eligible'] else ';'.join([name for name,ok in [('malignant_lt50',rec['malignant_cells']>=50),('macro_lt100',rec['macro_cells']>=100),('stromal_lt50',rec['stromal_pool_cells']>=50)] if not ok])
  rows.append(rec)
 d=pd.DataFrame(rows);tab(d,'MATCHED_SOURCE_RECEIVER_INPUTS.tsv');tab(allsource,'SOURCE_PROGRAM_COMPARISON.tsv')
 endpoints=[]
 for src in ['malignant','macro']:
  for y in ['CAF_activation','TGF_response','ECM_organization']:endpoints.append(('TYMP_receiver',src+'_TYMP','stromal_pool_'+y,src+'_cells'))
 endpoints += [('output_receiver','macro_TGFB1','stromal_pool_TGF_response','macro_cells'),('output_receiver','macro_SPP1','stromal_pool_ECM_organization','macro_cells'),('output_receiver','macro_CTSB','stromal_pool_ECM_remodeling','macro_cells')]
 stats=[]
 use=d.eligible.to_numpy()
 for family,x,y,cell in endpoints:
  for model in ['study_adjusted','cell_count_sensitivity']:
   ex=np.column_stack([np.log1p(d.loc[use,cell]),np.log1p(d.loc[use,'stromal_pool_cells'])]) if model=='cell_count_sensitivity' else None
   s=base.cross(d.loc[use,x],d.loc[use,y],d.loc[use,'study'],extra=ex,seed=289,B=9999)
   stats.append(dict(family=family,predictor=x,endpoint=y,model=model,**s))
 out=pd.DataFrame(stats);out['FDR']=np.nan
 from scipy import stats as ss
 for _,g in out.groupby(['family','model']):
  good=g.P.notna()
  if good.any():out.loc[g.index[good],'FDR']=ss.false_discovery_control(g.loc[good,'P'].to_numpy(),method='bh')
 tab(out,'SOURCE_RECEIVER_ASSOCIATIONS.tsv')
 # Specific activated or transition state capacity/descriptive means only, not added association fishing.
 cap=[]
 for pop in ['activated','transition','conservative']:
  z=d[d[pop+'_cells']>=20]
  cap.append(dict(population=pop,matched_patients=len(z),studies=z.study.nunique(),cells=int(z[pop+'_cells'].sum()),interpretation='descriptive_capacity_no_new_tests'))
 tab(cap,'SPECIFIC_STATE_MATCHED_CAPACITY.tsv')
 sourcefiles=[R/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_INDEX.tsv',R/'phase157_ordered_bridge_reassessment_20260915_v1/02_inputs/GROUP_PSEUDOBULK.npz']
 dump(P/'qc/PROCESS_INPUT_MANIFEST.json',[dict(path=str(f),sha256=sha(f)) for f in sourcefiles])
 dump(P/'qc/PROCESS_LINKAGE_COMPLETE.json',dict(status='COMPUTED',original_matched_patients=len(d),eligible=int(use.sum()),models=len(out),new_mediation=False,client_experiment_reanalysis=False))
 print('PROCESS_LINKAGE_COMPLETE',int(use.sum()),flush=True)
if __name__=='__main__':main()
