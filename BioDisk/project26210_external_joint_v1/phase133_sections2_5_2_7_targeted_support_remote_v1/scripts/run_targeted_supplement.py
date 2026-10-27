from pathlib import Path
import json,hashlib,warnings
import numpy as np,pandas as pd
from scipy.stats import rankdata,wilcoxon
w=Path(__file__).resolve().parents[1];q=w/'qc';root=w.parent
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
p129=root/'phase129_sections2_6_2_7_redraw_remote_v1';p128=root/'phase128_section2_5_metaviper_redraw_remote_v1'
sources=[]
def tab(p,index=0):
 p=Path(p);sources.append(p);return pd.read_csv(p,sep='\t',index_col=index)
m=tab(b/'07_pyscenic_v1/input/cell_metadata_v1.tsv.gz')
a=tab(p128/'data/SELECTED_METAVIPER_ACTIVITY.tsv');m=m.loc[a.index]
raw=tab(p129/'data/ENZYME_RAW_COUNTS.tsv.gz').loc[a.index]
e=tab(p129/'data/ENZYME_LOG1P_CP10K.tsv.gz').loc[a.index]
lib=tab(p129/'qc/RNA_EXTRACTION_ROWS.tsv').loc[a.index]
rna=tab('/mnt/BioDisk/wangmingyuan/project26210_part1_metaviper_downstream_phase29_20260814_v1/m13_expression_v1/m13_tymp_rrm2b_cell_expression_v1.tsv.gz').loc[a.index]
tr=tab(b/'09_metabolism_trajectory/trajectory_v1/source_data/trajectory_all_methods_cell_scores_v1.tsv.gz').loc[a.index]
assert len(a)==11921 and m.curated_patient_id.nunique()==20
assert np.array_equal(raw.TYMP,rna.TYMP_raw_count)
assert np.allclose(e.TYMP,rna.TYMP_log1p_cp10k,atol=2e-6)
raw['RRM2B']=rna.RRM2B_raw_count;e['RRM2B']=rna.RRM2B_log1p_cp10k
assert np.allclose(raw,raw.round()) and (raw>=0).all().all() and e.notna().all().all()
assert np.allclose(np.log1p(raw.RRM2B/lib.total_retained_UMI*1e4),e.RRM2B,atol=2e-6)
protein_targets=['CASP9','CHMP1A','PRSS8','PXN'];genes=['TYMP','TK1','SLC29A1','SLC29A2','ODC1','SRM','SMS','MTAP','CHKA','PCYT1A','PCYT2','CEPT1','RRM2B']
methods=['Monocle3','Slingshot','scTour']
def rho(x,y):
 if len(x)<30 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
 return float(np.corrcoef(rankdata(x),rankdata(y))[0,1])
def partial_rho(x,y,cov):
 if len(x)<30 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
 xx=rankdata(x);yy=rankdata(y);cc=np.column_stack([np.ones(len(x))]+[rankdata(cov[:,i]) for i in range(cov.shape[1])])
 rx=xx-cc@np.linalg.lstsq(cc,xx,rcond=None)[0];ry=yy-cc@np.linalg.lstsq(cc,yy,rcond=None)[0]
 if np.std(rx)<1e-10 or np.std(ry)<1e-10:return np.nan
 return float(np.corrcoef(rx,ry)[0,1])
def wp(x):
 x=np.asarray(x,float);x=x[np.isfinite(x)];x[np.abs(x)<1e-12]=0
 if len(x)<2:return np.nan
 if not np.any(x):return 1.
 with warnings.catch_warnings():
  warnings.simplefilter('ignore');return float(wilcoxon(x,method='approx',correction=True,zero_method='wilcox').pvalue)
def bh(x):
 p=np.asarray(x,float);o=np.argsort(p);out=np.empty(len(p));out[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return out
rng=np.random.default_rng(262133)
def summarize(vals):
 vals=np.asarray(vals,float);vals=vals[np.isfinite(vals)];vals[abs(vals)<1e-12]=0
 if len(vals)==0:return {'n':0,'median_effect':np.nan,'positive':0,'negative':0,'zero':0,'P':np.nan,'median_CI_low':np.nan,'median_CI_high':np.nan}
 boot=np.median(vals[rng.integers(0,len(vals),size=(5000,len(vals)))],axis=1);lo,hi=np.quantile(boot,[.025,.975])
 return {'n':len(vals),'median_effect':float(np.median(vals)),'positive':int(np.sum(vals>0)),'negative':int(np.sum(vals<0)),'zero':int(np.sum(vals==0)),'P':wp(vals),'median_CI_low':float(lo),'median_CI_high':float(hi)}
def status(p,q,eff):
 if not np.isfinite(eff):return 'not_estimable'
 direction='positive' if eff>0 else 'negative' if eff<0 else 'zero'
 return direction+('_FDR' if np.isfinite(q) and q<.05 else '_nominal' if p<.05 else '_direction_only')
# A: directly link TYMP RNA to activity, with a fixed technical-covariate sensitivity.
A=[]
for pid,sub in m.groupby('curated_patient_id',sort=True):
 ix=sub.index;x=e.loc[ix,'TYMP'].values;cov=sub[['raw_library_size_v1','raw_detected_genes_v1']].values
 assert np.isfinite(cov).all()
 for target in protein_targets+['RRM2B_RNA']:
  y=e.loc[ix,'RRM2B'].values if target=='RRM2B_RNA' else a.loc[ix,target].values
  for model,eff in [('unadjusted',rho(x,y)),('technical_adjusted',partial_rho(x,y,cov))]:
   A.append({'patient':pid,'study':sub.study.iloc[0],'target':target,'type':'RNA' if target=='RRM2B_RNA' else 'metaVIPER_activity','model':model,'n_cells':len(ix),'effect':eff})
A=pd.DataFrame(A);AS=[]
for (target,typ,model),d in A.groupby(['target','type','model'],sort=True):AS.append({'target':target,'type':typ,'model':model,**summarize(d.effect)})
AS=pd.DataFrame(AS);AS['FDR']=np.nan
for model in ['unadjusted','technical_adjusted']:
 mask=(AS.type=='metaVIPER_activity')&(AS.model==model);AS.loc[mask,'FDR']=bh(AS.loc[mask,'P'])
AS['evidence']=[status(r.P,r.FDR,r.median_effect) for r in AS.itertuples()]
A.to_csv(w/'data/TYMP_ACTIVITY_PATIENT_EFFECTS.tsv',sep='\t',index=False);AS.to_csv(q/'TYMP_ACTIVITY_SUMMARY.tsv',sep='\t',index=False)
# B: independently normalized early/terminal raw-count pseudobulks.
B=[];group_rows=[];correlations=[]
for pid,sub in m.groupby('curated_patient_id',sort=True):
 ix=sub.index
 for method in methods:
  t=tr.loc[ix,method+'_pseudotime_patient_v1'].values;q20,q80=np.quantile(t,[.2,.8]);early=t<=q20;late=t>=q80
  eligible=len(ix)>=30 and q80>q20 and early.sum()>=5 and late.sum()>=5
  if not eligible:continue
  assert not np.any(early&late)
  elo=lib.loc[ix[early],'total_retained_UMI'].sum();ehi=lib.loc[ix[late],'total_retained_UMI'].sum()
  for gene in genes:
   lowcount=float(raw.loc[ix[early],gene].sum());highcount=float(raw.loc[ix[late],gene].sum())
   low=float(np.log2(1+lowcount/elo*1e6));high=float(np.log2(1+highcount/ehi*1e6));effect=high-low
   B.append({'patient':pid,'study':sub.study.iloc[0],'method':method,'gene':gene,'n_cells':len(ix),'early_cells':int(early.sum()),'late_cells':int(late.sum()),'early_log2_CPM1':low,'late_log2_CPM1':high,'effect':effect})
   group_rows.append({'patient':pid,'method':method,'gene':gene,'early_gene_UMI':lowcount,'late_gene_UMI':highcount,'early_total_UMI':float(elo),'late_total_UMI':float(ehi),'q20':float(q20),'q80':float(q80)})
   correlations.append({'patient':pid,'method':method,'gene':gene,'rho':rho(t,e.loc[ix,gene].values)})
B=pd.DataFrame(B);BS=[]
for (method,gene),d in B.groupby(['method','gene'],sort=True):BS.append({'method':method,'gene':gene,**summarize(d.effect)})
BS=pd.DataFrame(BS);assert len(BS)==39;BS['FDR']=bh(BS.P);BS['evidence']=[status(r.P,r.FDR,r.median_effect) for r in BS.itertuples()]
B.to_csv(w/'data/TERMINAL_ENZYME_PATIENT_EFFECTS.tsv',sep='\t',index=False);pd.DataFrame(group_rows).to_csv(w/'data/TERMINAL_PSEUDOBULK_COUNTS.tsv',sep='\t',index=False)
pd.DataFrame(correlations).to_csv(w/'data/TERMINAL_ENZYME_CONTINUOUS_RHO_DESCRIPTIVE.tsv',sep='\t',index=False);BS.to_csv(q/'TERMINAL_ENZYME_SUMMARY.tsv',sep='\t',index=False)
# Leave-one-patient-out effect stability is reported without selecting a subset.
loo=[]
for label,frame,keys in [('TYMP_activity',A,['target','type','model']),('terminal_RNA',B,['method','gene'])]:
 for key,d in frame.groupby(keys,sort=True):
  d=d[np.isfinite(d.effect)];vals=d.effect.values
  for j,pid in enumerate(d.patient):loo.append({'analysis':label,'endpoint':'|'.join(key),'omitted_patient':pid,'median_effect':float(np.median(np.delete(vals,j)))})
pd.DataFrame(loo).to_csv(q/'LEAVE_ONE_PATIENT_OUT_EFFECTS.tsv',sep='\t',index=False)
# Carry forward all original strict/broad axis endpoints, not only the favorable broad TYMP signal.
oldaxis=root/'local_records_migration_20260908_150455/project/phase30_part1_enzyme_function_axis_local_v1/results_local_v1/enzyme_axis_endpoint_summary_v1.tsv'
axis=tab(oldaxis,index=None);axis.to_csv(q/'PRIOR_ENZYME_AXIS_COMPLETE_SUMMARY.tsv',sep='\t',index=False)
validation={'status':'COMPUTED_PENDING_INDEPENDENT_REVIEW','cells':len(m),'patients':20,'activity_patient_effects':len(A),'activity_summary_rows':len(AS),'terminal_patient_effects':len(B),'terminal_summary_rows':len(BS),'terminal_BH_family':39,'protein_BH_family_per_model':4,'new_cell_level_inferential_tests':False,'wetlab_directions':'pending user; all new tests two-sided'}
(q/'COMPUTATION_STATUS.json').write_text(json.dumps(validation,indent=2));(q/'SOURCE_PATHS.json').write_text(json.dumps([str(p) for p in sources],indent=2))
print('A_SUMMARY\n'+AS.to_string(index=False),flush=True)
print('B_SUMMARY\n'+BS.to_string(index=False),flush=True)
