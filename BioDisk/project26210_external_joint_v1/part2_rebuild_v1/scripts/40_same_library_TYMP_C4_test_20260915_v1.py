from pathlib import Path
import json,datetime,hashlib,subprocess,sys
import pandas as pd,numpy as np,h5py
from scipy import sparse
from scipy.stats import spearmanr,rankdata,false_discovery_control
from anndata.io import read_elem
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1'
O=P/'phase154_same_library_TYMP_C4_test_20260915_v1';O.mkdir(exist_ok=False)
H=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
M=R/'08_regulons/input_v1/CELL_METADATA.tsv.gz';F=P/'phase152_TYMP_SPP1_direct_bridge_plan_20260915_v1/SAME_LIBRARY_FEASIBILITY.tsv'
def gate(stage):
 with (O/'RESOURCE_GATE.log').open('a') as f:
  f.write(stage+'\n');f.flush();subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,stderr=subprocess.STDOUT,check=True)
def state(stage,**kw):
 d=dict(status='RUNNING',stage=stage,time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw)
 (O/'STATUS.json').write_text(json.dumps(d,indent=2));print(stage,flush=True)
plan=dict(primary_focus='Same-library malignant TYMP pseudobulk expression versus frozen C4 fraction in mature macrophages',endpoints=['C4_fraction','SPP1_pseudobulk','SPP1_detection_fraction'],population='Previously frozen 16 tumour libraries/16 patients; malignant eligibility flag; mature macrophages exclude candidate1 and Monocyte',minimums='20 malignant and 50 mature macrophages in same library; retain zero C4',normalization='raw counts pseudobulk; log2(CPM+1); no binary exposure threshold',tests='Spearman across patients; within-study rank residual association with 9999 study-preserving permutations; three endpoint BH separately by analysis type; no significance or effect-size pass/fail gate',robustness='study strata, leave one patient/study; alternate mean per-cell log1p CP10K TYMP exposure; technical library-size adjustment diagnostic',interpretation='Exploratory continuous effects; nominal P permitted with full FDR retained; no causal claims',seed=262154,figures=False,created_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'ANALYSIS_PLAN.json').write_text(json.dumps(plan,indent=2))
state('input_metadata');gate('input_metadata')
elig=pd.read_csv(F,sep='\t');elig=elig[elig.eligible_same_library==True].copy()
assert len(elig)==16 and elig.sample_id.nunique()==16 and elig.curated_patient_id.nunique()==16
mm=pd.read_csv(M,sep='\t',index_col=0,dtype={'candidate':str},low_memory=False)
mm=mm[(~mm.candidate.isin(['1','Monocyte']))&(mm.source_tissue.str.lower()=='tumor')]
with h5py.File(H,'r') as f:
 ids=f['obs']['cell_id'].asstr()[:]
 fields=['study','curated_patient_id','sample_id','source_tissue','total_counts_raw','eligible_for_subclone_selection_v1']
 obs=pd.DataFrame({k:np.asarray(read_elem(f['obs'][k])) for k in fields},index=ids)
 genes=f['var']['gene'].asstr()[:];idx=[int(np.flatnonzero(genes==gene)[0]) for gene in ['TYMP','SPP1']]
 assert all((genes==gene).sum()==1 for gene in ['TYMP','SPP1'])
mal=obs[(obs.eligible_for_subclone_selection_v1==True)&(obs.source_tissue.astype(str).str.lower()=='tumor')&obs.sample_id.isin(elig.sample_id)]
mac=mm[mm.sample_id.isin(elig.sample_id)]
assert set(mac.index)<=set(obs.index);assert not set(mac.index)&set(mal.index)
selected=obs.index.isin(set(mac.index)|set(mal.index));positions=np.flatnonzero(selected)
wanted=obs.loc[selected].copy();wanted['lineage']=np.where(wanted.index.isin(mal.index),'Malignant','Macrophage');wanted['candidate']=mac.candidate.reindex(wanted.index).fillna('malignant').to_numpy()
for _,row in elig.iterrows():
 z=wanted[wanted.sample_id==row.sample_id]
 assert set(z.curated_patient_id)=={row.curated_patient_id} and set(z.study)=={row.study}
 assert (z.lineage=='Malignant').sum()==row.malignant_cells and (z.lineage=='Macrophage').sum()==row.macrophage_cells
wanted['TYMP_counts']=0.;wanted['SPP1_counts']=0.;wanted['recomputed_library_counts']=0.
state('extract_two_genes_raw_counts',selected_cells=len(wanted));gate('extract_two_genes_raw_counts')
with h5py.File(H,'r',rdcc_nbytes=128*1024*1024) as f:
 x=f['layers/counts'];ptr=x['indptr'][:];n=len(ids);block=8192
 for start in range(0,n,block):
  end=min(start+block,n);mask=selected[start:end]
  if not mask.any():continue
  lo,hi=int(ptr[start]),int(ptr[end]);ind=x['indices'][lo:hi];data=x['data'][lo:hi]
  mat=sparse.csr_matrix((data,ind,ptr[start:end+1]-lo),shape=(end-start,len(genes)))[mask]
  assert np.isfinite(mat.data).all() and (mat.data>=0).all() and np.equal(mat.data,np.floor(mat.data)).all()
  lib=np.asarray(mat.sum(axis=1,dtype=np.float64)).ravel();vals=mat[:,idx].toarray()
  cells=ids[start:end][mask];wanted.loc[cells,['TYMP_counts','SPP1_counts']]=vals
  wanted.loc[cells,'recomputed_library_counts']=lib
assert np.allclose(wanted.total_counts_raw,wanted.recomputed_library_counts,rtol=0,atol=1e-5)
assert (wanted.recomputed_library_counts>0).all()
wanted.to_csv(O/'SELECTED_CELL_TWO_GENE_COUNTS.tsv.gz',sep='\t',index_label='cell_id')
state('aggregate_same_library');gate('aggregate_same_library')
rows=[]
for _,e in elig.iterrows():
 z=wanted[wanted.sample_id==e.sample_id];a=z[z.lineage=='Malignant'];b=z[z.lineage=='Macrophage']
 tm=float(a.TYMP_counts.sum());ml=float(a.recomputed_library_counts.sum());sp=float(b.SPP1_counts.sum());bl=float(b.recomputed_library_counts.sum())
 rows.append(dict(study=e.study,patient=e.curated_patient_id,sample_id=e.sample_id,malignant_cells=len(a),macrophage_cells=len(b),C4_cells=int((b.candidate=='4').sum()),malignant_total_counts=ml,macrophage_total_counts=bl,malignant_TYMP_counts=tm,macrophage_SPP1_counts=sp,TYMP_pseudobulk=float(np.log2(tm/ml*1e6+1)),TYMP_mean_log1pCP10K=float(np.log1p(a.TYMP_counts/a.recomputed_library_counts*10000).mean()),C4_fraction=float((b.candidate=='4').mean()),SPP1_pseudobulk=float(np.log2(sp/bl*1e6+1)),SPP1_detection_fraction=float((b.SPP1_counts>0).mean())))
d=pd.DataFrame(rows);d.to_csv(O/'SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t',index=False)
state('patient_statistics');gate('patient_statistics')
responses=['C4_fraction','SPP1_pseudobulk','SPP1_detection_fraction'];out=[];rng=np.random.default_rng(262154);nperm=9999
for exposure in ['TYMP_pseudobulk','TYMP_mean_log1pCP10K']:
 x=rankdata(d[exposure]);Y=np.column_stack([rankdata(d[t]) for t in responses]);xr=x.copy();yr=Y.copy();groups=[np.asarray(v) for v in d.groupby('study').groups.values()]
 for ix in groups:xr[ix]-=xr[ix].mean();yr[ix]-=yr[ix].mean(axis=0)
 den=np.sqrt((xr*xr).sum()*(yr*yr).sum(axis=0));obs=xr@yr/den;hits=np.zeros(3,int)
 for _ in range(nperm):
  xp=xr.copy()
  for ix in groups:xp[ix]=rng.permutation(xp[ix])
  r=xp@yr/den;hits+=(np.abs(r)>=np.abs(obs)-1e-14)
 for j,target in enumerate(responses):
  rr,pp=spearmanr(d[exposure],d[target]);out.append(dict(exposure=exposure,target=target,adjustment='raw',n_patients=len(d),rho=rr,P=pp))
  out.append(dict(exposure=exposure,target=target,adjustment='within_study_rank_permutation',n_patients=len(d),rho=obs[j],P=(hits[j]+1)/(nperm+1)))
out=pd.DataFrame(out);out['FDR_three_endpoints']=np.nan
for _,ix in out.groupby(['exposure','adjustment']).groups.items():out.loc[ix,'FDR_three_endpoints']=false_discovery_control(out.loc[ix,'P'])
out.to_csv(O/'TYMP_C4_AND_SPP1_RESULTS.tsv',sep='\t',index=False)
study=[]
for name,z in d.groupby('study'):
 for target in responses:
  rr=pp=np.nan
  if len(z)>=3 and z.TYMP_pseudobulk.nunique()>1 and z[target].nunique()>1:rr,pp=spearmanr(z.TYMP_pseudobulk,z[target])
  study.append(dict(study=name,target=target,n_patients=len(z),rho=rr,P=pp,reporting_tier='exploratory' if len(z)<10 else 'patient_association'))
st=pd.DataFrame(study);st['FDR_three_endpoints']=np.nan
for _,ix in st.groupby('study').groups.items():
 p=st.loc[ix,'P'].dropna()
 if len(p):st.loc[p.index,'FDR_three_endpoints']=false_discovery_control(p.to_numpy())
st.to_csv(O/'STUDY_STRATIFIED_RESULTS.tsv',sep='\t',index=False)
loo=[]
for typ,col in [('patient','patient'),('study','study')]:
 for value in d[col].unique():
  z=d[d[col]!=value]
  for target in responses:
   rr=float(spearmanr(z.TYMP_pseudobulk,z[target]).statistic)
   loo.append(dict(omission_type=typ,omitted=value,target=target,n_patients=len(z),rho=rr))
pd.DataFrame(loo).to_csv(O/'LEAVE_ONE_PATIENT_STUDY_DIRECTIONS.tsv',sep='\t',index=False)
# Independent summary checks from per-cell grouped counts.
group=wanted.groupby(['sample_id','lineage'],observed=True).agg(TYMP=('TYMP_counts','sum'),SPP1=('SPP1_counts','sum'),library=('recomputed_library_counts','sum'),cells=('study','size'))
checks=[]
for _,row in d.iterrows():
 ma=group.loc[(row.sample_id,'Malignant')];mc=group.loc[(row.sample_id,'Macrophage')]
 checks.append(dict(check=row.sample_id+'_pseudobulk',passed=bool(np.isclose(np.log2(ma.TYMP/ma.library*1e6+1),row.TYMP_pseudobulk) and np.isclose(np.log2(mc.SPP1/mc.library*1e6+1),row.SPP1_pseudobulk))))
 checks.append(dict(check=row.sample_id+'_C4_fraction',passed=bool(np.isclose(row.C4_cells/row.macrophage_cells,row.C4_fraction))))
assert all(c['passed'] for c in checks)
(O/'AGGREGATION_VALIDATION.json').write_text(json.dumps(dict(status='PASS',checks=checks,raw_library_count_identity_pass=True),indent=2))
sources=[]
for path in [M,F]:
 sources.append(dict(path=str(path),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
# Global H5AD SHA was previously validated; retain actual location/stat, no false claim of newly hashing 1.8GB here.
sources.append(dict(path=str(H),bytes=H.stat().st_size,mtime_ns=H.stat().st_mtime_ns,prior_validated_sha256='e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7',new_full_hash=False))
(O/'SOURCE_MANIFEST.json').write_text(json.dumps(sources,indent=2))
status=dict(status='COMPUTATION_COMPLETE_AWAITING_REVIEW',patients=16,libraries=16,selected_cells=len(wanted),C4_zero_patients=int((d.C4_cells==0).sum()),no_statistical_significance_gate=True,figures_created=False,model_retrained=False,time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print('PRIMARY_AND_SENSITIVITY');print(out.to_string(index=False));print('STUDIES');print(st.to_string(index=False))
print('LOPO_C4');print(pd.DataFrame(loo).query('target == "C4_fraction"').to_string(index=False))
print(json.dumps(status))

