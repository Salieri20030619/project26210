from pathlib import Path
import json,hashlib,datetime,subprocess,importlib.metadata
import pandas as pd,numpy as np,h5py
from scipy import sparse
H=Path.home();R=H/'project26210_external_joint_v1';O=R/'phase213_immune_restart_execution_20260919_v1';D=O/'I1_patient_inputs';D.mkdir(exist_ok=False);F=R/'phase211_immune_receiver_refreeze_20260919_v1/freeze_v2';P=R/'phase211_immune_receiver_refreeze_20260919_v1';SRC=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
def save(df,name):df.to_csv(D/name,sep='\t',index=False)
with (D/'RESOURCE_GATE.log').open('w') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,stderr=f,check=True)
seal=json.loads((F/'FREEZE_SHA256.json').read_text());assert all(hashlib.sha256((F/k).read_bytes()).hexdigest()==v for k,v in seal.items())
(O/'AUTHORIZATION_AND_SCOPE.md').write_text('2026-09-19 用户授权：按I1-I5顺序逐个展开，用户离开期间无需汇报。执行与结果仅记服务器，常规问题自主修复。无新图、空间或临床模型；前两部分和C4不变。依据Phase212计划执行全部五阶段，不在中间阶段等待确认。\n')
meta=pd.read_csv(F/'DISCOVERY_RECEIVER_CELL_FREEZE.tsv.gz',sep='\t',index_col='cell_id');cache=np.load(P/'RAW_PANEL_CACHE.npz');assert np.array_equal(cache['cell_id'],meta.index);genes=list(cache['genes']);raw=cache['counts'];tot=cache['total_counts'];assert np.isfinite(raw).all() and (raw>=0).all() and np.allclose(tot,meta.total_counts_raw,rtol=0,atol=1e-5)
macsource=R/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz';mac=pd.read_csv(macsource,sep='\t',index_col=0,dtype={'candidate':str});mac=mac[~mac.candidate.isin(['1','Monocyte'])].copy();assert len(mac)==22444 and (mac.candidate=='4').sum()==2674
mac['SPP1_counts']=0.;mac['recomputed_total']=0.
with h5py.File(SRC,'r',rdcc_nbytes=128*1024*1024) as f:
 ids=f['obs/cell_id'].asstr()[:];gi=f['var/gene'].asstr()[:];idx=int(np.flatnonzero(gi=='SPP1')[0]);assert (gi=='SPP1').sum()==1;sel=pd.Index(ids).isin(mac.index);assert sel.sum()==len(mac);x=f['layers/counts'];ptr=x['indptr'][:]
 for start in range(0,len(ids),4096):
  end=min(start+4096,len(ids));mask=sel[start:end]
  if not mask.any():continue
  lo,hi=int(ptr[start]),int(ptr[end]);a=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(gi)))[mask];assert np.isfinite(a.data).all() and (a.data>=0).all() and np.equal(a.data,np.floor(a.data)).all()
  cells=ids[start:end][mask];mac.loc[cells,'SPP1_counts']=a[:,idx].toarray().ravel();mac.loc[cells,'recomputed_total']=np.asarray(a.sum(axis=1,dtype=np.float64)).ravel()
assert np.allclose(mac.recomputed_total,mac.total_counts_raw,rtol=0,atol=1e-5);save(mac.reset_index(),'MACROPHAGE_SPP1_CELL_COUNTS.tsv.gz')
keys=['study','curated_patient_id','sample_id','source_tissue'];macstats={}
for key,g in mac.groupby(keys,observed=True):
 c=g[g.candidate=='4'];other=g[g.candidate!='4'];z=dict(n_mature=len(g),n_C4=len(c),SPP1_counts=float(g.SPP1_counts.sum()),mature_total_counts=float(g.recomputed_total.sum()),C4_SPP1_counts=float(c.SPP1_counts.sum()),C4_total_counts=float(c.recomputed_total.sum()),mature_SPP1_pseudobulk=float(np.log2(1+g.SPP1_counts.sum()/g.recomputed_total.sum()*1e6)),C4_fraction=len(c)/len(g),C4_SPP1_pseudobulk=float(np.log2(1+c.SPP1_counts.sum()/c.recomputed_total.sum()*1e6)) if len(c) else np.nan,other_SPP1_pseudobulk=float(np.log2(1+other.SPP1_counts.sum()/other.recomputed_total.sum()*1e6)) if len(other) else np.nan)
 macstats[key]=z
programs=json.loads((F/'FREEZE_DEFINITION.json').read_text())['programs'];act=['FOS','JUN','JUNB','NR4A1','CD69','IL2RA'];assert all(g in genes for g in act)
log=np.log1p(raw/tot[:,None]*10000);cell_scores={ep:log[:,[genes.index(g) for g in gs]].mean(axis=1) for ep,gs in programs.items()};cell_scores['activation']=log[:,[genes.index(g) for g in act]].mean(axis=1)
E=pd.read_csv(F/'DISCOVERY_ENDPOINT_LIBRARY_FREEZE.tsv',sep='\t');allrows=[];gene_rows=[]
for row in E.itertuples():
 k=(row.study,row.patient,row.library,row.tissue);mm=macstats.get(k,{});assert mm.get('n_mature',0)==row.n_mature and mm.get('n_C4',0)==row.n_C4
 base=(meta.study==row.study)&(meta.curated_patient_id==row.patient)&(meta.sample_id==row.library)&(meta.source_tissue==row.tissue)
 lab=['CD4_Tcon','CD4_Treg'] if row.endpoint=='Treg_CD4' else (['CD8_T'] if row.endpoint.startswith('CD8') else ['NK']);base&=meta.receiver_label_v1.isin(lab)
 for variant in ['primary','exclude_concentrated','high_confidence_only']:
  m=base.copy()
  if variant=='exclude_concentrated':m &= ~meta.patient_concentrated_cluster
  if variant=='high_confidence_only':m &= meta.identity_confidence_v1.eq('high')
  pos=np.flatnonzero(m.to_numpy());n=len(pos);nT=int((meta.iloc[pos].receiver_label_v1=='CD4_Treg').sum());det={};cov=True
  if row.endpoint in programs:
   gs=programs[row.endpoint];ix=[genes.index(g) for g in gs];dc=(raw[np.ix_(pos,ix)]>0).sum(axis=0);det=dict(zip(gs,dc.tolist()));cov=(dc>=5).sum()>=3;value=float(cell_scores[row.endpoint][pos].mean()) if n else np.nan
   if variant=='primary':
    for gene,jj,nd in zip(gs,ix,dc):gene_rows.append(dict(study=row.study,patient=row.patient,library=row.library,endpoint=row.endpoint,gene=gene,receiver_cells=n,detected_cells=int(nd),mean_log1pCP10K=float(log[pos,jj].mean()) if n else np.nan))
  else:value=nT/n if n else np.nan
  ownelig=n>=(50 if row.endpoint=='Treg_CD4' else 20) and cov
  eligible=bool(row.eligible and ownelig);c4eligible=bool(row.C4_intrinsic_eligible and ownelig)
  if variant=='primary':assert n==row.n_receiver and eligible==row.eligible and c4eligible==row.C4_intrinsic_eligible
  z=dict(study=row.study,patient=row.patient,library=row.library,tissue=row.tissue,endpoint=row.endpoint,variant=variant,outcome=value,receiver_cells=n,Treg_cells=nT,receiver_mean_n_genes=float(meta.iloc[pos].n_genes_raw.mean()) if n else np.nan,receiver_mean_total_counts=float(tot[pos].mean()) if n else np.nan,activation=float(cell_scores['activation'][pos].mean()) if n else np.nan,eligible=eligible,C4_intrinsic_eligible=c4eligible,original_unresolved_fraction=row.unresolved_fraction,coverage_counts=json.dumps(det),exclusion_reasons=str(row.exclusion_reasons) if not row.eligible else ('derived_receiver_count_or_coverage_failure' if not ownelig else ''))
  z.update(mm);allrows.append(z)
lib=pd.DataFrame(allrows);save(lib,'LIBRARY_METRICS_ALL_VARIANTS.tsv');save(pd.DataFrame(gene_rows),'PROGRAM_GENE_COVERAGE.tsv')
reg=pd.read_csv(R/'phase212_immune_restart_plan_20260919_v1/HYPOTHESIS_REGISTER.tsv',sep='\t');patientrows=[]
for rr in reg.itertuples():
 for variant in ['primary','exclude_concentrated','high_confidence_only','common_C4_libraries']:
  b=lib[(lib.endpoint==rr.endpoint)&(lib.variant==('primary' if variant=='common_C4_libraries' else variant))];b=b[b.C4_intrinsic_eligible if (rr.predictor=='C4_SPP1_pseudobulk' or variant=='common_C4_libraries') else b.eligible]
  for (study,pt),g in b.groupby(['study','patient']):
   z=dict(hypothesis_id=rr.hypothesis_id,predictor=rr.predictor,endpoint=rr.endpoint,variant=variant,study=study,patient=pt,n_libraries=len(g),libraries=';'.join(sorted(g.library)),x=float(g[rr.predictor].mean()),y=float(g.outcome.mean()),receiver_mean_n_genes=float(g.receiver_mean_n_genes.mean()),receiver_mean_total_counts=float(g.receiver_mean_total_counts.mean()),activation=float(g.activation.mean()),mean_receiver_cells=float(g.receiver_cells.mean()),mean_mature_cells=float(g.n_mature.mean()),mean_C4_cells=float(g.n_C4.mean()))
   assert np.isfinite(z['x']) and np.isfinite(z['y']);patientrows.append(z)
patient=pd.DataFrame(patientrows);save(patient,'PATIENT_HYPOTHESIS_INPUTS.tsv')
checks={'freeze_sha256_all_match':True,'mature_total_count_identity':True,'C4_membership_count':True,'receiver_counts_and_eligibility_match':True,'no_duplicate_patient_hypothesis_variant':not patient.duplicated(['hypothesis_id','variant','study','patient']).any()}
fp=pd.read_csv(F/'DISCOVERY_PATIENT_ENDPOINT_FREEZE.tsv',sep='\t')
for rr in reg.itertuples():
 a=patient[(patient.hypothesis_id==rr.hypothesis_id)&(patient.variant=='primary')];assert a.patient.nunique()==rr.frozen_patient_n
 expected=fp[fp.endpoint==rr.endpoint];expected=expected[expected.C4_intrinsic_libraries.notna()] if rr.predictor=='C4_SPP1_pseudobulk' else expected[expected.eligible]
 assert set(a.patient)==set(expected.patient)
 for z in a.itertuples():
  f=expected[expected.patient==z.patient].iloc[0];assert set(z.libraries.split(';'))==set((f.C4_intrinsic_libraries if rr.predictor=='C4_SPP1_pseudobulk' else f.eligible_libraries).split(';'))
checks['all12_patient_sets_and_library_sets_exact']=True
# Independent scalar aggregation comparison across every primary patient input.
for z in patient[patient.variant=='primary'].itertuples():
 b=lib[(lib.variant=='primary')&(lib.endpoint==z.endpoint)&lib.library.isin(z.libraries.split(';'))&(lib.patient==z.patient)];assert abs(sum(b[z.predictor])/len(b)-z.x)<1e-12 and abs(sum(b.outcome)/len(b)-z.y)<1e-12
checks['patient_equal_library_aggregation']=True;assert all(checks.values());(D/'VALIDATION.json').write_text(json.dumps(checks,indent=2))
(D/'DATA_DICTIONARY.md').write_text('Patient is the independent biological unit. x and y use exactly the same endpoint-qualified tumour libraries, averaged with equal library weights. x=mature_SPP1 or C4_SPP1 log2(CPM+1) pseudobulk; C4_fraction=C4/mature. Functional y is cell-equal then gene-equal mean log1p(CP10K); Treg y=Treg/(CD4_Tcon+CD4_Treg). Variants are derived checks and never replace primary freeze. common_C4_libraries compares all predictors on identical C4-intrinsic-qualified libraries. No P/FDR calculated in I1.\n')
(D/'SOFTWARE.json').write_text(json.dumps({k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','h5py']},indent=2));(D/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_VALIDATED','time':datetime.datetime.now().astimezone().isoformat(),'checks':checks,'association_tests_run':False,'primary_patient_counts':patient[patient.variant=='primary'].groupby('hypothesis_id').patient.nunique().to_dict()},indent=2));print((D/'STATUS.json').read_text(),flush=True)
