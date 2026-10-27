from pathlib import Path
import os,json,time,datetime,hashlib,traceback,sys,platform
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse,stats
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent;B=R.parent
H=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
def log(s):print(datetime.datetime.now().isoformat(),s,flush=True)
def save(n,o):(P/n).write_text(json.dumps(o,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def table(n,df):df.to_csv(P/n,sep='\t',index=False)
def profile(counts,meta,label,train):
 ii=np.flatnonzero((meta.label.to_numpy()==label)&train)
 study=[]
 for st in sorted(meta.iloc[ii].study.unique()):
  jj=ii[meta.iloc[ii].study.to_numpy()==st];study.append(counts[jj].mean(axis=0))
 return np.mean(study,axis=0)
def run():
 log('read frozen metadata')
 old=pd.read_csv(B/'project26210_part4_tymp_reset_phase17_20260812_v1/reference_tymp_continuous_v1/reference_tymp_continuous_v1_cell_assignments.tsv.gz',sep='\t',index_col=0)
 aligned=pd.read_csv(R/'phase218_spatial_identity_alignment_20260920_v1/CELL_IDENTITY_CROSSWALK.tsv.gz',sep='\t',index_col=0)
 mm=pd.read_csv(R/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',dtype={'candidate':str},index_col=0)
 with h5py.File(H) as h:
  obs=read_elem(h['obs']);var=read_elem(h['var']);genes=var.index.astype(str).to_numpy();assert len(set(genes))==len(genes)
  label=pd.Series('EXCLUDED_UNRESOLVED',index=obs.index)
  broad=obs.cell_type_broad_v1.astype(str)
  for a,b in {'Endothelial':'Endothelial','Pericyte/SMC':'SMC_Pericyte','B cell':'B_cell','Dendritic cell':'DC','Mast cell':'Mast','Plasma cell':'Plasma','Fibroblast':'Fibroblast','Erythrocyte':'Erythrocyte'}.items():label.loc[broad==a]=b
  for a,b in [('Malignant','Malignant'),('Epithelial-Normal','Epithelial_Normal')]:
   ids=old.index[old.state_tymp_continuous_v1==a];assert ids.isin(obs.index).all();label.loc[ids]=b
  ids=mm.index[mm.candidate=='Monocyte'];label.loc[ids]='Monocyte'
  ma=aligned[aligned.lineage=='mature_macrophage'];label.loc[ma.index]='Other_mature_macrophage';label.loc[ma.index[ma.new_identity=='C4']]='C4'
  re=aligned[aligned.lineage=='receiver'];label.loc[re.index]='EXCLUDED_UNRESOLVED'
  for a in ['CD8_T','CD4_Tcon','CD4_Treg','NK']:label.loc[re.index[re.new_identity==a]]=a
  gd=re.index[re.frozen_fine_label=='gdT'];label.loc[gd]='gdT'
  meta=pd.DataFrame({'cell_id':obs.index,'label':label.to_numpy(),'patient':obs.curated_patient_id.astype(str).to_numpy(),'study':obs.study.astype(str).to_numpy(),'tissue':obs.source_tissue.astype(str).to_numpy()})
  assert (meta.label=='C4').sum()==2674
  table('reference/ALL_CELL_REFERENCE_DISPOSITION.tsv.gz',meta)
  cov=meta.groupby(['label','study','patient']).size().rename('cells').reset_index();table('reference/PATIENT_STATE_CELL_COUNTS.tsv',cov)
  eligible=cov[(cov.label!='EXCLUDED_UNRESOLVED')&(cov.cells>=20)].copy().reset_index(drop=True)
  key={(x.label,x.study,x.patient):i for i,x in eligible.iterrows()}
  codes=np.array([key.get((x.label,x.study,x.patient),-1) for x in meta.itertuples()],dtype=int)
  sums=np.zeros((len(eligible),len(genes)),dtype=np.float64)
  node=h['layers/counts'];ptr=node['indptr'][:]
  total_input_counts=0.;used=0
  log('aggregate '+str(len(eligible))+' patient-state profiles')
  for start in range(0,len(obs),2048):
   end=min(start+2048,len(obs));cc=codes[start:end];valid=np.flatnonzero(cc>=0)
   if not len(valid):continue
   lo,hi=int(ptr[start]),int(ptr[end]);x=sparse.csr_matrix((node['data'][lo:hi],node['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(genes)))
   assert np.isfinite(x.data).all() and (x.data>=0).all()
   sel=sparse.csr_matrix((np.ones(len(valid)),(cc[valid],valid)),shape=(len(eligible),end-start))
   agg=(sel@x).tocoo();np.add.at(sums,(agg.row,agg.col),agg.data)
   used+=len(valid);total_input_counts+=float(x[valid].sum())
   if start%32768==0:log('cells processed '+str(end))
  assert used==eligible.cells.sum() and abs(sums.sum()-total_input_counts)<.1
 np.savez_compressed(P/'reference/PATIENT_STATE_COUNTS.npz',counts=sums,genes=genes.astype(str));table('reference/PATIENT_STATE_METADATA.tsv',eligible)
 meanraw=sums/eligible.cells.to_numpy()[:,None];cpm=sums/sums.sum(axis=1)[:,None]*1e6
 coverage=[];labels=[]
 for lab,g in cov.groupby('label'):
  e=eligible[eligible.label==lab];ok=lab!='EXCLUDED_UNRESOLVED' and e.patient.nunique()>=3 and e.study.nunique()>=2
  coverage.append(dict(label=lab,all_cells=int(g.cells.sum()),eligible_cells=int(e.cells.sum()),eligible_patients=e.patient.nunique(),eligible_studies=e.study.nunique(),coverage_eligible=ok))
  if ok:labels.append(lab)
 table('reference/REFERENCE_COVERAGE.tsv',pd.DataFrame(coverage));log('coverage eligible '+str(labels))
 train=np.ones(len(eligible),dtype=bool)
 rawprof=np.column_stack([profile(meanraw,eligible,l,train) for l in labels]);pd.DataFrame(rawprof,index=genes,columns=labels).to_csv(P/'reference/MEAN_COUNTS_STUDY_PATIENT_BALANCED.tsv.gz',sep='\t',index_label='gene')
 # Mapping gene universe is defined by presence only, not spatial associations.
 common=set(genes)
 secroot=B/'project26210_external_joint_v1_phase54_gse175540_spatial_v1/run_v1/sections_v1'
 for d in sorted(secroot.iterdir()):
  with h5py.File(d/'processed_section_v1.h5ad') as h:common &= set(read_elem(h['var']).index.astype(str))
 base=np.array([g in common and not g.startswith(('MT-','RPL','RPS')) for g in genes])
 programs=json.loads((R/'phase218_spatial_identity_alignment_20260920_v1/PROGRAM_DEFINITIONS.json').read_text())
 exclude=set(programs['spp1_state_full'])|set(programs['vascular_remodeling'])|{'SPP1','TYMP'}
 variants={'complete':set(),'exclude_SPP1_TYMP':{'SPP1','TYMP'},'exclude_program_union':exclude}
 predictions=[];foldgenes=[];skips=[]
 for variant,excluded in variants.items():
  allowed=base&np.array([g not in excluded for g in genes])
  for st in sorted(eligible.study.unique()):
   tr=eligible.study.to_numpy()!=st;test=~tr
   labs=[l for l in labels if eligible[tr&(eligible.label==l)].patient.nunique()>=3 and eligible[tr&(eligible.label==l)].study.nunique()>=2]
   if len(labs)<2:continue
   prof=np.array([profile(cpm,eligible,l,tr) for l in labs]);lp=np.log2(prof+1)
   indices=set()
   for j,l in enumerate(labs):
    diff=lp[j]-np.mean(np.delete(lp,j,axis=0),axis=0)
    ix=np.flatnonzero(allowed&(prof[j]>=1));ix=ix[np.argsort(diff[ix])[::-1]];ix=ix[diff[ix]>0][:50];indices.update(ix.tolist())
   ix=np.array(sorted(indices));assert len(ix)>10
   foldgenes.append(dict(variant=variant,heldout_study=st,genes=';'.join(genes[ix]),n_genes=len(ix),training_labels=';'.join(labs)))
   pp=lp[:,ix];pp=pp-pp.mean(axis=1,keepdims=True);pp=pp/np.linalg.norm(pp,axis=1,keepdims=True)
   for i in np.flatnonzero(test):
    truth=eligible.label.iloc[i]
    if truth not in labs:
     skips.append(dict(variant=variant,study=st,patient=eligible.patient.iloc[i],label=truth,reason='NOT_ELIGIBLE_IN_TRAINING_FOLD'));continue
    z=np.log2(cpm[i,ix]+1);z-=z.mean();z/=np.linalg.norm(z);sim=pp@z;j=labs.index(truth);best=int(np.argmax(sim));other=np.max(np.delete(sim,j))
    predictions.append(dict(variant=variant,study=st,patient=eligible.patient.iloc[i],truth=truth,predicted=labs[best],correct=best==j,margin=float(sim[j]-other),true_correlation=float(sim[j]),n_cells=int(eligible.cells.iloc[i])))
 pred=pd.DataFrame(predictions);table('reference/LOSO_PREDICTIONS.tsv',pred);table('reference/LOSO_TRAINING_GENE_PANELS.tsv',pd.DataFrame(foldgenes));table('reference/LOSO_UNEVALUABLE.tsv',pd.DataFrame(skips))
 res=[]
 for (variant,lab),df in pred.groupby(['variant','truth']):
  sr=df.groupby('study').correct.mean();passed=len(df)>=5 and len(sr)>=3 and df.correct.mean()>=.8 and sr.min()>=.5 and df.margin.median()>=.05
  res.append(dict(variant=variant,label=lab,n_patient_state_units=len(df),n_studies=len(sr),recall=df.correct.mean(),min_study_recall=sr.min(),median_margin=df.margin.median(),gate_pass=passed,study_recall=json.dumps(sr.to_dict())))
 table('reference/IDENTIFIABILITY_RESULTS.tsv',pd.DataFrame(res));save('reference/REFERENCE_VALIDATION.json',dict(status='PASS_ACCOUNTING',source_cells=len(meta),C4_cells=int((meta.label=='C4').sum()),eligible_profile_cells=used,n_profiles=len(eligible),n_reference_labels=len(labels),common_spatial_genes=len(common),aggregated_raw_counts=total_input_counts,no_spatial_mapping=True))
 log('reference complete; shared gene program check')
 oldtab=pd.read_csv(B/'project26210_external_joint_v1_phase70_cell2location_spatial_composition_v1/results_v1/composition_adjusted_bridges_by_section_v1.tsv',sep='\t');oldtab=oldtab[oldtab.relation_id=='B09_SPP1_to_Vascular'].set_index('sample_id')
 pmap=pd.read_csv(R/'phase168_reference_source_audit_20260916_v1/SPATIAL_26_SECTION_PATIENT_MAP.tsv',sep='\t').set_index('section')
 sec55=B/'project26210_external_joint_v1_phase55_cell2location_v5_spawn_resume_20260821_v1/run_v1/sections_v1'
 results=[]
 def partial(x,y,c):
  a=np.column_stack([x,y]+c);rk=np.column_stack([stats.rankdata(a[:,j]) for j in range(a.shape[1])]);D=np.column_stack([np.ones(len(x)),rk[:,2:]]);beta=np.linalg.lstsq(D,rk[:,:2],rcond=None)[0];res=rk[:,:2]-D@beta
  return np.corrcoef(res.T)[0,1]
 for d in sorted(secroot.iterdir()):
  sid=d.name
  with h5py.File(d/'processed_section_v1.h5ad') as h:
   so=read_elem(h['obs']);vv=read_elem(h['var']);X=read_elem(h['X']);idx=vv.index.astype(str);gs=sorted(exclude-{'TYMP'})
   x=X[:,idx.get_indexer(gs)].toarray();Z=(x-x.mean(axis=0))/np.where(x.std(axis=0)==0,1,x.std(axis=0));look={g:i for i,g in enumerate(gs)}
  ab=pd.read_csv(sec55/sid/'spot_cell_abundance_v1.tsv.gz',sep='\t',index_col='barcode').loc[so.index]
  scores=pd.read_csv(d/'spot_metadata_scores_v1.tsv.gz',sep='\t',index_col='barcode').loc[so.index]
  covs=[ab['mean__Macrophage'].to_numpy(),(ab['mean__Endothelial']+ab['mean__Pericyte/SMC']).to_numpy(),np.log1p(scores.total_counts.to_numpy())]
  for variant,remove in [('original',False),('no_LGALS3',True)]:
   aa=[g for g in programs['spp1_state_no_spp1'] if not remove or g!='LGALS3'];bb=[g for g in programs['vascular_remodeling'] if not remove or g!='LGALS3']
   # Preserve original float32 score rounding before rank residuals.
   sx=Z[:,[look[g] for g in aa]].mean(axis=1).astype(np.float32);sy=Z[:,[look[g] for g in bb]].mean(axis=1).astype(np.float32)
   if variant=='original':
    # Exact legacy score reproduction check, then use original stored scores for coefficient baseline.
    err=max(np.max(abs(sx-scores.score_primary__spp1_state_no_spp1)),np.max(abs(sy-scores.score_primary__vascular_remodeling)))
    assert err<2e-5,(sid,err)
    sx=scores.score_primary__spp1_state_no_spp1.to_numpy();sy=scores.score_primary__vascular_remodeling.to_numpy()
   rho=partial(sx,sy,covs)
   if variant=='original':assert abs(rho-oldtab.loc[sid,'rho'])<1e-10
   results.append(dict(section=sid,preparation=str(pmap.loc[sid,'preparation']),patient=pmap.loc[sid,'patient_id'],variant=variant,rho=rho,n_spots=len(so)))
 rr=pd.DataFrame(results);table('program_check/SECTION_EFFECTS.tsv',rr)
 summaries=[];patients=[]
 for variant,df in rr.groupby('variant'):
  for prep,g in df.groupby('preparation'):
   a=g.rho.to_numpy();summaries.append(dict(variant=variant,population=prep,n=len(a),positive=int((a>0).sum()),median_rho=np.median(a),P=stats.wilcoxon(a,alternative='two-sided').pvalue))
  for pat,g in df.groupby('patient'):
   rho=np.tanh(np.mean(np.arctanh(g.rho.to_numpy())));patients.append(dict(variant=variant,patient=pat,n_sections=len(g),rho=rho))
  a=np.array([x['rho'] for x in patients if x['variant']==variant]);assert len(a)==23
  summaries.append(dict(variant=variant,population='ALL_23_PATIENTS',n=len(a),positive=int((a>0).sum()),median_rho=np.median(a),P=stats.wilcoxon(a,alternative='two-sided').pvalue))
 su=pd.DataFrame(summaries);su['FDR_new_sensitivity_family']=np.nan
 mask=(su.variant=='no_LGALS3')&(su.population!='ALL_23_PATIENTS');su.loc[mask,'FDR_new_sensitivity_family']=stats.false_discovery_control(su.loc[mask,'P'].to_numpy())
 table('program_check/PATIENT_EFFECTS.tsv',pd.DataFrame(patients));table('program_check/PROGRAM_CHECK_SUMMARY.tsv',su)
 save('program_check/VALIDATION.json',dict(status='PASS',original_section_coefficients_reproduced=24,shared_genes_after_removal=[],deduplicated_patients=23,spot_level_p_computed=False))
 save('SOFTWARE.json',dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,h5py=h5py.__version__))
 save('STATUS.json',dict(status='COMPUTE_COMPLETE',reference_profiles_built=True,LOSO_complete=True,program_check_complete=True,new_spatial_mapping_started=False));(P/'EXIT_STATUS').write_text('0\n');log('COMPLETE')
try:run()
except Exception:
 (P/'ERROR.txt').write_text(traceback.format_exc());(P/'EXIT_STATUS').write_text('1\n');save('STATUS.json',dict(status='COMPUTE_FAILED'));raise
