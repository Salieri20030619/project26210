from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,h5py
from scipy import stats
from anndata.io import read_elem
P=Path(__file__).resolve().parents[1];R=P.parent;B=R.parent
assert (P/'EXIT_STATUS').read_text().strip()=='0'
meta=pd.read_csv(P/'reference/PATIENT_STATE_METADATA.tsv',sep='\t');z=np.load(P/'reference/PATIENT_STATE_COUNTS.npz');counts=z['counts'];genes=z['genes'];meanraw=counts/meta.cells.to_numpy()[:,None]
prof=pd.read_csv(P/'reference/MEAN_COUNTS_STUDY_PATIENT_BALANCED.tsv.gz',sep='\t',index_col=0)
assert list(prof.index)==list(genes) and np.isfinite(prof.to_numpy()).all() and (prof.to_numpy()>=0).all()
errors=[]
for label in prof:
 mask=meta.label==label;n_studies=meta.loc[mask,'study'].nunique();weights=np.zeros(len(meta))
 for st in meta.loc[mask,'study'].unique():
  idx=np.flatnonzero(mask&(meta.study==st));weights[idx]=1/(n_studies*len(idx))
 estimate=weights@meanraw;errors.append(float(np.max(abs(estimate-prof[label].to_numpy()))))
assert max(errors)<1e-8
allc=pd.read_csv(P/'reference/ALL_CELL_REFERENCE_DISPOSITION.tsv.gz',sep='\t',index_col='cell_id')
freeze=pd.read_csv(R/'phase218_spatial_identity_alignment_20260920_v1/CELL_IDENTITY_CROSSWALK.tsv.gz',sep='\t',index_col='cell_id')
assert (allc.label=='C4').sum()==2674 and set(allc.index[allc.label=='C4'])==set(freeze.index[freeze.new_identity=='C4'])
for label in ['CD8_T','CD4_Tcon','CD4_Treg','NK']:
 assert set(allc.index[allc.label==label])==set(freeze.index[freeze.new_identity==label])
assert set(allc.index[allc.label=='gdT_enriched'])==set(freeze.index[freeze.frozen_fine_label=='gdT_enriched'])
sourcehash=json.loads((R/'phase218_spatial_identity_alignment_20260920_v1/SOURCE_SHA256.json').read_text())
assert all(hashlib.sha256(Path(k).read_bytes()).hexdigest()==v for k,v in sourcehash.items())
programs=json.loads((R/'phase218_spatial_identity_alignment_20260920_v1/PROGRAM_DEFINITIONS.json').read_text());ex=set(programs['spp1_state_full'])|set(programs['vascular_remodeling'])|{'SPP1','TYMP'}
fold=pd.read_csv(P/'reference/LOSO_TRAINING_GENE_PANELS.tsv',sep='\t')
for row in fold.itertuples():
 geneset=set(row.genes.split(';'))
 if row.variant=='exclude_SPP1_TYMP':assert not geneset&{'SPP1','TYMP'}
 if row.variant=='exclude_program_union':assert not geneset&ex
 for label in row.training_labels.split(';'):
  tr=meta[(meta.study!=row.heldout_study)&(meta.label==label)];assert tr.patient.nunique()>=3 and tr.study.nunique()>=2
sec=pd.read_csv(P/'program_check/SECTION_EFFECTS.tsv',sep='\t');checks=[]
for sid in ['GSM5924030','GSM5924042']:
 d=B/'project26210_external_joint_v1_phase54_gse175540_spatial_v1/run_v1/sections_v1'/sid
 with h5py.File(d/'processed_section_v1.h5ad') as h:
  obs=read_elem(h['obs']);var=read_elem(h['var']);xx=read_elem(h['X']);gs=sorted((set(programs['spp1_state_no_spp1'])|set(programs['vascular_remodeling']))-{'LGALS3'});a=xx[:,var.index.get_indexer(gs)].toarray().astype(float)
  a=stats.zscore(a,axis=0,ddof=0);look={g:i for i,g in enumerate(gs)}
 x=a[:,[look[g] for g in programs['spp1_state_no_spp1'] if g!='LGALS3']].mean(axis=1).astype(np.float32)
 y=a[:,[look[g] for g in programs['vascular_remodeling'] if g!='LGALS3']].mean(axis=1).astype(np.float32)
 ab=pd.read_csv(B/'project26210_external_joint_v1_phase55_cell2location_v5_spawn_resume_20260821_v1/run_v1/sections_v1'/sid/'spot_cell_abundance_v1.tsv.gz',sep='\t',index_col='barcode').loc[obs.index]
 sc=pd.read_csv(d/'spot_metadata_scores_v1.tsv.gz',sep='\t',index_col='barcode').loc[obs.index]
 rank=np.column_stack([stats.rankdata(t) for t in [x,y,ab['mean__Macrophage'],ab['mean__Endothelial']+ab['mean__Pericyte/SMC'],np.log1p(sc.total_counts)]])
 precision=np.linalg.inv(np.corrcoef(rank.T));rho=-precision[0,1]/np.sqrt(precision[0,0]*precision[1,1]);expected=sec[(sec.section==sid)&(sec.variant=='no_LGALS3')].rho.iloc[0]
 assert abs(rho-expected)<1e-10;checks.append({'section':sid,'error':abs(rho-expected)})
summ=pd.read_csv(P/'program_check/PROGRAM_CHECK_SUMMARY.tsv',sep='\t');xx=summ[(summ.variant=='no_LGALS3')&(summ.population!='ALL_23_PATIENTS')];pv=xx.P.to_numpy();order=np.argsort(pv);q=np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1];assert np.allclose(q,xx.FDR_new_sensitivity_family.to_numpy()[order])
(P/'FINAL_VALIDATION.json').write_text(json.dumps({'status':'PASS','frozen_cells_and_sources_unchanged':True,'reference_balance_max_error':max(errors),'reference_labels':len(prof.columns),'gene_exclusion_folds_checked':len(fold),'independent_partial_correlations':checks,'BH_independently_recomputed':True,'raw_count_conservation':json.loads((P/'reference/REFERENCE_VALIDATION.json').read_text()),'computation_exit':0,'new_mapping_started':False},indent=2));print('FINAL_VALIDATION_PASS',len(prof.columns),'labels',len(fold),'folds')
