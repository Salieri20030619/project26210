from pathlib import Path
import json,hashlib,subprocess,datetime,importlib.metadata
import numpy as np,pandas as pd,h5py
from scipy import sparse
from anndata.io import read_elem
H=Path.home();R=H/'project26210_external_joint_v1';O=R/'phase210_immune_identity_pairing_check_20260919_v1'
P=Path('/mnt/BioDisk/ludingyu/project26210');SRC=P/'results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad';OLD=P/'object/01_annotation_adata_annotated.h5ad';MAC=R/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz';C4=R/'part2_rebuild_v1/09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz'
def save(d,n):d.to_csv(O/n,sep='\t',index=False)
def msg(s):print(s,flush=True)
plan={'scope':'Identity and pairing audit only; no new clustering, model, differential test or association. Legacy labels are provisional candidates, not certified identities.','created':datetime.datetime.now().astimezone().isoformat(),'unit':'study-qualified patient; same tumour library pairing before patient aggregation','mature_macrophage':'Frozen regulon metadata excluding candidate 1 and Monocyte, inherited Phase154 rule','minimums_for_readiness_screen':{'mature_macrophage':50,'CD8_or_NK_receiver':20,'CD4_denominator':50,'C4_expression':20},'threshold_comparison':'receiver counts 10/20/50; no outcome used to choose threshold; all are readiness screens','legacy_mapping':'Exact cell_id only; no barcode truncation; majority label is contextual, not independent replication','marker_evidence':'T backbone >=2 of CD3D CD3E TRAC; CD8A or CD8B; CD4 or IL7R with inherited helper candidate; Treg FOXP3 plus IL2RA or CTLA4 with T backbone; NK >=2 of KLRD1 NCR1 KLRF1 NCAM1 FCGR3A with T backbone <2. Detection absence alone not proof of absence.','mixed_flags':'>=2 epithelial EPCAM KRT8 KRT18 or myeloid LST1 LYZ FCER1G TYROBP genes AND mean log1pCP10K exceeds lineage backbone; review flag, not proven doublet','program_readiness':'At least 3 program genes detected in >=5 cells per library, at least 20 receiver cells; gene coverage only, no function score or inference','CD4_denominator':'Provisional inherited CD4 helper/cytotoxic + Treg labels; unresolved/gamma-delta/MAIT/NKT excluded and reported separately','C4_zero':'Keep zero C4 in fraction readiness if mature macrophage denominator qualifies; C4-intrinsic SPP1 requires >=20 C4','limitations':'Marker-positive subset is a strict check and cannot be treated as unbiased full population; old cytotoxic naming may preselect state. No clinical or causal claims.'}
(O/'PLAN_BEFORE_CHECK.json').write_text(json.dumps(plan,indent=2))
with (O/'MEMORY_GATE.log').open('w') as f:subprocess.run([str(P/'check_memory_gate.sh')],stdout=f,stderr=f,check=True)
fields=['study','curated_patient_id','sample_id','source_tissue','cell_type_broad_v1','annotation_confidence_v1','annotation_ambiguous_v1','total_counts_raw','n_genes_raw','scrublet_predicted_doublet','low_quality']
with h5py.File(SRC) as f:
 ids=f['obs']['cell_id'].asstr()[:];obs=pd.DataFrame({k:np.asarray(read_elem(f['obs'][k])) for k in fields},index=pd.Index(ids,name='cell_id'));genes=f['var']['gene'].asstr()[:]
assert obs.index.is_unique
with h5py.File(OLD) as f:
 old=read_elem(f['obs'])
assert old.index.is_unique
selected=obs.cell_type_broad_v1.isin(['T cell','NK cell']);d=obs[selected].copy();assert d.index.isin(old.index).all()
for k in ['cell_type','celltypist_label_raw','celltypist_label_majority','marker_label_v2','annotation_confidence','annotation_disagreement']:
 d['legacy_'+k]=old[k].reindex(d.index).astype(str)
assert (d.study.astype(str)==old.study.reindex(d.index).astype(str)).all()
assert (d.sample_id.astype(str)==old.sample_id.reindex(d.index).astype(str)).all()
maprows=[]
def typ(label):
 if label in ['Regulatory T cells','Treg(diff)']:return 'Treg'
 if label in ['CD8a/a','CD8a/b(entry)','Tcm/Naive cytotoxic T cells','Tem/Temra cytotoxic T cells','Tem/Trm cytotoxic T cells','Trm cytotoxic T cells']:return 'CD8_candidate'
 if label in ['Tcm/Naive helper T cells','Tem/Effector helper T cells','Tem/Effector helper T cells PD1+','Memory CD4+ cytotoxic T cells','Follicular helper T cells','Type 1 helper T cells','Type 17 helper T cells']:return 'CD4_candidate'
 if label in ['NK cells','CD16+ NK cells','CD16- NK cells','Cycling NK cells','Transitional NK']:return 'NK_candidate'
 if label in ['gamma-delta T cells','CRTAM+ gamma-delta T cells','MAIT cells','NKT cells']:return 'unconventional_T'
 if label=='Cycling T cells':return 'cycling_T_unresolved'
 return 'unresolved'
d['legacy_candidate']=d.legacy_celltypist_label_raw.map(typ);d['majority_candidate']=d.legacy_celltypist_label_majority.map(typ)
for label in sorted(d.legacy_celltypist_label_raw.unique()):maprows.append({'legacy_label':label,'candidate':typ(label)})
save(pd.DataFrame(maprows),'LEGACY_LABEL_MAPPING.tsv')
panels={'T_backbone':['CD3D','CD3E','TRAC'],'CD8_identity':['CD8A','CD8B'],'CD4_support':['CD4','IL7R'],'Treg_identity':['FOXP3','IL2RA','CTLA4'],'NK_identity':['KLRD1','NCR1','KLRF1','NCAM1','FCGR3A'],'myeloid':['LST1','LYZ','FCER1G','TYROBP'],'epithelial':['EPCAM','KRT8','KRT18'],'cytotoxicity':['NKG7','GNLY','GZMB','PRF1','CTSW'],'dysfunction':['PDCD1','CTLA4','LAG3','HAVCR2','TIGIT','TOX','CXCL13','ENTPD1']}
panel=sorted(set(sum(panels.values(),[])));assert all((genes==g).sum()==1 for g in panel);idx=[np.flatnonzero(genes==g)[0] for g in panel];counts=np.zeros((len(d),len(panel)),dtype=np.float64);totals=np.zeros(len(d));offset=0
msg('extracting selected raw counts: '+str(len(d)))
with h5py.File(SRC,'r',rdcc_nbytes=128*1024*1024) as f:
 x=f['layers/counts'];assert x.attrs['encoding-type']=='csr_matrix';ptr=x['indptr'][:]
 for start in range(0,len(obs),4096):
  end=min(start+4096,len(obs));mask=selected.iloc[start:end].to_numpy()
  if not mask.any():continue
  lo,hi=int(ptr[start]),int(ptr[end]);a=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(genes)))[mask]
  assert np.isfinite(a.data).all() and (a.data>=0).all() and np.equal(a.data,np.floor(a.data)).all()
  n=a.shape[0];counts[offset:offset+n]=a[:,idx].toarray();totals[offset:offset+n]=np.asarray(a.sum(axis=1)).ravel();offset+=n
assert offset==len(d) and np.allclose(totals,d.total_counts_raw,rtol=0,atol=1e-5)
log=np.log1p(counts/totals[:,None]*10000);raw=pd.DataFrame(counts,index=d.index,columns=panel)
for name,gs in panels.items():
 ix=[panel.index(g) for g in gs];d[name+'_detected']=(counts[:,ix]>0).sum(axis=1)
 if name not in ['cytotoxicity','dysfunction']:d[name+'_meanlog']=log[:,ix].mean(axis=1)
T=d.T_backbone_detected>=2;NK=d.NK_identity_detected>=2
specific=np.maximum(d.T_backbone_meanlog,d.NK_identity_meanlog)
d['mixed_lineage_review']=((d.myeloid_detected>=2)&(d.myeloid_meanlog>specific))|((d.epithelial_detected>=2)&(d.epithelial_meanlog>specific))
d['marker_supported']=False
for label,mask in [('CD8_candidate',T&(d.CD8_identity_detected>=1)),('CD4_candidate',T&(d.CD4_support_detected>=1)),('Treg',T&(raw.FOXP3>0)&((raw.IL2RA>0)|(raw.CTLA4>0))),('NK_candidate',NK&~T)]:
 d.loc[d.legacy_candidate==label,'marker_supported']=mask[d.legacy_candidate==label]
d['broad_compatible']=((d.cell_type_broad_v1=='T cell')&d.legacy_candidate.isin(['CD8_candidate','CD4_candidate','Treg','unconventional_T','cycling_T_unresolved']))|((d.cell_type_broad_v1=='NK cell')&(d.legacy_candidate=='NK_candidate'))
d['strict_identity_supported']=d.marker_supported&d.broad_compatible&~d.mixed_lineage_review
# No overwriting source labels; high dropout / mixed cells remain in audit.
save(d.reset_index(),'CELL_IDENTITY_AUDIT.tsv.gz');save(raw.reset_index(),'CELL_PANEL_RAW_COUNTS.tsv.gz')
rows=[]
for keys,g in d.groupby(['source_tissue','cell_type_broad_v1','legacy_candidate'],observed=True):
 rows.append(dict(tissue=keys[0],broad=keys[1],candidate=keys[2],cells=len(g),patients=g.curated_patient_id.nunique(),libraries=g.sample_id.nunique(),marker_supported=int(g.marker_supported.sum()),strict_supported=int(g.strict_identity_supported.sum()),mixed_review=int(g.mixed_lineage_review.sum()),raw_majority_agree=int((g.legacy_candidate==g.majority_candidate).sum())))
save(pd.DataFrame(rows),'IDENTITY_SUMMARY.tsv')
rows=[]
for keys,g in d.groupby(['source_tissue','legacy_candidate'],observed=True):
 a=raw.loc[g.index]
 for gene in panel:rows.append(dict(tissue=keys[0],candidate=keys[1],gene=gene,cells=len(g),detected=int((a[gene]>0).sum()),fraction_detected=float((a[gene]>0).mean()),mean_log1pCP10K=float(np.log1p(a[gene]/g.total_counts_raw*10000).mean())))
save(pd.DataFrame(rows),'MARKER_EVIDENCE_BY_CANDIDATE.tsv')
mac=pd.read_csv(MAC,sep='\t',index_col=0,dtype={'candidate':str});mature=mac[~mac.candidate.isin(['1','Monocyte'])].copy();c4=pd.read_csv(C4,sep='\t').cell_id;assert len(c4)==2674 and set(c4)==set(mature.index[mature.candidate=='4']);assert set(mature.index)<=set(obs.index)
assert not set(mature.index)&set(d.index)
keys=['study','curated_patient_id','sample_id','source_tissue'];tab=obs.groupby(keys,observed=True).size().rename('all_cells').to_frame();tab['mature_macrophage']=mature.groupby(keys,observed=True).size();tab['C4']=mature[mature.candidate=='4'].groupby(keys,observed=True).size()
for lab in ['CD8_candidate','CD4_candidate','Treg','NK_candidate','unconventional_T','cycling_T_unresolved','unresolved']:
 a=d[d.legacy_candidate==lab];tab[lab]=a.groupby(keys,observed=True).size();tab[lab+'_strict']=a[a.strict_identity_supported].groupby(keys,observed=True).size()
tab=tab.fillna(0).astype(int).reset_index();tab['CD4_denominator_candidate']=tab.CD4_candidate+tab.Treg;tab['CD4_denominator_strict']=tab.CD4_candidate_strict+tab.Treg_strict
save(tab,'ALL_LIBRARY_IDENTITY_COUNTS.tsv')
tum=tab[tab.source_tissue=='Tumor'].copy();qual=[]
for _,row in tum.iterrows():
 for endpoint,cc in [('CD8','CD8_candidate'),('NK','NK_candidate'),('Treg_CD4','CD4_denominator_candidate')]:
  for threshold in ([10,20,50] if endpoint!='Treg_CD4' else [20,50]):
   strictcol={'CD8':'CD8_candidate_strict','NK':'NK_candidate_strict','Treg_CD4':'CD4_denominator_strict'}[endpoint]
   qual.append(dict(study=row.study,patient=row.curated_patient_id,library=row.sample_id,endpoint=endpoint,receiver_threshold=threshold,n_mature=row.mature_macrophage,n_C4=row.C4,n_receiver_candidate=row[cc],n_receiver_strict=row[strictcol],candidate_count_ready=bool(row.mature_macrophage>=50 and row[cc]>=threshold),strict_count_ready=bool(row.mature_macrophage>=50 and row[strictcol]>=threshold),C4_expression_count_ready=bool(row.C4>=20 and row[cc]>=threshold),zero_C4=bool(row.C4==0)))
q=pd.DataFrame(qual);save(q,'SAME_TUMOR_LIBRARY_READINESS.tsv')
summary=[]
for keys,g in q.groupby(['endpoint','receiver_threshold']):
 for version,col in [('candidate','candidate_count_ready'),('strict_check','strict_count_ready')]:
  a=g[g[col]];summary.append(dict(endpoint=keys[0],threshold=keys[1],identity_version=version,libraries=len(a),patients=a.patient.nunique(),zero_C4_libraries=int(a.zero_C4.sum()),studies=a.study.nunique()))
save(pd.DataFrame(summary),'READINESS_SUMMARY.tsv')
coverage=[]
for (study,patient,lib),g in d[d.source_tissue=='Tumor'].groupby(['study','curated_patient_id','sample_id']):
 for label in ['CD8_candidate','NK_candidate']:
  a=g[g.legacy_candidate==label]
  for program in (['cytotoxicity','dysfunction'] if label=='CD8_candidate' else ['cytotoxicity']):
   gs=panels[program];n_detect=(raw.loc[a.index,gs]>0).sum(axis=0)
   coverage.append(dict(study=study,patient=patient,library=lib,candidate=label,program=program,receiver_cells=len(a),genes=len(gs),genes_detected_5cells=int((n_detect>=5).sum()),coverage_ready=bool(len(a)>=20 and (n_detect>=5).sum()>=3),gene_detection_counts=json.dumps(n_detect.to_dict())))
save(pd.DataFrame(coverage),'PROGRAM_EXPRESSION_COVERAGE.tsv')
# Multiple libraries are not independent patients: record the required pairing rule.
patient_rows=[]
for keys,g in q.groupby(['study','patient','endpoint','receiver_threshold']):
 a=g[g.candidate_count_ready];patient_rows.append(dict(study=keys[0],patient=keys[1],endpoint=keys[2],threshold=keys[3],eligible_library_count=len(a),eligible_libraries=';'.join(a.library),proposed_aggregation='Equal weight eligible tumour libraries within patient; do not add nonpaired libraries; not executed'))
save(pd.DataFrame(patient_rows),'PATIENT_PAIRING_REGISTER.tsv')
# Identity provenance and input integrity.
inputs=[]
for p in [SRC,OLD,MAC,C4]:inputs.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest() if p.stat().st_size<20000000 else 'large_file_not_rehashed'))
(O/'INPUT_MANIFEST.json').write_text(json.dumps(inputs,indent=2));(O/'SOFTWARE.json').write_text(json.dumps({k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','anndata','h5py']},indent=2))
res={'status':'COMPLETE_IDENTITY_AND_READINESS_AUDIT_NOT_FINAL_ANNOTATION','TNK_cells':len(d),'exact_legacy_ID_matches':len(d),'study_library_matches':True,'raw_count_totals_match':True,'mature_macrophage_cells':len(mature),'C4_matches':len(c4),'no_new_association':True,'identity_status':'PROVISIONAL_INHERITED_CANDIDATES_WITH_MARKER_AUDIT'}
(O/'STATUS.json').write_text(json.dumps(res,indent=2));msg(json.dumps(res));msg(pd.DataFrame(summary).to_string(index=False))
