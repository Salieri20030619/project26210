from pathlib import Path
import json,hashlib,datetime,shutil
import pandas as pd,numpy as np
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1'
O=P/'phase149_part2_scope_initial_comparison_20260915_v1';O.mkdir(exist_ok=False)
E=O/'evidence';E.mkdir();sources=[]
def read(path):
 path=Path(path);sources.append(path);return pd.read_csv(path,sep='\t',low_memory=False)
def dump(df,name):df.to_csv(E/name,sep='\t',index=False)
def js(path):
 path=Path(path);sources.append(path);return json.loads(path.read_text())
docs=js(P/'phase145_macrophage_provenance_trace_remote_v1/DOCX_TEXT_EXTRACTIONS.json')
for i,d in enumerate(docs[:3]):
 path=Path(d['text']);sources.extend([path,Path(d['path'])]);text=path.read_text()
 if i<2:
  start=text.index('二、');end=text.index('三、',start)
  (E/('ORIGINAL_PLAN_PART2_'+str(i+1)+'.txt')).write_text(text[start:end],encoding='utf-8')
 else:(E/'LEGACY_REPORT_TEXT.txt').write_text(text,encoding='utf-8')
 assert hashlib.sha256(Path(d['path']).read_bytes()).hexdigest()==d['sha256']
m=read(R/'07_trajectory/input_v1/CELL_METADATA.tsv.gz').set_index('cell_id');m.candidate=m.candidate.astype(str)
c4=read(R/'09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz').cell_id.astype(str)
assert set(c4)==set(m.index[m.candidate=='4'])
legacy=read(P/'phase145_macrophage_provenance_trace_remote_v1/legacy_object/metadata.tsv.gz').set_index('cell_id')
subset=legacy.reindex(c4)
lab=subset.macro_subtype.fillna('not_in_legacy_final_object')
counts=lab.value_counts().rename_axis('legacy_subtype').reset_index(name='new_C4_cells')
counts['fraction_of_C4']=counts.new_C4_cells/len(c4);dump(counts,'NEW_C4_VS_JULY_SEVEN_SUBTYPES.tsv')
study=m.loc[c4].groupby('study').size().rename('C4_cells').reset_index();study['fraction']=study.C4_cells/len(c4);dump(study,'C4_STUDY_COMPOSITION.tsv')
oldsp=legacy.index[legacy.macro_subtype=='SPP1+Macro'];hit=len(set(oldsp)&set(c4))
overlap={'old_SPP1_cells':len(oldsp),'new_C4_cells':len(c4),'intersection':hit,'new_C4_fraction_in_old_SPP1':hit/len(c4),'old_SPP1_fraction_in_new_C4':hit/len(oldsp),'jaccard':hit/len(set(oldsp)|set(c4))}
(E/'JULY_SPP1_OVERLAP.json').write_text(json.dumps(overlap,indent=2))
pq=read(R/'07_trajectory/monocle3_v3/PARTITION_QC.tsv')
cov=pq.groupby('mode',as_index=False).agg(cells=('cells','sum'),C4_cells=('C4_cells','sum'),C4_with_pseudotime=('C4_finite','sum'),all_with_pseudotime=('finite_primary','sum'))
cov['C4_coverage']=cov.C4_with_pseudotime/cov.C4_cells;dump(cov,'TRAJECTORY_COVERAGE.tsv')
dump(pq,'PARTITION_QC.tsv')
pt=read(R/'07_trajectory/dynamics_v2/PSEUDOTIME_SUMMARY.tsv');dump(pt,'NEW_PSEUDOTIME_ALL_SENSITIVITIES.tsv')
dt=read(R/'07_trajectory/dynamics_v2/DYNAMIC_TREND_SUMMARY.tsv');dump(dt,'NEW_DYNAMIC_ALL_SENSITIVITIES.tsv')
sr=R/'08_regulons/statistics_v2_20260915'
rc=read(sr/'AUC_CONTRAST_SUMMARY.tsv');rt=read(sr/'AUC_TRAJECTORY_SUMMARY.tsv')
selected=['JUN','JUNB','CEBPD','PPARG','THRA','TFEB','RXRB','CEBPA']
dump(rc,'NEW_ALL_REGULON_CONTRASTS.tsv');dump(rt[rt.TF.isin(selected)],'SELECTED_TF_ALL_TRAJECTORY_SENSITIVITIES.tsv')
lopo=read(sr/'LEAVE_ONE_PATIENT_CONTRAST_STABILITY.tsv');dump(lopo[lopo.regulon.str.split('(').str[0].isin(selected)],'SELECTED_TF_LOPO.tsv')
patient=read(sr/'PATIENT_AUC_CONTRASTS.tsv');patient=patient[patient.TF.isin(selected)]
study_tf=patient.groupby(['scope','endpoint','study','regulon'],as_index=False).agg(n_patients=('patient','nunique'),median_effect=('effect','median'),positive=('effect',lambda s:int((s>0).sum())),negative=('effect',lambda s:int((s<0).sum())))
dump(study_tf,'SELECTED_TF_STUDY_DIRECTIONS_DESCRIPTIVE.tsv')
oldtf=read(P/'phase144_parts2_3_4_audit_remote_v1/evidence/Y03.tsv');dump(oldtf,'INITIAL_RULE_STATE_TF_RESULTS.tsv')
oldpt=read(P/'phase144_parts2_3_4_audit_remote_v1/evidence/Y02.tsv');dump(oldpt,'INITIAL_RULE_STATE_TRAJECTORY.tsv')
oldreg=read(P/'phase34_part2_y03_myeloid_pyscenic_v1/pyscenic_core_v1/regulon_registry_v1.tsv')
newreg=read(R/'08_regulons/network_v2_20260915/REGULON_REGISTRY.tsv')
network=[]
for _,row in newreg[newreg.TF.isin(selected)].iterrows():
 match=oldreg[oldreg.regulon==row.regulon]
 a=set(row.targets.split(';'));b=set(match.iloc[0].targets.split(';')) if len(match) else set()
 network.append(dict(regulon=row.regulon,old_targets=len(b),new_targets=len(a),shared_targets=len(a&b),target_Jaccard=len(a&b)/len(a|b)))
dump(pd.DataFrame(network),'OLD_NEW_REGULON_TARGET_COMPARISON.tsv')
pairs=[]
for tf in ['JUN','JUNB','CEBPD']:
 for old_ep,new_ep in [('SPP1_minus_Monocyte','C4_minus_Monocyte'),('SPP1_minus_C1QC','C4_minus_other_macrophages')]:
  a=oldtf[(oldtf.tf==tf)&(oldtf.endpoint==old_ep)].iloc[0]
  b=rc[(rc.TF==tf)&(rc.scope=='all_tissues')&(rc.endpoint==new_ep)].iloc[0]
  pairs.append(dict(TF=tf,comparison_type='state_contrast',old_endpoint=old_ep,new_endpoint=new_ep,old_n=a.n_eligible_patients,new_n=b.n_patients,old_effect=a.median_patient_effect,new_effect=b.median_effect,old_P=a.wilcoxon_p_two_sided,new_P=b.P,old_FDR_3TF=a.BH_FDR_within_endpoint3TF,new_FDR_3TF=b.FDR_prespecified_3TF,new_FDR_all=b.FDR_all_regulons,comparability='different_state_and_network; C1QC_vs_other_comparator_not_equivalent'))
 a=oldtf[(oldtf.tf==tf)&(oldtf.endpoint=='Monocle3_continuous')].iloc[0]
 b=rt[(rt.TF==tf)&(rt.scope=='all_tissues')&(rt['mode']=='primary')&(rt.pt_version=='pt_primary')&(rt.adjustment=='rho')].iloc[0]
 pairs.append(dict(TF=tf,comparison_type='trajectory_correlation',old_endpoint='Monocle3_continuous',new_endpoint='primary_partition_monocle3',old_n=a.n_eligible_patients,new_n=b.n_patients,old_effect=a.median_patient_effect,new_effect=b.median_effect,old_P=a.wilcoxon_p_two_sided,new_P=b.P,old_FDR_3TF=a.BH_FDR_within_endpoint3TF,new_FDR_3TF=b.FDR_prespecified_3TF,new_FDR_all=b.FDR_all_regulons,comparability='different_network_cells_graph_and_coverage'))
dump(pd.DataFrame(pairs),'INITIAL_VS_CURRENT_TF_COMPARISON.tsv')
dump(read(sr/'REGULON_TARGET_MARKER_OVERLAP.tsv'),'NEW_REGULON_MARKER_OVERLAP.tsv')
dump(read(R/'09_downstream_crosswalk/remaining_analysis_v2_20260915/C4_LEGACY_OVERLAP_METRICS.tsv'),'NEW_C4_VS_AUGUST_RULE_STATES.tsv')
for fn in ['C4_CONTENT_REVIEW_KEY_RESULTS.tsv','C4_ALL_PROGRAMS.tsv']:
 dump(read(R/'06_functional_programs/run_v1'/fn),fn)
for path in [R/'STATUS.json',R/'provenance/pre_scenic_acceptance_20260915_v2/STATUS.json',R/'provenance/remaining_analysis_validation_v2_20260915/STATUS.json',R/'08_regulons/network_v2_20260915/STATUS.json',R/'08_regulons/network_v2_20260915/REGISTRY_STATUS.json',R/'00_plan/remaining_analysis_v1/AUTHORIZATION_AND_SCOPE_CN.md',P/'phase33_part2_y02_myeloid_monocle3_v1/provenance/02_run_y02_monocle3_v1.R',R/'scripts/18b_run_partitioned_monocle3_header_fixed.R',P/'phase145_macrophage_provenance_trace_remote_v1/MACROPHAGE_SPP1_PROVENANCE_REPORT_CN.md']:
 sources.append(path)
manifest=[]
for p in sorted(set(sources)):
 manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(manifest).to_csv(O/'SOURCE_MANIFEST.tsv',sep='\t',index=False)
(O/'EVIDENCE_COLLECTION_STATUS.json').write_text(json.dumps(dict(status='PASS',source_files=len(manifest),C4_IDs_frozen=True,original_document_hashes_match=True,model_reruns=False,figures_created=False,time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2))
print('COVERAGE');print(cov.to_string(index=False));print('C4_STUDY');print(study.to_string(index=False));print('JULY_OVERLAP',overlap);print(counts.to_string(index=False));print('TARGETS');print(pd.DataFrame(network).to_string(index=False))
print('TF_PRIMARY');print(rt[(rt.TF.isin(selected))&(rt.scope=='all_tissues')&(rt['mode']=='primary')&(rt.pt_version=='pt_primary')].to_string(index=False))
print('LOPO');print(lopo[(lopo.regulon.str.split('(').str[0].isin(selected))&(lopo.scope=='all_tissues')].to_string(index=False))
print('STUDY');print(study_tf[(study_tf.scope=='all_tissues')&(study_tf.endpoint=='C4_minus_other_macrophages')].to_string(index=False))

