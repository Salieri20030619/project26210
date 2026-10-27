from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
H=Path.home();R=H/'project26210_external_joint_v1';O=R/'phase211_immune_receiver_refreeze_20260919_v1';F=O/'freeze_v2';F.mkdir(exist_ok=False)
d=pd.read_csv(O/'RECEIVER_CELL_METADATA_WITH_CLUSTERS.tsv.gz',sep='\t',index_col=0);cache=np.load(O/'RAW_PANEL_CACHE.npz');genes=list(cache['genes']);raw=pd.DataFrame(cache['counts'],index=d.index,columns=genes);assert np.array_equal(cache['cell_id'],d.index)
def save(x,name):x.to_csv(F/name,sep='\t',index=False)
ann={0:('CD4_Tcon','moderate','T backbone + IL7R; low CD8A/B; helper legacy support; depth and dominant patient flagged'),1:('CD8_T','moderate','T backbone with CD8A/B and little CD4; NK-like phenotype does not override T identity'),2:('gdT_enriched','moderate','T backbone and coherent TRDC/TRGC1/TRGC2; keep separate from conventional CD8 and NK'),3:('T_unresolved','unresolved','Mixed CD4/CD8 parent; only marker-coherent intersections 3.0 and 3.8 assigned'),4:('CD8_T','high','T backbone with CD8A/B and little CD4; gamma-gene signal recorded, no functional gate'),5:('CD4_Tcon','high','T backbone CD4 IL7R with low CD8A/B, low FOXP3'),6:('T_unresolved','unresolved','CD4 and CD8 programs mixed without separate coherent boundary'),7:('T_unresolved','unresolved','Low-complexity T population with unresolved CD4/CD8'),8:('ILC_candidate','unresolved','IL7R high with little TCR/CD3; not conventional CD4 or NK'),9:('T_NK_unresolved','unresolved','Small unstable mixed-label population'),10:('CD8_T','high','Strong T backbone CD8A/B and little CD4'),11:('T_unresolved','unresolved','Mixed helper/cytotoxic legacy and insufficient coherent subtype'),12:('NK','high','KLRD1/NK profile with low TCR/CD3 and CD8B'),13:('CD8_T','moderate','T backbone and CD8A/B in cycling/mixed-state group; state not used for lineage exclusion'),14:('CD4_Treg','high','Coherent T backbone CD4 FOXP3 IKZF2 IL2RA CTLA4; keep all members, no single-cell FOXP3 gate'),15:('NK','high','NK receptor program; low TCR/CD3/CD8B'),16:('T_NK_unresolved','unresolved','Low-depth mixed parent remains unresolved; seed fragments alone insufficient'),17:('non_TNK_unresolved','unresolved','No coherent T/NK backbone; legacy pDC-like, no specific new pDC claim'),18:('CD8_T','moderate','CD8A/B with T backbone and little CD4, low-depth flag retained'),19:('T_unresolved','unresolved','Patient-concentrated mixed CD4/CD8/NK-like evidence')}
rows=[dict(unit=str(cl),label=z[0],confidence=z[1],reason=z[2]) for cl,z in ann.items()]
rows += [dict(unit='3.0',label='CD4_Tcon',confidence='moderate',reason='Mixed parent resolved by seed42 intersection: 2518 cells, IL7R76.4%, CD8A/B3.4/3.3%, T backbone coherent'),dict(unit='3.8',label='CD8_T',confidence='moderate',reason='Mixed parent resolved by seed42 intersection: 2629 cells, CD8A/B48.3/50.2%, CD4 2.3%, T backbone coherent')]
mapdf=pd.DataFrame(rows);save(mapdf,'DISCOVERY_UNIT_ANNOTATION.tsv');lookup=mapdf.set_index('unit')
d['annotation_unit']=d.receiver_cluster.astype(str)
for seed in [0,8]:d.loc[(d.receiver_cluster==3)&(d.receiver_cluster_seed42==seed),'annotation_unit']=f'3.{seed}'
d['receiver_label_v1']=d.annotation_unit.map(lookup.label);d['identity_confidence_v1']=d.annotation_unit.map(lookup.confidence);d['identity_reason_v1']=d.annotation_unit.map(lookup.reason)
review=pd.read_csv(O/'CLUSTER_IDENTITY_REVIEW.tsv',sep='\t').set_index('cluster');d['patient_concentrated_cluster']=d.receiver_cluster.map(review.max_patient_fraction)>=0.5;d['cluster_partition_Jaccard']=d.receiver_cluster.map(review.seed42_Jaccard);d['lineage_resolved']=d.receiver_label_v1.isin(['CD4_Tcon','CD4_Treg','CD8_T','NK','gdT_enriched']);d['primary_receiver_type']=d.receiver_label_v1.where(d.receiver_label_v1.isin(['CD4_Tcon','CD4_Treg','CD8_T','NK']),'excluded_unresolved_or_nonconventional')
assert d.receiver_label_v1.notna().all();save(d.reset_index(),'DISCOVERY_RECEIVER_CELL_FREEZE.tsv.gz')
# Raw expression probe counts used for eligibility only, never for selecting subtype.
programs={'CD8_cytotoxicity':['NKG7','GNLY','GZMB','PRF1','CTSW'],'CD8_dysfunction':['PDCD1','CTLA4','LAG3','HAVCR2','TIGIT','TOX','CXCL13','ENTPD1'],'NK_cytotoxicity':['NKG7','GNLY','GZMB','PRF1','CTSW']}
mac=pd.read_csv(R/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str});mac=mac[~mac.candidate.isin(['1','Monocyte'])];c4=pd.read_csv(R/'part2_rebuild_v1/09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz',sep='\t').cell_id
assert len(c4)==2674 and set(c4)==set(mac.index[mac.candidate=='4']) and not set(d.index)&set(mac.index)
keys=['study','curated_patient_id','sample_id','source_tissue'];mcounts=mac.groupby(keys,observed=True).size().to_dict();c4counts=mac[mac.candidate=='4'].groupby(keys,observed=True).size().to_dict();librows=[];ep=[]
for k,g in d.groupby(keys,observed=True):
 ct=g.receiver_label_v1.value_counts();unknown=int((~g.lineage_resolved).sum());frac=unknown/len(g);nm=int(mcounts.get(k,0));nc=int(c4counts.get(k,0));base=dict(study=k[0],patient=k[1],library=k[2],tissue=k[3],TNK_cells=len(g),unresolved_cells=unknown,unresolved_fraction=frac,n_mature=nm,n_C4=nc)
 for label in sorted(d.receiver_label_v1.unique()):base[label]=int(ct.get(label,0))
 base['CD4_denominator']=int(ct.get('CD4_Tcon',0)+ct.get('CD4_Treg',0));librows.append(base)
 for endpoint in list(programs)+['Treg_CD4']:
  receiver='CD8_T' if endpoint.startswith('CD8') else 'NK';selected=g[g.receiver_label_v1==receiver] if endpoint!='Treg_CD4' else g[g.receiver_label_v1.isin(['CD4_Tcon','CD4_Treg'])];nr=len(selected);threshold=50 if endpoint=='Treg_CD4' else 20
  detected=(raw.loc[selected.index,programs[endpoint]]>0).sum() if endpoint in programs else pd.Series(dtype=int);covered=int((detected>=5).sum()) if endpoint in programs else 0
  reasons=[]
  if k[3]!='Tumor':reasons.append('non_tumour')
  if nm<50:reasons.append('mature_macrophage_lt50')
  if nr<threshold:reasons.append('receiver_below_threshold')
  if frac>0.2:reasons.append('unresolved_identity_gt20pct')
  if endpoint in programs and covered<3:reasons.append('program_coverage_lt3genes_5cells')
  ep.append(dict(dataset='discovery75',study=k[0],patient=k[1],library=k[2],tissue=k[3],endpoint=endpoint,n_mature=nm,n_C4=nc,n_receiver=nr,n_Treg=int(ct.get('CD4_Treg',0)) if endpoint=='Treg_CD4' else np.nan,unresolved_fraction=frac,covered_genes=covered,eligible=len(reasons)==0,exclusion_reasons=';'.join(reasons),C4_intrinsic_eligible=(len(reasons)==0 and nc>=20),gene_detection_counts=json.dumps(detected.to_dict())))
lib=pd.DataFrame(librows);save(lib,'DISCOVERY_LIBRARY_COUNTS.tsv');e=pd.DataFrame(ep);save(e,'DISCOVERY_ENDPOINT_LIBRARY_FREEZE.tsv')
patient=[]
for (study,pt,endpoint),g in e.groupby(['study','patient','endpoint']):
 good=g[g.eligible];c4good=g[g.C4_intrinsic_eligible];patient.append(dict(study=study,patient=pt,endpoint=endpoint,eligible=len(good)>0,eligible_libraries=';'.join(good.library),n_libraries=len(good),C4_intrinsic_libraries=';'.join(c4good.library),aggregation='equal_library_weights_within_patient_same_pairing_set',excluded_libraries_with_reasons='|'.join(g.loc[~g.eligible,'library']+':'+g.loc[~g.eligible,'exclusion_reasons'])))
pf=pd.DataFrame(patient);save(pf,'DISCOVERY_PATIENT_ENDPOINT_FREEZE.tsv')
summary=[]
for endpoint,g in e[e.eligible].groupby('endpoint'):summary.append(dict(dataset='discovery75',endpoint=endpoint,libraries=len(g),patients=g.patient.nunique(),studies=g.study.nunique(),C4_intrinsic_patients=g[g.C4_intrinsic_eligible].patient.nunique()))
save(pd.DataFrame(summary),'DISCOVERY_ENDPOINT_SUMMARY.tsv')
# Audit sensitivity register, never substitute primary IDs automatically.
sens=[]
for endpoint,g in e.groupby('endpoint'):
 for cap in [0.1,0.2,0.3]:
  gate=(g.tissue=='Tumor')&(g.n_mature>=50)&(g.n_receiver>=(50 if endpoint=='Treg_CD4' else 20))&(g.unresolved_fraction<=cap)
  if endpoint!='Treg_CD4':gate&=g.covered_genes>=3
  z=g[gate];sens.append(dict(endpoint=endpoint,unresolved_cap=cap,primary=cap==0.2,patients=z.patient.nunique(),libraries=len(z)))
save(pd.DataFrame(sens),'IDENTITY_MISSINGNESS_THRESHOLD_REGISTER.tsv')
# External source clusters: conserve ambiguity rather than perpetuate NK misannotation.
extmaps={'GSE264586':{0:'T_NK_unresolved',1:'CD4_Tcon',4:'T_unresolved',7:'CD4_Tcon',15:'CD8_T',21:'CD8_T',23:'T_NK_unresolved'},'GSE294109':{0:'CD8_T',1:'CD8_T',2:'CD4_Tcon',3:'T_NK_unresolved',5:'CD4_Treg',7:'CD8_T',13:'gdT_NKlike_unresolved',14:'CD8_T',23:'CD8_T'}}
extrows=[];extp=[];extann=[]
for ds,mapping in extmaps.items():
 x=pd.read_csv(O/f'external/{ds}_CELL_MAPPING.tsv.gz',sep='\t',index_col='cell_id');rr=pd.read_csv(O/f'external/{ds}_TNK_PANEL_RAW_COUNTS.tsv.gz',sep='\t',index_col='cell_id')
 t=x.loc[rr.index].copy();t['receiver_label_v1']=t.rds_phase7_cluster_res0_8.map(mapping);assert t.receiver_label_v1.notna().all();t['identity_confidence_v1']=t.rds_phase7_annotation_confidence;t['resolved_receiver']=t.receiver_label_v1.isin(['CD8_T','CD4_Tcon','CD4_Treg','NK']);t['patient_concentrated_cluster']=(ds=='GSE294109')&(t.rds_phase7_cluster_res0_8==7)
 save(t.reset_index(),f'{ds}_RECEIVER_CELL_FREEZE.tsv.gz')
 for cl,label in mapping.items():extann.append(dict(dataset=ds,source_cluster=cl,new_label=label,reason='Raw TCR/CD3 plus gamma/NK signals conflict with old NK; excluded from conventional NK/CD8 primary' if 'NK' in label else ('TOX state does not settle CD4/CD8 identity' if label=='T_unresolved' else 'Existing identity retained after source-cluster and raw-marker check; cluster7 dominance separate from CD8 identity')))
 for k,g in t.groupby(['patient_key','sample_id','timepoint','treatment'],dropna=False):
  allunit=x[(x.patient_key==k[0])&(x.sample_id==k[1])&(x.timepoint==k[2])];nm=int(allunit.phase7_fine_annotation.isin(['Macrophage_C1QC','Macrophage_SPP1']).sum());uf=float((~g.resolved_receiver).mean());tissue_values=';'.join(sorted(set(g.rds_source_tissue_status.astype(str))));tissue_ok=set(g.rds_source_tissue_status.astype(str))=={'Tumor'}
  for endpoint in list(programs)+['Treg_CD4']:
   lab='CD8_T' if endpoint.startswith('CD8') else 'NK';sel=g[g.receiver_label_v1==lab] if endpoint!='Treg_CD4' else g[g.receiver_label_v1.isin(['CD4_Tcon','CD4_Treg'])];nr=len(sel);det=(rr.loc[sel.index,programs[endpoint]]>0).sum() if endpoint in programs else pd.Series(dtype=int);ncovered=int((det>=5).sum()) if endpoint in programs else 0;reasons=[]
   if not tissue_ok:reasons.append('tumour_status_not_confirmed_in_source')
   if nm<50:reasons.append('mature_macrophage_lt50')
   if nr<(50 if endpoint=='Treg_CD4' else 20):reasons.append('receiver_below_threshold')
   if uf>0.2:reasons.append('unresolved_identity_gt20pct')
   if endpoint in programs and ncovered<3:reasons.append('program_coverage_lt3genes_5cells')
   if endpoint=='Treg_CD4' and ds=='GSE264586':reasons.append('no_resolved_Treg_annotation')
   extrows.append(dict(dataset=ds,patient=k[0],sample_id=k[1],timepoint=k[2],treatment=k[3],source_tissue_status=tissue_values,endpoint=endpoint,n_mature=nm,n_receiver=nr,n_Treg=int((sel.receiver_label_v1=='CD4_Treg').sum()) if endpoint=='Treg_CD4' else np.nan,unresolved_fraction=uf,covered_genes=ncovered,eligible=len(reasons)==0,exclusion_reasons=';'.join(reasons),C4_available=False))
save(pd.DataFrame(extann),'EXTERNAL_CLUSTER_ANNOTATION_FREEZE.tsv');ee=pd.DataFrame(extrows);save(ee,'EXTERNAL_ENDPOINT_SAMPLE_FREEZE.tsv')
for k,g in ee.groupby(['dataset','patient','timepoint','treatment','endpoint'],dropna=False):
 z=g[g.eligible];extp.append(dict(dataset=k[0],patient=k[1],timepoint=k[2],treatment=k[3],endpoint=k[4],eligible=len(z)>0,n_samples=len(z),eligible_samples=';'.join(z.sample_id),aggregation='equal_paired_samples_within_patient_within_timepoint_and_treatment_only'))
save(pd.DataFrame(extp),'EXTERNAL_PATIENT_ENDPOINT_FREEZE.tsv')
extsum=[]
for k,g in ee[ee.eligible].groupby(['dataset','timepoint','treatment','endpoint']):extsum.append(dict(dataset=k[0],timepoint=k[1],treatment=k[2],endpoint=k[3],patients=g.patient.nunique(),samples=len(g)))
save(pd.DataFrame(extsum),'EXTERNAL_ENDPOINT_STRATUM_SUMMARY.tsv')
# Source and freeze definitions.
(O/'ANNOTATION_DECISIONS.json').write_text(json.dumps({'discovery':rows,'external':extann,'mixed_parent_policy':'Only primary3 seed42 intersections0/8 are resolved from coherent opposing marker patterns. Other fragments remain unresolved; seed42 is an aid to subdivision, not independent validation.','no_clustering_change_for_C4':True,'external_NK':'No conventional NK freeze from these old mixed/NKlike external clusters; not biological absence.','external_GSE264586':'Formal tumour-only validation not eligible because source tumour status Unknown; organs are sampling sites, not confirmation of tumour. Existing results preserved as historical context.'},indent=2))
(F/'FREEZE_DEFINITION.json').write_text(json.dumps({'version':'phase211_receiver_freeze_v2','created':datetime.datetime.now().astimezone().isoformat(),'cell_membership':'Exact ID tables authoritative; no changing markers or thresholds based on downstream results','unresolved_included_in_audit':True,'unresolved_excluded_from_type_specific_denominators':True,'unknown_fraction_primary_max':0.2,'minimum_mature_macrophages':50,'minimum_CD8_NK':20,'minimum_CD4_denominator':50,'minimum_C4_intrinsic':20,'minimum_program_genes':3,'minimum_detected_cells_per_gene':5,'external_strata_separate':True,'C4_membership_sha256':hashlib.sha256((R/'part2_rebuild_v1/09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz').read_bytes()).hexdigest(),'programs':programs,'association_tests_run':False},indent=2))
(F/'STATUS.json').write_text(json.dumps({'status':'FROZEN_PENDING_INDEPENDENT_VALIDATION','discovery_cells':len(d),'discovery_identity_counts':d.receiver_label_v1.value_counts().to_dict(),'external_cells':{ds:len(pd.read_csv(F/f'{ds}_RECEIVER_CELL_FREEZE.tsv.gz',sep='\t',usecols=['cell_id'])) for ds in extmaps},'C4_cells_unchanged':2674},indent=2))
print('FROZEN',d.receiver_label_v1.value_counts().to_dict());print('DISCOVERY',pd.DataFrame(summary).to_string(index=False));print('EXTERNAL',pd.DataFrame(extsum).to_string(index=False))
