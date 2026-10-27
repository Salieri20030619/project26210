from pathlib import Path
import json,hashlib,datetime,csv
import pandas as pd
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');P=R/'phase281_macrophage_version_audit_20260929_v1';P.mkdir(exist_ok=False)
for d in ['code','delivery','qc']:(P/d).mkdir()
paths={
 'original_recluster':R/'part2_rebuild_v1/03_macrophage_clustering/run_v4_exploratory/CELL_CLUSTER_ASSIGNMENTS.tsv.gz',
 'new_regulatory_metadata':R/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz',
 'phase200_mature_macrophages':R/'phase200_union_C4_C10_C15_reassessment_20260918_v1/inputs/MACRO_METADATA.tsv.gz',
 'phase179_STAT3_cache':R/'phase179_STAT3_external_activity_validation_20260916_v1/run_v2/SC_PRIMARY_SCORE_CACHE.tsv.gz',
 'phase280_current_cache':R/'phase280_TYMP_source_niche_function_20260928_v1/inputs/CELL_CACHE.tsv.gz'}
frames={}
for k,f in paths.items():
 frames[k]=pd.read_csv(f,sep='\t',dtype={'candidate':str,'cluster':str})
a=frames['original_recluster'].set_index('cell_id');b=frames['new_regulatory_metadata'].set_index('cell_id');c=frames['phase200_mature_macrophages'].set_index('cell_id')
tum=c[c.source_tissue.str.lower().eq('tumor')]
rows=[]
for name,z in [('recluster_C4',a[a.cluster=='4']),('regulatory_C4',b[b.candidate=='4']),('phase200_C4',c[c.candidate=='4'])]:
 ids=set(z.index);ref=set(a.index[a.cluster=='4'])
 rows.append(dict(check=name,n_cells=len(ids),reference_cells=len(ref),missing_from_reference=len(ref-ids),extra_vs_reference=len(ids-ref),passed=ids==ref))
for name,z in [('phase200_all_mature',c),('phase200_tumor_mature',tum)]:
 common=z.index.intersection(a.index)
 mismatch=int((z.loc[common,'candidate'].astype(str)!=a.loc[common,'cluster'].astype(str)).sum())
 rows.append(dict(check=name,n_cells=len(z),reference_cells=len(common),missing_from_reference=0,extra_vs_reference=len(z)-len(common),label_mismatches=mismatch,passed=(len(common)==len(z) and mismatch==0)))
for k in ['phase179_STAT3_cache','phase280_current_cache']:
 d=frames[k].drop_duplicates(['cell_id','candidate'])
 assert not d.cell_id.duplicated().any()
 ids=set(d.cell_id);ref=set(tum.index);common=d.cell_id.isin(tum.index)
 mismatches=sum(str(v)!=str(tum.at[i,'candidate']) for i,v in zip(d.loc[common,'cell_id'],d.loc[common,'candidate']))
 rows.append(dict(check=k,n_cells=len(ids),reference_cells=len(ref),missing_from_reference=len(ref-ids),extra_vs_reference=len(ids-ref),label_mismatches=mismatches,passed=not(ids-ref) and mismatches==0,interpretation='subset permitted by original patient eligibility'))
# Old STAT3 source eligibility actually uses the same mature definition.
expected=b[(~b.candidate.isin(['1','Monocyte']))&b.source_tissue.str.lower().eq('tumor')]
assert set(expected.index)==set(tum.index)
rows.append(dict(check='phase179_source_filter_equals_current_tumor_mature',n_cells=len(expected),reference_cells=len(tum),passed=True))
# AUC cell coverage is checked without reading 122 MB of scores into memory.
auc=R/'part2_rebuild_v1/08_regulons/network_v2_20260915/AUC_ALL_CELLS.csv'
ac=pd.read_csv(auc,usecols=[0]).iloc[:,0].astype(str)
rows.append(dict(check='new_regulatory_AUC_cell_ids',n_cells=len(ac),reference_cells=len(b),missing_from_reference=len(set(b.index)-set(ac)),extra_vs_reference=len(set(ac)-set(b.index)),passed=set(ac)==set(b.index) and ac.is_unique))
assert all(r['passed'] for r in rows)
pd.DataFrame(rows).to_csv(P/'delivery/CELL_ID_VERSION_AUDIT.tsv',sep='\t',index=False)
counts=[]
for source,d,label in [('original_recluster',a,'cluster'),('new_regulatory_input',b,'candidate'),('current_mature',c,'candidate'),('current_tumor_mature',tum,'candidate')]:
 for state,n in d[label].value_counts().items():counts.append(dict(source=source,state=state,cells=int(n)))
pd.DataFrame(counts).to_csv(P/'delivery/STATE_CELL_COUNTS.tsv',sep='\t',index=False)
manifest=[dict(source=k,path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for k,f in paths.items()]
(P/'qc/INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
validation={'status':'PASS','at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':rows,'current_C4_cells':2674,'current_mature_cells':22444,'regulatory_input_cells':33098,'new_statistics':False}
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(validation,indent=2))
print(json.dumps(validation,ensure_ascii=False))
