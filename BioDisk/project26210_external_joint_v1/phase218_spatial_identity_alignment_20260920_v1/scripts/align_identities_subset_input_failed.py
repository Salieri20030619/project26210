from pathlib import Path
import csv,gzip,json,collections,hashlib,datetime
P=Path(__file__).resolve().parents[1];R=P.parent;B=R.parent
sources={'old':B/'project26210_part4_tymp_reset_phase17_20260812_v1/reference_tymp_continuous_v1/reference_tymp_continuous_v1_cell_assignments.tsv.gz','c4':R/'part2_rebuild_v1/09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz','macro':R/'phase154_same_library_TYMP_C4_test_20260915_v1/SELECTED_CELL_TWO_GENE_COUNTS.tsv.gz','receiver':R/'phase211_immune_receiver_refreeze_20260919_v1/freeze_v2/DISCOVERY_RECEIVER_CELL_FREEZE.tsv.gz','genes':B/'project26210_external_joint_v1_phase54_gse175540_spatial_v1/scripts_v1/frozen_signature_genes_v1.tsv'}
def read(f):
 with (gzip.open(f,'rt') if f.suffix=='.gz' else f.open()) as h:return list(csv.DictReader(h,delimiter='\t'))
def write(n,rows,fields):
 with (gzip.open(P/n,'wt') if n.endswith('.gz') else (P/n).open('w')) as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
oldrows=read(sources['old']);old={x['exact_cell_id']:x for x in oldrows};assert len(old)==len(oldrows)
c4=read(sources['c4']);assert len(c4)==2674;ids={x['cell_id'] for x in c4}
mac=[x for x in read(sources['macro']) if x['lineage']=='Macrophage'];rec=read(sources['receiver'])
joined=[]
for kind,rows in [('mature_macrophage',mac),('receiver',rec)]:
 for x in rows:
  cid=x['cell_id'];o=old.get(cid,{})
  new=('C4' if cid in ids else 'Macrophage_C'+x['candidate']) if kind=='mature_macrophage' else x['primary_receiver_type']
  joined.append({'cell_id':cid,'lineage':kind,'new_identity':new,'patient':x['curated_patient_id'],'study':x['study'],'sample':x['sample_id'],'tissue':x['source_tissue'],'old_state':o.get('state_tymp_continuous_v1','NOT_IN_OLD_REFERENCE'),'old_patient':o.get('patient',''),'old_study':o.get('study_v1',''),'in_old_reference':bool(o),'is_C4':cid in ids})
assert len({x['cell_id'] for x in joined})==len(joined)
assert sum(x['is_C4'] for x in joined)==2674
write('CELL_IDENTITY_CROSSWALK.tsv.gz',joined,list(joined[0]))
counts=collections.Counter((x['lineage'],x['new_identity'],x['old_state']) for x in joined)
write('IDENTITY_CROSSTAB.tsv',[dict(lineage=k[0],new_identity=k[1],old_state=k[2],cells=v) for k,v in sorted(counts.items())],['lineage','new_identity','old_state','cells'])
ct=collections.Counter(x['old_state'] for x in joined if x['is_C4']);oldSP={k for k,v in old.items() if v['state_tymp_continuous_v1']=='SPP1-state Macrophage'}
inter=ids&oldSP
coverage=[]
for new in sorted({x['new_identity'] for x in joined}):
 rows=[x for x in joined if x['new_identity']==new]
 for scope in ['all','tumor']:
  rr=[x for x in rows if scope=='all' or x['tissue']=='Tumor']
  c=collections.Counter(x['study'] for x in rr);pc=collections.Counter(x['patient'] for x in rr)
  coverage.append({'identity':new,'scope':scope,'cells':len(rr),'patients_present':len(pc),'studies':len(c),'largest_patient_fraction':max(pc.values())/len(rr) if rr else None,'matched_old':sum(x['in_old_reference'] for x in rr),'study_cells':json.dumps(c)})
write('REFERENCE_IDENTITY_COVERAGE.tsv',coverage,list(coverage[0]))
groups=collections.defaultdict(list)
for x in read(sources['genes']):groups[x['signature_id']].append(x['gene'])
pairs=[]
for a in groups:
 if 'SPP1' in a.upper():
  for b in groups:
   if 'vascular' in b.lower() or 'dysfunction' in b.lower():
    aa,bb=set(groups[a]),set(groups[b]);pairs.append(dict(left=a,right=b,left_n=len(aa),right_n=len(bb),shared=sorted(aa&bb),jaccard=len(aa&bb)/len(aa|bb)))
(P/'PROGRAM_DEFINITIONS.json').write_text(json.dumps(groups,indent=2));(P/'PROGRAM_GENE_OVERLAPS.json').write_text(json.dumps(pairs,indent=2))
summary={'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'old_reference_cells':len(old),'old_state_counts':dict(collections.Counter(x['state_tymp_continuous_v1'] for x in oldrows)),'mature_macro_cells':len(mac),'receiver_cells':len(rec),'C4_cells':2674,'C4_old_state_counts':dict(ct),'old_SPP1_state_cells':len(oldSP),'C4_oldSPP1_overlap':len(inter),'fraction_C4_in_oldSPP1':len(inter)/len(ids),'fraction_oldSPP1_in_C4':len(inter)/len(oldSP),'C4_oldSPP1_jaccard':len(inter)/len(ids|oldSP),'patient_label_disagreements':sum(x['in_old_reference'] and x['patient']!=x['old_patient'] for x in joined),'study_label_disagreements':sum(x['in_old_reference'] and x['study']!=x['old_study'] for x in joined),'scope':'EXACT_CELL_ID_ALIGNMENT_AND_REFERENCE_COVERAGE; NO_MAPPING_OR_NEW_SPATIAL_ASSOCIATIONS'}
(P/'ALIGNMENT_SUMMARY.json').write_text(json.dumps(summary,indent=2));(P/'SOURCE_SHA256.json').write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources.values()},indent=2))
print(json.dumps(summary));print('OVERLAPS',json.dumps(pairs));print('PROGRAMS',json.dumps({k:v for k,v in groups.items() if 'spp1' in k.lower() or 'vascular' in k.lower()}))
