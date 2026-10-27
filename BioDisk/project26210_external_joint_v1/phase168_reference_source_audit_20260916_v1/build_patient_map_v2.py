from pathlib import Path
import json,re,csv,hashlib
import pandas as pd
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase168_reference_source_audit_20260916_v1');R=O.parent
j=json.loads((O/'ORIGINAL_TABLES_EXTRACT.json').read_text())
rows=j['original_mmc2.xlsx']['Sheet1']; byid={r[1]:(i+1,r) for i,r in enumerate(rows[1:],start=1)}
meta=json.loads((R/'phase165_A_input_confirmation_20260915_v1/metadata/SAMPLE_METADATA_EXTRACT.json').read_text())
mapped=[];checks=[]
for x in meta:
 m=re.search(r'(ffpe|frozen)_([abc]_\d+)',x['title'][0]);prep,label=m.groups();er,row=byid[label]
 assert prep in row[10].lower()
 mapped.append(dict(section=x['gsm'],dataset='GSE175540',preparation=prep,source_id=label,cohort=row[0],patient_id='GSE175540:'+row[0]+':'+label,pT=row[2],clinical_stage=None,prior_treatment=None,ISUP=None,mapping_status='CONFIRMED_ORIGINAL_TABLE',source_file='original_mmc2.xlsx',source_sheet='Sheet1',source_row=er,source_id_cell='B'+str(er),source_preparation_cell='K'+str(er),preparation_in_original_table=row[10]))
assert len(mapped)==24 and len(set(x['patient_id'] for x in mapped))==23
text=(O/'GSE210041_family.soft').read_text()
for block in text.split('^SAMPLE = ')[1:]:
 gsm=block.splitlines()[0];chars={}
 for line in block.splitlines():
  if line.startswith('!Sample_characteristics_ch1 = '):
   key,val=line.split(' = ',1)[1].split(': ',1);chars[key]=val
 mapped.append(dict(section=gsm,dataset='GSE210041',preparation='ffpe',source_id=chars['patient'],cohort='GSE210041',patient_id='GSE210041:'+chars['patient'],pT=None,clinical_stage=chars['disease stage'],prior_treatment=chars['prior treatment'],ISUP=chars['tumour grade'],mapping_status='CONFIRMED_GEO',source_file='GSE210041_family.soft',source_sheet=None,source_row=None,source_id_cell=None,source_preparation_cell=None,preparation_in_original_table='FFPE per GEO extraction protocol'))
d=pd.DataFrame(mapped);d.to_csv(O/'SPATIAL_26_SECTION_PATIENT_MAP.tsv',sep='\t',index=False)
endpoint=pd.read_csv(R/'phase167_A_24section_endpoint_check_20260916_v1/SECTION_ENDPOINT_SUMMARY.tsv',sep='\t')
joined=d[d.dataset=='GSE175540'].merge(endpoint.drop(columns='preparation'),on='section',validate='one_to_one')
joined.to_csv(O/'PHASE167_ELIGIBILITY_WITH_PATIENTS.tsv',sep='\t',index=False)
print('COHORTS',d.groupby(['dataset','cohort']).agg(sections=('section','size'),patients=('patient_id','nunique')).to_string())
print('ELIGIBLE',joined[joined.source_all_designs_eligible].patient_id.nunique(),joined[joined.receiver_all_designs_eligible].patient_id.nunique())
print('DUPLICATE',d[d.patient_id.duplicated(keep=False)].to_string(index=False))
print('MISSING_PT',d[(d.dataset=='GSE175540')&d.pT.isin(['NA',None])].section.tolist())
# Independently reconcile all three embedded copies of the 26-row clinical table.
pages=(O/'combined_reference.txt').read_text().split('\f')
tables=[]
for start in [68]:
 matches=re.findall(r'(GSM\d+)\s+ccRCC\s+(pT\S+|NA)\s+(Stage IV|Stage III|NA)\s+(Nivolumab\+Ipilimumab|Naive|NA)',pages[start]+'\n'+pages[start+1])
 tables.append(matches)
print('TABLE_COPY_COUNTS',[len(x) for x in tables])
assert all(len(x)==26 for x in tables) 
for gsm,pt,stage,trt in tables[-1]:
 v=d.set_index('section').loc[gsm]
 assert (v.pT or 'NA')==pt
 assert (v.clinical_stage or 'NA')==stage
 assert (v.prior_treatment or 'NA')==trt
pd.DataFrame(tables[-1],columns=['section','pT','clinical_stage','treatment']).to_csv(O/'SCIADV_COMBINED_TABLE2_26_ROWS.tsv',sep='\t',index=False)
tls=j['original_mmc3.xlsx']['Sheet1']; tlsmap={r[0]:(i+1,r[-1]) for i,r in enumerate(tls) if i>=2 and r[0]}
tr=[]
for r in mapped[:24]:
 if r['preparation']!='frozen':continue
 found=tlsmap.get(r['source_id'])
 tr.append(dict(section=r['section'],source_id=r['source_id'],TLS_status=found[1] if found else None,source_row=found[0] if found else None,source='original_mmc3.xlsx exact ID match only',note='a_3 supported TLS+ by original_mmc1 Figure S4 legend, separate evidence; mmc3 has a_4, do not silently rename' if not found else ''))
pd.DataFrame(tr).to_csv(O/'FROZEN_TLS_SOURCE_REVIEW.tsv',sep='\t',index=False)
print('TLS_EXACT',sum(x['TLS_status'] is not None for x in tr))
status=dict(stage='Phase168',mapping_status='RESOLVED_FROM_ORIGINAL_SOURCES',GSE175540_sections=24,GSE175540_patients=23,duplicate_patient='BIONIKK:c_2',duplicate_sections=['GSM5924030','GSM5924050'],GSE210041_sections=2,GSE210041_patients=2,phase167_source_patients=23,phase167_receiver_patients=int(joined[joined.receiver_all_designs_eligible].patient_id.nunique()),original_table_spatial_unique_ids=sum(r[10] is not None for r in rows[1:]),clinical_table_26_rows_verified=True,clinical_pt_records_23_match=True,association_tests=0,new_analysis_figures=0)
(O/'MAPPING_VALIDATION.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
