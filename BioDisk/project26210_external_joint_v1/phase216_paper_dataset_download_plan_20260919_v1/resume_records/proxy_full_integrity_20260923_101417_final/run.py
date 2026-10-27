import pathlib,json,csv,os,select,datetime,sys,traceback,hashlib,shutil
C=pathlib.Path(__file__).resolve().parent
P=C.parents[1];ROOT=P.parent
V=C.with_name('proxy_full_integrity_20260923_101417_verification')
R=C.with_name('proxy_full_integrity_20260923_101417')/'ncbi_recovery'
def now():return datetime.datetime.now().astimezone().isoformat()
def save(path,obj):
 tmp=path.with_name(path.name+'.tmp');tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False));tmp.replace(path)
def emit(event,**kw):print(json.dumps(dict(event=event,time=now(),**kw)),flush=True)
def wait_pid(file):
 pid=int(file.read_text())
 try:fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
code=1
try:
 emit('WAITING_FOR_VERIFICATION_AND_RECOVERY')
 wait_pid(V/'PID');wait_pid(R/'RECOVERY_PID')
 assert (V/'VERIFICATION_FINISHED').exists(),'Full verification did not finish'
 v=json.loads((V/'QA.json').read_text());r=json.loads((R/'RECOVERY_QA.json').read_text())
 assert v['verified_files']==119 and v['failed_files']==0
 assert (R/'RECOVERY_EXIT_STATUS').read_text().strip()=='0' and r['status']=='PASS'
 assert r['archive_md5_match'] and r['vdb_validate_exit_code']==0 and r['FASTQ_structure_pass'] and r['R1_full_text_matches_ENA'] and r['R2_partial_prefix_matches'] and r['gzip_exit_code']==0
 manifest=list(csv.DictReader((P/'metadata/PUBLIC_FASTQ_MANIFEST.tsv').open(),delimiter='\t'))
 rows=[];integrity=[]
 for sample in manifest:
  for u,n,h in zip(sample['fastq_urls'].split(';'),sample['fastq_bytes'].split(';'),sample['fastq_md5'].split(';')):
   name=u.rsplit('/',1)[1];q=P/'fastq'/name
   if name=='SRR32389763_2.fastq.gz':
    q=pathlib.Path(r.get('canonical_file',r['recovered_file']));size=r['recovered_bytes'];observed=r['recovered_md5'];exact=r['exact_ENA_compressed_bytes_match'];evidence=str(R/'RECOVERY_QA.json');source='NCBI_SRA_exact_ENA_reconstruction' if exact else 'NCBI_SRA_validated_FASTQ_reconstruction'
   else:
    x=json.loads((V/'files'/(name+'.json')).read_text());assert x['status']=='PASS' and x['gzip_exit_code']==0 and x['md5']==x['observed_md5'] and x['unchanged_during_verification'];size=x['observed_bytes'];observed=x['observed_md5'];exact=True;evidence=str(V/'files'/(name+'.json'));source='ENA_original_FASTQ'
   assert q.is_file() and q.stat().st_size==size
   row=dict(run=sample['run_accession'],sample_alias=sample['sample_alias'],group=sample['group'],filename=name,resolved_path=str(q),source=source,observed_bytes=size,observed_md5=observed,ENA_expected_bytes=int(n),ENA_expected_md5=h,exact_ENA_compressed_match=exact,integrity_status='PASS',gzip_integrity='PASS',evidence=evidence)
   rows.append(row)
 assert len(rows)==120 and len({x['filename'] for x in rows})==120
 runs={}
 for x in rows:runs.setdefault(x['run'],[]).append(x)
 assert len(runs)==60 and all(len(x)==2 for x in runs.values())
 with (C/'RESOLVED_FASTQ_MANIFEST.tsv').open('w') as f:
  w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 groups={}
 for x in rows:
  g=groups.setdefault(x['group'],dict(files=0,bytes=0));g['files']+=1;g['bytes']+=x['observed_bytes']
 exact=sum(x['exact_ENA_compressed_match'] for x in rows)
 qa=dict(status='COMPLETE_VERIFIED' if exact==120 else 'COMPLETE_VERIFIED_WITH_NCBI_RECOVERY',finished_at=now(),expected_files=120,verified_files=120,complete_runs=60,groups=groups,verified_bytes=sum(x['observed_bytes'] for x in rows),exact_ENA_size_MD5_files=exact,NCBI_archive_recovered_files=1,all_resolved_files_gzip_integrity_pass=True,failed_files=0,missing_resolved_files=0,original_ENA_endpoint_unavailable='SRR32389763_2.fastq.gz',recovery_evidence=str(R/'RECOVERY_QA.json'),full_integrity_evidence=str(V/'QA.json'),resolved_manifest=str(C/'RESOLVED_FASTQ_MANIFEST.tsv'),source_identity_and_pairing_check=True,treatment_pairing_confirmed=False,biological_analysis_started=False)
 save(C/'QA.json',qa)
 (C/'backups').mkdir()
 for name in ['STATUS.json','DOWNLOAD_SUMMARY.json','EXIT_STATUS']:
  if (P/name).exists():shutil.copy2(P/name,C/'backups'/name)
 save(P/'STATUS.json',qa);save(P/'DOWNLOAD_SUMMARY.json',qa)
 save(P/'FINAL_FULL_INTEGRITY_QA.json',qa)
 shutil.copy2(C/'RESOLVED_FASTQ_MANIFEST.tsv',P/'RESOLVED_FASTQ_MANIFEST.tsv')
 if exact==120:
  full=[]
  for x in rows:
   o=dict(run=x['run'],sample_alias=x['sample_alias'],group=x['group'],filename=x['filename'],bytes=x['observed_bytes'],md5=x['observed_md5'],observed_bytes=x['observed_bytes'],observed_md5=x['observed_md5'],status='VERIFIED',completed_at=now(),full_integrity_evidence=x['evidence'],gzip_exit_code=0)
   save(P/'verification'/(x['filename']+'.json'),o);full.append(o)
  fields=['run','sample_alias','group','filename','bytes','md5','observed_bytes','observed_md5','status','completed_at']
  shutil.copy2(P/'VERIFIED_FILES.tsv',C/'backups/VERIFIED_FILES.tsv')
  with (P/'VERIFIED_FILES.tsv').open('w') as f:
   w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(full)
  (P/'COMPLETE_VERIFIED').write_text(now()+'\n'+str(C/'QA.json')+'\n')
 else:
  (P/'COMPLETE_WITH_NCBI_RECOVERY').write_text(now()+'\n'+str(C/'QA.json')+'\nUse RESOLVED_FASTQ_MANIFEST.tsv; original ENA compressed MD5 not claimed for reconstructed R2.\n')
 (P/'EXIT_STATUS').write_text('0\n');(C/'COMPLETE_VERIFIED').write_text(now()+'\n')
 note='\n\n## '+now()+' Phase216 download and full integrity complete\nAll 120 resolved FASTQs / 60 paired runs are complete and verified. '+str(exact)+'/120 match official ENA compressed size/MD5. The unavailable SRR32389763 R2 was recovered from the same NCBI SRA archive (official archive MD5 verified; vdb-validate passed; complete exported R1 matches ENA text; R2 matches all available partial text; paired FASTQ structure and gzip integrity passed). Exact ENA compressed match for recovered R2: '+str(r['exact_ENA_compressed_bytes_match'])+'.\nUse '+str(C/'RESOLVED_FASTQ_MANIFEST.tsv')+' for actual input paths and verification provenance. Full evidence: '+str(C/'QA.json')+'. Original data, partials and prior failed attempts preserved. This is download/integrity completion only; no biological analysis; treatment identity remains unconfirmed.\n'
 for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
  q=ROOT/name;q.write_text(note+q.read_text())
 with (ROOT/'RESULTS_MANIFEST.tsv').open('a') as f:f.write('\n'+'\t'.join([C.name,'phase216_download_integrity','json',str(C/'QA.json'),now(),'PASS','120 files; 60 paired runs; all integrity checks passed','NCBI recovery provenance; use resolved manifest'])+'\n')
 save(C/'ARTIFACT_SHA256.json',{name:hashlib.sha256((C/name).read_bytes()).hexdigest() for name in ['QA.json','RESOLVED_FASTQ_MANIFEST.tsv','run.py']})
 emit('COMPLETE',**qa);code=0
except BaseException:
 err=traceback.format_exc();(C/'ERROR.txt').write_text(err);emit('FINALIZATION_FAILED',error=err)
finally:(C/'EXIT_STATUS').write_text(str(code)+'\n')
sys.exit(code)
