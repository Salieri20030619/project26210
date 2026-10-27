import pathlib,subprocess,json,os,select,hashlib,gzip,zlib,datetime,sys,traceback,shutil
C=pathlib.Path(__file__).resolve().parent
P=C.parents[2]
B=C/'sratoolkit.3.2.1-ubuntu64/bin'
def now():return datetime.datetime.now().astimezone().isoformat()
def save(name,obj):(C/name).write_text(json.dumps(obj,indent=2))
def emit(event,**kw):print(json.dumps(dict(event=event,time=now(),**kw)),flush=True)
def digest(q):
 h=hashlib.md5()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def run(args,name):
 emit('START',step=name)
 with (C/(name+'.log')).open('w') as log:cp=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT)
 (C/(name+'.exit_code')).write_text(str(cp.returncode)+'\n')
 assert cp.returncode==0,name+' failed: '+str(cp.returncode)
 emit('PASS',step=name)
def main():
 pid=int((C/'DOWNLOAD_PID').read_text())
 try:
  fd=os.pidfd_open(pid);emit('WAITING_FOR_ARCHIVE',pid=pid);select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
 qa=json.loads((C/'ARCHIVE_QA.json').read_text());assert qa['status']=='PASS'
 assert (C/'DOWNLOAD_EXIT_STATUS').read_text().strip()=='0'
 emit('ARCHIVE_MD5_VERIFIED',md5=qa['observed_md5'],bytes=qa['observed_bytes'])
 run([str(B/'vdb-validate'),str(C/'SRR32389763.sra')],'vdb_validate')
 out=C/'export';out.mkdir(exist_ok=True);tmp=C/'temp';tmp.mkdir(exist_ok=True)
 assert shutil.disk_usage(C).free>40*10**9
 args=[str(B/'fasterq-dump'),str(C/'SRR32389763.sra'),'--split-files','--threads','4','--mem','1024MB','--temp',str(tmp),'--outdir',str(out),'--seq-defline','@SRR32389763.$si $sn/$ri','--qual-defline','+']
 save('EXPORT_COMMAND.json',args);run(args,'fasterq_export')
 r1=out/'SRR32389763_1.fastq';r2=out/'SRR32389763_2.fastq'
 emit('PAIR_AND_STRUCTURE_CHECK_STARTED')
 n=0;h1=hashlib.md5();h2=hashlib.md5()
 with r1.open('rb') as a,r2.open('rb') as b:
  while True:
   x=[a.readline() for _ in range(4)];y=[b.readline() for _ in range(4)]
   if not x[0] and not y[0]:break
   assert all(x) and all(y),'Truncated or unequal records'
   assert x[0].rstrip().endswith(b'/1') and y[0].rstrip().endswith(b'/2')
   assert x[0].rstrip()[:-2]==y[0].rstrip()[:-2],'Pair identity mismatch'
   assert x[0].startswith(b'@') and x[2].strip()==b'+' and y[2].strip()==b'+'
   assert len(x[1].rstrip(b'\r\n'))==len(x[3].rstrip(b'\r\n')) and len(y[1].rstrip(b'\r\n'))==len(y[3].rstrip(b'\r\n'))
   for z in x:h1.update(z)
   for z in y:h2.update(z)
   n+=1
 emit('PAIR_AND_STRUCTURE_CHECK_PASS',paired_reads=n)
 h=hashlib.md5()
 with gzip.open(P/'fastq/SRR32389763_1.fastq.gz','rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 assert h.hexdigest()==h1.hexdigest(),'Full R1 text differs from ENA'
 emit('FULL_R1_TEXT_MATCH',md5=h.hexdigest())
 dec=zlib.decompressobj(31);prefix=0
 with (P/'fastq/SRR32389763_2.fastq.gz.part').open('rb') as f,r2.open('rb') as g:
  for b in iter(lambda:f.read(1048576),b''):
   raw=dec.decompress(b);assert g.read(len(raw))==raw,'R2 partial text mismatch';prefix+=len(raw)
 assert prefix>0
 emit('R2_PARTIAL_TEXT_MATCH',uncompressed_prefix_bytes=prefix)
 recovered=C/'SRR32389763_2.recovered.fastq.gz'
 emit('COMPRESSING_R2')
 with recovered.open('wb') as outgz:cp=subprocess.run(['gzip','-n','-6','-c',str(r2)],stdout=outgz,stderr=subprocess.PIPE)
 assert cp.returncode==0,cp.stderr.decode(errors='replace')
 run(['gzip','-t',str(recovered)],'recovered_gzip_test')
 observed=digest(recovered)
 expected='d5f9be6d4c699cc12bb697b1f04233f5'
 result=dict(status='PASS',finished_at=now(),archive_official_md5=qa['expected_md5'],archive_md5_match=True,vdb_validate_exit_code=0,export_exit_code=0,paired_reads=n,FASTQ_structure_pass=True,R1_full_text_matches_ENA=True,R1_uncompressed_md5=h1.hexdigest(),R2_uncompressed_md5=h2.hexdigest(),R2_partial_prefix_matches=True,R2_prefix_uncompressed_bytes=prefix,recovered_file=str(recovered),recovered_bytes=recovered.stat().st_size,recovered_md5=observed,expected_ENA_md5=expected,exact_ENA_compressed_bytes_match=observed==expected,gzip_exit_code=0)
 if observed==expected:
  assert recovered.stat().st_size==1095895219
  dest=P/'fastq/SRR32389763_2.fastq.gz';assert not dest.exists();os.link(recovered,dest)
  result['canonical_file']=str(dest)
  save('EXACT_ENA_RESTORED.json',result)
 save('RECOVERY_QA.json',result);emit('RECOVERY_COMPLETE',**result)
code=1
try:main();code=0
except BaseException:
 error=traceback.format_exc();(C/'RECOVERY_ERROR.txt').write_text(error);emit('RECOVERY_FAILED',error=error)
finally:(C/'RECOVERY_EXIT_STATUS').write_text(str(code)+'\n')
sys.exit(code)
