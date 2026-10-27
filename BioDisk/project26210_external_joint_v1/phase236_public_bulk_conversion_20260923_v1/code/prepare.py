import pathlib,subprocess,concurrent.futures,json,gzip,hashlib,zipfile,os,sys,traceback,datetime,shutil
P=pathlib.Path(__file__).resolve().parents[1];S=P/'software';R=P/'reference'
def run(cmd,name):
 with (P/'logs'/(name+'.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 (P/'logs'/(name+'.exit_code')).write_text(str(r.returncode)+'\n');assert r.returncode==0,name
 print(name+' COMPLETE',flush=True)
def sha(q):
 h=hashlib.sha256()
 with q.open('rb') as f:
  for x in iter(lambda:f.read(8388608),b''):h.update(x)
 return h.hexdigest()
def fetch(url,q):
 run(['curl','-q','--fail','--location','--silent','--show-error','--connect-timeout','60','--retry','10','--retry-delay','15','--retry-all-errors','--continue-at','-','--proxy','http://127.0.0.1:7897','--noproxy','','--output',str(q),url],'download_'+q.name)
 return dict(url=url,path=str(q),bytes=q.stat().st_size,sha256=sha(q))
def envsetup():
 base='/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python'
 run([base,'-m','venv',str(S/'venv')],'venv')
 run([str(S/'venv/bin/pip'),'install','--proxy','http://127.0.0.1:7897','HTSeq==2.0.5','RSeQC==5.0.4','pandas==2.2.3'],'python_packages')
 run([str(S/'venv/bin/pip'),'freeze'],'python_packages_frozen')
code=1
try:
 urls=[('https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download',S/'hisat2-2.1.0.zip'),('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz',R/'GRCh38.primary_assembly.genome.fa.gz'),('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz',R/'gencode.v44.primary_assembly.annotation.gtf.gz'),('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/MD5SUMS',R/'MD5SUMS')]
 with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
  futures=[pool.submit(fetch,*x) for x in urls];env=pool.submit(envsetup)
  result=[f.result() for f in futures];env.result()
 (P/'reference/SOURCE_DOWNLOADS.json').write_text(json.dumps(result,indent=2))
 hashes={x.split()[-1].lstrip('./'):x.split()[0] for x in (R/'MD5SUMS').read_text().splitlines() if len(x.split())==2}
 for name in ['GRCh38.primary_assembly.genome.fa.gz','gencode.v44.primary_assembly.annotation.gtf.gz']:
  q=R/name;assert name in hashes,name+' missing from official MD5SUMS'
  h=hashlib.md5()
  with q.open('rb') as f:
   for x in iter(lambda:f.read(8388608),b''):h.update(x)
  assert h.hexdigest()==hashes[name],'Reference MD5 mismatch'
  with gzip.open(q,'rb') as src,q.with_suffix('').open('wb') as dest:shutil.copyfileobj(src,dest,8388608)
 with zipfile.ZipFile(S/'hisat2-2.1.0.zip') as z:
  assert z.testzip() is None
  for info in z.infolist():
   q=(S/info.filename).resolve();assert S.resolve() in q.parents or q==S.resolve()
   z.extract(info,S)
   perm=(info.external_attr>>16)&0o777
   if perm:q.chmod(perm)
 hisat=S/'hisat2-2.1.0'
 run([str(hisat/'hisat2'),'--version'],'hisat2_version')
 run(['/usr/bin/samtools','faidx',str(R/'GRCh38.primary_assembly.genome.fa')],'fasta_index')
 gtf=R/'gencode.v44.primary_assembly.annotation.gtf'
 with (R/'splice_sites.tsv').open('w') as f:
  r=subprocess.run([str(S/'venv/bin/python'),str(hisat/'hisat2_extract_splice_sites.py'),str(gtf)],stdout=f,stderr=subprocess.PIPE)
 assert r.returncode==0,r.stderr.decode(errors='replace')
 ref=dict(species='Homo sapiens',assembly='GRCh38',annotation='GENCODE v44 primary assembly',fasta=str(R/'GRCh38.primary_assembly.genome.fa'),gtf=str(gtf),fasta_sha256=sha(R/'GRCh38.primary_assembly.genome.fa'),gtf_sha256=sha(gtf),splice_sites_sha256=sha(R/'splice_sites.tsv'),index_prefix=str(R/'hisat2_GRCh38'),official_md5_verified=True,HISAT2_version='2.1.0',HTSeq_version='2.0.5',paper_HTSeq_version='0.11.2',version_difference_documented=True)
 (R/'REFERENCE_FROZEN.json').write_text(json.dumps(ref,indent=2));(P/'SOFTWARE_READY').write_text(datetime.datetime.now().isoformat()+'\n')
 run([str(hisat/'hisat2-build'),'-p','8',ref['fasta'],ref['index_prefix']],'hisat2_build')
 ix=list(R.glob('hisat2_GRCh38.*.ht2*'));assert len(ix)==8
 ref['index_sha256']={q.name:sha(q) for q in ix};(R/'REFERENCE_FROZEN.json').write_text(json.dumps(ref,indent=2));(P/'REFERENCE_READY').write_text(datetime.datetime.now().isoformat()+'\n');code=0
except BaseException:(P/'logs/prepare_ERROR.txt').write_text(traceback.format_exc());print(traceback.format_exc(),flush=True)
finally:(P/'PREPARE_EXIT_STATUS').write_text(str(code)+'\n')
sys.exit(code)
