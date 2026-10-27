from common import *
import concurrent.futures,select,traceback,shutil,re,fcntl,csv
samples=json.loads((P/'SAMPLES.json').read_text());states={s['run']:'PENDING' for s in samples};outcomes={}

def state(run_id,value):
 with LOCK:
  states[run_id]=value
  save(P/'STATUS.json',dict(status='RUNNING',updated_at=now(),samples=60,states=states,completed=sum(x=='COMPLETE' for x in states.values()),failed=sum(x=='FAILED' for x in states.values()),config=str(P/'CONFIG.json')))
 emit('SAMPLE_STAGE',run=run_id,stage=value)
def diskcheck():
 if shutil.disk_usage(P).free<100*10**9:raise RuntimeError('Disk reserve below predeclared 100 GB; original data retained')
def qc(s):
 d=P/'samples'/s['run'];d.mkdir(exist_ok=True);state(s['run'],'FASTP_RUNNING');diskcheck()
 args=['/usr/bin/fastp','-i',s['R1'],'-I',s['R2'],'-o',str(d/'clean_R1.fastq.gz'),'-O',str(d/'clean_R2.fastq.gz'),'--thread','4','--average_qual','20','--n_base_limit','5','--length_required','76','--unqualified_percent_limit','100','--detect_adapter_for_pe','--json',str(d/'fastp.json'),'--html',str(d/'fastp.html')]
 save(d/'fastp_command.json',args);run(args,d/'fastp.log')
 j=json.loads((d/'fastp.json').read_text());before=j['summary']['before_filtering']['total_reads'];after=j['summary']['after_filtering']['total_reads'];assert before%2==0 and after%2==0 and 0<after<=before
 save(d/'QC.json',dict(raw_pairs=before//2,clean_pairs=after//2,retained_fraction=after/before,fastp_version='0.20.1',status='PASS',finished_at=now()))
 state(s['run'],'QC_COMPLETE_WAITING_REFERENCE');return s

def prepare_wait():
 pid=int((P/'PREPARE_PID').read_text())
 try:fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
 assert (P/'PREPARE_EXIT_STATUS').read_text().strip()=='0' and (P/'REFERENCE_READY').exists(),'Reference preparation failed; see logs/prepare_ERROR.txt'
 py=P/'software/venv/bin/python'
 run([str(P/'software/venv/bin/pip'),'install','--proxy','http://127.0.0.1:7897','numpy==1.26.4'],P/'logs/numpy_compatibility.log')
 run([str(py),'-c','import HTSeq,pysam,pandas,numpy; print(HTSeq.__version__,pysam.__version__,pandas.__version__,numpy.__version__)'],P/'logs/import_validation.log')
 run([str(P/'software/venv/bin/pip'),'freeze'],P/'logs/final_software_freeze.log')
 run([str(py),str(P/'code/annotation.py')],P/'logs/annotation.log')
 save(P/'SOFTWARE_VERSIONS.json',dict(fastp='0.20.1',HISAT2='2.1.0',HTSeq='2.0.5',RSeQC='5.0.4',samtools='1.13',numpy='1.26.4',pandas='2.2.3',freeze=str(P/'logs/final_software_freeze.log'),note='HTSeq differs from paper v0.11.2; current union counting choices fixed in CONFIG.json'))
 run([str(py),str(P/'code/counting_smoke.py')],P/'logs/counting_smoke.log')
 (P/'CONVERSION_ENV_READY').write_text(now()+'\n');return json.loads((P/'reference/REFERENCE_FROZEN.json').read_text())

def align_count(s,ready):
 ref=ready.result();d=P/'samples'/s['run'];diskcheck();state(s['run'],'HISAT2_RUNNING')
 hisat=P/'software/hisat2-2.1.0/hisat2';bam=d/'aligned.bam';part=d/'aligned.bam.part'
 cmd=[str(hisat),'-p','8','-x',ref['index_prefix'],'--known-splicesite-infile',str(P/'reference/splice_sites.tsv'),'-1',str(d/'clean_R1.fastq.gz'),'-2',str(d/'clean_R2.fastq.gz'),'--rg-id',s['run'],'--rg','SM:'+s['run'],'--summary-file',str(d/'hisat2_summary.txt'),'--new-summary']
 save(d/'hisat2_command.json',cmd)
 with (d/'hisat2.log').open('w') as log,(d/'samtools_sort.log').open('w') as sortlog:
  align=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log)
  sort=subprocess.Popen(['/usr/bin/samtools','sort','-@','4','-m','1G','-O','BAM','-T',str(d/'sort_tmp'),'-o',str(part),'-'],stdin=align.stdout,stderr=sortlog);align.stdout.close();srt=sort.wait();aln=align.wait()
 save(d/'alignment_exit_codes.json',dict(hisat2=aln,samtools_sort=srt));assert aln==srt==0,'Alignment/sort failed'
 run(['/usr/bin/samtools','quickcheck','-v',str(part)],d/'bam_quickcheck.log');part.rename(bam)
 run(['/usr/bin/samtools','index','-@','2',str(bam)],d/'bam_index.log')
 with (d/'flagstat.txt').open('w') as f:run(['/usr/bin/samtools','flagstat','-@','2',str(bam)],d/'flagstat.log',stdout=f)
 def count_flags(required):
  r=subprocess.run(['/usr/bin/samtools','view','-c','-f',str(required),'-F','2304',str(bam)],capture_output=True,text=True);assert r.returncode==0,r.stderr;return int(r.stdout.strip())
 total=count_flags(64);both_unmapped=count_flags(76);mapped=total-both_unmapped
 quality=json.loads((d/'QC.json').read_text());assert total==quality['clean_pairs'],'Primary paired-fragment conservation failed'
 state(s['run'],'STRAND_INFERENCE')
 run([str(P/'software/venv/bin/infer_experiment.py'),'-i',str(bam),'-r',str(P/'reference/transcripts.bed12'),'-s','200000','-q','30'],d/'strand.log')
 text=(d/'strand.log').read_text();a=re.search(r'Fraction of reads explained by "1\+\+,1--,2\+-,2-\+":\s*([0-9.]+)',text);b=re.search(r'Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s*([0-9.]+)',text)
 assert a and b,'Cannot parse paired library strand inference'
 forward=float(a.group(1));reverse=float(b.group(1));assert forward+reverse>0,'No informative strand observations'
 f=forward/(forward+reverse);r=1-f
 mode='yes' if f>=0.8 else 'reverse' if r>=0.8 else 'no' if max(f,r)<=0.6 else None
 save(d/'STRAND_QA.json',dict(forward_fraction=forward,reverse_fraction=reverse,normalized_forward=f,normalized_reverse=r,HTSeq_stranded=mode,technical_rule='CONFIG.json',status='PASS' if mode else 'STRAND_REVIEW_REQUIRED'))
 assert mode is not None,'STRAND_REVIEW_REQUIRED: BAM retained; no guessed count orientation'
 state(s['run'],'GENE_COUNTING')
 namebam=d/'primary_name_sorted.bam'
 with (d/'primary_filter.log').open('w') as log,(d/'name_sort.log').open('w') as slog:
  view=subprocess.Popen(['/usr/bin/samtools','view','-b','-F','2304',str(bam)],stdout=subprocess.PIPE,stderr=log)
  sort=subprocess.Popen(['/usr/bin/samtools','sort','-n','-@','4','-m','1G','-O','BAM','-T',str(d/'name_sort_tmp'),'-o',str(namebam),'-'],stdin=view.stdout,stderr=slog);view.stdout.close();sr=sort.wait();vr=view.wait()
 assert sr==vr==0,'Primary/name sort failed'
 cmd=[str(P/'software/venv/bin/htseq-count'),'-f','bam','-r','name','-s',mode,'-a','10','-t','exon','-i','gene_id','-m','union','--nonunique','none',str(namebam),ref['gtf']]
 save(d/'htseq_command.json',cmd)
 with (d/'gene_counts.tsv').open('w') as f:run(cmd,d/'htseq.log',stdout=f)
 assigned=0;special={};n=0
 with (d/'gene_counts.tsv').open() as f:
  for line in f:
   gene,value=line.rstrip().split('\t');value=int(value);assert value>=0
   if gene.startswith('__'):special[gene]=value
   else:assigned+=value;n+=1
 assert n>10000 and 0<assigned<=mapped,'Invalid gene counts or fragment unit'
 q=dict(status='PASS',finished_at=now(),raw_pairs=quality['raw_pairs'],clean_pairs=quality['clean_pairs'],primary_fragments=total,mapped_fragments=mapped,assigned_fragments=assigned,assigned_fraction_of_mapped=assigned/mapped,HTSeq_stranded=mode,strand_forward=forward,strand_reverse=reverse,genes_counted=n,HTSeq_special=special,BAM=str(bam),count_sha256=sha(d/'gene_counts.tsv'))
 save(d/'QA.json',q);(d/'COMPLETE').write_text(now()+'\n')
 assert namebam.parent.resolve()==d.resolve() and namebam.name=='primary_name_sorted.bam';namebam.unlink()
 state(s['run'],'COMPLETE');return q

def failed(s,e):
 d=P/'samples'/s['run'];d.mkdir(exist_ok=True);save(d/'FAILURE.json',dict(status='FAILED',error=str(e),traceback=traceback.format_exc(),time=now()));state(s['run'],'FAILED')

def main():
 lock=(P/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 with concurrent.futures.ThreadPoolExecutor(max_workers=1) as rp,concurrent.futures.ThreadPoolExecutor(max_workers=2) as qp,concurrent.futures.ThreadPoolExecutor(max_workers=2) as ap:
  ready=rp.submit(prepare_wait);qfs={qp.submit(qc,s):s for s in samples};afs={}
  for f in concurrent.futures.as_completed(qfs):
   s=qfs[f]
   try:f.result();afs[ap.submit(align_count,s,ready)]=s
   except BaseException as e:failed(s,e)
  for f in concurrent.futures.as_completed(afs):
   s=afs[f]
   try:outcomes[s['run']]=f.result()
   except BaseException as e:failed(s,e)
 failures=[x for x,v in states.items() if v!='COMPLETE'];save(P/'tables/job_outcomes.json',dict(states=states,outcomes=outcomes,failures=failures))
 if failures:raise RuntimeError('Incomplete conversion: '+','.join(failures))
 run([str(P/'software/venv/bin/python'),str(P/'code/aggregate.py')],P/'logs/aggregate.log')
 qa=json.loads((P/'MATRIX_QA.json').read_text());assert qa['status']=='PASS' and qa['sample_libraries']==60
 save(P/'STATUS.json',dict(status='COMPLETE_VERIFIED',finished_at=now(),sample_libraries=60,matrices=str(P/'matrices'),QA=str(P/'MATRIX_QA.json')));(P/'COMPLETE_VERIFIED').write_text(now()+'\n')
code=1
try:main();code=0
except BaseException:
 err=traceback.format_exc();(P/'logs/PIPELINE_ERROR.txt').write_text(err);save(P/'STATUS.json',dict(status='INCOMPLETE_REVIEW_REQUIRED',states=states,error=err,updated_at=now()));emit('PIPELINE_FAILED',error=err)
finally:
 (P/'EXIT_STATUS').write_text(str(code)+'\n');(P/'RUN_FINISHED').write_text(now()+'\n');emit('PIPELINE_EXIT',exit_code=code)
