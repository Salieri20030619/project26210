from pathlib import Path
import os,select,json,csv,datetime,hashlib,sys,fcntl
P=Path(__file__).resolve().parents[1];R=P.parent
lock=open(P/'logs/completion_watcher.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
pid=int((P/'logs/LAUNCH_PID').read_text())
try:
 fd=os.pidfd_open(pid);poll=select.poll();poll.register(fd,select.POLLIN);poll.poll();os.close(fd)
except ProcessLookupError:pass
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
expected=list(csv.DictReader((P/'metadata/PUBLIC_FASTQ_MANIFEST.tsv').open(),delimiter='\t'))
errors=[];count=0;size=0
for row in expected:
 for url,nbytes,md5 in zip(row['fastq_urls'].split(';'),row['fastq_bytes'].split(';'),row['fastq_md5'].split(';')):
  name=url.rsplit('/',1)[1];f=P/'fastq'/name;j=P/'verification'/(name+'.json')
  try:
   e=json.loads(j.read_text())
   assert e['status']=='VERIFIED' and e['observed_md5']==md5 and e['observed_bytes']==int(nbytes) and f.stat().st_size==int(nbytes)
   count+=1;size+=int(nbytes)
  except Exception as exc:errors.append({'file':name,'error':str(exc)})
exitcode=(P/'EXIT_STATUS').read_text().strip() if (P/'EXIT_STATUS').exists() else 'MISSING'
passed=count==120 and exitcode=='0' and (P/'COMPLETE_VERIFIED').exists() and not errors
v=dict(checked_at=now,download_exit_status=exitcode,verified_files=count,verified_bytes=size,completion_marker_exists=(P/'COMPLETE_VERIFIED').exists(),status='PASS' if passed else 'INCOMPLETE',errors=errors,verification_method='Check all per-file measured MD5 records, official expected MD5, current size, exit status and completion marker')
(P/'FINAL_DOWNLOAD_VALIDATION.json').write_text(json.dumps(v,indent=2))
state='下载完成并通过120文件验收' if passed else '下载未完整通过，须查看失败或中断记录'
body='\n\n## '+now+' Phase216下载结束自动验收\n'+state+'；已验证'+str(count)+'/120文件、'+str(round(size/1e9,3))+'GB，退出状态'+exitcode+'。详见phase216_paper_dataset_download_plan_20260919_v1/FINAL_DOWNLOAD_VALIDATION.json。样本处理/配对冲突仍未解决，分析计划已完成但D2–D5未启动。无自动生物学分析。\n'
with (R/'phase216_record_update.lock').open('w') as rec:
 fcntl.flock(rec,fcntl.LOCK_EX)
 for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
  f=R/n;f.write_text(body+f.read_text())
 with (R/'RESULTS_MANIFEST.tsv').open('a') as f:
  f.write('\t'.join(['phase216_download_final_validation','phase216','validation',str(P/'FINAL_DOWNLOAD_VALIDATION.json'),now,v['status'],'EXIT_STATUS; COMPLETE_VERIFIED; file size and measured MD5 records','Downloaded '+str(count)+'/120 FASTQs; pairing unresolved'])+'\n')
(P/'COMPLETION_WATCHER_EXIT').write_text('0\n' if passed else '1\n')
print(json.dumps(v),flush=True)
