from pathlib import Path
import os,signal,select,json,datetime,shutil,csv,fcntl
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
P=R/'phase216_paper_dataset_download_plan_20260919_v1'
stamp=datetime.datetime.now().astimezone();when=stamp.isoformat();D=P/'pause_records'/stamp.strftime('%Y%m%d_%H%M%S');D.mkdir(parents=True,exist_ok=False)
for name in ['STATUS.json','RESUME_STATUS.json','DOWNLOAD_SUMMARY.json','EXIT_STATUS','COMPLETE_VERIFIED']:
 f=P/name
 if f.exists():shutil.copy2(f,D/name)
before=[]
def cmdline(pid):
 try:return Path('/proc')/str(pid)/'cmdline'
 except Exception:return None
def inspect(pid):
 p=Path('/proc')/str(pid)
 try:
  args=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode()
  state=(p/'stat').read_text().split(') ',1)[1].split()[0]
  return dict(pid=pid,pgid=os.getpgid(pid),cmd=args,state=state)
 except (FileNotFoundError,ProcessLookupError):return None
def finish(pid,group=False):
 info=inspect(pid)
 if info is None:return
 assert str(P/'scripts') in info['cmd'],info
 if group:assert info['pgid']==pid
 fd=os.pidfd_open(pid)
 os.killpg(pid,signal.SIGTERM) if group else os.kill(pid,signal.SIGTERM)
 poll=select.poll();poll.register(fd,select.POLLIN)
 done=bool(poll.poll(5000))
 if not done:
  os.killpg(pid,signal.SIGKILL) if group else os.kill(pid,signal.SIGKILL)
  done=bool(poll.poll(5000))
 os.close(fd);assert done,info
 return dict(**info,termination='SIGTERM' if done else 'UNKNOWN')
root=inspect(3002854);watcher=inspect(3002855)
for v in [root,watcher]:
 if v:before.append(v)
# Stop completion watcher first, so intentional pause is not mistaken for completion.
if watcher:
 assert 'wait_and_finalize.py' in watcher['cmd']
 finish(watcher['pid'])
if root:
 assert 'download_public_fastq.py' in root['cmd']
 for f in Path('/proc').iterdir():
  if f.name.isdigit():
   v=inspect(int(f.name))
   if v and v['pgid']==root['pgid']:
    assert str(P) in v['cmd'],v
    if v['pid']!=root['pid']:before.append(v)
 finish(root['pid'],group=True)
# Verify every downloader/worker in this project has stopped, including any surviving child.
remaining=[]
for f in Path('/proc').iterdir():
 if f.name.isdigit():
  v=inspect(int(f.name))
  if v and v['state']!='Z' and str(P) in v['cmd'] and any(t in v['cmd'] for t in ['download_public_fastq.py','wait_and_finalize.py','curl ','wget ','aria2']):
   remaining.append(v)
assert not remaining,remaining
expected=list(csv.DictReader((P/'metadata/PUBLIC_FASTQ_MANIFEST.tsv').open(),delimiter='\t'))
verified=[];partials=[];other=[]
for row in expected:
 for url,size,md5 in zip(row['fastq_urls'].split(';'),row['fastq_bytes'].split(';'),row['fastq_md5'].split(';')):
  name=url.rsplit('/',1)[1];f=P/'fastq'/name;e=P/'verification'/(name+'.json')
  if f.exists() and e.exists():
   ev=json.loads(e.read_text())
   if ev.get('status')=='VERIFIED' and ev.get('observed_md5')==md5 and ev.get('observed_bytes')==int(size) and f.stat().st_size==int(size):
    verified.append(dict(file=name,bytes=int(size),md5=md5))
  p=(P/'fastq'/name).with_suffix('.gz.part')
  if p.exists():partials.append(dict(file=p.name,bytes=p.stat().st_size,expected_bytes=int(size)))
allfiles=[dict(file=f.name,bytes=f.stat().st_size) for f in (P/'fastq').iterdir() if f.is_file()]
state=dict(status='PAUSED_BY_USER',paused_at=when,reason='User explicitly requested pause background downloads',expected_files=120,verified_files=len(verified),verified_bytes=sum(x['bytes'] for x in verified),partial_files=len(partials),partial_bytes=sum(x['bytes'] for x in partials),total_retained_bytes=sum(x['bytes'] for x in allfiles),download_processes_remaining=0,completion_watcher_stopped=True,source_files_deleted=0,automatic_resume=False,biological_analysis_started=False,pause_record=str(D),completion=False)
(D/'PROCESS_STOP_RECORD.json').write_text(json.dumps(dict(before=before,remaining=remaining,action='stop watcher then own download process group with SIGTERM; pidfd wait; preserve all data'),indent=2))
(D/'RETAINED_FILES.json').write_text(json.dumps(dict(verified=verified,partial=partials,all_files=allfiles),indent=2))
(D/'PAUSE_VALIDATION.json').write_text(json.dumps(dict(**state,validation_status='PASS_PAUSED'),indent=2))
for name in ['STATUS.json','RESUME_STATUS.json']:(P/name).write_text(json.dumps(state,indent=2))
(P/'PAUSED_PROGRESS_CN.md').write_text(f'# 下载已按用户要求暂停\n\n暂停时间：{when}\n\n下载主进程、4个下载工作进程及完成监听均已停止；断点文件保留，不自动恢复。\n已验证文件：{len(verified)}/120，{state["verified_bytes"]/1e9:.3f}GB；片段：{len(partials)}个，{state["partial_bytes"]/1e9:.3f}GB。\n本记录是主动暂停，不是下载完成验收。\n记录：{D}\n')
body=f'## {when} 用户要求暂停后台下载\n\nPhase216 FASTQ下载主进程3002854及其工作进程、完成监听3002855均已停止，已核对无残留下载进程。保留{len(verified)}/120个已验证文件（{state["verified_bytes"]/1e9:.3f}GB）和{len(partials)}个片段（{state["partial_bytes"]/1e9:.3f}GB），不自动恢复。状态PAUSED_BY_USER，断点记录{D}。Phase235试点结果保持，本次只暂停下载，不宣称下载完成。\n\n'
with (R/'phase216_record_update.lock').open('a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
  f=R/name;shutil.copy2(f,D/name);f.write_text(body+f.read_text())
 with (R/'RESULTS_MANIFEST.tsv').open('a',newline='') as f:
  csv.writer(f,delimiter='\t').writerow(['phase216_pause_'+stamp.strftime('%Y%m%d_%H%M%S'),'phase216','pause_record',str(D/'PAUSE_VALIDATION.json'),when,'PASS_PAUSED','pidfd process termination; no residual download workers; retained file inventory','Paused by explicit user request; no automatic resume'])
print(json.dumps(state,ensure_ascii=False),flush=True)
