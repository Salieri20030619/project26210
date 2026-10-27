from pathlib import Path
import os,sys,time,json,datetime,subprocess,shutil,traceback
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');P=R.parent
O=R/'provenance/automatic_continuation_20260915_v1';S=R/'scripts/automatic_continuation_20260915_v1';N=R/'08_regulons/network_v2_20260915'
PY='/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python';PYS='/mnt/BioDisk/ludingyu/project26210/.venv-pyscenic/bin/python'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def identity(pid):
 try:return Path('/proc/'+str(pid)+'/stat').read_text().rsplit(')',1)[1].split()[19]
 except FileNotFoundError:return None
def state(status,stage,**kw):
 d=dict(time_utc=now(),status=status,stage=stage,pid=os.getpid(),interval_seconds=900,figures_created=False,**kw)
 tmp=O/'STATUS.tmp';tmp.write_text(json.dumps(d,indent=2));tmp.replace(O/'STATUS.json')
 with (O/'EVENTS.jsonl').open('a') as f:f.write(json.dumps(d)+'\n')
 print(json.dumps(d),flush=True)
 note='## 2026-09-15 Phase148 ??????\n??????????900?????????????? '+status+'??? '+stage+'?????part2_rebuild_v1/provenance/automatic_continuation_20260915_v1/STATUS.json?SCENIC??????????????????????????????????????????????????????????????????????????????\n\n'
 for fn in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
  p=P/fn;p.write_text(note+p.read_text())
 p=R/'STATUS.json';d=json.loads(p.read_text());d['remaining_analysis'].update(status=status,current_step=stage,supervisor_pid=os.getpid(),supervisor_status=str(O/'STATUS.json'));p.write_text(json.dumps(d,indent=2,ensure_ascii=False))
 with (P/'RESULTS_MANIFEST.tsv').open('a') as f:f.write('\t'.join(['phase148_auto_'+stage,'phase148','supervision',str(O/'STATUS.json'),now(),status,'900 second checks; fail-stop','No figures'])+'\n')
if '--launch' in sys.argv:
 for p in S.glob('*.py'):compile(p.read_text(),str(p),'exec')
 O.mkdir(exist_ok=False);b=O/'state_before';b.mkdir()
 for p in [P/'CODEX_HANDOFF.md',P/'PROJECT_STATUS.md',P/'DECISIONS.md',P/'RESULTS_MANIFEST.tsv',R/'STATUS.json']:shutil.copy2(p,b/p.name)
 with (O/'supervisor.log').open('ab') as f:
  c=subprocess.Popen([PY,str(Path(__file__).resolve()),'--run'],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
 (O/'PID.txt').write_text(str(c.pid));print('SUPERVISOR_LAUNCHED',c.pid);sys.exit(0)
if '--wait' in sys.argv:
 pid=int((O/'PID.txt').read_text());token=identity(pid)
 print('PERSISTENT_CONNECTION_ACTIVE interval=900 pid='+str(pid),flush=True)
 while token is not None and identity(pid)==token:
  if (O/'STATUS.json').exists():
   d=json.loads((O/'STATUS.json').read_text())
   if d['status'] in ['FAILED','CONTENT_DRAFT_READY']:print(json.dumps(d),flush=True);sys.exit(0)
  time.sleep(900)
 print((O/'STATUS.json').read_text(),flush=True);sys.exit(0)
stage='waiting_scenic'
try:
 state('WAITING',stage);pid=int((N/'PID.txt').read_text());token=identity(pid)
 while True:
  n=json.loads((N/'STATUS.json').read_text());alive=token is not None and identity(pid)==token
  with (O/'MONITOR.jsonl').open('a') as f:f.write(json.dumps(dict(time_utc=now(),scenic=n,alive=alive))+'\n')
  if n['status']=='COMPLETE_SCENIC_CORE_VERIFIED':
   assert (N/'RUN_COMPLETE').exists()
   for v in ['GRNBoost2','motif_pruning','AUCell']:assert json.loads((N/(v+'_VERIFIED.json')).read_text())['status']=='PASS'
   break
  if n['status']=='FAILED' or not alive:raise RuntimeError('SCENIC_NOT_VERIFIED '+str(n))
  time.sleep(900)
 env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUNBUFFERED='1')
 for stage,python,script in [('registry',PYS,'24_extract_new_regulon_targets.py'),('patient_statistics',PY,'25_new_regulon_statistics_v2.py'),('crosswalk',PY,'26_crosswalk_and_dependencies.py'),('validation',PY,'29_validate_remaining_results_v2.py'),('technical_draft',PY,'30_build_no_figure_content_review.py')]:
  state('RUNNING',stage)
  with (O/'RESOURCE_GATE.log').open('a') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,stderr=subprocess.STDOUT,check=True)
  cmd=[python,str(S/script)]
  with (O/(stage+'.log')).open('wb') as f:
   c=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,stdin=subprocess.DEVNULL)
   (O/'CHILD_PID.txt').write_text(str(c.pid));code=c.wait()
  (O/(stage+'_EXIT.json')).write_text(json.dumps(dict(returncode=code,time_utc=now(),command=cmd)))
  if code!=0:raise RuntimeError(stage+' failed returncode='+str(code))
 state('CONTENT_DRAFT_READY','await_biological_review')
except BaseException as e:
 state('FAILED',stage,error=repr(e));traceback.print_exc();sys.exit(1)

