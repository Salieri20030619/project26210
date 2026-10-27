from pathlib import Path
import os,time,json,datetime,subprocess
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'provenance/low_frequency_supervision_v1';O.mkdir(exist_ok=True)
PIPELINE_PID=1855310
INTERVAL=900

def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def procstart(pid):
 try:return Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):return None
start=procstart(PIPELINE_PID);previous={};idle=0
while True:
 now=stamp();exitfile=R/'logs/phase148_sequence_v2_exit.txt'
 if exitfile.exists():
  text=exitfile.read_text().strip();state={'time_utc':now,'status':'COMPUTATION_EXITED','exit_record':text,'next_action':'Review results and run remaining independent validation; no figures.','interval_seconds':INTERVAL}
  (O/'STATUS.json').write_text(json.dumps(state,indent=2));print(json.dumps(state),flush=True);break
 if start is None or procstart(PIPELINE_PID)!=start:
  state={'time_utc':now,'status':'ALERT_PIPELINE_PROCESS_EXITED_WITHOUT_MARKER','pipeline_pid':PIPELINE_PID};(O/'STATUS.json').write_text(json.dumps(state,indent=2));print(json.dumps(state),flush=True);break
 lines=subprocess.check_output(['ps','-eo','pid,ppid,stat,etime,pcpu,rss,args'],text=True).splitlines()[1:];rows={}
 for line in lines:
  b=line.split(None,6)
  if len(b)==7:rows[int(b[0])]={'pid':int(b[0]),'ppid':int(b[1]),'state':b[2],'elapsed':b[3],'cpu_lifetime_pct':float(b[4]),'rss_kb':int(b[5]),'command':b[6]}
 ids={PIPELINE_PID}
 for _ in range(10):
  nxt=ids|{pid for pid,v in rows.items() if v['ppid'] in ids}
  if nxt==ids:break
  ids=nxt
 procs=[rows[i] for i in sorted(ids) if i in rows];ticks={}
 for v in procs:
  try:
   parts=Path(f"/proc/{v['pid']}/stat").read_text().rsplit(')',1)[1].split();ticks[str(v['pid'])]=int(parts[11])+int(parts[12])
  except (FileNotFoundError,ProcessLookupError):pass
 delta=sum(max(0,val-previous.get(pid,val)) for pid,val in ticks.items());idle=idle+1 if previous and delta==0 else 0
 logs=[]
 for name in ['phase148_pyscenic_v1.log','phase148_regulon_statistics_v1.log','phase148_crosswalk_v1.log']:
  p=R/'logs'/name
  if p.exists():logs.append({'name':name,'bytes':p.stat().st_size,'modified_epoch':p.stat().st_mtime})
 state={'time_utc':now,'status':'ALERT_NO_CPU_PROGRESS_30MIN' if idle>=2 else 'RUNNING','pipeline_pid':PIPELINE_PID,'interval_seconds':INTERVAL,'cpu_ticks_since_previous_check':delta if previous else None,'descendant_processes':len(procs),'processes':procs,'logs':logs,'automatic_termination':False,'figures_authorized':False}
 (O/'STATUS.json').write_text(json.dumps(state,indent=2))
 with (O/'CHECKS.jsonl').open('a') as f:f.write(json.dumps(state)+'\n')
 previous=ticks
 time.sleep(INTERVAL)
