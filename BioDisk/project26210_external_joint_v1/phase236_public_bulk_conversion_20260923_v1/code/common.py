import pathlib,json,datetime,subprocess,hashlib,os,threading
P=pathlib.Path(__file__).resolve().parents[1]
LOCK=threading.RLock()
def now():return datetime.datetime.now().astimezone().isoformat()
def save(q,obj):
 with LOCK:
  t=q.with_name(q.name+'.tmp');t.write_text(json.dumps(obj,indent=2,ensure_ascii=False));t.replace(q)
def sha(q):
 h=hashlib.sha256()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def run(args,log,stdout=None):
 with log.open('w') as err:
  r=subprocess.run([str(x) for x in args],stdout=stdout if stdout is not None else err,stderr=err)
 log.with_name(log.name+'.exit_code').write_text(str(r.returncode)+'\n')
 if r.returncode:raise RuntimeError(str(args[0])+' exit '+str(r.returncode)+'; log '+str(log))
def emit(event,**kw):
 with LOCK:
  row=dict(time=now(),event=event,**kw)
  with (P/'logs/events.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
  print(json.dumps(row),flush=True)
