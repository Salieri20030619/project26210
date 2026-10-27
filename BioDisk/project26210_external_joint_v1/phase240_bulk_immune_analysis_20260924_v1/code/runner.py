import pathlib,subprocess,sys,json,datetime,traceback
P=pathlib.Path(__file__).resolve().parents[1]
def now():return datetime.datetime.now().astimezone().isoformat()
try:
 for name in ['preflight','pipeline']:
  with open(P/'logs'/f'{name}.log','w') as f:
   r=subprocess.run([sys.executable,'-u',str(P/'code'/f'{name}.py')],stdout=f,stderr=subprocess.STDOUT)
  (P/'logs'/f'{name}.exit_code').write_text(str(r.returncode)+'\n')
  if r.returncode:raise RuntimeError(name+' failed; see logs')
except Exception:
 if not (P/'EXIT_STATUS').exists():(P/'EXIT_STATUS').write_text('1\n')
 (P/'RUN_FINISHED').write_text(now())
 (P/'STATUS.json').write_text(json.dumps({'status':'FAILED','at':now(),'error':traceback.format_exc()},indent=2))
 raise
