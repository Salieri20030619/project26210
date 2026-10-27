from pathlib import Path
import os,select,subprocess,sys,json,traceback,datetime
P=Path(__file__).resolve().parents[1]
rc=1
try:
 f=P/'logs/acquire_pilot_v2.exit'
 if not f.exists():
  pid=int((P/'logs/acquire_pilot_v2.pid').read_text())
  try:
   fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
  except ProcessLookupError:pass
 assert f.exists() and f.read_text().strip()=='0','acquisition not successful'
 assert (P/'qc/PILOT_INPUT_VALIDATION.json').exists()
 exe='/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python'
 with (P/'logs/reconstruct_pilot.log').open('w') as log:r=subprocess.run([exe,'-u',str(P/'code/reconstruct_pilot.py')],stdout=log,stderr=subprocess.STDOUT)
 (P/'logs/reconstruct_pilot.exit').write_text(str(r.returncode)+'\n');assert r.returncode==0
 rc=0
except Exception:traceback.print_exc()
finally:(P/'logs/controller.exit').write_text(str(rc)+'\n')
print('CONTROLLER_EXIT',rc,datetime.datetime.now().astimezone().isoformat(),flush=True)
raise SystemExit(rc)
