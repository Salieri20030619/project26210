from pathlib import Path
import subprocess,os,select,traceback
P=Path(__file__).resolve().parents[1];rc=1
try:
 f=P/'logs/proseg.exit'
 if not f.exists():
  try:
   fd=os.pidfd_open(int((P/'logs/proseg.pid').read_text()));select.select([fd],[],[]);os.close(fd)
  except ProcessLookupError:pass
 assert f.exists() and f.read_text().strip()=='0','Proseg failed or no recorded exit'
 env=os.environ.copy();env.update({'OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMBA_NUM_THREADS':'4'})
 exe='/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python'
 with (P/'logs/cluster_identity.log').open('w') as log:r=subprocess.run([exe,'-u',str(P/'code/cluster_identity.py')],stdout=log,stderr=subprocess.STDOUT,env=env)
 (P/'logs/cluster_identity.exit').write_text(str(r.returncode)+'\n');rc=r.returncode
except Exception:traceback.print_exc()
finally:(P/'logs/identity_controller.exit').write_text(str(rc)+'\n')
raise SystemExit(rc)
