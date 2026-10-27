from pathlib import Path
import os,select,subprocess,json,traceback
P=Path(__file__).resolve().parents[1];rc=1
try:
 try:
  fd=os.pidfd_open(int((P/'logs/acquire_region1.pid').read_text()));select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
 assert (P/'logs/acquire_region1.exit').read_text().strip()=='0','download failed'
 assert json.loads((P/'qc/REGION1_DOWNLOAD_VALIDATION.json').read_text())['status']=='PASS'
 env=os.environ.copy();env.update({'OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMBA_NUM_THREADS':'4'})
 exe='/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python'
 for stage in ['reconstruct_region1','cluster_region1']:
  r=subprocess.run([exe,'-u',str(P/'code/run_stage.py'),stage],env=env);assert r.returncode==0,stage
 rc=0
except Exception:traceback.print_exc()
finally:(P/'logs/region1_controller.exit').write_text(str(rc)+'\n')
raise SystemExit(rc)
