from pathlib import Path
import subprocess,os,select,json,traceback
P=Path(__file__).resolve().parents[1];rc=1
try:
 f=P/'logs/validate_segmentation.exit'
 if not f.exists():
  try:
   fd=os.pidfd_open(int((P/'logs/validate_segmentation.pid').read_text()));select.select([fd],[],[]);os.close(fd)
  except ProcessLookupError:pass
 assert f.exists() and f.read_text().strip()=='0','segmentation validation failed'
 assert json.loads((P/'qc/PROSEG_IDENTITY_VALIDATION.json').read_text())['status']=='PASS'
 exe='/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python';env=os.environ.copy();env.update({'OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
 for stage in ['run_regional','validate_regional']:
  r=subprocess.run([exe,'-u',str(P/'code/run_stage.py'),stage],env=env);assert r.returncode==0,stage
 rc=0
except Exception:traceback.print_exc()
finally:(P/'logs/regional_controller.exit').write_text(str(rc)+'\n')
raise SystemExit(rc)
