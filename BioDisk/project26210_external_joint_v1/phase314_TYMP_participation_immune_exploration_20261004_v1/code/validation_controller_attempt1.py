from pathlib import Path
import os,select,subprocess
P=Path(__file__).resolve().parents[1]
pid=int((P/'logs/analyze.pid').read_text())
try:
 fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
except ProcessLookupError:pass
code=1
try:
 assert (P/'logs/analyze.exit').read_text().strip()=='0'
 q=subprocess.run(['/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python',str(P/'code/run_stage.py'),'validate'])
 code=q.returncode
finally:(P/'logs/validation_controller.exit').write_text(str(code))
raise SystemExit(code)
