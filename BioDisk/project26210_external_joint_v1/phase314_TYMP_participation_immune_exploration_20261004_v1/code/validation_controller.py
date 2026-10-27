from pathlib import Path
import subprocess
P=Path(__file__).resolve().parents[1]
code=1
try:
 assert (P/'logs/analyze.exit').read_text().strip()=='0'
 result=subprocess.run(['/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python',str(P/'code/run_stage.py'),'validate'])
 code=result.returncode
finally:
 (P/'logs/validation_controller.exit').write_text(str(code))
raise SystemExit(code)
