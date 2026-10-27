from pathlib import Path
import os,ctypes,select,json,subprocess,datetime
P=Path(__file__).resolve().parents[1];SRC=P.parent/'phase249_independent_spatial_sources_20260924_v1'
pid=int((SRC/'DOWNLOADER_PID').read_text());proc=Path('/proc')/str(pid)/'cmdline'
(P/'STATUS.json').write_text(json.dumps(dict(status='WAITING_FOR_VERIFIED_DOWNLOAD',at=datetime.datetime.now().astimezone().isoformat(),download_pid=pid),indent=2))
if proc.exists() and b'download_emt.py' in proc.read_bytes():
 fd=ctypes.CDLL(None,use_errno=True).syscall(434,pid,0)
 if fd>=0:select.select([fd],[],[]);os.close(fd)
assert (SRC/'DOWNLOAD_EXIT_STATUS').exists() and (SRC/'DOWNLOAD_EXIT_STATUS').read_text().strip()=='0','download failed; do not start analysis'
with open(P/'logs/analysis.log','w') as f:
 result=subprocess.run(['/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python','-u',str(P/'code/analyze_mainline.py')],stdout=f,stderr=subprocess.STDOUT)
(P/'RUNNER_EXIT_STATUS').write_text(str(result.returncode)+'\n')
raise SystemExit(result.returncode)

