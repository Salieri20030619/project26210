from pathlib import Path
import subprocess,zipfile,io
r=Path(__file__).resolve().parents[1];out=Path('C:/Users/31937/Downloads/26210/01_单细胞注释总览_v3');out.mkdir(parents=True,exist_ok=True)
code=(r/'scripts/03_remote_redraw_v3.py').read_text(encoding='utf-8-sig').encode()
p=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10-codex','/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python -'],input=code,capture_output=True,timeout=600)
(r/'qc/remote_stderr_v3.txt').write_bytes(p.stderr);p.check_returncode();z=zipfile.ZipFile(io.BytesIO(p.stdout))
assert all('/' not in n and '\\' not in n for n in z.namelist())
for n in z.namelist():
 dest=out/('核查记录' if n.endswith('.json') else '')/n;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(n))
print('Received',len(z.namelist()),'figure and aggregate files')
