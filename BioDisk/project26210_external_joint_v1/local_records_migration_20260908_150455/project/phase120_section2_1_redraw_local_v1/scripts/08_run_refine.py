from pathlib import Path
import subprocess,zipfile,io
r=Path(__file__).resolve().parents[1];d=Path('C:/Users/31937/Downloads/26210/02_2.1恶性上皮亚克隆_v1');d.mkdir(parents=True,exist_ok=True)
p=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10-codex','/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python -'],input=(r/'scripts/07_remote_refine.py').read_text(encoding='utf-8-sig').encode(),capture_output=True,timeout=600)
(r/'qc/remote_refine_stderr.txt').write_bytes(p.stderr);p.check_returncode();z=zipfile.ZipFile(io.BytesIO(p.stdout));assert all('/' not in n and '\\' not in n for n in z.namelist())
for n in z.namelist():
 f=d/('核查记录' if n.endswith(('.json','.csv')) else '')/n;f.parent.mkdir(parents=True,exist_ok=True);f.write_bytes(z.read(n))
print('Received',len(z.namelist()),'figures and aggregates')
