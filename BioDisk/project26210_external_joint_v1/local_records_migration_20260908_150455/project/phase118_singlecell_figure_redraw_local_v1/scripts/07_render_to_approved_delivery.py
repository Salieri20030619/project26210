from pathlib import Path
import subprocess,zipfile,io,json,shutil
root=Path(__file__).resolve().parents[1]
out=Path('C:/Users/31937/Downloads/26210/01_单细胞注释总览_v1');out.mkdir(parents=True,exist_ok=True)
remote_code=(root/'scripts/06_remote_render.py').read_text(encoding='utf-8-sig').encode()
r=subprocess.run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10-codex','/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python -'],input=remote_code,capture_output=True,timeout=900)
(root/'qc/remote_render_stderr.txt').write_bytes(r.stderr)
r.check_returncode();z=zipfile.ZipFile(io.BytesIO(r.stdout))
assert all('/' not in n and '\\' not in n for n in z.namelist())
for n in z.namelist():
 sub='单个Marker' if n.startswith('Marker_') else ('核查记录' if n.endswith('.json') else '')
 dest=out/sub/n;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(n))
for p in (root/'figures').glob('*'):shutil.copy2(p,out/p.name)
print(json.dumps({'delivery':str(out),'files_received':len(z.namelist())},ensure_ascii=False))
