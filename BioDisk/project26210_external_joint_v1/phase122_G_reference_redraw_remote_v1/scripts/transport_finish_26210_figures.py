from pathlib import Path
import subprocess,json,hashlib,io,tarfile,shlex,time
ssh=['ssh','-T','-F','C:/Users/31937/.ssh/config','-o','BatchMode=yes','-o','ClearAllForwardings=yes','-o','UserKnownHostsFile=C:/Users/31937/.ssh/known_hosts','-i','C:/Users/31937/.ssh/id_ed25519','sgs10-codex']
root='/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1'
w=root+'/phase122_G_reference_redraw_remote_v1'
def remote(code):
 for attempt in range(4):
  p=subprocess.run(ssh+['python3 -'],input=code.encode(),capture_output=True)
  if p.returncode!=255:break
  time.sleep(2)
 print(p.stdout.decode(errors='replace'),p.stderr.decode(errors='replace'),flush=True);p.check_returncode();return p.stdout
src=Path('D:/project26210_external_joint_v1/phase121_G_I_reference_revision_local_v1/scripts/02_I_remote.py').read_text(encoding='utf-8-sig')
src=src.replace("['violin_points','box_points']","['violin_points']").replace("bb.set_edgecolor(col);bb.set_alpha(.4)","bb.set_edgecolor('#333333');bb.set_alpha(.85)")
src=src.replace("z.close();sys.stdout.buffer.write(buf.getvalue())","z.close();dest=Path("+repr(w)+")/'figures'; archive=zipfile.ZipFile(io.BytesIO(buf.getvalue()));\nfor name in archive.namelist():\n target=(dest.parent/'qc' if name.endswith('.json') else dest)/name.replace('_v2.','_v3.');target.write_bytes(archive.read(name))\nprint('I_RENDER_COMPLETE')")
remote('from pathlib import Path\nPath('+repr(w+'/scripts/draw_I.py')+').write_text('+repr(src)+')')
p=subprocess.run(ssh+['test -f '+w+'/figures/I_marker_violin_points_v3.png || /mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python '+w+'/scripts/draw_I.py > '+w+'/qc/render_I.log 2>&1']);p.check_returncode()
remote('''from pathlib import Path
import json,hashlib
w=Path('''+repr(w)+''');rows=[]
for f in sorted((w/'figures').iterdir()):
 rows.append({'file':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(w/'qc'/'FIGURE_MANIFEST.json').write_text(json.dumps(rows,indent=2))
print(json.dumps(rows))
''')
local=Path('C:/Users/31937/Downloads/26210/02_2.1_G_I参考图重绘_v3');local.mkdir(exist_ok=True)
# Receive figures/readme only. Keep all audit records and cell data on server.
p=subprocess.run(ssh+['tar -czf - -C '+shlex.quote(w+'/figures')+' .'],capture_output=True);p.check_returncode()
with tarfile.open(fileobj=io.BytesIO(p.stdout),mode='r:gz') as tar:
 for m in tar.getmembers():
  if not m.isfile():continue
  name=Path(m.name).name;assert m.name in ['./'+name,name]
  (local/name).write_bytes(tar.extractfile(m).read())
p=subprocess.run(ssh+['cat '+w+'/qc/FIGURE_MANIFEST.json'],capture_output=True);p.check_returncode();rows=json.loads(p.stdout)
for row in rows:assert hashlib.sha256((local/row['file']).read_bytes()).hexdigest()==row['sha256']
remote('from pathlib import Path\nimport json\nPath('+repr(w+'/qc/DELIVERY_VALIDATION.json')+').write_text(json.dumps('+repr({'status':'PASS','files':len(rows),'all_sha256_match':True,'destination':str(local)})+',indent=2))')
print('DELIVERED '+str(local))
