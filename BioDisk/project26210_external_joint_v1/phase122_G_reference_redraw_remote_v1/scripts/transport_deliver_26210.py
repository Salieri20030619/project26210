from pathlib import Path
import subprocess,io,tarfile,json,hashlib,time
ssh=['ssh','-T','-F','C:/Users/31937/.ssh/config','-o','BatchMode=yes','-o','ClearAllForwardings=yes','-o','UserKnownHostsFile=C:/Users/31937/.ssh/known_hosts','-i','C:/Users/31937/.ssh/id_ed25519','sgs10-codex']
w='/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase122_G_reference_redraw_remote_v1'
def call(code):
 for i in range(4):
  p=subprocess.run(ssh+['python3 -'],input=code.encode(),capture_output=True)
  if p.returncode!=255:break
  time.sleep(3)
 if p.returncode:print(p.stderr.decode(errors='replace'))
 p.check_returncode();return p.stdout
code='''from pathlib import Path
import json,hashlib,tarfile,io,sys
w=Path('''+repr(w)+''');rows=[]
for f in sorted((w/'figures').iterdir()):
 rows.append({'file':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(w/'qc'/'FIGURE_MANIFEST.json').write_text(json.dumps(rows,indent=2))
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as tar:
 for f in sorted((w/'figures').iterdir()):tar.add(f,arcname=f.name)
 tar.add(w/'qc'/'FIGURE_MANIFEST.json',arcname='manifest.json')
'''
b=call(code);local=Path('C:/Users/31937/Downloads/26210/02_2.1_G_I参考图重绘_v3');local.mkdir(exist_ok=True)
with tarfile.open(fileobj=io.BytesIO(b),mode='r:gz') as tar:
 rows=json.load(tar.extractfile('manifest.json'))
 for row in rows:
  content=tar.extractfile(row['file']).read();assert hashlib.sha256(content).hexdigest()==row['sha256'];(local/row['file']).write_bytes(content)
receipt={'status':'PASS','files':len(rows),'all_sha256_match':True,'destination':str(local)}
time.sleep(3)
call('from pathlib import Path\nimport json\nPath('+repr(w+'/qc/DELIVERY_VALIDATION.json')+').write_text(json.dumps('+repr(receipt)+',indent=2))')
print(json.dumps(receipt,ensure_ascii=False))
