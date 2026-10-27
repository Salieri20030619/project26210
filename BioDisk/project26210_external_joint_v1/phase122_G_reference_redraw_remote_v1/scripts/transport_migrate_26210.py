from pathlib import Path
import subprocess,tarfile,io,json,hashlib,datetime,shlex
root=Path('D:/project26210_external_joint_v1')
remote='/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1'
stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
dest=remote+'/local_records_migration_'+stamp
ssh=['ssh','-T','-F','C:/Users/31937/.ssh/config','-o','BatchMode=yes','-o','ClearAllForwardings=yes','-o','UserKnownHostsFile=C:/Users/31937/.ssh/known_hosts','-i','C:/Users/31937/.ssh/id_ed25519','sgs10-codex']
def run(code):
 p=subprocess.run(ssh+['python3 -'],input=code.encode(),capture_output=True);print(p.stdout.decode(errors='replace'),p.stderr.decode(errors='replace'),flush=True);p.check_returncode();return p.stdout
run('from pathlib import Path\np=Path('+repr(dest)+');p.mkdir()\n')
record_extensions={'.md','.txt','.json','.jsonl','.tsv','.csv','.log','.py','.r','.sh','.ps1','.rmd','.yml','.yaml','.html','.gz','.rproj','.renviron','.gitignore'}
records=[];skipped=[]
for f in root.rglob('*'):
 if not f.is_file():continue
 rel=f.relative_to(root).as_posix()
 if f.suffix.lower() not in record_extensions or any(s in rel for s in ['/site-packages/','/__pycache__/']):
  skipped.append({'path':rel,'bytes':f.stat().st_size,'reason':'binary data, figure/report artifact, runtime or incomplete cache; outside textual project records migration'});continue
 records.append((f,'project/'+rel))
for base,prefix in [(Path('C:/Users/31937/project26210_figure_redraw_20260908_v1'),'references')]:
 for f in base.rglob('*'):
  if f.is_file():records.append((f,prefix+'/'+f.relative_to(base).as_posix()))
print('Transferring',len(records),'record/artifact files;',sum(f.stat().st_size for f,_ in records),'bytes',flush=True)
p=subprocess.Popen(ssh+['tar -xzf - -C '+shlex.quote(dest)],stdin=subprocess.PIPE)
manifest=[]
with tarfile.open(fileobj=p.stdin,mode='w|gz',compresslevel=9) as tar:
 for f,rel in records:
  h=hashlib.sha256()
  with f.open('rb') as rd:
   for b in iter(lambda:rd.read(1024*1024),b''):h.update(b)
  tar.add(f,arcname=rel,recursive=False)
  manifest.append({'source':str(f),'path':rel,'bytes':f.stat().st_size,'sha256':h.hexdigest()})
 for name,obj in [('MIGRATION_MANIFEST.json',manifest),('EXCLUDED_DATA_RUNTIME_INVENTORY.json',skipped)]:
  b=json.dumps(obj,ensure_ascii=False,indent=2).encode();i=tarfile.TarInfo(name);i.size=len(b);tar.addfile(i,io.BytesIO(b))
p.stdin.close();assert p.wait()==0
verify='''from pathlib import Path
import json,hashlib,shutil,datetime
p=Path(DEST);root=p.parent
rows=json.loads((p/'MIGRATION_MANIFEST.json').read_text());bad=[]
for r in rows:
 f=p/r['path'];h=hashlib.sha256()
 with f.open('rb') as rd:
  for b in iter(lambda:rd.read(1048576),b''):h.update(b)
 if h.hexdigest()!=r['sha256'] or f.stat().st_size!=r['bytes']:bad.append(r['path'])
result={'status':'PASS' if not bad else 'FAIL','verified_files':len(rows),'mismatches':bad,'bytes':sum(r['bytes'] for r in rows),'excluded_inventory':'EXCLUDED_DATA_RUNTIME_INVENTORY.json'}
(p/'MIGRATION_VALIDATION.json').write_text(json.dumps(result,indent=2));assert not bad
backup=p/'previous_server_state';backup.mkdir()
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']:
 if (root/name).exists():shutil.copy2(root/name,backup/name)
 if name!='AGENTS.md':shutil.copy2(p/'project'/name,root/name)
note=''' + repr('## 2026-09-08 服务器唯一权威记录入口\n\n用户要求：本地项目记录全部迁移到服务器，之后不再在本地记录。\n当前唯一权威根目录：`/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1`。以下历史段落中 D 盘权威入口已失效；Windows 历史路径通过迁移清单映射到服务器快照。\n当前任务：Phase121 未完成的 G/I 参考格式重绘；用户最新反馈 G 仍不满意，要求按参考图整体重做。已接受的单细胞总览不重做。仅图件及非识别性汇总交付 Downloads/26210，分析与日志留服务器。\n') + '''
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=root/name;f.write_text(note+'\\n完整历史记录及脚本/图件快照：`'+str(p)+'`。哈希核查见 MIGRATION_VALIDATION.json。\\n\\n'+f.read_text(encoding='utf-8-sig'),encoding='utf8')
(root/'AGENTS.md').write_text('# Project26210 server operating rules\\n\\n'+note+'\\nPreserve source data and previous outputs. New work uses versioned server directories. Record scripts, input identity, software versions, statistical units and validation. Update PROJECT_STATUS.md, DECISIONS.md, RESULTS_MANIFEST.tsv and CODEX_HANDOFF.md on completion. No new local progress records. Historical Windows scripts are archived provenance, not directly executable server scripts.\\n')
with (root/'RESULTS_MANIFEST.tsv').open('a') as f:f.write('migration_'+p.name+'\\tserver_records_migration\\tvalidated_archive\\t'+str(p)+'\\t'+datetime.datetime.now().isoformat()+'\\tPASS\\tSHA256 verified '+str(len(rows))+' files\\tServer is sole future record authority; raw binary objects listed separately\\n')
print(json.dumps({'destination':str(p),**result}))
'''.replace('DEST',repr(dest))
run(verify)
