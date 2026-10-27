from pathlib import Path
import subprocess,json,tarfile,io,hashlib,time
ssh=['ssh','-T','-F','C:/Users/31937/.ssh/config','-o','BatchMode=yes','-o','ClearAllForwardings=yes','-o','UserKnownHostsFile=C:/Users/31937/.ssh/known_hosts','-i','C:/Users/31937/.ssh/id_ed25519','sgs10-codex']
root='/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1';w=root+'/phase122_G_reference_redraw_remote_v1'
def call(code):
 for i in range(4):
  p=subprocess.run(ssh+['python3 -'],input=code.encode(),capture_output=True)
  if p.returncode!=255:break
  time.sleep(3)
 print(p.stdout.decode(errors='replace'),p.stderr.decode(errors='replace'),flush=True);p.check_returncode();return p.stdout
helpers={f.name:f.read_text(encoding='utf8') for f in Path('C:/codex_test').glob('*26210*.py')}
supp=[];b=io.BytesIO()
with tarfile.open(fileobj=b,mode='w:gz') as tar:
 for f in Path('D:/project26210_external_joint_v1').rglob('*'):
  if not f.is_file() or f.suffix or any(s in f.parts for s in ['site-packages','__pycache__']):continue
  raw=f.read_bytes()
  if b'\0' in raw:continue
  try:raw.decode('utf8')
  except UnicodeDecodeError:continue
  rel=f.relative_to('D:/project26210_external_joint_v1').as_posix();tar.add(f,arcname=rel);supp.append({'path':rel,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
import base64
code='''from pathlib import Path
import json,hashlib,datetime,tarfile,io,base64
root=Path(ROOT);w=Path(W)
snapshots=[p for p in sorted(root.glob('local_records_migration_*')) if (p/'MIGRATION_VALIDATION.json').exists()]
assert snapshots,'MIGRATION_STILL_PENDING'
m=snapshots[-1];mv=json.loads((m/'MIGRATION_VALIDATION.json').read_text());assert mv['status']=='PASS'
for name,body in HELPERS.items():(w/'scripts'/('transport_'+name)).write_text(body)
sdir=m/'extensionless_status_records';sdir.mkdir(exist_ok=True)
with tarfile.open(fileobj=io.BytesIO(base64.b64decode(SUPP_BYTES)),mode='r:gz') as tar:
 for item in tar.getmembers():
  target=sdir/item.name;assert target.resolve().is_relative_to(sdir.resolve());target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(tar.extractfile(item).read())
for row in SUPP_ROWS:assert hashlib.sha256((sdir/row['path']).read_bytes()).hexdigest()==row['sha256']
(m/'EXTENSIONLESS_RECORDS_VALIDATION.json').write_text(json.dumps({'status':'PASS','files':len(SUPP_ROWS),'manifest':SUPP_ROWS},indent=2))
v=json.loads((w/'qc/G_VALIDATION.json').read_text());iv=json.loads((w/'qc/I_aggregate_validation.json').read_text());delivery=json.loads((w/'qc/DELIVERY_VALIDATION.json').read_text())
assert v['matched_every_selected_cell'] and delivery['all_sha256_match']
review={'status':'PASS_VISUAL_REVIEW','G_composite':'all four panels inspected; framed compact layout, common scale, unclipped labels and true continuous distributions','G_R‌CC4_single':'inspected','I':'inspected; three marker panels and patient-median points readable','user_review':'PENDING','pdf':'matplotlib exports from same figure objects; separate PDF parser unavailable in local environment'}
(w/'qc/VISUAL_REVIEW.json').write_text(json.dumps(review,ensure_ascii=False,indent=2))
sources=[]
for pt in ['GSE207493__RCC104','GSE207493__RCC86','GSE156632__RCC4','GSE156632__RCC5']:
 p=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/infercnv_v2/patients')/pt
 for f in [p/'cell_id_mapping_v2.tsv',p/'infercnv_run/preliminary.infercnv_obj']:
  h=hashlib.sha256()
  with f.open('rb') as rd:
   for chunk in iter(lambda:rd.read(1048576),b''):h.update(chunk)
  sources.append({'path':str(f),'bytes':f.stat().st_size,'sha256':h.hexdigest()})
(w/'qc/SOURCE_MANIFEST.json').write_text(json.dumps(sources,indent=2))
note=NOTE
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=root/name;old=f.read_text(encoding='utf8')
 if '## Phase122 G/I reference redraw complete' not in old:f.write_text(note+'\\n\\n'+old,encoding='utf8')
rows=json.loads((w/'qc/FIGURE_MANIFEST.json').read_text());mf=root/'RESULTS_MANIFEST.tsv';existing=mf.read_text()
with mf.open('a') as out:
 for i,row in enumerate(rows):
  key='phase122_review_figure_'+str(i+1)
  if key+'\\t' not in existing:out.write(key+'\\tphase122_G_I_redraw\\treview_deliverable\\t'+str(w/'figures'/row['file'])+'\\t'+datetime.datetime.now().isoformat()+'\\tPASS_REVIEW_CHECKS\\tSHA256 '+row['sha256']+'\\tSource labels frozen; new continuous descriptive CNV score explicitly disclosed\\n')
(w/'qc/COMPLETION.json').write_text(json.dumps({'status':'PASS_COMPLETE_PENDING_USER_VISUAL_FEEDBACK','migration':str(m),'migration_files':mv['verified_files'],'supplement_status_records':len(SUPP_ROWS),'delivered_files':len(rows),'authority':str(root),'future_local_records':False},indent=2))
print((w/'qc/COMPLETION.json').read_text())
'''
note='''## Phase122 G/I reference redraw complete (2026-09-08)

- 当前唯一权威记录位置：服务器 `/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1`。用户明确要求全部本地记录迁到服务器，后续不在本地记录。旧本地文件为静态历史副本；本地 Downloads/26210 只用于图件交付。
- 已恢复此前漏登的要求：G/I 按第三页小提琴和点图重绘；用户最新反馈 G 不满意，并明确选择“按参考图整体重做”。Phase121 旧 G v2 被本次 G v3 替代。
- G v3：四位原有患者、6,034 个合格细胞、正式亚克隆标签不变；紧凑四分面、统一纵轴、细轮廓、柔和实色、实际观测点和中位数/IQR。每患者 C01 等标签仅在该患者内有意义，不表示跨患者同一克隆。
- 指标更换已公开说明：旧图是 HMM 组级 altered-gene burden 相对阈值；新版为 preliminary inferCNV 相对信号的逐细胞 mean(abs(log2(signal)))。新版提供真实连续分布，不与旧数值直接比较，不改变恶性判定或克隆，不增加推断性检验。不同患者有独立处理/参考，统一轴用于展示，不据此作跨患者恶性程度推断。
- I v3：NDUFA4L2/CA9/MET 三基因，保留原始 CNV 分组，细胞表达分布加患者中位数散点；未新增 P/FDR。
- 正式服务器目录：`phase122_G_reference_redraw_remote_v1`。绘图脚本、逐细胞中间表、来源哈希、审核和交付校验均在服务器。
- 已交付 `C:/Users/31937/Downloads/26210/02_2.1_G_I参考图重绘_v3`，共 6 PNG、6 PDF、1 阅读说明，13/13 哈希一致。G 总图、G 单图和 I 视觉检查完成。当前断点：用户审阅 G v3；不自动推进下一节。
'''
code=code.replace('ROOT',repr(root)).replace('W)',repr(w)+')').replace('SUPP_BYTES',repr(base64.b64encode(b.getvalue()).decode())).replace('SUPP_ROWS',repr(supp)).replace('NOTE',repr(note)).replace('HELPERS',repr(helpers))
call(code)
