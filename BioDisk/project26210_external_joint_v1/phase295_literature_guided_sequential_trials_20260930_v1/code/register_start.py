from pathlib import Path
import json,hashlib,datetime,shutil,os,csv,io
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
now=datetime.datetime.now().astimezone().isoformat()
status=dict(status='PLAN_FROZEN_STEP1_PREFLIGHT',at=now,plan_sha256=sha(P/'CURRENT_STATE_AND_SEQUENTIAL_PLAN.md'),step1='PREFLIGHT',step2='NOT_STARTED',step3='NOT_STARTED',new_statistics=False)
(P/'qc/STATUS.json').write_text(json.dumps(status,indent=2))
entry="""## Phase295（2026-09-30）：当前状态已记录，按文献方法顺序尝试
用户授权先记录、制定计划并依次尝试。已冻结phase295_literature_guided_sequential_trials_20260930_v1/CURRENT_STATE_AND_SEQUENTIAL_PLAN.md，保留Phase280–294全部来源定义、正负结果与访问限制。
顺序：1接收端靶基因驱动的NicheNet候选网络；2免疫冷/浸润/排斥表型原方法资格；3独立纤维化状态/发表程序参考对齐。每步终态后再进入下一步。原始权重/网络/参考不到位不得自造替代评分；细胞不补患者n，不要求阳性，不等待客户新实验。
当前第一步预检开始，第二/三步未启动。无新统计/图，旧C4/H5AD/边界/文章不变，记录仅服务器。

"""
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
 t=R/name;old=t.read_bytes()
 with (P/'root_record_backups'/name).open('xb') as f:f.write(old)
 if name.endswith('.md'):new=entry.encode()+old
 else:
  fields=old.decode('utf-8-sig').splitlines()[0].split('\t');out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
  w.writerow(dict(result_id='phase295_sequential_plan',stage='Phase295',artifact_type='protocol',absolute_path=str(P/'CURRENT_STATE_AND_SEQUENTIAL_PLAN.md'),generated_at=now,validation_status='FROZEN',validation_evidence=str(P/'qc/STATUS.json'),notes='user authorized sequential attempts;SHA256='+status['plan_sha256']))
  new=old+(b'' if old.endswith(b'\n') else b'\n')+out.getvalue().encode()
 assert t.read_bytes()==old
 tmp=R/(name+'.phase295start.tmp')
 with tmp.open('xb') as f:f.write(new)
 shutil.copymode(t,tmp);os.replace(tmp,t)
print(json.dumps(status))
