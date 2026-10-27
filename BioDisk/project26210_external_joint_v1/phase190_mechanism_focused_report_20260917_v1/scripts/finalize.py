from pathlib import Path
import json,shutil,hashlib,datetime,re
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');old=r/'phase189_integrated_TYMP_MTA_macrophage_report_20260917_v1';d=r/'phase190_mechanism_focused_report_20260917_v1'
a=json.loads((old/'report/report_materialized.json').read_text());b=json.loads((d/'report/report_materialized.json').read_text())
assert [x for x in a if x['type']=='table']==[x for x in b if x['type']=='table']
same=[]
for i in range(5,17):
 assert (old/'qc/pages_final'/f'page-{i:02}.png').read_bytes()==(d/'qc/pages_final'/f'page-{i:02}.png').read_bytes(),i
 same.append(i)
build=json.loads((d/'qc/BUILD.json').read_text());assert build['pages']==33 and build['tables']==31 and build['figures']==7
(d/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'unchanged_tables':31,'figures':7,'pages':33,'visually_reviewed_pages':[1,2,3,4]+list(range(17,34)),'pixel_identical_previously_reviewed_pages':same,'no_new_statistics':True,'main_boundaries_concentrated':'Discussion3.5','unrelated_exploration':'separate supplement','layout_pass':True},indent=2))
delivery=d/'delivery';delivery.mkdir(exist_ok=True)
name='26210_TYMP_MTA_SPP1_机制讨论修订版_v2';stem='26210_TYMP_MTA_SPP1_macrophage_integrated_report_v2'
pdf=d/'report'/f'{stem}.pdf';shutil.copy2(pdf,delivery/(name+'.pdf'))
txt=(d/'report'/f'{stem}.md').read_text()
for pl in json.loads((d/'qc/PLACEMENTS.json').read_text()):
 txt=re.sub(r'图'+str(pl['figure'])+r' \[既有图件\]\([^\n]+\)',f"图{pl['figure']}（图件见同名PDF第{pl['page']}页）",txt)
(delivery/(name+'_正文与表格.md')).write_text(txt)
shutil.copy2(d/'supplement/MECHANISM_EXPLORATION_SUPPLEMENT_CN.md',delivery/'26210_机制探索补充材料_v2.md')
now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
note=f"""## Phase190 正文机制讨论与补充材料分工修订完成（{now}）
按用户本次要求，未证实且无关主线的分支从正文移除；直接相关边界集中Discussion3.5；探索过程另列MECHANISM_EXPLORATION_SUPPLEMENT_CN.md。重写Discussion3.1-3.4，围绕TYMP核苷及蛋氨酸/多胺背景、MTA生成清除和局部暴露、SPP1脂质/溶酶体状态及工作模型。
新版PDF：{pdf}
33页、7原矢量图、31表；31表结构及内容逐字保持，无新统计/新图，主要阴性保留。正文不再列IL10/LGALS3/IL6探索失败，完整材料保存于单独补充及旧分析档案。主要阅读变化为19-24页。
验证：变化页视觉检查，5-16页像素与已审版完全一致，其余页均检查；无排版问题。旧v1保留。正式交付PDF、可编辑MD与探索补充MD。无活跃作业。
"""
backup=d/'project_record_backups';backup.mkdir(exist_ok=True)
for n in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']:
 p=r/n;shutil.copy2(p,backup/n)
 if n=='RESULTS_MANIFEST.tsv':
  with p.open('a') as f:f.write('\n'+'\t'.join(['phase190_mechanism_report_v2','phase190','pdf',str(pdf),now,'PASS','31 tables identical; 33 pages visual coverage','mainline discussion; exploration separate'])+'\n')
 else:p.write_text(note+'\n'+p.read_text())
(d/'STATUS.json').write_text(json.dumps({'time':now,'complete':True,'pages':33,'tables':31,'figures':7,'active_jobs':[]},indent=2))
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in delivery.iterdir()}
(delivery/'SHA256.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2))
print(json.dumps(hashes,ensure_ascii=False))
