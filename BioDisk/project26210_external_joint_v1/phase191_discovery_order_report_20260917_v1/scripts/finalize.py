from pathlib import Path
import json,shutil,hashlib,datetime,re
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase191_discovery_order_report_20260917_v1';old=r/'phase190_mechanism_focused_report_20260917_v1'
b=json.loads((d/'report/report_materialized.json').read_text());cv=json.loads((d/'qc/CONTENT_VALIDATION.json').read_text())
assert cv['tables_preserved']==31 and cv['figures_preserved']==7
titles=[x['text'] for x in b if x['type']=='h2']
assert titles.index(next(t for t in titles if t.startswith('2.9 无监督')))<titles.index(next(t for t in titles if t.startswith('2.10 围绕')))
cmp=json.loads((d/'qc/PAGE_COMPARISON.json').read_text())
(d/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'pages':33,'tables':31,'figures':7,'visual_review_changed_pages':cmp['changed'],'unchanged_pages_reuse':cmp['unchanged'],'all_pass':True,'no_new_statistics':True,'sequence_review_pass':True},indent=2))
delivery=d/'delivery';delivery.mkdir(exist_ok=True);name='26210_TYMP_MTA_SPP1_主线顺序修订版_v3';stem='26210_TYMP_MTA_SPP1_macrophage_integrated_report_v3'
pdf=d/'report'/f'{stem}.pdf';shutil.copy2(pdf,delivery/(name+'.pdf'))
txt=(d/'report'/f'{stem}.md').read_text()
for pl in json.loads((d/'qc/PLACEMENTS.json').read_text()):
 txt=re.sub(r'图'+str(pl['figure'])+r' \[既有图件\]\([^\n]+\)',f"图{pl['figure']}（图件见同名PDF第{pl['page']}页）",txt)
(delivery/(name+'_正文与表格.md')).write_text(txt)
shutil.copy2(old/'supplement/MECHANISM_EXPLORATION_SUPPLEMENT_CN.md',d/'supplement/MECHANISM_EXPLORATION_SUPPLEMENT_CN.md')
shutil.copy2(d/'supplement/MECHANISM_EXPLORATION_SUPPLEMENT_CN.md',delivery/'26210_机制探索补充材料_沿用v2.md')
now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
note=f"""## Phase191 代谢铺垫与SPP1先发现后关联的正文顺序修订（{now}）
用户指出2.8不能出现第一部分等编排用语，MTA须先铺垫，2.9/2.10应先发现SPP1群再研究联系。已落实：
2.8从Cys/Met、Arg/Pro、ODC1/SRM观察出发，经SAM和SRM/SMS反应引出完整名称MTA，再用配对实测检验；
2.9独立无监督髓系鉴定发现SPP1富集C4；2.10围绕该发现检验MTA-SPP1组织联系及直接髓系依据。
摘要、1.1研究设计、3.4工作模型、3.7结论同步；表8/9交换，所有原表头数据保持，31表/7图、33页，无新统计/图。正文已无第一部分字样。补充探索沿用v2。
PDF：{pdf}
主要修订第19-22页；变化页面全视觉检查，未变页面与已审版像素一致；无排版问题。旧v1/v2保留，交付PDF及可编辑MD。无活跃任务。
"""
backup=d/'project_record_backups';backup.mkdir(exist_ok=True)
for n in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']:
 p=r/n;shutil.copy2(p,backup/n)
 if n=='RESULTS_MANIFEST.tsv':
  with p.open('a')as f:f.write('\n'+'\t'.join(['phase191_discovery_order_v3','phase191','pdf',str(pdf),now,'PASS','sequence and table reference checks; 33 pages visual coverage','metabolic lead-in; SPP1 discovery before association'])+'\n')
 else:p.write_text(note+'\n'+p.read_text())
(d/'STATUS.json').write_text(json.dumps({'time':now,'complete':True,'pages':33,'active_jobs':[]},indent=2))
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in delivery.iterdir()}
(delivery/'SHA256.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2));print(json.dumps(hashes,ensure_ascii=False))
