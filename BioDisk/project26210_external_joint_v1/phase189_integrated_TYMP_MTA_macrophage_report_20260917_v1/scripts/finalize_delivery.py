from pathlib import Path
import json,hashlib,shutil,datetime,re,platform
from pypdf import PdfReader
import reportlab
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase189_integrated_TYMP_MTA_macrophage_report_20260917_v1'
stem='26210_TYMP_MTA_SPP1_macrophage_integrated_report_v1'
pdf=d/'report'/f'{stem}.pdf';md=d/'report'/f'{stem}.md'
original=json.loads((r/'phase143_sections25_27_candidate_report_remote_v1/report/report_materialized.json').read_text())
new=json.loads((d/'report/report_materialized.json').read_text())
ot=[x for x in original if x['type']=='table'];nt=[x for x in new if x['type']=='table']
assert len(nt)==31 and all(any(t['headers']==u['headers'] and t['rows']==u['rows'] for u in nt) for t in ot)
assert len(PdfReader(pdf).pages)==33
pl=json.loads((d/'qc/PLACEMENTS.json').read_text());assert len(pl)==7 and [p['figure'] for p in pl]==list(range(1,8))
vis={'pages_rendered':33,'visual_coverage':'all pages reviewed; unchanged pages matched earlier inspected raster; remaining final pages inspected','unresolved_clipping_or_overlap':False,'pagination':'short tables kept together; long appendices split with repeated headers','new_figures':0,'render_builder':'render_report_v3.py'}
(d/'qc/VISUAL_REVIEW.json').write_text(json.dumps(vis,ensure_ascii=False,indent=2))
manifest=json.loads((d/'qc/SOURCE_MANIFEST.json').read_text())
manifest+=json.loads((r/'phase188_mainline_reorganization_20260917_v1/SOURCE_MANIFEST.json').read_text())
for p in (d/'figures').glob('*.pdf'):manifest.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(d/'qc/SOURCE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
delivery=d/'delivery';delivery.mkdir(exist_ok=True)
outname='26210_TYMP_MTA_SPP1_巨噬状态整合报告_v1'
shutil.copy2(pdf,delivery/(outname+'.pdf'))
portable=md.read_text()
for p in pl:
 portable=re.sub(r'图'+str(p['figure'])+r' \[既有图件\]\([^\n]+\)',f"图{p['figure']}（图件见同名PDF第{p['page']}页）",portable)
(delivery/(outname+'_正文与表格.md')).write_text(portable,encoding='utf-8')
now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
status={'time':now,'pages':33,'original_figures_preserved':7,'original_tables_preserved':25,'total_tables':31,'new_statistics':False,'new_figures':False,'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'python':platform.python_version(),'reportlab':reportlab.Version,'scope':'first report expanded through macrophage state localization','active_jobs':[]}
(d/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
(d/'CHANGE_SUMMARY_CN.md').write_text('第一部分正文为基础统一扩展。新增方法1.9–1.11、结果2.8–2.12、讨论3.2–3.3及阶段结论、附表S11与参考文献。原表7实验对应重编号为表12并移到讨论3.6，原25张表的表头及数据逐项保持。保留7幅原矢量图，无新增图或统计。新增结果第19–21页；全报告33页。')
backup=d/'project_record_backups';backup.mkdir(exist_ok=True)
note=f"""## Phase189 第一部分扩展至巨噬状态定位的整合报告完成（{now}）
用户要求在第一部分报告基础上展开至定位巨噬状态。正式PDF：{pdf}
统一33页，原7矢量图、25原统计表均保留，新增6表（总31）。新增方法1.9-1.11、结果2.8-2.12：MTA候选生化及实测分段、MTA—组织SPP1名义联系、C4多基因身份、脂质/溶酶体程序与SPP1贡献、髓系状态排序。
新增结果位于19-21页。实验对应移至统一讨论3.6，原表7重编号12；主要阴性及不同患者分母保留，未确认完整MTA诱导机制。STAT3等仅讨论辅助，不展开IL6旁支。免疫空间只交代后续问题，未新增下游模型。
最终render_report_v3.py；33页视觉覆盖完成，修正附表留白/小标题分离及短表孤行。25表数据保持校验、7图和关键数字检查通过，无新图/新统计。原PDF只读，源正文/字体/矢量图/脚本/校验在本目录。
交付：delivery内PDF及可编辑正文MD，后续先据本稿审核表达；原免疫/空间资格任务仍为后续，不在本轮执行。无活跃任务。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']:
 p=r/n;shutil.copy2(p,backup/n)
 if n=='RESULTS_MANIFEST.tsv':
  with p.open('a')as f:f.write('\n'+'\t'.join(['phase189_integrated_report','phase189','pdf',str(pdf),now,'PASS','33 pages reviewed; 25 original tables preserved; 7 vector figures','scope ends macrophage state; no new analyses'])+'\n')
 else:p.write_text(note+'\n'+p.read_text(),encoding='utf-8')
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in delivery.iterdir()}
(delivery/'SHA256.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2))
print(json.dumps(status,ensure_ascii=False));print(json.dumps(hashes,ensure_ascii=False))
