from pathlib import Path
import json,shutil,hashlib
P=Path(__file__).resolve().parents[1];R=P.parent;O=R/'phase331_manuscript_immune_figure_integration_20261006_v4'
b=json.loads((O/'report/report_materialized.json').read_text());old=json.loads(json.dumps(b));changes=[]
b[3]['text']=b[3]['text'].replace('肾癌细胞的代谢状态与巨噬细胞异质性及组织分布之间的联系尚需在患者层面系统解析。','肾癌TYMP相关代谢状态与SPP1巨噬细胞及免疫屏障之间的联系，是理解肿瘤免疫微环境的重要问题。')
b[8]['text']=b[8]['text'].replace('提供功能对应。','提供功能对应。SPP1相关状态还与成纤维细胞及基质程序正相关；GSE175540中19/20名合格患者的肿瘤表达区内侧CD8信号较低。公开CRC模型中，Spp1缺失后CAF减少、CD8数量和增殖升高。')
b[9]['text']=b[9]['text'].replace('空间分析进一步显示，SPP1相关状态与成纤维细胞和细胞外基质程序相联系。','这些发现从代谢、巨噬细胞免疫状态和基质环境三个层面，为研究TYMP参与SPP1相关免疫屏障提供了机制线索。')
for i,(a,c) in enumerate(zip(old,b)):
 if a!=c:changes.append({'block':i,'before':a['text'],'after':c['text']})
assert [x['block'] for x in changes]==[3,8,9]
assert b[10:]==old[10:]
(P/'report/report_materialized.json').write_text(json.dumps(b,ensure_ascii=False,indent=2))
for f in (O/'figures').iterdir():shutil.copy2(f,P/'figures'/f.name)
s=(O/'code/render_report.py').read_text().replace(O.name,P.name).replace('26210_manuscript_immune_figure_integrated_v4','26210_manuscript_immune_figure_integrated_v5');(P/'code/render_report.py').write_text(s)
(P/'qc/ABSTRACT_CHANGES.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2))
(P/'delivery/摘要修改说明.txt').write_text('仅微调摘要三个位置，正文、统计表及图件保持v4。\n\n'+'\n\n'.join('【修改后】\n'+x['after'] for x in changes),encoding='utf-8-sig')
print('ABSTRACT_REVISED',[(x['block'],len(x['before']),len(x['after'])) for x in changes])
