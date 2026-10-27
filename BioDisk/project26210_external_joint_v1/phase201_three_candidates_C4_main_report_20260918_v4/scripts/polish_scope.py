from pathlib import Path
import json
d=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase201_three_candidates_C4_main_report_20260918_v4')
p=d/'report/report_materialized.json';b=json.loads(p.read_text())
for x in b:
 if 'text' in x:
  x['text']=x['text'].replace('这一组合将本项目的SPP1信号定位于脂质处理和细胞内降解相关的巨噬状态。','这一组合将C4中的SPP1富集定位于脂质处理和细胞内降解相关的巨噬状态。')
  x['text']=x['text'].replace('单细胞分析独立发现SPP1富集C4后，','单细胞候选比较确定C4为主分析状态后，')
p.write_text(json.dumps(b,ensure_ascii=False,indent=2))
