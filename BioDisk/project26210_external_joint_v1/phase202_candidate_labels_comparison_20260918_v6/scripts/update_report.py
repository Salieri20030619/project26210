from pathlib import Path
import json,shutil
r=Path.cwd();d=r/'phase201_three_candidates_C4_main_report_20260918_v4';n=r/'phase202_candidate_labels_comparison_20260918_v6'
p=d/'report/report_materialized.json';b=json.loads(p.read_text())
idx=next(i for i,v in enumerate(b) if v['type']=='h2' and v.get('text','').startswith('2.10 '))
add=[{'type':'figure'},{'type':'caption','text':'图9 三个候选群、三群联合与整体巨噬细胞的 SPP1 表达和来源贡献','keep_with_next':True},{'type':'p','text':'A–B：完整发现队列中，按患者汇总 SPP1 伪批量表达（log2[CPM+1]）及阳性细胞比例，展示患者中位数；候选群和联合群要求每例至少20个目标细胞及20个其他巨噬细胞，整体要求每例至少20个巨噬细胞。有效患者数分别为 C4 26例、C10 2例、C15 9例、联合群28例和整体38例。C：原16例同文库匹配队列中各群 SPP1 原始 RNA 计数占整体巨噬细胞 SPP1 计数的比例，缺失该群的患者按0计入。各组存在嵌套，患者中位数不能相加；图中均为描述性比较，不作组间推断检验。C10 的较高表达仅由2例有效患者支持，不能替代跨患者覆盖和多基因状态证据。'}]
b[idx:idx]=add
for v in b:
 if v.get('type')=='p' and '59.63%' in v.get('text','') and '图9' not in v['text']:v['text']+='三个候选群、联合群及整体巨噬细胞的表达与来源贡献比较见图9。'
p.write_text(json.dumps(b,ensure_ascii=False,indent=2))
shutil.copy2(n/'figures/Candidate_states_vs_whole_macrophages_v6.pdf',d/'figures/Figure9_candidate_states_vs_whole.pdf')
f=d/'scripts/render_report.py';s=f.read_text().replace("'Figure8_three_SPP1_candidates.pdf']","'Figure8_three_SPP1_candidates.pdf','Figure9_candidate_states_vs_whole.pdf']").replace('if figi==8:','if figi in [8,9]:');f.write_text(s)
f=n/'scripts/assemble3.py';s=f.read_text().replace("functional_candidate_locations.pdf',[0,0,1,1]","functional_candidate_locations.pdf',[0,.065,1,1]");f.write_text(s)
