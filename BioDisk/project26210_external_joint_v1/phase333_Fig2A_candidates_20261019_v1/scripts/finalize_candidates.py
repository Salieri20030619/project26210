from pathlib import Path
from datetime import datetime
import json,hashlib,csv,io
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase333_Fig2A_candidates_20261019_v1'
v=json.loads((d/'qc/VALIDATION.json').read_text(encoding='utf-8'))
v.update(status='PASS',visual_review='Both final PDF renders reviewed; horizontal labels legible, no clipping or label-strip overlap.',user_approval='PENDING',integration='NOT_PERFORMED')
(d/'qc/VALIDATION.json').write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8')
caption='''26210 Fig2A 候选单图（独立预览，待用户选择）

候选①：上皮细胞恶性身份
同一冻结 UMAP 中展示全部38,670个肿瘤组织上皮细胞；每点代表一个细胞，颜色为原CNV分类。双方法共同支持恶性10,844个、仅inferCNV支持恶性1,106个、共同支持非恶性498个、方法不一致22,216个、未定4,006个。前两类构成既有11,950个恶性细胞分析集。比例分母为38,670，正常组织上皮细胞不在本图分类范围内。原坐标、分类不变，散点绘制顺序不能用于推断类别丰度。

候选③：cNMF程序代表基因
列为F1–F7，行为既有固定35个代表基因，左右分块连续展示F1–F4和F5–F7各自的代表基因。每个程序原先按载荷降序、排除此前程序已选基因后取5个；未因本次绘图重新挑选。颜色为每个基因跨7个程序的载荷除以该基因最大载荷，范围0–1；右侧色条为原选择程序。载荷不是差异表达效应，也不是P值或FDR；本图无新增推断检验。

两图均190×116 mm，提供PDF、300 dpi PNG及可编辑SVG。数据数值与PDF渲染检查通过。当前为独立候选图，未经用户选定，不替换Fig2A、不重拼Fig2、不改完整稿。Fig2其他标签问题待选图后处理。
'''
(d/'delivery/候选单图说明.txt').write_text(caption,encoding='utf-8')
files=[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted((d/'delivery').iterdir()) if p.is_file()]
(d/'qc/DELIVERY_MANIFEST.json').write_text(json.dumps(files,ensure_ascii=False,indent=2),encoding='utf-8')
note='''## Phase333（2026-10-19）：Fig2A候选①及③独立预览完成，待用户选择
用户指出Fig1A与Fig2A重复、Fig2多处标签倾斜。已核实二者均源于细胞图谱。用户明确先报告候选、生成独立单图，经其通过后再加入；本轮指定生成候选①和③。
候选①使用原CNV分类和冻结UMAP，保留38670肿瘤组织上皮细胞五类别；候选③使用Phase123固定35基因×7程序相对载荷，数值与原载荷重建最大差9.75e-17。两图190×116mm，PDF/PNG/SVG及中文说明已交付Downloads/26210/Phase333_Fig2A候选单图。脚本、源数据及验收在phase333_Fig2A_candidates_20261019_v1；qc/VALIDATION.json及DELIVERY_MANIFEST.json，数据及最终PDF视觉检查PASS。
当前等待用户选图/修改意见；未替换Fig2A、未重拼Fig2、未改完整稿。正式稿仍Phase332 v5；其他标签修订在选图后继续。不因候选出图自动视作用户通过。

'''
backup=d/'root_record_backups';backup.mkdir(exist_ok=True)
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
    p=r/name
    if not (backup/name).exists():(backup/name).write_bytes(p.read_bytes())
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
    p=r/name;s=p.read_text(encoding='utf-8')
    if '## Phase333（2026-10-19）' not in s:p.write_text(note+s,encoding='utf-8')
p=r/'RESULTS_MANIFEST.tsv';s=p.read_text(encoding='utf-8')
if 'phase333_candidate1' not in s:
    buf=io.StringIO();w=csv.writer(buf,delimiter='\t',lineterminator='\n')
    for num,name in [(1,'Candidate_1_epithelial_malignancy'),(3,'Candidate_3_cNMF_representative_genes')]:
        w.writerow([f'phase333_candidate{num}','Phase333','figure',str(d/'delivery'/f'{name}.pdf'),datetime.now().isoformat(),'PASS',str(d/'qc/VALIDATION.json'),'Standalone preview; user selection pending; manuscript unchanged'])
    p.write_text(s+('' if s.endswith('\n') else '\n')+buf.getvalue(),encoding='utf-8')
print('COMPLETE: two previews; data and visual checks PASS; awaiting user selection')
