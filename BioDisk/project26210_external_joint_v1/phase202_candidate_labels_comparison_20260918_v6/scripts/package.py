from pathlib import Path
import subprocess,shutil,json,hashlib,re
r=Path.cwd();d=r/'phase202_candidate_labels_comparison_20260918_v6';rep=r/'phase201_three_candidates_C4_main_report_20260918_v4';out=d/'delivery'
for p in (d/'figures').glob('Main_Figure_*.pdf'):
 subprocess.run(['pdftocairo','-svg',str(p),str(p.with_suffix('.svg'))],check=True,stderr=subprocess.PIPE)
 subprocess.run(['pdftoppm','-r','600','-singlefile','-png',str(p),str(p.with_suffix(''))],check=True,stderr=subprocess.PIPE)
for p in (d/'figures').iterdir():
 if p.suffix in ['.pdf','.svg','.png']:shutil.copy2(p,out/p.name)
p=rep/'report/26210_TYMP_MTA_SPP1_macrophage_C4_primary_report_v4.pdf';shutil.copy2(p,out/p.name)
md=p.with_suffix('.md').read_text();md=re.sub(r'\[既有图件\]\(\.\./figures/[^)]+\)','（见 PDF 正文）',md);(out/p.with_suffix('.md').name).write_text(md)
(out/'补充记录').mkdir(exist_ok=True)
for p in (rep/'supplement').glob('*.md'):shutil.copy2(p,out/'补充记录'/p.name)
(out/'说明.md').write_text('''# 26210 本次更新

- Fig2D 与 Fig3B：保留 C4 紫色标记，新增 C10 蓝色、C15 绿色位置标记。轮廓来自原冻结群的细胞密度，不表示表达峰值或专属表达。
- 独立比较图：C4、C10、C15、三群联合、整体巨噬细胞，比较 SPP1 表达、阳性率与 RNA 来源贡献。
- 比较图 A/B：完整发现队列，患者等权中位数，有效患者数 26/2/9/28/38。群或联合要求目标与其余各至少20个细胞；整体要求至少20个巨噬细胞。
- 比较图 C：原16例同文库匹配患者，RNA贡献统一以各例全部巨噬细胞SPP1计数为分母，缺失群记0。C10中位贡献0不表示每位患者均无贡献。
- C10只有2例满足表达汇总条件，因此只作描述，不进行推断检验；各群患者构成不相同，柱高不构成配对效应量。群之间存在包含关系，中位数不可相加。
- 报告 v4：原报告基础新增2.9节及图8、图9，第20–21页；后续主分析继续以C4为对象。34页、31张原表不变、7张原图保持，新增2张说明图。
- 主图 Fig2/Fig3 属于单独投稿组图版本，与原报告图号是两套编号；Fig1沿用暂定稿。
- 提供三张图的矢量PDF/SVG及600 dpi PNG。报告MD保留正文和表格，完整图件见PDF。
- 脚本、患者级数据与审计记录均保留在sgs10的phase201、phase202目录；旧交付物未覆盖。
''')
checks=json.loads((d/'qc/FINAL_CHECKS.json').read_text());checks['visual_review']='All 34 report pages, both new report figure pages, main Fig2/Fig3 and comparison inspected; no clipping or overlapping labels.';(d/'qc/FINAL_CHECKS.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
files=[p for p in out.rglob('*') if p.is_file()];manifest=[{'path':str(p.relative_to(out)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files];(out/'FILE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print('PACKAGE_READY',len(files),sum(x['bytes'] for x in manifest),flush=True)
