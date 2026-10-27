from pathlib import Path
import shutil,json,hashlib,datetime
r=Path.cwd();d=r/'phase203_compact_comparison_matrix_20260918_v7';out=d/'delivery'
for p in (d/'figures').iterdir():shutil.copy2(p,out/p.name)
(out/'图注与使用说明.md').write_text('''# SPP1候选群比较：主图面板版

列依次为C4、C10、C15、三群联合及整体巨噬细胞。上层色块表示患者级SPP1伪批量表达中位数，统一色标8–12 log2(CPM+1)；中层气泡面积表示患者级SPP1阳性细胞比例中位数；下层圆环弧长表示患者级SPP1原始RNA贡献比例中位数，灰环为剩余部分。各项数值与v6条形图完全一致。

上两行使用完整发现队列：候选群/联合要求每例目标与其他成熟巨噬细胞各至少20个，整体要求至少20个巨噬细胞，患者数为26、2、9、28、38。下行使用原16例同文库匹配队列，以该患者全部成熟巨噬细胞SPP1原始RNA为分母，缺失群按0计入。C10的0%是16例患者贡献的中位数，不代表每例均为0。

图中为描述性患者中位数，未进行新的推断检验。各组患者构成不同且存在包含关系，不能把数值之差当作配对效应，不能相加。C10仅2例有效表达患者，C15为9例；RNA贡献行的16例不得用于上两行样本量。

图件198.1 × 102.9 mm，推荐约160–180 mm宽作为主图单个面板，详细口径放入主图图注。PDF/SVG为矢量，PNG为600 dpi。本次独立交付新版比较面板，尚未替换整张Fig3或报告，以便先审阅新的图形样式。C4主分析口径不变。
''')
files=[p for p in out.iterdir() if p.is_file()];(out/'FILE_MANIFEST.json').write_text(json.dumps([dict(path=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in files],ensure_ascii=False,indent=2))
print('PACKAGE_READY')
