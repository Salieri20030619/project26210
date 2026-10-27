from pathlib import Path
import shutil,json,hashlib
r=Path.cwd();d=r/'phase204_concentric_comparison_20260918_v8';out=d/'delivery'
for p in (d/'figures').iterdir():shutil.copy2(p,out/p.name)
text=(r/'phase203_compact_comparison_matrix_20260918_v7/delivery/图注与使用说明.md').read_text()
start=text.index('上两行使用完整发现队列：');end=text.index('图件198.1')
text='# SPP1同心双环比较图 v8\n\n五组横排，依次为C4、C10、C15、三群联合和整体成熟巨噬细胞。中心填写SPP1伪批量表达中位数（log2[CPM+1]）；浅色内环弧长表示SPP1阳性细胞比例中位数；深色外环弧长表示SPP1原始RNA贡献比例中位数。各环统一完整圆为100%，从12点顺时针增加，浅灰表示剩余比例。每组下方依次给出内环/外环的百分比。数据与v6、v7完全一致。\n\n'+text[start:end].replace('上两行','中心表达及内环').replace('下行','外环').replace('上两行样本量','表达和阳性率样本量')+'\n图件208.3 × 78.7 mm，适合作为主图中的横向面板；详细统计口径可纳入主图图注。PDF/SVG为矢量，PNG为600 dpi。本次按要求合并三行图形，尚未重排整张Fig3或报告。\n'
(out/'图注与使用说明.md').write_text(text)
files=list(out.iterdir());(out/'FILE_MANIFEST.json').write_text(json.dumps([dict(path=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files],ensure_ascii=False,indent=2))
print('PACKAGE_READY')
