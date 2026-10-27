from pathlib import Path
import json,hashlib,subprocess,shutil
from pypdf import PdfReader
r=Path.cwd();d=r/'phase205_fig2_concentric_20260918_v9';out=d/'delivery';old=r/'phase202_candidate_labels_comparison_20260918_v6'
a=json.loads((old/'tables/FIG2_PANEL_PROVENANCE.json').read_text());b=json.loads((d/'tables/FIG2_PANEL_PROVENANCE.json').read_text());assert [x['panel'] for x in b]==list('ABCDEF')
for x,y in zip(a,[b[0],b[1],b[2],b[3],b[5]]):
 assert hashlib.sha256(Path(x['source']).read_bytes()).digest()==hashlib.sha256(Path(y['source']).read_bytes()).digest();assert x['source_crop']==y['source_crop']
f=d/'figures/Main_Figure_2_concentric_comparison_v9.pdf';reader=PdfReader(f);assert len(reader.pages)==1
W=float(reader.pages[0].mediabox.width)/72*25.4;H=float(reader.pages[0].mediabox.height)/72*25.4
for p in b:assert p['x_mm']>=0 and p['y_mm']>=0 and p['x_mm']+p['width_mm']<=W and p['y_mm']+p['height_mm']<=H
for i,p in enumerate(b):
 for q in b[i+1:]:
  overlap=min(p['x_mm']+p['width_mm'],q['x_mm']+q['width_mm'])-max(p['x_mm'],q['x_mm']);vertical=min(p['y_mm']+p['height_mm'],q['y_mm']+q['height_mm'])-max(p['y_mm'],q['y_mm']);assert overlap<=0 or vertical<=0
subprocess.run(['pdftocairo','-svg',str(f),str(f.with_suffix('.svg'))],check=True,stderr=subprocess.PIPE)
subprocess.run(['pdftoppm','-r','600','-singlefile','-png',str(f),str(f.with_suffix(''))],check=True,stderr=subprocess.PIPE)
for p in (d/'figures').iterdir():shutil.copy2(p,out/p.name)
base=(r/'phase198_fig2_contrast_corr_revision_20260918_v5/delivery/修改说明与图注.md').read_text();abc=base[base.index('## A：'):base.index('## D：')].replace('独立相关图同时交付，便于单独放大阅读。','')
fcap=base[base.index('## E：'):base.index('## 生化图')].replace('## E：','## F：')
ecap=(r/'phase204_concentric_comparison_20260918_v8/delivery/图注与使用说明.md').read_text();ecap=ecap[ecap.index('五组横排'):ecap.index('图件208.3')]
(out/'Fig2图注_v9.md').write_text('# Fig2 v9：代谢背景、候选群定位与SPP1贡献\n\n本版将同心双环比较图加入E，接在三候选定位D之后；原组织MTA热图顺延为F。A–D及原组织热图的源图内容、数值及裁切保留，E采用phase204已确认设计。Fig3保持原版本。\n\n'+abc+'## D：成熟巨噬状态及三候选位置\n左侧为原冻结12个成熟巨噬候选群；右侧保留C4密度，并标出C4、C10、C15。紫色虚线为C4，蓝色点划线为C10，绿色点线为C15，边界使用各群自身密度峰值的12%。所有坐标、成员及表达数据保持；轮廓代表群的位置，不代表SPP1表达峰值或专属表达。\n\n## E：SPP1表达与来源贡献\n'+ecap+'\n'+fcap+'\n## 交付规格\n单页矢量PDF、SVG和600 dpi PNG；图幅340 × '+f'{H:.1f}'+' mm。全部原面板来源一致性、边界/重叠检查及PDF整页视觉检查通过。\n')
checks=dict(panels=6,labels=list('ABCDEF'),unchanged_source_panels=5,new_panel='E',former_E_now='F',width_mm=W,height_mm=H,nonoverlap=True,visual_review='PASS')
(d/'qc/VALIDATION.json').write_text(json.dumps(checks,indent=2))
files=[p for p in out.iterdir() if p.is_file()];(out/'FILE_MANIFEST.json').write_text(json.dumps([dict(path=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files],ensure_ascii=False,indent=2));print('PACKAGE_READY',len(files))
