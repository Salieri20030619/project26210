from pathlib import Path
import json,re,hashlib,shutil
import fitz,numpy as np
r=Path.cwd();d=r/'phase207_manuscript_figures_crossrefs_20260918_v5';out=d/'delivery';figdir=out/'figures';figdir.mkdir(exist_ok=True)
p=d/'report/26210_TYMP_MTA_SPP1_four_figures_crossrefs_v5.pdf';shutil.copy2(p,out/p.name);md=p.with_suffix('.md').read_text();valid={1:'ABC',2:'ABCDEFG',3:'ABCDEF',4:'ABCDE'}
for n in valid:
 anchor='\n'.join(f'<a id="fig{n}{a}"></a>' for a in valid[n]);md=md.replace(f'<a id="fig{n}"></a>',f'<a id="fig{n}"></a>\n'+anchor)
lines=[]
for line in md.splitlines():
 if not line.startswith(('![','<a ')):line=re.sub(r'图\s*([1-4])([A-G]?)',lambda m:f'[图{m[1]}{m[2]}](#fig{m[1]}{m[2]})',line)
 lines.append(line)
(out/p.with_suffix('.md').name).write_text('\n'.join(lines))
for n in [1,2,3,4]:
 doc=fitz.open(d/'figures'/f'Figure{n}.pdf');doc[0].get_pixmap(matrix=fitz.Matrix(200/72,200/72),alpha=False).save(str(figdir/f'Figure{n}.png'))
# Validate each clickable panel destination, including the vertical position.
doc=fitz.open(p);pl=json.loads((d/'qc/PLACEMENTS.json').read_text());fp={x['figure']:x['page']-1 for x in pl};pan=json.loads((r/'phase206_all_mainfigures_A4_20260918_v10/tables/PANEL_PROVENANCE.json').read_text());expected={f'1{k}':(10,(y-5)) for k,y in [('A',24),('B',85),('C',170)]}
for x in pan:expected[str(x['figure']+1)+x['panel']]=(x['x'],x['y']-5)
matched=0
for pg in doc:
 for link in pg.get_links():
  if link['kind']!=1:continue
  tx=pg.get_textbox(link['from']+(-1,-3,1,3));m=re.search(r'图\s*([1-4])([A-G])',tx)
  if m:
   key=m[1]+m[2];assert link['page']==fp[int(m[1])];assert np.allclose(link['to'],np.array(expected[key])*72/25.4,atol=.02);matched+=1
assert matched>=40
qa=json.loads((d/'qc/VALIDATION.json').read_text());qa.update(panel_destination_coordinates_verified=matched,visual_review='All pages reviewed; index table kept on final page');(d/'qc/VALIDATION.json').write_text(json.dumps(qa,indent=2))
(out/'更新说明.md').write_text(f'''# 四图整合文章 v5

- 图1：保留原单细胞预处理与注释图，第4页。
- 图2：新版恶性上皮主图，第9页。
- 图3：新版代谢与巨噬来源主图，第15页。
- 图4：新版C4功能与髓系状态主图，第18页。
- 共{qa['pages']}页，四张图均独立占页、原尺寸嵌入，图注放在图后。
- PDF内{qa['clickable_internal_links']}处图/面板引用可点击跳转，包含{matched}处面板坐标核验；MD保留对应锚点。
- 原31张统计表数据完全保留；新增附表S12用于图2D的P01–P31表达谱身份对照。31个profile来自19名患者，未改变统计样本单位。
- 新版GSEA方法与既有结果已补入1.11和2.12节；吞噬体排序富集与原程序评分检验的不同结果均保留。
- PDF正文是排版版本；Markdown正文和figures/图片用于后续编辑。四张原矢量图和脚本保存在服务器本版目录。
- 旧报告和之前交付版本保留，未覆盖。
''')
files=[f for f in out.rglob('*') if f.is_file() and f.name!='FILE_MANIFEST.json'];manifest=[dict(path=str(f.relative_to(out)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files];(out/'FILE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));print('PACKAGE_READY',len(files),'FILES',matched,'PANEL_LINKS',flush=True)
