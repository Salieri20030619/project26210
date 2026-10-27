from pathlib import Path
import json,re,hashlib,shutil
import fitz,numpy as np
r=Path.cwd();d=r/'phase208_client_wording_SPP1_narrative_20260918_v6';out=d/'delivery';figdir=out/'figures';figdir.mkdir(exist_ok=True)
p=d/'report/26210_TYMP_MTA_SPP1_four_figures_narrative_v6.pdf';shutil.copy2(p,out/p.name);md=p.with_suffix('.md').read_text();valid={1:'ABC',2:'ABCDEFG',3:'ABCDEF',4:'ABCDE'}
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
(out/'更新说明.md').write_text("""# 正文措辞与SPP1发现叙事修订 v6

- 全文31页，四张图及正文21个面板引用保留。
- 将笼统的分析标签改为实际执行的距离阈值比较、文库量及检出基因数调整、排除应激和增殖细胞后的模型等具体说明。
- 正文、表题、表内模型名称、图注及图4曲线图例同步统一；图例为Model 1和Model 2，图注明确各自模型及患者数。
- 2.9按“髓系无监督分群—标志基因比较识别SPP1—三个群的来源比较—选择C4为主分析”重写，并调整2.10开头衔接。
- 数值、P、FDR、样本量及不同模型的结果均保持。四图页码仍为4、9、15、18。
- PDF用于正式阅读；Markdown及figures目录用于编辑。旧版本保留。
""")
for item in [out/p.with_suffix('.md').name,out/'更新说明.md']:
 assert not re.search(r'敏感|Sensitivity|稳健性|鲁棒性|审计|审核',item.read_text(),re.I)
files=[f for f in out.rglob('*') if f.is_file() and f.name!='FILE_MANIFEST.json'];manifest=[dict(path=str(f.relative_to(out)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files];(out/'FILE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));print('PACKAGE_READY',len(files),'FILES',matched,'PANEL_LINKS',flush=True)
