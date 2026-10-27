from pathlib import Path
import json,hashlib,zipfile,shutil,math,csv,io,gc
import pymupdf as fitz
from PIL import Image
w=Path(__file__).resolve().parent;q=w/'qc';package=w/'package';package.mkdir(exist_ok=True)
for name in ['26210_Main_Manuscript_v7.pdf','26210_Supplement_v3.pdf']:shutil.copy2(w/'delivery'/name,package/name)
sources={}
for n in range(1,7):sources[f'Main/Figure{n}']=w/'inputs'/f'Figure{n}.pdf'
for n in range(1,12):
 if n<=2:p=w/'inputs'/f'Supplementary_Figure_S{n}.pdf'
 elif n<=7:p=Path(r'C:\Users\31937\tmp\26210_supplement\figures')/f'S{n}.pdf'
 else:p=Path(r'C:\Users\31937\tmp\26210_supplement_full\figures')/f'S{n}.pdf'
 sources[f'Supplement/FigureS{n}']=p
records=[]
for key,p in sources.items():
 target=package/'600DPI_PNG'/f'{key}.png';target.parent.mkdir(parents=True,exist_ok=True)
 vector=package/'Vector_PDF'/f'{key}.pdf';vector.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,vector)
 d=fitz.open(p);assert len(d)==1
 if not target.exists():
  pix=d[0].get_pixmap(dpi=600,alpha=False);pix.save(target);del pix
 with Image.open(target) as im:
  im.verify()
 with Image.open(target) as im:
  dpi=im.info.get('dpi');assert dpi and all(abs(x-600)<.1 for x in dpi)
  expected=(d[0].rect.width*600/72,d[0].rect.height*600/72)
  assert all(abs(a-b)<=1.01 for a,b in zip(im.size,expected)),(im.size,expected)
  records.append({'figure':key,'pixels':list(im.size),'dpi':list(dpi),'png_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'source_pdf_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'width_mm':d[0].rect.width*25.4/72,'height_mm':d[0].rect.height*25.4/72})
 d.close();gc.collect();print('600DPI',key,records[-1]['pixels'],flush=True)
stats=package/'Statistics';stats.mkdir(exist_ok=True)
prior=Path(r'C:\Users\31937\Downloads\26210\Phase338_扩充补充材料_v2\26210_Expanded_Statistics_v2.zip')
with zipfile.ZipFile(prior) as src,zipfile.ZipFile(stats/'26210_Statistics_v3.zip','w',zipfile.ZIP_DEFLATED) as dst:
 for name in src.namelist():
  data=src.read(name)
  if name=='MODULE_INDEX.json':
   a=json.loads(data.decode('utf-8'))
   for m in a:m.pop('boundary',None)
   data=json.dumps(a,ensure_ascii=False,indent=2).encode('utf-8')
  if name=='阅读说明.txt':data='统计数据包与正文v7、补充材料v3配套。统计数据值及来源清单保持原值；模块索引已删除解释与适用范围段落。BASE保留早期配套统计来源。\n'.encode('utf-8')
  dst.writestr(name,data)
 moved=json.loads((q/'MOVED_APPENDIX_BLOCKS.json').read_text(encoding='utf-8'))
 tables=[x for x in moved if x['type']=='table'];assert len(tables)==28
 for j,t in enumerate(tables,1):
  b=io.StringIO();cw=csv.writer(b,delimiter='\t',lineterminator='\n');cw.writerow(t['headers']);cw.writerows(t['rows']);dst.writestr(f'Migrated_Supplement_Tables/Table_{j:02}.tsv',b.getvalue().encode('utf-8'))
 dst.writestr('Migrated_Supplement_Tables/INDEX.json',json.dumps([{'file':f'Table_{j:02}.tsv','caption':t['caption']} for j,t in enumerate(tables,1)],ensure_ascii=False,indent=2).encode('utf-8'))
counts=json.loads((q/'CONTENT_VALIDATION.json').read_text(encoding='utf-8'))
readme=f'''26210 正文、补充材料与600 DPI图件交付包

正文：26210_Main_Manuscript_v7.pdf（{counts['main_pages']}页）
补充材料：26210_Supplement_v3.pdf（{counts['supplement_pages']}页）

原正文第4部分的28张统计附表及说明、补充图S1/S2和图注已迁入补充册。
补充册27个“结果解释与适用范围”小节及其内容已删除。
主图1–6、补充图S1–S11保留原编号，原表S1–S20保留原编号及统计数值。

600DPI_PNG：6张主图+11张补充图，共17张PNG；按源PDF实际尺寸以600 DPI输出。
Vector_PDF：对应17张原始矢量/混合内容PDF，供排版和编辑。
Statistics：完整配套统计数据包，另含迁入的28张附表TSV及索引。

请先完整解压，再打开两份PDF。正文与补充册放在同一目录，以保留跨文件引用跳转。
旧正文与旧补充版本保留；本次未重算统计或修改原数据。
'''
(package/'阅读说明.txt').write_text(readme,encoding='utf-8')
(q/'FIGURE_600DPI_VALIDATION.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
(package/'图件规格与校验.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
files=[{'file':str(p.relative_to(package)).replace('\\','/'),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(package.rglob('*')) if p.is_file()]
(package/'FILE_MANIFEST.json').write_text(json.dumps(files,ensure_ascii=False,indent=2),encoding='utf-8')
target=Path(r'C:\Users\31937\Downloads')/'26210_正文_补充材料_600DPI图件_20261020.zip'
with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
 for p in sorted(package.rglob('*')):
  if p.is_file():z.write(p,'26210_delivery/'+str(p.relative_to(package)).replace('\\','/'))
with zipfile.ZipFile(target) as z:
 assert z.testzip() is None
 for x in files:assert hashlib.sha256(z.read('26210_delivery/'+x['file'])).hexdigest()==x['sha256']
(q/'PACKAGE_VALIDATION.json').write_text(json.dumps({'status':'PASS','zip':str(target),'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'png_600dpi_count':len(records),'vector_pdf_count':len(records),'files':len(files)+1,'crc_and_hash_verified':True},ensure_ascii=False,indent=2),encoding='utf-8')
print('PACKAGE_PASS',target,target.stat().st_size,flush=True)
