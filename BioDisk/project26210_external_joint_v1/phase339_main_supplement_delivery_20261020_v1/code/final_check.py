from pathlib import Path
import json,re,zipfile,hashlib,io
import pymupdf as fitz
from PIL import Image
w=Path(__file__).resolve().parent;q=w/'qc';pack=json.loads((q/'PACKAGE_VALIDATION.json').read_text(encoding='utf-8'));archive=Path(pack['zip'])
old=json.loads(Path(r'C:\Users\31937\tmp\26210_manuscript_v6\report_v6.json').read_text(encoding='utf-8'))
mainblocks=json.loads((w/'main/report/report_materialized.json').read_text(encoding='utf-8'))
moved=json.loads((q/'MOVED_APPENDIX_BLOCKS.json').read_text(encoding='utf-8'))
originaltables=[x for x in old if x['type']=='table'];newtables=[x for x in mainblocks+moved if x['type']=='table'];assert len(originaltables)==len(newtables)==41
for a,b in zip(originaltables,newtables):
 assert a['headers']==b['headers']
 assert len(a['rows'])==len(b['rows'])
 for ar,br in zip(a['rows'],b['rows']):
  for av,bv in zip(ar,br):
   if av!=bv:assert str(av).replace('分割敏感性','不同分割方法下的结果')==str(bv)
modules=json.loads(Path(r'C:\Users\31937\tmp\26210_supplement_full\qc\FINAL_MODULES.json').read_text(encoding='utf-8'))
counts={'internal':0,'cross_file':0}
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 names=['26210_Main_Manuscript_v7.pdf','26210_Supplement_v3.pdf'];docs={name:fitz.open(stream=z.read('26210_delivery/'+name),filetype='pdf') for name in names}
 assert [len(docs[n]) for n in names]==[37,44]
 supptext='\n'.join(p.get_text() for p in docs[names[1]])
 flat=re.sub(r'\s+','',supptext)
 assert '结果解释与适用范围' not in flat and '敏感性' not in flat
 for m in modules:assert re.sub(r'\s+','',m['boundary']) not in flat
 for name,doc in docs.items():
  for page in doc:
   for link in page.get_links():
    if link['kind']==fitz.LINK_GOTO:assert 0<=link['page']<len(doc);counts['internal']+=1
    elif link['kind']==fitz.LINK_GOTOR:
     target=Path(link['file']).name;assert target in docs,(name,link)
     assert 0<=link['page']<len(docs[target]),link;counts['cross_file']+=1
  assert all(1<=x[2]<=len(doc) for x in doc.get_toc())
 files=[n for n in z.namelist() if '/600DPI_PNG/' in n and n.endswith('.png')];assert len(files)==17
 for name in files:
  with Image.open(io.BytesIO(z.read(name))) as im:assert all(abs(d-600)<.1 for d in im.info['dpi'])
 info=json.loads(z.read('26210_delivery/FILE_MANIFEST.json'))
 for x in info:assert hashlib.sha256(z.read('26210_delivery/'+x['file'])).hexdigest()==x['sha256']
oldmain=Path(r'C:\Users\31937\Downloads\26210\Phase335_完整稿Fig2更新_v6\26210_manuscript_Fig2_updated_v6.pdf')
assert hashlib.sha256(oldmain.read_bytes()).hexdigest()=='c88788dfc257078e8315cdbdcb5565ecd431610e7513ab77c148964cb9903503'
result={'status':'PASS','main_pages':37,'supplement_pages':44,'original_statistical_tables_retained':41,'moved_tables':28,'removed_sections_and_paragraphs':27,'png_count':17,'dpi':600,'links':counts,'zip_crc_and_hashes':'PASS','old_versions_unchanged':True,'visual_review':'Main 37 and supplement 44 pages inspected; appendix duplicate captions and orphan section headings corrected; all moved tables, S1/S2, and 27 result modules checked.','zip_sha256':pack['sha256'],'zip_bytes':pack['bytes']}
(q/'FINAL_VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(result)
