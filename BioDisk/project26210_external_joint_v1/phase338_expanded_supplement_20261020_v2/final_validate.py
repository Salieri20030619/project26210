from pathlib import Path
import json,zipfile,hashlib,gzip,codecs
import pymupdf as fitz
w=Path(__file__).resolve().parent;q=w/'qc';o=w/'delivery'
sel=json.loads((q/'FINAL_TABLE_SELECTION.json').read_text(encoding='utf-8'))
enc={};errors=[]
with zipfile.ZipFile(o/'26210_Expanded_Statistics_v2.zip') as z:
 assert z.testzip() is None
 for x in sel['included']:
  b=z.read(x['file']);assert hashlib.sha256(b).hexdigest()==x['sha256']
  data=gzip.decompress(b) if x['file'].endswith('.gz') else b
  encoding='utf-16' if data.startswith((codecs.BOM_UTF16_LE,codecs.BOM_UTF16_BE)) else 'utf-8-sig' if data.startswith(codecs.BOM_UTF8) else 'utf-8'
  text=data.decode(encoding);assert '\ufffd' not in text,x['file'];enc[encoding]=enc.get(encoding,0)+1
  assert '\t' in text.splitlines()[0],x['file']
mods=json.loads((w/'MODULES.json').read_text(encoding='utf-8'));mods=[x for x in mods if x['id'] not in ['E12','E29','E30']]
(q/'FINAL_MODULES.json').write_text(json.dumps(mods,ensure_ascii=False,indent=2),encoding='utf-8')
doc=fitz.open(o/'26210_Expanded_Supplement_v2.pdf');assert len(doc)==26
text='\n'.join(p.get_text() for p in doc);assert '敏感性' not in text and '\ufffd' not in text
assert all(x['id'] in text for x in mods)
assert all(1<=row[2]<=len(doc) for row in doc.get_toc())
for p in doc:
 assert abs(p.rect.width*25.4/72-210)<.02 and abs(p.rect.height*25.4/72-297)<.02
 for block in p.get_text('dict')['blocks']:
  for line in block.get('lines',[]):
   for span in line['spans']:
    if span['text'].strip() and (span['bbox'][0]<-1 or span['bbox'][1]<-1 or span['bbox'][2]>p.rect.width+1 or span['bbox'][3]>p.rect.height+1):errors.append((p.number+1,span['text']))
assert not errors,errors
main=Path(r'C:\Users\31937\Downloads\26210\Phase335_完整稿Fig2更新_v6\26210_manuscript_Fig2_updated_v6.pdf')
assert hashlib.sha256(main.read_bytes()).hexdigest()=='c88788dfc257078e8315cdbdcb5565ecd431610e7513ab77c148964cb9903503'
v=json.loads((q/'BUILD_VALIDATION.json').read_text(encoding='utf-8'));v.update(status='PASS',visual_review='All 26 pages inspected. S8 point labels separated and four new figure-page placements corrected; final four pages rechecked.',statistical_archive_hashes_verified=len(sel['included']),text_encodings=enc,main_manuscript_unchanged=True,content_scope='Directory-level registration for all indexed stages; result/source review for included modules, not a re-audit of all original analyses.')
(q/'FINAL_VALIDATION.json').write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8')
files=[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in o.iterdir() if p.is_file()]
(q/'DELIVERY_MANIFEST.json').write_text(json.dumps(files,ensure_ascii=False,indent=2),encoding='utf-8')
print('PASS',len(doc),'pages;',len(mods),'modules;',len(sel['included']),'tables;',enc)
