from pathlib import Path
import os,json,datetime,zipfile,xml.etree.ElementTree as ET,hashlib,csv
q=Path('/mnt/BioDisk/ludingyu/project26210');o=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase145_macrophage_provenance_trace_remote_v1')
rows=[];hits=[];docs=[];scanned=0
skip={'.git','.venv-pyscenic','envs','external','data_public','__pycache__','.codex','.agents'}
needles=['02_macrophage_final_subtypes','SPP1+Macro','macro_subtype','FindSubCluster','macrophage_subcluster']
for root,dirs,files in os.walk(q):
 dirs[:]=[d for d in dirs if d not in skip and not d.startswith('renv')]
 for name in files:
  f=Path(root)/name
  if f.suffix.lower() in ['.r','.rmd','.py','.ipynb','.sh','.md','.log','.txt'] and f.stat().st_size<30000000:
   scanned+=1;txt=f.read_text(errors='replace');hs=[(i+1,l[:1800]) for i,l in enumerate(txt.splitlines()) if any(n in l for n in needles)]
   if hs:
    raw=f.read_bytes(); rows.append({'path':str(f),'bytes':len(raw),'mtime_utc':datetime.datetime.fromtimestamp(f.stat().st_mtime,datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(raw).hexdigest(),'hits':len(hs)})
    for i,l in hs:hits.append({'path':str(f),'line':i,'text':l})
  if f.suffix.lower()=='.docx' and (f.parent==q or '/docs/' in str(f)):
   try:
    with zipfile.ZipFile(f) as z:
     ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'};r=ET.fromstring(z.read('word/document.xml'))
     ps=[''.join(t.text or '' for t in p.findall('.//w:t',ns)) for p in r.findall('.//w:p',ns)]
    if any('巨噬' in t or 'SPP1' in t for t in ps):
     dest=o/'evidence'/(f.stem+'_text.txt');dest.write_text('\n'.join(ps));docs.append({'path':str(f),'text':str(dest),'paragraphs':len(ps),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
   except Exception as e:docs.append({'path':str(f),'error':str(e)})
(o/'SOURCE_SEARCH_HITS.json').write_text(json.dumps(hits,ensure_ascii=False,indent=2));(o/'SOURCE_SEARCH_FILES.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2));(o/'DOCX_TEXT_EXTRACTIONS.json').write_text(json.dumps(docs,ensure_ascii=False,indent=2))
print('SCANNED',scanned,'MATCHED_FILES',len(rows),'HITS',len(hits));print('FILES',json.dumps(rows,ensure_ascii=False));print('DOCS',json.dumps(docs,ensure_ascii=False))
print('GENERATORS',json.dumps([h for h in hits if h['path'].endswith(('.ipynb','.Rmd')) or 'saveRDS' in h['text']],ensure_ascii=False))
