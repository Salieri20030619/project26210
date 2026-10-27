from pathlib import Path
import json,hashlib,zipfile
import pandas as pd
import pymupdf as fitz
w=Path(__file__).resolve().parent;i=w/'inputs';q=w/'qc';out=w/'delivery'
man=json.loads((q/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
with zipfile.ZipFile(out/'26210_Supplementary_Data_D01-D06_v1.zip') as z:
 assert z.testzip() is None
 for a in man:
  if a['name'].startswith('D'):assert hashlib.sha256(z.read(a['name'])).hexdigest()==a['sha256']
de=pd.read_csv(i/'D03_clone_DE.tsv.gz',sep='\t');sig=de[de.edgeR_FDR<.05]
assert (sig.edgeR_logFC>0).sum()==233 and (sig.edgeR_logFC<0).sum()==527
assert len(pd.read_csv(i/'D04_FORMAL_PATIENT_EFFECTS.tsv',sep='\t'))==836
assert len(pd.read_csv(i/'D04_LEAVE_ONE_GENE_SUMMARIES.tsv',sep='\t'))==3036
old=json.loads(Path(r'C:\Users\31937\tmp\26210_manuscript_v6\report_materialized.json').read_text(encoding='utf-8'))
new=json.loads(Path(r'C:\Users\31937\tmp\26210_manuscript_v6\report_v6.json').read_text(encoding='utf-8'))
key=pd.read_csv(r'C:\Users\31937\tmp\26210_fig2_rebuild\qc\PROFILE_LABEL_KEY.tsv',sep='\t')
table=next(x for x in new if x['type']=='table' and '附表 S12 ' in x['caption'])
assert len(table['rows'])==len(key)==31
for (k,label),row in zip(table['rows'],key.itertuples(index=False,name=None)):
 assert k==row[-1] and label.replace(' ','')==row[0].replace('__','').replace(' ','')
doc=fitz.open(out/'26210_Supplementary_Materials_S3-S7_v1.pdf');assert len(doc)==7
assert [t[2] for t in doc.get_toc()]==list(range(1,8))
for p in doc:
 assert abs(p.rect.width*25.4/72-210)<.01
 assert '\ufffd' not in p.get_text()
v=json.loads((q/'VALIDATION.json').read_text(encoding='utf-8'));v.update(status='PASS',visual_review='All 7 final pages reviewed; S6 colourbar overlap corrected and rechecked; remaining pages unchanged in content.',archive_integrity='all copied statistical files SHA256 identical',main_profile_keys_checked=31)
(q/'VALIDATION.json').write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8')
manifest=[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()]
(q/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print('PASS: seven pages; complete archive intact; 31 main-figure profile keys; counts and labels verified')
