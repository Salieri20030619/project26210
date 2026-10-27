from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import pandas as pd,numpy as np,json,hashlib,datetime,csv,io,os,shutil
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
doc=P/'delivery/26210_免疫微环境与基质屏障探索性分析_客户交付版.docx'
data=json.loads((P/'inputs/DOCUMENT_CONTENT.json').read_text())
audit=[x for x in data if x['type']=='table']
with ZipFile(doc) as z:
 assert z.testzip() is None
 tree=ET.fromstring(z.read('word/document.xml'))
 styles=ET.fromstring(z.read('word/styles.xml'))
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
tables=tree.findall('.//w:tbl',ns);assert len(tables)==len(audit)==8
count=0
for t,a in zip(tables,audit):
 rows=t.findall('w:tr',ns);expected=[a['headers']]+a['rows'];assert len(rows)==len(expected)
 for rr,ee in zip(rows,expected):
  cells=rr.findall('w:tc',ns)
  for cell,value in zip(cells,ee):
   text=''.join(x.text or '' for x in cell.findall('.//w:t',ns))
   assert text==str(value).replace('\n',''),(text,value)
   count+=1
# Verify every numerical record copied to the main patient table and 48-row spatial appendix.
sec=pd.read_csv(R/'phase298_matrix_to_immune_entry_20260930_v1/tables/SECTION_RESULTS.tsv',sep='\t')
co=pd.read_csv(R/'phase298_matrix_to_immune_entry_20260930_v1/tables/COHORT_RESULTS.tsv',sep='\t')
z=sec[(sec.dataset=='TJRCC')&(sec.radius==200)&(sec.model=='vascular')&(sec.program=='ECM_organization')&(sec.endpoint=='T')&sec.reason.eq('OK')]
assert len(z)==9 and z.effect.lt(0).all()
for row,r in zip(audit[2]['rows'],z.itertuples()):
 assert row==[r.patient.replace('TJRCC:',''),f'{r.effect:.4f}',f'{r.spatial_P:.3f}',f'{r.spatial_FDR:.3f}']
for tab,ds in zip(audit[-3:],['TJRCC','GSE175540','GLMF_GSE210041']):
 d=co[co.dataset.eq(ds)].sort_values(['radius','model','program','endpoint']);assert len(d)==len(tab['rows'])==16
 for row,r in zip(tab['rows'],d.itertuples()):
  assert int(row[4])==r.n and row[5]==f'{r.negative}/{r.n}' and row[6]==f'{r.median_rho:.4f}'
  assert row[7]==('—' if pd.isna(r.P) else f'{r.P:.5f}')
  assert row[8]==('—' if pd.isna(r.FDR) else f'{r.FDR:.5f}')
src=json.loads((P/'inputs/SOURCE_MANIFEST.json').read_text());assert all(sha(Path(p))==h for p,h in src.items())
texts=[s for s in (P/'qc/render_text_final.txt').read_text().split('\f') if s.strip()]
assert len(texts)==10 and all(len(s)>500 for s in texts)
for i in range(1,10):assert sha(P/f'render_final/page-{i}.png')==sha(P/f'render_v2/page-{i}.png')
alltext=''.join(x.text or '' for x in tree.findall('.//w:t',ns))
assert 'Treg_CD4' not in alltext and 'Phase29' not in alltext and '/mnt/' not in alltext
assert '0.057' in alltext and '0.04248' in alltext and '9位患者全部' in alltext
assert all(e.attrib.get('{'+ns['w']+'}val')=='000000' for e in styles.findall('.//w:color',ns))
now=datetime.datetime.now().astimezone().isoformat()
validation=dict(status='PASS',at=now,document_sha256=sha(doc),pages=10,tables=8,table_cells_checked=count,patient_rows_checked=9,spatial_summary_rows_checked=48,immune_primary_rows=12,source_files_unchanged=True,visual_review=dict(pages=list(range(1,11)),status='PASS',notes='Every final page visually inspected; pages1-9 match reviewed v2 hashes; page10 re-inspected after final wording. No clipping, blank pages, orphan paragraph, missing glyphs or split table rows.'),first_render_issues_resolved=['method overflow','blank page','Treg label'],delivery='DOCX only;PDF and PNG internal QA')
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2))
entry=f"""## Phase299（2026-09-30）：独立客户免疫与基质屏障文档完成
用户明确探索性可用，要求从免疫分析到全部同向信号降低整理为单独客户文档。
交付：{P.name}/delivery/{doc.name}。10页Word，含摘要、方法、结果、讨论、文献、12项免疫功能完整附表、9患者逐例空间P/FDR和48行空间队列结果、补充方法。
重点呈现TJ-RCC 200μm血管调整9/9负向，并同时给出主模型7/9及400μm敏感性。准确解释为ECM组织与T细胞相对内侧信号差的负向关联；不是所有免疫细胞绝对减少或已证明迁移下降。n9无患者P，严格空间检验未通过及其他队列反向集中保留。
8表数值/源哈希/10页渲染视觉验收通过；首轮分页与术语已修正留档。未新分析/绘图，未改原报告、对象或前两部分。进度服务器保存，客户交付文件下载至本机Downloads/26210，未向客户发送。

"""
names=['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv'];old={n:(R/n).read_bytes() for n in names}
for n,b in old.items():
 with (P/'root_record_backups'/n).open('xb') as f:f.write(b)
updates={}
for n,b in old.items():
 if n.endswith('.md'):updates[n]=entry.encode()+b
 else:
  fields=b.decode('utf-8-sig').splitlines()[0].split('\t');s=io.StringIO();w=csv.DictWriter(s,fieldnames=fields,delimiter='\t',lineterminator='\n')
  for rid,typ,path in [('phase299_client_docx','document',doc),('phase299_validation','validation',P/'qc/FINAL_VALIDATION.json')]:
   w.writerow(dict(result_id=rid,stage='Phase299',artifact_type=typ,absolute_path=str(path),generated_at=now,validation_status='PASS',validation_evidence=str(P/'qc/FINAL_VALIDATION.json'),notes='SHA256='+sha(path)))
  updates[n]=b+(b'' if b.endswith(b'\n') else b'\n')+s.getvalue().encode()
for n in names:assert (R/n).read_bytes()==old[n]
for n,b in updates.items():
 tmp=R/(n+'.phase299.tmp')
 with tmp.open('xb') as f:f.write(b)
 shutil.copymode(R/n,tmp);os.replace(tmp,R/n)
assert all((R/n).read_bytes()==updates[n] for n in names)
(P/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE',document=str(doc),pages=10,sha256=sha(doc),background_jobs=0,new_statistics=False),indent=2))
print(json.dumps(validation,ensure_ascii=False))
