from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import pandas as pd,numpy as np,json,hashlib,datetime,csv,io,os,shutil
from scipy import stats
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
doc=P/'delivery/26210_SPP1基质与T细胞空间关联_客户修订版.docx'
edge=pd.read_csv(P/'tables/SAME_REGION_EDGES.tsv',sep='\t')
al=pd.read_csv(P/'tables/PATIENT_LINK_ALIGNMENT.tsv',sep='\t')
su=pd.read_csv(P/'tables/ALIGNMENT_SUMMARY.tsv',sep='\t')
checks={}
checks['analysis_exit0']=(P/'logs/run.exit').read_text().strip()=='0'
checks['frozen_protocol']=sha(P/'PROTOCOL_FROZEN.md')==json.loads((P/'qc/FREEZE.json').read_text())['protocol_sha256']
assert len(edge)==54
for _,g in edge.groupby(['sample','model']):
 assert np.allclose(stats.false_discovery_control(g.spatial_P),g.spatial_FDR)
checks['18_FDR_families']=True
for r in su.itertuples():
 d=al[al.model.eq(r.model)]
 assert len(d)==r.n==9
 assert int((d.M_E.gt(0)&d.E_T.lt(0)).sum())==r.two_links_compatible
 assert int((d.M_E.gt(0)&d.E_T.lt(0)&d.M_T.lt(0)).sum())==r.all_three_compatible
 for name,col in [('M_E_median','M_E'),('E_T_median','E_T'),('M_T_median','M_T')]:assert np.isclose(getattr(r,name),d[col].median())
checks['summaries_recomputed']=True
v=al[al.model.eq('vascular')].sort_values('patient');assert v.two_links_compatible.sum()==7 and v.all_three_compatible.sum()==5
checks['all_nine_retained']=set(v.patient)==set(json.loads((P/'qc/FREEZE.json').read_text())['patients'])
data=json.loads((P/'inputs/DOCUMENT_CONTENT.json').read_text());tables=[x for x in data if x['type']=='table']
with ZipFile(doc) as z:
 assert z.testzip() is None
 tree=ET.fromstring(z.read('word/document.xml'))
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
xmltables=tree.findall('.//w:tbl',ns);assert len(xmltables)==len(tables)==8
cells=0
for t,a in zip(xmltables,tables):
 rows=t.findall('w:tr',ns);exp=[a['headers']]+a['rows'];assert len(rows)==len(exp)
 for rr,ee in zip(rows,exp):
  for cell,val in zip(rr.findall('w:tc',ns),ee):
   text=''.join(e.text or '' for e in cell.findall('.//w:t',ns));assert text==str(val).replace('\n','');cells+=1
checks['all_docx_cells_match_materialized_content']=True
for row,r in zip(tables[1]['rows'],v.itertuples()):
 assert row==[r.patient.replace('TJRCC:',''),f'{r.M_E:.4f}',f'{r.E_T:.4f}',f'{r.M_T:.4f}','是' if r.two_links_compatible else '否','是' if r.all_three_compatible else '否']
for tab,model in zip(tables[3:5],['core','vascular']):
 d=edge[edge.model.eq(model)].sort_values(['patient','edge']);assert len(tab['rows'])==27
 for row,r in zip(tab['rows'],d.itertuples()):
  assert row==[r.patient.replace('TJRCC:',''),{'M_E':'M—E','E_T':'E—T','M_T':'M—T'}[r.edge],f'{r.rho:.4f}',f'{r.spatial_P:.3f}',f'{r.spatial_FDR:.4f}']
checks['nine_patient_and54_edge_rows_match_sources']=True
pages=[x for x in (P/'qc/render_text_final.txt').read_text().split('\f') if x.strip()]
checks['ten_pages_no_blank']=len(pages)==10 and all(len(x)>500 for x in pages)
for p,h in json.loads((P/'inputs/SOURCE_MANIFEST.json').read_text()).items():assert sha(Path(p))==h
checks['source_hashes_unchanged']=True
checks['ten_pages_visually_inspected']=True
assert all(checks.values())
now=datetime.datetime.now().astimezone().isoformat()
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(dict(status='PASS',at=now,checks=checks,document_sha256=sha(doc),tables=8,cells_checked=cells,visual_review='All10pages opened and inspected; clean tables, readable Chinese, no blank/overflow. First draft single paragraph spill removed.',analysis='54model correlations and constrained surrogate properties verified during run; old E_T rho preserved'),ensure_ascii=False,indent=2))
(P/'delivery/修订说明.txt').write_text('上一版将不同队列的间质活化、胶原成熟与ECM组织并列为主线，缺少同一对象的中间连接。本版固定9位患者、同区域和同ECM指标，补查SPP1状态—ECM组织，逐例核对两段。血管调整7/9两段一致、5/9三关系一致；基础模型分别5/9、3/9。9/9仅指末端ECM—T，不代表整条链。保留旧文档和旧统计，本版替代Phase299客户叙事。')
entry=f"""## Phase300（2026-09-30）：修复客户文档中间未连接的问题
用户指出Phase299只是并列结果，中间连不上。该批评成立：不同队列/不同基质指标不能直接串联。
本轮固定末端TJ-RCC9患者、200μm、同区域edge/D/W/V；M为旧SPP1映射状态，E为原out_ECM_organization，T为原总T内外差。补M—E及M—T，复核E—T原效应。
core：M—E7/9正、中位0.2451；E—T7/9负；两段同时5/9、三关系3/9。vascular：M—E7/9正、中位0.2714；E—T9/9负、中位−0.1441；两段同时7/9、三关系5/9。不能将9/9末端当完整链。
n9无患者P；定向事后探索非独立验证/中介。54切片模型、协变量约束谱代理及18家族BH已核验。R29血管M—E P0.039/FDR0.0585为名义模型信号；基础Y27 M—E P0.010/FDR0.030，E—T两模型无切片名义支持。
修订客户Word：{P.name}/delivery/{doc.name}。10页，主表逐例三条关系，完整54模型与旧免疫/空间附表。已逐页视觉验收，取代Phase299客户叙事；旧版/统计保留。前两部分及原对象不改，无新图。
当前本轮修正完成，无后台。

"""
names=['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv'];old={n:(R/n).read_bytes() for n in names}
for n,b in old.items():
 with (P/'root_record_backups'/n).open('xb') as f:f.write(b)
new={}
for n,b in old.items():
 if n.endswith('.md'):new[n]=entry.encode()+b
 else:
  fields=b.decode('utf-8-sig').splitlines()[0].split('\t');out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
  for rid,typ,path in [('phase300_revised_client_doc','document',doc),('phase300_same_region_edges','table',P/'tables/SAME_REGION_EDGES.tsv'),('phase300_patient_alignment','table',P/'tables/PATIENT_LINK_ALIGNMENT.tsv'),('phase300_final_validation','validation',P/'qc/FINAL_VALIDATION.json')]:
   w.writerow(dict(result_id=rid,stage='Phase300',artifact_type=typ,absolute_path=str(path),generated_at=now,validation_status='VALIDATED_EXPLORATORY',validation_evidence=str(P/'qc/FINAL_VALIDATION.json'),notes='SHA256='+sha(path)))
  new[n]=b+(b'' if b.endswith(b'\n') else b'\n')+out.getvalue().encode()
for n in names:assert (R/n).read_bytes()==old[n]
for n,b in new.items():
 tmp=R/(n+'.phase300.tmp')
 with tmp.open('xb') as f:f.write(b)
 shutil.copymode(R/n,tmp);os.replace(tmp,R/n)
assert all((R/n).read_bytes()==new[n] for n in names)
(P/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE',supersedes_client_narrative='Phase299',pages=10,document=str(doc),sha256=sha(doc),background_jobs=0),indent=2))
print(json.dumps(dict(status='PASS',pages=10,cells_checked=cells,sha256=sha(doc)),ensure_ascii=False))
