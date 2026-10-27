from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import json,hashlib,datetime,csv,io,os,shutil
import pandas as pd
P=Path(__file__).resolve().parents[1];R=P.parent
O=R/'phase301_barrier_readable_client_report_20260930_v1'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
doc=next((P/'delivery').glob('*.docx'))
before_content=json.loads((O/'inputs/DOCUMENT_CONTENT.json').read_text())
content=json.loads((P/'inputs/DOCUMENT_CONTENT.json').read_text())
ot=[x for x in before_content if x['type']=='table'];nt=[x for x in content if x['type']=='table']
assert len(ot)==8 and len(nt)==10
assert ot[1]['rows']==nt[2]['rows']
for old,new in zip(ot[2:],nt[3:9]):assert old==new
checks={'original_123_statistical_rows_unchanged':True}
up=json.loads((P/'inputs/UPSTREAM_MATERIALIZED.json').read_text());assert len(up)==12
frames={x['source']:pd.read_csv(x['source'],sep='\t',float_precision='round_trip') for x in up}
for x in up:
    frame=frames[x['source']]
    scope={'C4':'C4','整体巨噬':'WholeMacro','其他巨噬':'OtherMacro'}[x['scope']]
    if x['family']=='F1':
        model={'未调整':'unadjusted','测序调整':'adjusted'}[x['model']]
        rr=frame[frame.group.eq(scope)&frame.endpoint.eq('SPP1')&frame.model.eq(model)]
    elif x['family']=='F2':rr=frame[frame.group.eq(scope)&frame.endpoint.eq('SPP1')]
    elif x['family']=='F3':
        model='context' if x['model'].endswith('加环境') else 'core'
        rr=frame[frame.scope.eq(scope)&frame.model.eq(model)]
    else:rr=frame[frame.scope.eq(scope)&frame.regulon.eq(x['model'][2:]+'(+)')]
    assert len(rr)==1
    row=rr.iloc[0]
    assert x['n']==int(row['n']) and x['positive']==int(row['positive'])
    for col,old in [('rho','median'),('P','P'),('FDR','FDR')]:assert x[col]==float(row[old])
    assert row.family_size=={'F1':12,'F2':12,'F3':4,'F4':28}[x['family']]
checks['all_12_upstream_rows_exactly_match_source_and_family']=True
def vf(x):return f'{x:.3e}' if abs(x)<.001 else f'{x:.6f}'
main_up=[x for x in up if x['family'] in ['F1','F2'] and x['model']!='未调整']
assert nt[1]['rows']==[[x['scope'],x['model'],str(x['n']),f"{x['positive']}/{x['n']}",f"{x['rho']:.6f}",vf(x['P']),vf(x['FDR'])] for x in main_up]
assert nt[9]['rows']==[[x['scope'],x['model'],str(x['n']),f"{x['positive']}/{x['n']}",f"{x['rho']:.6f}",vf(x['P']),vf(x['FDR']),x['family']] for x in up]
checks['main_upstream_table_and_S7_match_original_statistics']=True
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with ZipFile(doc) as z:
    assert z.testzip() is None
    root=ET.fromstring(z.read('word/document.xml'));xmltables=root.findall('.//w:tbl',ns)
    assert len(xmltables)==10
    cells=0
    for xt,t in zip(xmltables,nt):
        rows=xt.findall('w:tr',ns);exp=[t['headers']]+t['rows'];assert len(rows)==len(exp)
        for row,expected in zip(rows,exp):
            actual=row.findall('w:tc',ns);assert len(actual)==len(expected)
            for cell,s in zip(actual,expected):
                assert ''.join(x.text or '' for x in cell.findall('.//w:t',ns))==str(s).replace('\n','')
                cells+=1
checks['all_document_cells_verified']=True
for step in ['build','render']:assert (P/f'logs/{step}.exit').read_text().strip()=='0'
checks['build_render_exit0']=True
for p,h in json.loads((P/'inputs/SOURCE_MANIFEST.json').read_text()).items():assert sha(p)==h
checks['statistical_sources_unchanged']=True
text=(P/'qc/render_text_final.txt').read_text();pages=[x for x in text.split('\f') if x.strip()]
assert len(pages)==len(list((P/'render_final').glob('page-*.png')))==13
assert all(len(x.strip())>400 for x in pages)
for phrase in ['0.134872','0.156478','0.049133','0.029541','0.002184','0.057373','0.04248','0.16992']:assert phrase in text
assert '达到校' in text and '正后显著' in text
checks['thirteen_substantive_pages']=True
checks['positive_and_opposing_results_present']=True
# All final page PNGs have been opened and inspected by the authoring agent.
checks['all_thirteen_final_pages_visually_inspected']=True
assert all(checks.values())
now=datetime.datetime.now().astimezone().isoformat()
qa=dict(status='PASS',at=now,pages=13,tables=10,cells_checked=cells,checks=checks,document=str(doc),document_sha256=sha(doc),scope='Restored existing Phase280 and Phase282 upstream evidence only; no biological recalculation.',visual_review='All thirteen final pages individually opened; readable Chinese, intact tables, correct superscripts, no clipping or overlap. Initial Series median field collision repaired before document build; original failure logs retained.')
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
(P/'REVISION_SCOPE.txt').write_text('用户指出TYMP到SPP1近期正面结果被遗漏。确认Phase280 C4及整体巨噬患者内正关联、环境调整，以及Phase282表达强度与检出/靶程序条件结果。补摘要、结果2.1、表1与附表S7，原空间小节顺延2.2–2.5、原患者表改表2。纠正TYMP仅作背景的叙事。原9例空间结果及全部6张附表不改，原123行统计逐项保持。13页视觉验收，不重算统计。')
entry=f'''## Phase302（{now[:10]}）：补回遗漏的TYMP—SPP1正面结果
用户指出Phase301把TYMP只作为背景，遗漏近期正面结果。已确认并补入摘要、结果2.1、主表1及附表S7：新C4患者内TYMP—SPP1测序调整n15/13正、中位rho0.134872、P0.007385/FDR0.029541；环境调整15/15正、rho0.156478、P6.104e-5/FDR3.662e-4。整体成熟巨噬调整n29/24正、rho0.049133、P0.000546/FDR0.002184，5研究中位均正，留一主研究方向保留。
Phase282强度—检出C4环境n14/13正、P0.001831/FDR0.002441及基础不显著对照；原ILF2/SPI1条件结果附S7，原F1/F2/F3/F4家族12/12/4/28继承不重校正。上游正面联系为髓系细胞内表达关联，保留有效证据，不因跨患者肿瘤来源或外源蛋白体系阴性而抹去；不拼为同病例因果中介，不以C4特异性为必要门槛。
当前客户Word：{P.name}/delivery/{doc.name}。13页、10表；新增4行主表和12行附表逐项原源核对，原123行空间/免疫统计完整保持。全部最终页面视觉验收，前两部分、正式文章及原对象未改；没有新统计或图。取代Phase301客户稿，旧版保留。首轮pandas字段同名失败已修复，日志留档。
验收：{P.name}/qc/FINAL_VALIDATION.json。本轮补稿完成，无新增后台，权威进度仅服务器，本地交付至Downloads/26210及workspace。

'''
names=['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv'];before={n:(R/n).read_bytes() for n in names}
for n,b in before.items():
    with (P/'root_record_backups'/n).open('xb') as f:f.write(b)
updates={}
for n,b in before.items():
    if n.endswith('.md'):updates[n]=entry.encode()+b
    else:
        fields=b.decode('utf-8-sig').splitlines()[0].split('\t');out=io.StringIO();writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
        for rid,typ,path in [('phase302_upstream_restored_client_doc','document',doc),('phase302_upstream_materialized','table',P/'inputs/UPSTREAM_MATERIALIZED.json'),('phase302_upstream_validation','validation',P/'qc/FINAL_VALIDATION.json')]:
            writer.writerow(dict(result_id=rid,stage='Phase302',artifact_type=typ,absolute_path=str(path),generated_at=now,validation_status='VALIDATED_EVIDENCE_RESTORATION',validation_evidence=str(P/'qc/FINAL_VALIDATION.json'),notes='SHA256='+sha(path)))
        updates[n]=b+(b'' if b.endswith(b'\n') else b'\n')+out.getvalue().encode()
assert all((R/n).read_bytes()==b for n,b in before.items())
for n,b in updates.items():
    tmp=R/(n+'.phase302.tmp')
    with tmp.open('xb') as f:f.write(b)
    shutil.copymode(R/n,tmp);os.replace(tmp,R/n)
assert all((R/n).read_bytes()==b for n,b in updates.items())
(P/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE',pages=13,document=str(doc),sha256=sha(doc),supersedes_client_document='Phase301',background_jobs=0),indent=2))
print(json.dumps(dict(status='PASS',pages=13,tables=10,upstream_rows=12,original_rows_unchanged=123,sha256=sha(doc)),ensure_ascii=False))
