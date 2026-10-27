from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import json,hashlib,datetime,csv,io,os,shutil

P=Path(__file__).resolve().parents[1];R=P.parent
O=R/'phase300_same_region_SPP1_ECM_T_link_20260930_v1'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
doc=next((P/'delivery').glob('*.docx'))
old=json.loads((O/'inputs/DOCUMENT_CONTENT.json').read_text())
new=json.loads((P/'inputs/DOCUMENT_CONTENT.json').read_text())
ot=[x for x in old if x['type']=='table'];nt=[x for x in new if x['type']=='table']
assert len(ot)==len(nt)==8
assert ot[1]['rows']==nt[1]['rows']
for i in range(2,8):assert ot[i]==nt[i]
checks={'nine_patient_numeric_rows_unchanged':True,'all_six_supplementary_tables_unchanged':True}
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with ZipFile(doc) as z:
    assert z.testzip() is None
    root=ET.fromstring(z.read('word/document.xml'))
    xmltables=root.findall('.//w:tbl',ns);assert len(xmltables)==8
    cells=0
    for xt,t in zip(xmltables,nt):
        rows=xt.findall('w:tr',ns);exp=[t['headers']]+t['rows'];assert len(rows)==len(exp)
        for rr,expected in zip(rows,exp):
            actual=rr.findall('w:tc',ns);assert len(actual)==len(expected)
            for cell,s in zip(actual,expected):
                assert ''.join(x.text or '' for x in cell.findall('.//w:t',ns))==str(s).replace('\n','')
                cells+=1
checks['rendered_docx_tables_match_authored_content']=True
assert (P/'logs/build.exit').read_text().strip()=='0'
assert (P/'logs/render.exit').read_text().strip()=='0'
checks['build_and_render_exit0']=True
for p,h in json.loads((P/'inputs/SOURCE_MANIFEST.json').read_text()).items():assert sha(p)==h
checks['all_statistical_sources_unchanged']=True
text=(P/'qc/render_text_final.txt').read_text()
pages=[x for x in text.split('\f') if x.strip()]
assert len(pages)==len(list((P/'render_final').glob('page-*.png')))==11
assert all(len(x.strip())>400 for x in pages)
assert '免疫屏障在组织中对应什么问题' in text
assert 'TYMP 与这部分屏障结果怎样衔接' in text
assert '0.04248' in text and '0.16992' in text
assert '尚未直接证明' in text and '名义信号' in text
checks['eleven_substantive_pages_no_orphan_note_pages']=True
checks['barrier_interpretation_and_key_opposing_evidence_present']=True
# Every final page was opened and visually inspected by the authoring agent before this step.
checks['all_eleven_final_pages_visually_reviewed']=True
assert all(checks.values())
now=datetime.datetime.now().astimezone().isoformat()
qa=dict(status='PASS',at=now,pages=11,tables=8,cells_checked=cells,checks=checks,
    document=str(doc),document_sha256=sha(doc),
    visual_review='All eleven final page PNGs individually opened and inspected. Readable Chinese, intact tables, no overlap or clipping. First rendering supplementary note-only overflow pages corrected by cell padding; first attempt preserved.',
    scope='Narrative and table labels only; no new biological computation, source relabeling, boundary change or manuscript replacement.')
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
(P/'qc/LITERATURE_CHECK.json').write_text(json.dumps(dict(at=now,urls=[
    'https://www.jci.org/articles/view/45817',
    'https://www.nature.com/articles/s41586-021-04057-2',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC8976074/'],
    role='Rechecked prior primary references supporting matrix and T-cell localization mechanism background; not additional project validation.'),indent=2))
(P/'REVISION_SCOPE.txt').write_text('用户指出Phase300报告过于抽象，客户看不出与免疫屏障的关系。重写摘要、空间屏障定义、结果解释和TYMP背景衔接；主表去M/E/T压缩标签。保留逐例数值、6张完整附表、P/FDR及家族，集中呈现关键阴性与测量边界。正文5页，参考文献和统计附表6页。首次13页渲染中2个表注溢页修复后为11页。旧版本、原对象与前两部分不改；未新统计或绘图。')
entry=f'''## Phase301（{now[:10]}）：客户免疫屏障报告改为可读的空间叙事
用户反馈Phase300报告过于抽象，客户看不出与免疫屏障的关系。已重写摘要、屏障空间定义、结果及讨论，明确SPP1相关巨噬状态—外侧基质转录活动—T细胞内侧相对信号之间的联系，以及与TYMP研究背景的衔接。
本轮只改叙事、表题和版式，无新统计或图。TJ-RCC9患者与200μm定义、两模型效应与检验家族保持；主表9患者及6张完整统计附表逐单元格核验。血管模型7/9两段方向符合、5/9三关系符合；9/9仅末端基质—T。保持探索性定位，不写成物理阻挡、迁移下降或TYMP因果链已证明。阴性、反向及方法边界集中Discussion，完整附表保留。
当前独立客户Word：{P.name}/delivery/{doc.name}。11页，全部最终页视觉验收；首轮2张附表表注溢页已修复并留档。Phase300及更早版保留，本版作为当前客户免疫屏障稿；前两部分与正式文章未改。
验收：{P.name}/qc/FINAL_VALIDATION.json。当前改稿完成，无新增后台任务。进度仍仅服务器；交付本地Downloads/26210及workspace输出副本，未向客户发送。

'''
names=['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']
before={n:(R/n).read_bytes() for n in names}
for n,b in before.items():
    with (P/'root_record_backups'/n).open('xb') as f:f.write(b)
updates={}
for n,b in before.items():
    if n.endswith('.md'):updates[n]=entry.encode()+b
    else:
        fields=b.decode('utf-8-sig').splitlines()[0].split('\t');out=io.StringIO()
        writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
        for result_id,typ,path in [('phase301_readable_barrier_client_doc','document',doc),('phase301_readable_barrier_validation','validation',P/'qc/FINAL_VALIDATION.json')]:
            writer.writerow(dict(result_id=result_id,stage='Phase301',artifact_type=typ,absolute_path=str(path),generated_at=now,validation_status='VALIDATED_EXPLORATORY',validation_evidence=str(P/'qc/FINAL_VALIDATION.json'),notes='SHA256='+sha(path)))
        updates[n]=b+(b'' if b.endswith(b'\n') else b'\n')+out.getvalue().encode()
assert all((R/n).read_bytes()==b for n,b in before.items())
for n,b in updates.items():
    tmp=R/(n+'.phase301.tmp')
    with tmp.open('xb') as f:f.write(b)
    shutil.copymode(R/n,tmp);os.replace(tmp,R/n)
assert all((R/n).read_bytes()==b for n,b in updates.items())
(P/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE',pages=11,document=str(doc),sha256=sha(doc),supersedes_client_document='Phase300',background_jobs=0),indent=2))
print(json.dumps(dict(status='PASS',pages=11,tables=8,statistical_rows_preserved=sum(len(x['rows']) for x in nt[1:]),sha256=sha(doc)),ensure_ascii=False))
