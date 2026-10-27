from pathlib import Path
import json,re,copy,hashlib,sys,shutil,subprocess,platform
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
P=R/'phase263_spatial_reintegrated_report_20260926_v1'
OLD=R/'phase208_client_wording_SPP1_narrative_20260918_v6'
SP=R/'phase262_report_spatial_figure5_20260926_v3'
for name in ['code','report','figures','delivery','qc','logs']:(P/name).mkdir(parents=True,exist_ok=True)
rawpath=OLD/'report/report_materialized.json'
spacepath=SP/'text/26210_spatial_figure5_body_v3.md'
source=json.loads(rawpath.read_text());space=spacepath.read_text()
inputs=[rawpath,spacepath,OLD/'scripts/render_report.py',OLD/'report/26210_TYMP_MTA_SPP1_four_figures_narrative_v6.pdf']
fig5=R/'phase261_spatial_ternary_integrated_20260926_v5/delivery/26210_Spatial_A4_GSEA_Ternary_v5.pdf'
for n in range(1,6):
 f=OLD/'figures'/f'Figure{n}.pdf' if n<5 else fig5
 inputs.append(f);shutil.copy2(f,P/'figures'/f'Figure{n}.pdf')
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
input_hashes=[dict(path=str(f),sha256=sha(f)) for f in inputs]
groups={1:[],2:[],3:[],'refs':[]}
intro=space.split('## 1 研究设计与分析方法')[0].split('\n\n',1)[1].strip()
current=None;caption=False
for line in space.splitlines():
 line=line.strip()
 if not line:continue
 if line.startswith('# '):continue
 if line.startswith('## 1 '):current=1;continue
 if line.startswith('## 2 '):current=2;continue
 if line.startswith('## 3 '):current=3;caption=False;continue
 if line=='## 新增参考文献与数据来源':current='refs';continue
 if current is None:continue
 if line=='[[FIGURE5]]':
  groups[current].append(dict(type='figure',figure=5,text='Figure5.pdf'));continue
 if line=='## 图5图注':caption=True;continue
 if line=='[[CAPTION5]]':
  groups[current].append(dict(type='caption',text='图5 代谢背景与SPP1相关状态的组织分布及区域功能特征',keep_with_next=True));continue
 typ='h2' if line.startswith('### ') else 'note' if caption else 'p'
 groups[current].append(dict(type=typ,text=line[4:] if typ=='h2' else line))
changes=[]
new=[]
for i,old in enumerate(source):
 if old['type']=='h1' and old['text'].startswith('2 '):new.extend(groups[1])
 if old['type']=='h1' and old['text'].startswith('3 '):
  new.append(dict(type='p',text=intro));new.extend(groups[2])
 if old['type']=='h1' and old['text'].startswith('4 '):new.extend(groups[3])
 x=copy.deepcopy(old);x['_source_index']=i
 if i==1:
  x['text']='项目 26210　单细胞与空间转录组研究　2026-09-26'
 if x.get('text','').startswith('本文将已完成结果与四张图统一整合：'):
  x['text']=x['text'].replace('四张图','五张图').replace('图4展示C4标志、功能富集和髓系排序。','图4展示C4标志、功能富集和髓系排序，图5展示空间组织分布、区域功能富集及三元距离。')
 if x.get('text')!=old.get('text'):changes.append(dict(source_index=i,before=old.get('text'),after=x.get('text')))
 new.append(x)
 if old.get('text','').startswith('[10] '):new.extend(groups['refs'])
assert len(changes)==2
assert [x for x in new if x['type']=='table']==[dict(x,_source_index=i) for i,x in enumerate(source) if x['type']=='table']
for g,start,end in [(1,12,18),(2,14,21),(3,8,11)]:
 actual=[x['text'].split()[0] for x in groups[g] if x['type']=='h2']
 assert actual==[f'{g}.{n}' for n in range(start,end+1)],actual
assert len([x for x in new if x['type']=='figure'])==5
assert [int(re.match(r'\[(\d+)\]',x['text'])[1]) for x in new if re.match(r'^\[\d+\]',x.get('text',''))]==list(range(1,19))
for i,x in enumerate(source):
 y=next(z for z in new if z.get('_source_index')==i);y={k:v for k,v in y.items() if k!='_source_index'}
 if i not in [c['source_index'] for c in changes]:assert x==y,i
for g in [1,2,3,'refs']:
 for x in groups[g]:assert x in new
(P/'report/report_materialized.json').write_text(json.dumps(new,ensure_ascii=False,indent=2))
(P/'qc/INPUTS.json').write_text(json.dumps(input_hashes,ensure_ascii=False,indent=2))
(P/'qc/EDITORIAL_CHANGES.json').write_text(json.dumps(dict(changes=changes,original_blocks=len(source),spatial_blocks=sum(map(len,groups.values())),intro_moved_before_section='2.14',original_tables=32,new_tests=0),ensure_ascii=False,indent=2))
s=(OLD/'scripts/render_report.py').read_text()
s=s.replace('phase208_client_wording_SPP1_narrative_20260918_v6',P.name)
s=s.replace('26210_TYMP_MTA_SPP1_four_figures_narrative_v6','26210_integrated_spatial_report_v1')
s=s.replace('for n in [1,2,3,4]','for n in [1,2,3,4,5]')
s=s.replace("4:set('ABCDE')}","4:set('ABCDE'),5:set('ABCDEFGHI')}")
s=s.replace('[1-4][A-G]?','[1-5][A-I]?').replace('([1-4])([A-G]?)','([1-5])([A-I]?)')
s=s.replace("level=0,closed=False)","level=1,closed=False)")
s=s.replace("  if n==1:\n", """  if n==5:
   anchors=[('A',3,3),('B',75,3),('C',150,3),('D',3,83),('E',105,83),('F',3,130),('G',3,178),('H',3,220),('I',140,220)]
  elif n==1:
""")
s=s.replace(" return Paragraph(''.join(rich),S[sty])", """ p=Paragraph(''.join(rich),S[sty])
 if sty in ['h1','h2']:p.outline_title=t;p.outline_level=0 if sty=='h1' else 1
 return p""")
s=s.replace("doc=SimpleDocTemplate(str(tmp)", """class ReportDoc(SimpleDocTemplate):
 def afterFlowable(self,flow):
  if hasattr(flow,'outline_title'):
   self._outline_serial=getattr(self,'_outline_serial',0)+1
   key='section_'+str(self._outline_serial)
   self.canv.bookmarkPage(key)
   self.canv.addOutlineEntry(flow.outline_title,key,level=flow.outline_level,closed=False)
doc=ReportDoc(str(tmp)""")
s=s.replace("'四图整合与可点击面板交叉引用'","'空间方法结果讨论整合 五图及可点击交叉引用'")
s=s.replace('figures=4','figures=5')
# Figure5 panel destinations point to the figure page; exact coordinates are not guessed.
s=s.replace("c.bookmarkHorizontalAbsolute(f'fig{n}{letter}',A4[1]-(y-5)*72/25.4,left=x*72/25.4)","c.bookmarkPage(f'fig{n}{letter}') if n==5 else c.bookmarkHorizontalAbsolute(f'fig{n}{letter}',A4[1]-(y-5)*72/25.4,left=x*72/25.4)")
(P/'code/render_report.py').write_text(s)
print(json.dumps(dict(status='PREPARED',blocks=len(new),changes=changes,python=sys.executable,version=platform.python_version()),ensure_ascii=False))

