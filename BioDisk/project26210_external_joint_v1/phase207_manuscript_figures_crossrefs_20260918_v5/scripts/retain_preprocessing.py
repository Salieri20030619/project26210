from pathlib import Path
import json,re,copy,shutil
r=Path.cwd();d=r/'phase207_manuscript_figures_crossrefs_20260918_v5';old=r/'phase201_three_candidates_C4_main_report_20260918_v4'
b=json.loads((d/'report/report_materialized.json').read_text());orig=json.loads((old/'report/report_materialized.json').read_text())
# New main figures are 2-4; retain the original preprocessing/annotation composite as Figure1.
for n in [3,2,1]:shutil.copy2(r/'phase206_all_mainfigures_A4_20260918_v10/figures'/f'Main_Figure_{n}_A4_v10.pdf',d/'figures'/f'Figure{n+1}.pdf')
shutil.copy2(old/'figures/Figure1_annotation_composite_v1.pdf',d/'figures/Figure1.pdf')
for x in b:
 for key in ['text','caption']:
  if key in x:x[key]=re.sub(r'图\s*([123])',lambda m:'图'+str(int(m[1])+1),x[key])
 if x['type']=='figure':x['figure']+=1;x['text']=f'Figure{x["figure"]}.pdf'
 x.pop('keep_with_next',None)
idx=next(i for i,x in enumerate(b) if x['type']=='h2' and x.get('text','').startswith('1.9 '))
pre=[dict(type='figure',figure=1,text='Figure1.pdf')]+copy.deepcopy(orig[38:42]);b[idx:idx]=pre
for x in b:
 if x.get('text','').startswith('按患者分层选择 8,000 个高变基因'):x['text']+='单细胞预处理后的注释总览及标志基因表达见图1A、图1B、图1C。'
 if x.get('text','').startswith('单细胞整合确立了覆盖'):x['text']=x['text'].replace('总体细胞图谱见图2A。','完整细胞注释及标志表达见图1A、图1B、图1C，主线图谱摘要见图2A。')
 if x.get('text','').startswith('本文将已完成结果与三张主图统一整合'):x['text']=x['text'].replace('与三张主图统一整合：','与四张图统一整合：图1保留原单细胞预处理与注释总览，')
b[1]['text']='项目 26210　四图整合与交叉引用版　2026-09-18'
assert [x['figure'] for x in b if x['type']=='figure']==[1,2,3,4]
assert [x for x in orig if x['type']=='table']==[x for x in b if x['type']=='table'][:31]
(d/'report/report_materialized.json').write_text(json.dumps(b,ensure_ascii=False,indent=2))
refs=[]
for i,x in enumerate(b):
 if x['type'] in ['p','note']:
  for m in re.finditer(r'图\s*([1-4])([A-G]?)',x.get('text','')):refs.append(dict(block=i,figure=int(m[1]),panel=m[2],text=x['text']))
(d/'tables/CROSS_REFERENCES.json').write_text(json.dumps(refs,ensure_ascii=False,indent=2));(d/'qc/CONTENT_CHECKS.json').write_text(json.dumps(dict(original_tables_identical=31,new_index_table=1,total_figures=4,retained_preprocessing_figure=1,new_main_figures=[2,3,4],cross_reference_occurrences=len(refs),superseded_figure_blocks=8),indent=2));print('FOUR_FIGURES_READY',len(refs),'refs')
