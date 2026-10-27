from pathlib import Path
import json,re,copy,shutil
import pandas as pd
r=Path.cwd();d=r/'phase207_manuscript_figures_crossrefs_20260918_v5';old=r/'phase201_three_candidates_C4_main_report_20260918_v4';b=json.loads((old/'report/report_materialized.json').read_text());orig=copy.deepcopy(b)
for n in [1,2,3]:shutil.copy2(r/'phase206_all_mainfigures_A4_20260918_v10/figures'/f'Main_Figure_{n}_A4_v10.pdf',d/'figures'/f'Figure{n}.pdf')
# Remove superseded figure blocks and their legends; all original statistical tables remain.
remove=set()
for i,x in enumerate(b):
 if x['type']=='figure':
  remove.add(i);j=i+1
  while j<len(b) and b[j]['type'] in ['caption','note']:remove.add(j);j+=1
remove.add(154)
b[1]['text']='项目 26210　三张主图整合与交叉引用版　2026-09-18'
# Figure citations only point to panels that actually display the associated evidence.
extra={55:'总体细胞图谱见图1A。',56:'染色体层面的CNV代表性示例见图1B。',68:'主导程序、F1使用量及TYMP表达的空间分布见图1C。患者×程序表达谱的原相似性矩阵见图1D。',81:'患者平衡的调控子活性矩阵见图1E，三个候选TF的原预测靶基因网络见图1F。',94:'F1使用量与四种候选蛋白推断活性的原密度分布见图1G。',109:'各程序内TYMP阴性与阳性细胞的Top20通路展示见图2A，患者配对效应仍以表5为依据。',113:'F1与通路的患者内相关方向及跨患者一致性见图2C。',140:'按实测MTA排序的组织表达分布见图2F。',156:'C4与其他成熟巨噬状态的代表标志基因比较见图3A。',164:'组织实测MTA及SPP1的排序分布见图2F。',170:'原局部表达和程序分布见图3B，患者配对检验见表10。',173:'C4、C10、C15、三群联合及整体巨噬细胞的表达与贡献比较见图2E。',175:'五项研究的原髓系分区轨迹见图3D，SPP1及功能基因的原动态曲线见图3E。',189:'相应标志和功能分布见图3A、图3B，基因排序富集见图3C。',191:'三群位置与来源贡献分别见图2D、图2E。'}
for i,t in extra.items():b[i]['text']+=t
b[138]['text']=b[138]['text'].replace('图6','图2A、图2C')
b[146]['text']=b[146]['text'].replace('（图8A）','（图3A）')+'三群在冻结巨噬图谱中的位置见图2D。'
b[147]['text']=b[147]['text'].replace('图9','图2E')
b[148]['text']=b[148]['text'].replace('（图8B）','（有效配对患者数见图2E）')
b[178]['text']='可估计轨迹中的C4位于单核之后，与其他成熟巨噬的相对排序没有统一方向（图3D、图3E）。本节C4为冻结巨噬候选，与图1D中患者×程序表达谱的C1–C7展示簇采用不同命名体系。'
b[275]['text']='本文将已完成结果与三张主图统一整合：图1展示恶性上皮程序、调控与蛋白活性，图2展示代谢背景及巨噬候选来源，图3展示C4标志、功能富集和髓系排序。既有统计表、原始数值与完整分析档案保留，图内方法和统计说明集中于配套图注。'
b[281]['text']=b[281]['text'].replace('新正文复用既有结果，未新增检验或重新拟合。','本文复用已完成的GSEA、候选群来源比较及既有结果，未新增计算或重新拟合。')
method={'type':'p','text':'C4的预排序GSEA采用26名符合配对条件的患者，按C4与其他成熟巨噬细胞的基因表达配对t统计量对12,048个基因作连续排序；正值表示C4较高。沿用已完成的KEGG分析，在299条合格通路的同一假设家族内作BH校正，并展示预先指定的胆固醇代谢、溶酶体、吞噬体及抗原呈递四条功能通路（图3C）。该分析检验通路成员在基因排序中的富集，与患者配对程序评分检验分别报告。'}
g=pd.read_csv(r/'phase199_fig3_marker_GSEA_20260918_v5/tables/GSEA_FOUR_PROGRAMS_primary.tsv',sep='\t');names={'Cholesterol metabolism':'胆固醇代谢','Lysosome':'溶酶体','Phagosome':'吞噬体','Antigen processing and presentation':'抗原呈递'}
gstats='；'.join(f"{names[z['ID']]}：NES={z['NES']:+.3f}，P={z['pvalue']:.3g}，FDR={z['p.adjust']:.3g}" for z in g.to_dict('records'))+'。'
gresult={'type':'p','text':'基因排序GSEA显示，胆固醇代谢、溶酶体和吞噬体成员向C4上调端富集，抗原呈递的负向富集未达统计支持（图3C）。'+gstats+'其中，吞噬体的基因排序富集与患者配对程序评分未获支持（P=0.54216）是不同统计量的结果，两者分别保留。'}
capdoc=(r/'phase206_all_mainfigures_A4_20260918_v10/delivery/三张主图图注与排版说明_v10.md').read_text()
caps={}
for n in [1,2,3]:
 sect=capdoc.split(f'## Fig{n}：',1)[1].split('\n## ',1)[0];parts=sect.split('\n\n');caps[n]=[dict(type='caption',text=f'图{n} '+parts[0].strip())]
 for p in parts[1:]:
  if not p.strip():continue
  p=p.strip().replace('原效应量、n、P和FDR附于统计表','原效应量、n、P和FDR见表4及附表S3').replace('完整患者/研究/程序对应见索引表','完整患者/研究/程序对应见附表S12').replace('统计数据中','分析档案中').replace('完整非四舍五入统计见附表','完整非四舍五入统计见附表S4A').replace('完整NES、P和FDR见下方统计表','完整NES、P和FDR见正文2.12节').replace('统计附件','统计附表')
  caps[n].append(dict(type='note',text=p))
new=[]
for i,x in enumerate(b):
 if i in remove:continue
 if i in [108,169,180]:
  n={108:1,169:2,180:3}[i];new.append(dict(type='figure',figure=n,text=f'Figure{n}.pdf'));new.extend(caps[n])
 if i==52:new.append(method)
 if i==123:new.append(dict(type='p',text='恶性上皮细胞的原Monocle3轨迹用于展示状态排序背景（图2B）；三种轨迹方法的患者内关联和末端比较分别汇总于表6及附表S6–S8。'))
 new.append(x)
 if i==170:new.append(gresult)
keys=json.loads((r/'phase206_all_mainfigures_A4_20260918_v10/tables/PATIENT_LABEL_KEY.json').read_text())
new.extend([dict(type='h2',text='4.13 图1D患者×程序表达谱索引'),dict(type='p',text='图1D的P01–P31为31个患者×程序表达谱，涉及19名独立患者。短索引只用于图形标记，不改变原样本身份或统计单位。'),dict(type='table',caption='附表 S12 图1D表达谱短索引与原身份',headers=['图中索引','研究、患者及程序'],rows=[[v['display'],v['original']] for v in keys])])
assert [x for x in orig if x['type']=='table']==[x for x in new if x['type']=='table'][:31]
assert [x['figure'] for x in new if x['type']=='figure']==[1,2,3]
assert not any(re.search(r'图\s*[4-9](?!\d)',x.get('text','')) for x in new)
(d/'report/report_materialized.json').write_text(json.dumps(new,ensure_ascii=False,indent=2))
refs=[]
for i,x in enumerate(new):
 if x['type'] in ['p','note']:
  for m in re.finditer(r'图\s*([123])([A-G]?)',x.get('text','')):refs.append(dict(block=i,figure=int(m[1]),panel=m[2],text=x['text']))
(d/'tables/CROSS_REFERENCES.json').write_text(json.dumps(refs,ensure_ascii=False,indent=2));(d/'qc/CONTENT_CHECKS.json').write_text(json.dumps(dict(original_tables_identical=31,new_index_table=1,main_figures=3,cross_reference_occurrences=len(refs),old_figure_blocks_removed=9),indent=2));print('REVISED',len(new),'blocks',len(refs),'refs')
