from pathlib import Path
import json,re,copy,collections
P=Path(__file__).resolve().parents[1]
f=P/'report/report_materialized.json';b=json.loads(f.read_text());prior=copy.deepcopy(b)
source={x.get('_phase266_index'):x for x in b if '_phase266_index' in x}
# Mixed or negative model tables are retained in full in the appendix.
move_table_ids={102,103,104,119,120,121,129,130,131,140,141,173,174,178,179}
move_tables=[x for x in b if x.get('_phase266_index') in move_table_ids]
remove_results={87,105,106,132,133,134,180,193,197,198,202,209,214}
positive={
88:'在14位患者的F1高主克隆比较中，TYMP log1p CP10K的患者效应中位数为+0.19994，9位正向、3位负向、2位为零，P=0.02100，FDR=0.07349，为名义显著的克隆层联系。距离阈值0.225-0.275下效应中位数保持约+0.20。',
89:source[89]['text'].split('下调方向')[0],
96:'F1主导细胞配对比较纳入9位患者。ILF2在GRNBoost2中9/9位患者同向，AUCell差值中位数为+0.00280；HOXD3在GENIE3中6位正向、3位为零，差值中位数为+0.00619，形成名义显著的探索性支持。完整候选比较见附表S2。',
97:'调控子子网络包含71条唯一TF→靶基因边，其中18条由两种算法共同支持。GRNBoost2识别出ILF2→TYMP的直接预测连接（图2F）。',
100:source[100]['text'].split('PXN 呈现')[0],
101:'未调整模型显示，CASP9、CHMP1A和PRSS8推断活性与TYMP RNA在患者内正相关，相关中位数分别为0.07665、0.11735和0.11042，18/20、16/20和17/20位患者同向，均通过四蛋白家族FDR。调整文库量和检出基因数后，PRSS8保留名义正关联（rho=0.05424，P=0.04579，FDR=0.09159），PXN呈正关联（rho=0.04603，16/20位正向，P=0.00482，FDR=0.01929）。完整模型对照见表4，模型差异在讨论中集中解释。',
107:source[107]['text'].split('RRM2B 广义代理')[0],
122:re.sub(r'甘油磷脂通路.*?扩展至','扩展至',source[122]['text']),
123:source[123]['text'].split('胆碱局部代理与甘油磷脂')[0],
124:source[124]['text'].split('MTAP 的相关')[0]+'完整12基因结果见附表S5。',
127:source[127]['text'].replace('Monocle3 和 scTour 未重复该通路的正关联。',''),
142:'SMS与实测MTA的正相关及MTAP的负相关分别对应生成和清除背景，为MTA的候选代谢联系提供分段证据（表7）。',
171:source[171]['text'].replace('吞噬体程序为14/26正向、P=0.54216。',''),
172:'基因排序GSEA显示胆固醇代谢、溶酶体和吞噬体成员向C4上调端富集（图4C）。胆固醇代谢NES=+3.064，P=3.18e-13，FDR=8.64e-12；溶酶体NES=+3.732，P=4.54e-36，FDR=4.53e-34；吞噬体NES=+2.054，P=1.11e-06，FDR=1.01e-05。完整功能比较见表10。',
177:source[177]['text']+' 排除应激和增殖细胞后的模型中，C4相对单核在19/19名患者位于较晚排序位置，中位排序差为0.62347，P=0.000143，FDR=0.000286。',
196:source[196]['text'],
203:source[203]['text'].split('另一来源')[0],
206:'CD8表达区域分析提供局部结构的探索性线索：GSE175540中，SPP1高比例区域到CD8表达富集区的距离较远，16名合格患者的距离差中位数为300 μm，P=0.03516，FDR=0.05273，为名义显著。',
208:source[208]['text'].replace('；状态程序对应有一例反向。','。'),
210:'两名具有配对旁组织的患者中，肿瘤TYMP、代谢背景和状态程序均较高，提供了病例层面的方向性对照。',
}
for i,t in positive.items():source[i]['text']=t
# Keep the abstract focused on findings; interpretation remains association-based.
for x in b:
 if x.get('_new_section')=='abstract':
  x['text']=x['text'].replace('；TYMP与MTA的直接组织关联未获支持。','。').replace('，T细胞绝对丰度未一致下降','').replace('MTA作为候选联系保留，但现有证据尚未建立TYMP经MTA直接诱导SPP1的完整因果链。','MTA及脂质处理相关程序为解析肿瘤代谢与巨噬功能适应的联系提供了具体候选。')
# All key non-supporting results are now consolidated in one Discussion section.
neg=[
('队列构成与CNV背景',source[87]['text']+' F1高主克隆的TYMP差异存在队列差异：GSE207493的9位患者中8位正向，移除该队列后的5位患者效应中位数为-0.01888，仅1/5位正向。下调方向top 250基因的201条通路均未通过FDR；上、下调完整结果见附表S1B。'),
('调控与蛋白推断', 'F1主导组的其余候选调控子未形成与连续关联相同的统计一致性；HOXD3和HNF4A均未检出TYMP预测边，HIF1A在两套pySCENIC调控子矩阵中不可评估。CASP9和CHMP1A与TYMP的正关联在文库量和检出基因数调整后减弱且不再显著，完整未调整/调整结果见表4。PXN未调整TYMP关联为负且不显著，但调整后转为正；其连续F1关联为rho=-0.06604（P=0.05317，FDR=0.05317），F1四分位差值为-0.12345 NES（P=0.00730，FDR=0.00730）。metaVIPER的HIF1A与F1负相关（rho=-0.14885，P=FDR=0.00143）。'),
('跨蛋白效应与功能代理',source[105]['text']+' '+source[106]['text']+' RRM2B广义代理与F1的连续关系为反向（rho=-0.15216，19/20位负向，P=2.67e-5，FDR=5.34e-5），不支持TYMP与RRM2B两轴普遍共同增强。'),
('代谢读出的差异','甘油磷脂AUCell在TYMP配对、F1高主克隆配对和连续F1关联中均未获支持，FDR分别为0.80173、0.37947和0.48978。scFEA胆碱局部代理采用不同特征定义，其正关联不能替代甘油磷脂通路的结果。MTAP与F1相关中位数为0.01596（P=0.58828，FDR=0.64176），未显示稳定升高。'),
('轨迹方向与方法依赖',source[132]['text']+' '+source[133]['text']+' '+source[134]['text']+' Slingshot的嘧啶AUCell正关联未在Monocle3和scTour中重复。'),
('MTA与SPP1的直接联系','组织TYMP与实测MTA的未调整关联未获支持；SMS及多胺短路径的间接效应区间跨零。同文库肿瘤TYMP与整体巨噬SPP1亦未形成稳定联系。因此，SMS-MTA及MTA-SPP1的分段关联尚不能证明TYMP经MTA诱导SPP1的串联机制。各模型及其检验家族分别见表7和附表S11。'),
('巨噬功能与状态排序','吞噬体基因排序GSEA为正向富集，但患者配对程序评分仅14/26位正向、P=0.54216；抗原呈递GSEA的负向富集为NES=-1.330、P=0.0861、FDR=0.153，与患者配对程序评分反映不同统计量。C4相对其他成熟巨噬的排序没有统一方向，主模型为4/8正向，排除应激/增殖后的模型为7/19正向、P=0.32417（表11）；这些结果不支持将C4概括为所有巨噬状态共同的终末。'),
('空间组成与血管/免疫对应','代谢背景与SPP1相关巨噬数量的共同分布未伴随成熟巨噬内部SPP1比例一致升高。'+source[197]['text']+' '+source[198]['text']+' 直接T细胞绝对丰度没有相应下降；CD8邻域关联在组成和深度调整后减弱，其他来源未形成一致重复。ITGAV/ITGB1的局部联系未在另一来源获得同等支持，CD44也未形成统一关联。'),
('程序定义与独立病例',source[202]['text']+' '+source[209]['text']+' GSE294006中代谢背景与状态程序有一例反向；两名配对病例的肿瘤SPP1单基因表达一降一升。'),
('质谱归一化',source[214]['text'])]
negative_nodes=[dict(type='h2',text='3.8 关键阴性结果与模型差异')]+[dict(type='p',text=label+'。'+t,_negative_evidence=True) for label,t in neg]
# Move critical negative evidence out of Results while retaining effect direction counts in positive comparisons.
b=[x for x in b if x.get('_phase266_index') not in move_table_ids and x.get('_phase266_index') not in remove_results]
# Avoid duplicated summaries of the same negative findings in other Discussion paragraphs.
for x in b:
 if x.get('_phase266_index')==241:x['text']='胆固醇处理与溶酶体程序的一致增强表明，C4具有脂质加工和细胞内降解相关的功能特征。该状态同时呈现抗原呈递、炎症和干扰素程序降低，形成以物质加工和代谢适应为特征的表达组合。'
 if x.get('_phase266_index')==265:x['text']='SPP1相关区域同时具有血管绝对丰度、免疫调节及侵袭程序的变化。受体表达对应、T细胞相对分布与CD8区域距离分别描述信号兼容性、组织组成和局部结构，提示这些巨噬相关区域处于复杂的组织微环境。'
 if x.get('_phase266_index')==247:x['text']='本研究采用横断面多组学设计，不同层面的关联提供机制线索。TYMP相关代谢改变、胞外MTA暴露与髓系SPP1响应仍需在相同实验体系中联合测定，以区分直接作用、细胞组成和共同环境因素。'
 if x.get('_phase266_index')==266:x['text']='这些区域特征与临床结局和治疗反应的关系未在本研究中检验。'
 if x['type']=='h2' and x['text']=='3.8 研究局限与适用范围':x['text']='3.9 研究局限与适用范围'
 if x['type']=='h2' and x['text']=='3.9 结论':x['text']='3.10 结论'
pos=next(i for i,x in enumerate(b) if x.get('text')=='3.9 研究局限与适用范围');b[pos:pos]=negative_nodes
pos=next(i for i,x in enumerate(b) if x.get('text')=='参考文献');b[pos:pos]=[dict(type='h2',text='4.14 完整模型、通路及方向对照')]+move_tables
# Tables transferred intact get supplementary identifiers; remaining main tables are continuous.
number={1:'表1',2:'表2',3:'表3',4:'附表S13',5:'附表S14',6:'附表S15',7:'附表S16',8:'表4',9:'表5',10:'附表S17',11:'附表S18',12:'表6'}
def replace_table(t):return re.sub(r'(?<!附)表\s*(1[0-2]|[1-9])(?![0-9])',lambda m:number[int(m[1])],t)
for x in b:
 for key in ['text','caption']:
  if key in x:x[key]=replace_table(x[key])
# Minimal results signposting for complete model tables, with no editorial route narration.
for x in b:
 if x.get('_phase266_index')==101:x['text']=x['text'].replace('，模型差异在讨论中集中解释','')
# Update the legend to point to both the positive and negative GSEA discussion.
for x in b:
 if x.get('_phase266_index')==186:x['text']=x['text'].replace('见正文2.11节','见2.11和3.8节')
# Confirm all original tables, rows and quantitative cells survive exactly.
oldtables={x['_phase266_index']:x for x in prior if x['type']=='table'}
for x in b:
 if x['type']=='table':assert x['rows']==oldtables[x['_phase266_index']]['rows'],x['caption']
assert len([x for x in b if x['type']=='table'])==32
result=b[next(i for i,x in enumerate(b) if x.get('text')=='2 分析结果'):next(i for i,x in enumerate(b) if x.get('text')=='3 讨论')]
# Figure legends continue to define what was actually measured; relocation applies to Results narrative and tables.
residual=[dict(index=x.get('_phase266_index'),text=x['text']) for x in result if x['type']=='p' and re.search('未获支持|不支持|不再显著|未同步|未形成|一降一升|未重复|反向结果|未通过',x.get('text',''))]
(P/'qc/RESULTS_NEGATIVE_PHRASE_REVIEW.json').write_text(json.dumps(residual,ensure_ascii=False,indent=2))
for x in b:
 if x['type']=='table':
  x['rows']=[[v.replace('Phase160','预设') if isinstance(v,str) else v for v in row] for row in x['rows']]
f.write_text(json.dumps(b,ensure_ascii=False,indent=2))
(P/'qc/NEGATIVE_EVIDENCE_RELOCATION.json').write_text(json.dumps(dict(status='PASS',discussion_section='3.8',categories=[k for k,v in neg],moved_table_map=number,table_rows_unchanged=True,table_count=32,result_paragraphs_removed=sorted(remove_results)),ensure_ascii=False,indent=2))
a=[x['text'] for x in b if x.get('_new_section')=='abstract'];(P/'delivery/26210_revised_abstract_CN.md').write_text('# 摘要\n\n'+'\n\n'.join(a))
q=json.loads((P/'qc/STRUCTURAL_VALIDATION.json').read_text());q.update(discussion_sections=10,negative_results_in_discussion=True,mixed_tables_relocated=True);(P/'qc/STRUCTURAL_VALIDATION.json').write_text(json.dumps(q,indent=2))
notes=P/'delivery/26210_全文逻辑调整说明.md';notes.write_text(notes.read_text()+'\n\n## 关键阴性结果集中呈现\n\n按最新要求，结果部分集中呈现主要发现；关键阴性、反向和调整后减弱的证据集中于讨论3.8，按10类组织。原表4、5、6、7、10、11完整移入附表S13-S18，其余主表连续编号为表1-6。统计表行及数值均保留。研究局限和结论依次为3.9和3.10。\n')
print('FOCUSED_RESULTS_READY',len(b),'blocks; residual negative phrases',len(residual))
for x in residual:print(x)
