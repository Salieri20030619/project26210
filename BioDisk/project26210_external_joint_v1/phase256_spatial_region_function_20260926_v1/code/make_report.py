from pathlib import Path
import json,sys,hashlib,zipfile,datetime,shutil
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'code'))
from extract import gate,stamp
gate('report_generation')
assert (P/'retry_EXIT_STATUS').read_text().strip()=='0'
qc=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());assert qc['status']=='PASS'
d=pd.read_csv(P/'tables/ALL_FUNCTIONAL_EVIDENCE.tsv.gz',sep='\t',low_memory=False)
pb=pd.read_csv(P/'tables/PSEUDOBULK_FUNCTION_DIRECTION.tsv.gz',sep='\t')
d['name']=d.term.str.split('__').str[-1];pb['name']=pb.term.str.split('__').str[-1]
def row(ds,model,lib,name,version='exclude_definition_genes'):
 a=d[(d.dataset==ds)&(d.model==model)&(d.library==lib)&(d.name.str.lower()==name.lower())&(d.version==version)]
 assert len(a)==1,(ds,model,lib,name,len(a));return a.iloc[0]
def fmt(v):
 if pd.isna(v):return 'NA'
 if isinstance(v,str):return v
 return f'{v:.3g}'
def mdtable(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(v) for v in row)+' |' for row in rows])
selected=[('KEGG','Lysosome','溶酶体'),('KEGG','Cholesterol metabolism','胆固醇代谢'),('KEGG','Phagosome','吞噬体'),('Hallmark','Oxidative Phosphorylation','氧化磷酸化'),('KEGG','Cysteine and methionine metabolism','半胱氨酸/甲硫氨酸代谢')]
tables={}
for model in ['absolute','adjusted','macro_share']:
 out=[]
 for lib,term,cn in selected:
  for ds in ['GSE175540','TJRCC']:
   r=row(ds,model,lib,term)
   out.append([cn,ds,int(r.patient_n),fmt(r.median_Z_difference),f'{int(r.positive)}/{int(r.patient_n)}',fmt(r.P_patient_sign),fmt(r.FDR_patient_sign),fmt(r.NES),fmt(r.P_GSEA),fmt(r.FDR_GSEA)])
 tables[model]=mdtable(['功能','来源','患者n','中位Z效应','正向患者','患者P','患者FDR','NES','GSEA P','GSEA FDR'],out)
# Reference ten GO term display complete and unselected.
rf=d[(d.library=='GO_BP')&(d.model=='absolute')&(d.version=='exclude_definition_genes')&(~d.dataset.str.contains('GLMF'))]
focus=json.loads((P/'resources/FOCUS_QUESTIONS_BEFORE_RESULT_REVIEW.json').read_text())
rf=rf[rf.name.str.lower().isin([s.lower() for s in focus['GO_BP']])]
go_rows=[[r.name,r.dataset,f'{r.positive}/{r.patient_n}',fmt(r.median_Z_difference),fmt(r.P_patient_sign),fmt(r.FDR_patient_sign),fmt(r.NES),fmt(r.P_GSEA),fmt(r.FDR_GSEA)] for r in rf.itertuples(index=False)]
# itertuples row.name above is column named name, preserved as positional attribute.
go_table=mdtable(['参考GO条目','来源','正向/n','Z效应','患者P','患者FDR','NES','GSEA P','GSEA FDR'],go_rows)
pbrows=[]
for lib,term,cn in selected[:3]+[('GO_BP','negative regulation of immune response','免疫反应负调控')]:
 for ds in ['GLMF_GSE210041','GSE175540','TJRCC']:
  r=pb[(pb.dataset==ds)&(pb.model=='absolute')&(pb.version=='exclude_definition_genes')&(pb.library==lib)&pb.name.str.lower().eq(term.lower())].iloc[0]
  pbrows.append([cn,ds,fmt(r.median_mean_gene_log2CPM_difference),f'{r.positive}/{r.patient_n}'])
pbtable=mdtable(['功能','来源','伪bulk平均基因log2CPM差的患者中位数','正向/n'],pbrows)
retry=pd.read_csv(P/'qc/GSEA_NUMERICAL_RETRY_AUDIT.tsv',sep='\t')
counts=pd.read_csv(P/'tables/RESULT_COUNTS.tsv',sep='\t')
counts=counts[(counts.model=='absolute')&(counts.version=='exclude_definition_genes')]
counttable=mdtable(['来源','功能库','可评条目','患者FDR正向','GSEA FDR正向','GSEA FDR负向','GSEA不可估计/描述性'],[[r.dataset,r.library,r.terms,r.patient_FDR_positive,r.GSEA_FDR_positive,r.GSEA_FDR_negative,r.GSEA_NA] for r in counts.itertuples()])
go_min=rf.FDR_GSEA.min()
# Complete versus exclusion comparison for all priority terms and all models.
chosen=d.apply(lambda r:any(r.library==lib and r['name'].lower()==term.lower() for lib,term,_ in selected),axis=1)
d[chosen].to_csv(P/'tables/PRIORITY_PATHWAYS_ALL_MODELS_AND_VERSIONS.tsv',sep='\t',index=False)
text=f"""# 26210：SPP1相关空间区域的功能富集分析
Phase256｜{stamp()}｜参考论文Fig.2G的分析目的｜文字统计交付，无新增图

## 结论
SPP1相关高丰度区域显示明确的溶酶体与胆固醇处理特征；这两类功能在GSE175540（23患者）和TJ-RCC（12患者）的区域比较、组成调整及伪bulk方向检查中均获支持。吞噬体和氧化磷酸化的原始区域信号保留，但其平均功能差异在调整后明显减弱。
当前数据不支持直接制作“普遍免疫抑制功能富集”的Fig.2G式主图。参考论文十类免疫负调控条目在两来源主GSEA中均未达到完整GO家族FDR<0.05；这些条目的原始区域平均表达升高及调整后的反向结果均完整保留。
建议功能图优先展示溶酶体和胆固醇代谢，吞噬体/氧化磷酸化作为带有调整对照的补充。3G热点分析尚未开始。

## 1 方法
### 1.1 数据和独立样本
复用Phase234冻结输入：38张切片、37登记患者；GLMF/GSE210041为2患者2片，TJ-RCC为12患者12片，GSE175540为23患者24片。各来源独立分析；同患者多制备先算切片效果，再等权平均为一个患者效果。读取counts、逐spot输入和原cell2location状态映射，不重训模型。
权威输入为 inputs/INPUT_SECTION_REGISTER.tsv、tables/SECTION_IDENTITY.tsv、tables/PATIENT_MEMBERSHIP.tsv；逐文件身份与检查见qc/SOURCE_MANIFEST_BEFORE.json、AFTER.json。

### 1.2 区域定义和表达效果
主比较absolute采用每片旧SPP1-state q05丰度第80百分位及以上为高区，其余QC spot为背景。adjusted使用同一区域，对log1p_UMI、总体巨噬比例、恶性比例、内皮比例和周细胞比例调整。macro_share限总体巨噬比例>=0.05的spot，以巨噬内SPP1-state比例第80百分位分组。两区均要求>=50spot、高区比例10%–35%；三个模式全部38片通过，共114个切片比较。
表达值为log1p(counts/spot总counts×10000)，每基因在有效spot内Z标准化。absolute/macro_share为高区减背景的平均Z差；adjusted为同一基因Z经协变量调整后的区域系数。条件模型不计算spot层P。基因需在该来源共同测量，并在至少半数有效切片中检出于>=10spot。
辅助伪bulk把各区域整数counts相加，计算log2(CPM+0.5)之差，再对同患者制备等权平均；通路效果为所含基因log2CPM差的平均，仅报告患者方向和效应。

### 1.3 功能库和排序
GO_BP来自现有org.Hs.eg.db与GO.db的GOALL生物过程注释；Hallmark2020和KEGG2021沿用2026-09-14冻结GMT。映射、版本、来源和hash均在resources/。全部三个功能库分别分析，集合与测量宇宙交集须10–500基因且覆盖原有效集合>=20%。
主版本去除TYMP、SPP1、GPNMB、LGALS3、FABP5、CD9、LPL；完整基因为敏感性。从全部合格基因排序，不按单基因显著性筛选。n>=10采用患者Z效应mean/(SD/√n)作为排序量；GLMF仅按平均效果计算描述性ES，不作推断P/NES。
fgseaMultilevel参数：gseaParam=1、minSize10、maxSize500、eps=1e-6、nPermSimple1000、sampleSize101、nproc2、seed262256+job索引。初轮无法估计者单独增加至100000次简单抽样，seed362256+job索引，保留原排序、ES和家族；原结果留档于gsea/initial_1000perm。

### 1.4 患者统计与多重检验
患者通路效果为该患者通路内合格基因Z效果的等权平均；报告中位数、正负患者数。n>=10使用双侧精确符号检验。BH家族为来源×功能库×主分析/敏感性；主分析是absolute+去定义基因，其余五种模式/基因版本共同构成对应敏感性家族。GSEA以同样分家族方式校正，但与患者检验分开。
GSEA随机化P表示给定患者排序中的基因集相对富集；患者符号检验表示通路平均效果的患者一致性。两类P/FDR分别保留。不可估计条目显示NA，并在其已注册GSEA BH池中按P=1计入；eps下限结果按P≤1e-6理解。

## 2 结果
### 2.1 SPP1高值区域同时呈现多类功能升高，选择性信号以代谢和吞噬溶酶体为主
原始区域比较中，GSE175540全部50个Hallmark及TJ-RCC的48/50个Hallmark在患者平均效果层为正向FDR显著；这种广泛升高不能单独用来突出某一种免疫功能。GSEA进一步识别在完整基因排序中相对靠前的程序。完整库计数如下；GLMF的NA为预设描述性定位，不是计算失败。

{counttable}

预先列出的主线功能结果：

{tables['absolute']}

### 2.2 溶酶体和胆固醇代谢在组成调整后仍保留
组成/深度调整后，溶酶体在GSE175540为22/23患者正向、TJ-RCC为12/12正向；胆固醇代谢分别22/23与11/12正向，均达到相应患者层FDR<0.05，并有同向GSEA富集。效果较原始区域比较明显缩小。
吞噬体在调整后GSE175540平均差异转为负中位数、TJ-RCC接近零；其GSEA仍可正向富集，表明相对排序富集与通路整体均值并不等价。氧化磷酸化的调整后患者一致性也未在两来源同时保留。不能只挑GSEA曲线宣称这些平均功能均上升。

{tables['adjusted']}

### 2.3 巨噬内比例比较显示来源差异
macro_share比较中，TJ-RCC溶酶体/胆固醇代谢方向保留较好；GSE175540两者患者正向数分别15/23与16/23，未达对应患者FDR阈值。绝对丰度高区的功能特征不能直接等同于普遍的巨噬状态诱导。

{tables['macro_share']}

### 2.4 参考论文的十类免疫负调控条目未形成同类GSEA主图依据
十类条目在两个外部来源中主GSEA FDR最小值为{fmt(go_min)}，均未通过0.05。原始患者平均表达上升是真实的本轮区域对比结果，但不是相对特异富集的同义词。
例如“免疫反应负调控”在GSE175540原始20/23患者正向；调整后19/23患者负向，中位Z差约−0.0088，患者P约0.00260、FDR约0.00996；TJ-RCC调整后7/12负向，未获统计支持。主要反向结果必须与原始上升一起解释。

{go_table}

### 2.5 区域伪bulk支持溶酶体/胆固醇主线，免疫负调控和吞噬体存在案例差异
伪bulk采用区域counts加和后的表达比例，补充每spot log表达均值。GSE175540和TJ-RCC的溶酶体与胆固醇代谢多数患者方向一致；GLMF溶酶体2/2正向，胆固醇1/2正向，吞噬体2/2负向。下表是辅助描述效果，无新增P检验。

{pbtable}

### 2.6 定义基因排除与数值复核完整保存
上述正文数值均来自排除七个定义/目标基因后的主版本。完整基因与删除版本、三个模式、每患者效果、leading edge及全部阳性/阴性均保留；priority表中可逐项比较，不仅保留方向有利的版本。
初轮GSEA有{int(retry.initial_NA.sum())}条因正负排序量不平衡导致采样不足而无法估计；增加抽样后仍不可估计{int(retry.final_NA.sum())}条，完整保留NA，不当作阴性。到达eps下限的条目见qc统计，P不能理解为比1e-6更精确的值。

## 3 Discussion
本轮最适合的功能叙事是“SPP1相关高丰度空间区域呈现溶酶体和胆固醇处理特征”，用于连接已有巨噬代谢发现与组织空间背景。它补充了现有地图缺少的功能层信息，同时保留丰度、比例和条件模型的不同答案。
区域来自既有表达推断丰度，同数据定义与功能验证仍有依赖；去除七基因降低直接重用，不能消除整套反卷积的依赖。混合spot、组成与测序检出特性均可影响区域均值；本轮以调整和伪bulk结果共同限定解释。模型调整后的剩余差异仍是关联，不代表细胞特异因果效应。
旧空间SPP1-state不是新C4，GLMF两例不能承担人群推断。该分析没有实测局部MTA、没有确认肿瘤病理边界，也不据此宣称免疫屏障、CD8排斥或TYMP→MTA→SPP1因果链；前两部分既有证据按原定义保留。半胱氨酸/甲硫氨酸代谢未在两来源显示一致的调整后平均功能增高，应保持候选背景定位。

## 4 图件采用建议
1. 优先候选：Lysosome与Cholesterol metabolism的规范GSEA曲线，配患者效应及组成调整对照。标注原始区域/调整模型和完整家族FDR。
2. Phagosome和Oxidative phosphorylation可作为“原始区域富集、调整后减弱”的补充，不能只画有利版本。
3. 不制作预设为免疫抑制的十连图；保留全部十类参考GO条目的统计表。
4. 本轮没有生成新图。后续3G宜另行评估SPP1相关信号与上述功能区域的空间聚集；本轮高20%区域不改称热点。

## 5 复现、验收和文件
Python：batch_bench，实际解释器/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python。
R：R_4.5.2，实际解释器/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript。
执行顺序：bash code/run_extract.sh → bash code/run_analysis.sh → bash code/run_retry.sh；再python code/pseudobulk_auxiliary.py与code/make_report.py。必要内存门禁均保存logs/memory_gate.log；未安装或升级依赖。
分析及重试退出状态均为0；完整验收qc/FINAL_VALIDATION.json为PASS：{qc['functional_rows']}条患者功能汇总、{qc['gsea_rows']}条推断性GSEA记录、54组独立程序/ES检查。逐基因效应稠密复算、患者去重、P/BH重算及输入前后身份核对通过；源数据和旧RDS不变。
内部NPZ基因名格式错误的首次失败及修复记录保留；仅转换自建中间文件标签类型，数值不变。数字型格式修复与GSEA数值重试均不被隐藏。
所有分析和进度保留服务器。交付包括本报告、统计ZIP和精简阅读结果表；完整源码、排序、测量基因集和患者矩阵位于本Phase目录。新图0，3G任务0；现有报告与图件未修改。
"""
# Replace accidental pandas Series .name semantics only in explicit tuple tables handled above.
(P/'REPORT_CN.md').write_text(text,encoding='utf-8')
summary=f"""# Phase256：区域功能分析结论
已完成38切片/37患者的功能分析，GLMF2、GSE17554023、TJ-RCC12分开。独立数值核验通过。

1. 溶酶体和胆固醇代谢最适合进入功能图：两个外部来源的原始区域和组成调整均支持，伪bulk方向亦支持。
2. 吞噬体/氧化磷酸化保留原始区域富集，但调整后平均功能差异减弱，不能只展示GSEA阳性。
3. 参考Fig2G十类免疫负调控条目在两来源主GSEA均无完整GO家族FDR显著，不能照搬“免疫抑制”十连图。
4. GLMF：溶酶体两例同向；胆固醇伪bulk1正1负、吞噬体伪bulk2负，案例差异完整保留。
5. 已排除七个定义/目标基因并保存完整基因敏感性。新图与3G尚未启动。

最适合先审核的调整后结果：
{tables['adjusted']}

全部原始/调整/比例、阳性/阴性、患者P/FDR与GSEA P/FDR见统计附件。GSEA P属于给定排序中的基因集随机化；患者P属于独立患者符号检验。详细方法、输入、版本和未解决数值条目见REPORT_CN.md。
"""
(P/'SUMMARY_CN.md').write_text(summary,encoding='utf-8')
# Copy user-facing artifacts only; progress remains server-side.
dest=P/'delivery'
shutil.copyfile(P/'REPORT_CN.md',dest/'26210_区域功能富集_完整报告.md')
shutil.copyfile(P/'SUMMARY_CN.md',dest/'26210_区域功能富集_结论与关键统计.md')
shutil.copyfile(P/'tables/PRIORITY_PATHWAYS_ALL_MODELS_AND_VERSIONS.tsv',dest/'26210_重点通路_全部比较.tsv')
zp=dest/'26210_区域功能富集_完整统计与复现.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for folder in ['tables','qc','code','resources']:
  for f in sorted((P/folder).rglob('*')):
   if f.is_file() and '__pycache__' not in str(f):z.write(f,f.relative_to(P))
 for f in [P/'PLAN_CN.md',P/'inputs/INPUT_SECTION_REGISTER.tsv',P/'gsea/JOBS.tsv']:
  z.write(f,f.relative_to(P))
 for f in sorted((P/'effects').glob('PROGRAM_*.npz')):z.write(f,f.relative_to(P))
 for f in sorted((P/'effects').glob('PATIENT_*.npz')):z.write(f,f.relative_to(P))
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
manifest=[dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted(dest.iterdir()) if f.is_file()]
(P/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
(P/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE_VALIDATED',time=stamp(),new_figures=0,phase3G_started=False),indent=2))
print('REPORT_COMPLETE')
print(json.dumps(manifest,ensure_ascii=False,indent=2))
