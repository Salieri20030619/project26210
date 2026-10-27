from common import *
import shutil,zipfile
gate('report_generation')
assert (P/'validation_EXIT_STATUS').read_text().strip()=='0'
qc=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());assert qc['status']=='PASS'
sp=pd.read_csv(P/'tables/SPATIAL_COHORT_SUMMARY.tsv',sep='\t')
pair=pd.read_csv(P/'tables/PAIRED_ASSOCIATIONS.tsv',sep='\t')
coverage=pd.read_csv(P/'tables/PAIRED_PROGRAM_COVERAGE.tsv',sep='\t')
def fmt(v):
 if pd.isna(v):return 'NA'
 if isinstance(v,str):return v
 return f'{v:.4g}'
def md(headers,rows):
 return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(x) for x in r)+' |' for r in rows])
names={'tumor_TYMP':'肿瘤TYMP','macro_TYMP':'巨噬TYMP','macro_SPP1':'巨噬SPP1','tumor_lys':'肿瘤溶酶体程序','tumor_chol':'肿瘤胆固醇程序','macro_lys':'巨噬溶酶体程序','macro_chol':'巨噬胆固醇程序','TYMP_RNA':'TYMP RNA','SPP1_RNA':'SPP1 RNA','V1_SPP1_q05':'SPP1-state绝对丰度','V3_SPP1_macro_share':'巨噬内SPP1-state比例','Lysosome':'溶酶体程序','Cholesterol metabolism':'胆固醇程序'}
def pairtable(df):return md(['关系','患者n','调整后rho','P','FDR'],[[names[r.x]+'—'+names[r.y],r.patient_n,fmt(r.rho),fmt(r.P),fmt(r.FDR)] for r in df.itertuples()])
def sptable(df):return md(['来源','关系','有效n','正/负','中位rho','患者P','FDR'],[[r.dataset,names[r.x]+'—'+names[r.y],r.patient_n,f'{r.positive}/{r.negative}',fmt(r.median_rho),fmt(r.P),fmt(r.FDR)] for r in df.itertuples()])
primarysp=sp[sp.role=='primary'];primarypair=pair[pair.role=='primary']
pair_primary=pairtable(primarypair)
tymp_sp=sptable(primarysp[primarysp.family=='TYMP_function'])
spp_sp=sptable(primarysp[primarysp.family=='function_SPP1'])
direct_sp=sptable(primarysp[primarysp.family=='TYMP_SPP1'])
sens=sp[(sp.family=='TYMP_function')&(sp.role=='sensitivity')&~sp.dataset.str.contains('GLMF')]
sensitive_table=md(['来源','程序','尺度','调整','版本','n','正/负','中位rho','P','FDR'],[[r.dataset,names[r.y],r.scale,r.adjustment,r.version,r.patient_n,f'{r.positive}/{r.negative}',fmt(r.median_rho),fmt(r.P),fmt(r.FDR)] for r in sens.itertuples()])
ps=pair[pair.family.isin(['tumor_intracellular','macro_SPP1_function'])&(pair.role=='sensitivity')]
paired_sens=md(['关系','模型','版本','rho','P','FDR'],[[names[r.x]+'—'+names[r.y],r.adjustment,r.version,fmt(r.rho),fmt(r.P),fmt(r.FDR)] for r in ps.itertuples()])
references="""[1] Distinctive metabolic remodeling in TYMP deficiency beyond mitochondrial dysfunction. 2023. https://pubmed.ncbi.nlm.nih.gov/37603049/
[2] Lysosomal dysfunction and overload of nucleosides in thymidine phosphorylase deficiency of MNGIE. 2024. https://link.springer.com/article/10.1186/s12967-024-05275-8
[3] Macrophages in epididymal adipose tissue secrete osteopontin to regulate bone homeostasis. https://pmc.ncbi.nlm.nih.gov/articles/PMC8776868/
"""
section="""# 独立小节草稿：TYMP与溶酶体/胆固醇处理的机制联系
建议接在现有空间正文后；以下编号承接Phase253的1.15、2.19、3.10。本文件为独立待审草稿，尚未写入原报告。

## 1.16 TYMP与溶酶体/胆固醇程序的定向关联分析
在既有38张空间切片中，采用冻结KEGG溶酶体和胆固醇代谢基因集，排除TYMP、SPP1及五个旧SPP1状态定义基因后计算程序分数。分别计算TYMP、功能程序与SPP1相关读出的偏Spearman相关，控制测序深度及细胞组成；同患者多切片以Fisher-z平均汇总，各来源独立进行患者层双侧符号检验及预定义家族BH校正。另在16位同文库患者的肿瘤和巨噬细胞伪bulk中计算程序分数，采用研究固定效应及研究内残差置换检验。去除两通路共享的LIPA、NPC1、NPC2、SORT1，以及空间格块、附加协变量和留一研究分析作为预设敏感性检查。

## 2.20 TYMP与胆固醇处理程序相关，SPP1相关区域呈现相应功能特征
同文库细胞区室分析显示，肿瘤细胞TYMP与肿瘤内胆固醇程序呈正相关（n=16，研究调整rho=0.729，P=0.00905，FDR=0.0181）。排除与溶酶体集合共享的四个基因后，该关联仍保留（rho=0.704，P=0.0152，敏感性家族FDR=0.0345）；进一步调整肿瘤总counts所得方向相同。肿瘤TYMP与肿瘤内溶酶体程序的rho为0.552，但研究调整检验未获支持（P=0.0635，FDR=0.0635）。在巨噬区室内，SPP1与胆固醇程序呈名义正相关（rho=0.647，P=0.0282，FDR=0.0563），与溶酶体程序的关联未达到统计支持（rho=0.601，P=0.0723，FDR=0.0723）。

空间分析表明，控制深度与组成后，TYMP与溶酶体程序在GSE175540和TJ-RCC中的患者中位rho分别为0.0298和0.0303，18/23和11/12患者正向（P=0.0106/0.00635，FDR=0.0212/0.00635）；与胆固醇程序的中位rho分别为0.0180和0.0201，17/23和12/12患者正向（P=0.0347/0.000488，FDR=0.0347/0.000977）。这些效应较弱，去共享基因后方向保留。相较之下，两功能程序与旧SPP1-state的绝对丰度和巨噬内比例关联更一致：溶酶体程序与状态比例的中位rho为0.189/0.147，胆固醇程序为0.118/0.104，两来源均为全部23/23和12/12患者正向（对应患者P=2.38×10⁻⁷/0.000488，六端点家族FDR=3.58×10⁻⁷/0.000732）。

四乘四格块的调整分析未复现TYMP—功能联系的患者层统计支持。同文库跨区室比较也未支持肿瘤TYMP与巨噬溶酶体、巨噬胆固醇程序或巨噬SPP1表达的稳定联系（rho分别−0.027、0.245、−0.108；P分别0.937、0.451、0.795；FDR分别0.937、0.902、0.795）。因此，本轮支持肿瘤TYMP相关胆固醇表达程序，以及与SPP1状态相关的空间功能背景；当前数据不支持将两者写成已建立的跨细胞调控链。

## 3.11 TYMP相关代谢稳态与SPP1状态的功能衔接
已有TYMP缺陷研究显示，恢复胸苷磷酸化酶活性能够挽救SREBP相关胆固醇/脂肪酸代谢改变[1]；另一项研究显示TYMP不足伴随溶酶体核苷积累、酸化及降解功能受损[2]。这些细胞内实验为本项目TYMP—代谢程序关联提供机制背景，但其证据主要来自MNGIE及相关细胞模型。SPP1编码的骨桥蛋白也参与巨噬细胞脂质摄取及溶酶体相关脂解调节，且其作用并非单向增强脂解[3]。因此，溶酶体与胆固醇程序可作为连接既有肿瘤代谢发现和SPP1相关巨噬表型的功能层描述。
本轮是基于既有结果的定向验证，空间读出来自混合spot及同数据反卷积；旧SPP1-state不等同新C4。16位配对患者中10位来自GSE207493，删除该研究后仅剩6位，研究调整模型自由度不足，尚不能称为独立多研究重复。转录程序不代表实际脂质通量或溶酶体酸化功能，同片分段关联也不构成中介证据。后续若需要确立肿瘤TYMP对巨噬功能的方向性作用，应依赖区室明确的干预证据；当前保留原有MTA研究主线，不用本小节替换该主线。

## 本小节参考文献（整合时接续原编号）
"""+references
(P/'SECTION_DRAFT_CN.md').write_text(section,encoding='utf-8')
text=f"""# Phase257：TYMP—溶酶体/胆固醇—SPP1针对性验证
完成时间：{stamp()}。本轮为新增独立小节提供数据依据，无新图、无热点分析、无模型重训。

## 结论
可以新增一个机制联系小节，但应采用分层关联叙事。最明确的新结果是肿瘤细胞TYMP与肿瘤内胆固醇程序正相关（n16，研究调整rho0.729，P0.00905，FDR0.0181），去共享基因及总counts调整后保留。空间中TYMP与两程序为弱而方向一致的联系，程序与旧SPP1-state丰度/比例的联系更一致。当前数据不支持肿瘤TYMP—巨噬功能或巨噬SPP1的稳定跨区室关联；不能写成已确认的TYMP经溶酶体/脂质处理诱导SPP1机制。
可直接审阅的正文草稿：SECTION_DRAFT_CN.md，包含方法1.16、结果2.20、讨论3.11；既有报告尚未修改。

## 1 Methods
### 1.1 输入与问题冻结
本轮已知Phase256结果后，冻结两条KEGG程序Lysosome与Cholesterol metabolism。空间输入沿用Phase234/256：38片37患者，GLMF2、GSE17554023和TJRCC12分别分析。同文库输入为Phase154固定16患者清单及Phase157分区室伪bulk，包含GSE207493 10、GSE156632 4、GSE159115 1、GSE178481 1患者。肿瘤counts和巨噬(C4+OtherMacrophage)counts逐样本与原表总量、细胞数、TYMP/SPP1整数counts精确核对。
排除TYMP/SPP1/GPNMB/LGALS3/FABP5/CD9/LPL；另以删除两集合共享LIPA/NPC1/NPC2/SORT1的版本作敏感性。每集合至少10可测基因且覆盖>=20%，成员完整保存。

### 1.2 空间效应
程序为log1p(CP10k)逐片基因Z均值，精确复核其高低区域差与Phase256对应效果。全片调整偏Spearman为主模型，控制log1pUMI、总巨噬比例、恶性比例、内皮比例和周细胞比例；当目标是SPP1绝对q05时用排除SPP1的其他巨噬比例替代总巨噬比例。未调整、冻结4x4格块中位数的未调/调整、删除共享基因后的全片调整为敏感性。
全片至少50spot，格块至少8区，残差自由度>=3、设计满秩且条件数<30；不计算spot或格块P。同患者多制备rho以Fisher-z等权平均为一个患者效果，n>=10使用双侧符号检验及10000次患者bootstrap中位rho区间；n<10仅方向与效果。BH按来源×关系家族×主/敏感池进行。三个家族分别TYMP—功能2端点、功能—SPP1三读出6端点、TYMP—SPP1三读出3端点；所有预定敏感版本在各自敏感池内共同校正。
对TYMP—两程序逐次删除原4x4单个格块，记录方向稳定性，不新增P。

### 1.3 同文库区室关联
各区室伪bulk为log2(CPM+1)，基因需CPM>=1于至少半患者且有变异；逐基因跨患者Z后平均为程序。主模型为研究固定效应偏Spearman，研究内残差置换19999次；P=(hits+1)/20000，实际seed由端点标识的SHA256前8位与262257异或，逐行保存。研究内14患者可交换，两个单例研究患者不参与交换，登记n仍为16。置换空间若不超过20000则全排列精确计算。
敏感性为未调Spearman、附加C4比例(巨噬端点)或肿瘤总counts(肿瘤内)、删除共享基因的研究调整。五个预定义家族分别为肿瘤内TYMP—两程序、肿瘤TYMP—巨噬两程序、巨噬SPP1—两程序、巨噬TYMP—两程序、肿瘤TYMP—巨噬SPP1；主与敏感池分别BH。留一研究只报告rho，不对小样本追加推断。
来源路径、代码、冻结方案、家族大小、设计条件数、程序覆盖及实际置换规模均保存。

## 2 Results
### 2.1 肿瘤TYMP与肿瘤胆固醇程序联系最明确
{pair_primary}

肿瘤内胆固醇关联在研究调整、额外控制肿瘤总counts及去共享基因后均保留。肿瘤内溶酶体关联原始模型有信号，但研究调整P0.06345未获支持；巨噬SPP1—胆固醇为名义显著，FDR0.05630，不能称校正后显著。巨噬TYMP与两程序均未获支持。

### 2.2 TYMP的空间功能联系方向重复，但全片调整后效应较弱
{tymp_sp}

GLMF两例中TYMP—溶酶体均正向，胆固醇一正一负；不作P检验。外部两来源调整后rho约0.02–0.03，方向可重复，但不能称强相关。

### 2.3 两功能与SPP1-state丰度及比例的联系更一致
{spp_sp}

与SPP1单基因RNA的结果较弱且来源不同，不能将映射状态、单基因RNA和新C4互换。各端点按实际变量解释。

### 2.4 TYMP直接对应SPP1绝对丰度，未对应一致的状态比例或单基因表达
{direct_sp}

这表明绝对量和状态占比回答不同问题。上游TYMP不能仅凭绝对丰度相关被解释为SPP1状态诱导。

### 2.5 删除共享基因保留方向，格块调整与研究集中度限制外推
{sensitive_table}

去掉LIPA/NPC1/NPC2/SORT1后，全片TYMP—两功能方向仍保留。4x4格块调整后的患者检验均未获支持：GSE175540有效21患者，TJ-RCC12患者。资格不足原因逐片保存，没有放宽阈值。留一格块方向汇总作为描述性稳定性检查，不把删块次数当样本数。

同文库敏感性：
{paired_sens}

留一研究删除GSE156632/两个单例研究后，肿瘤TYMP—肿瘤胆固醇rho分别约0.790/0.754/0.769。删除占10/16的GSE207493后剩6例，研究调整自由度不足，记不可估计；不能表述为四独立队列均验证。跨区室阴性以及巨噬SPP1—胆固醇的名义层级保持不变。

## 3 Discussion
最合适的机制联系是：TYMP相关肿瘤代谢背景涉及胆固醇处理，SPP1相关巨噬状态及其空间分布也具有相关功能程序。项目证据是细胞区室内和组织空间关联，未建立肿瘤TYMP到巨噬功能的稳定跨区室连接。
文献提供两类机制背景：TYMP缺陷及酶活恢复影响SREBP相关脂质代谢；TYMP缺乏引起核苷积累和溶酶体酸化/降解异常[1,2]。SPP1/OPN参与脂质摄取及溶酶体相关脂解调节，但不同环节作用方向不同[3]。这些研究支持功能层面的生物学合理性，不代表在本项目肾癌细胞中已经完成干预复现。
本次定向检验基于Phase256既有结果，不是新独立发现；两个空间程序同批RNA推算、SPP1状态同数据反卷积存在依赖。程序富集不等于实测通量/酶活，GLMF两例、区室队列研究集中及格块敏感性结果均应保留。原MTA主线不被此小节取代，也不将分段关联写成因果中介。

## 4 Validation与交付
FINAL_VALIDATION PASS。分析和验收退出0。核对{qc['spatial_correlations']}项空间相关、{qc['spatial_cohort_rows']}项队列汇总、{qc['paired_associations']}项配对关联；{qc['spatial_non_estimable']}项空间相关不可估计，完整保留原因。76项程序复算与Phase256原区域效果一致，1904次患者合并核对，35项配对矩阵偏相关复算，精确符号P、BH及置换前50个逐个重算通过。121项输入身份前后保持不变。
首次严格表达核对发现float32与既有float64结果约1.3e-7精度差；已改为原counts先转float64后计算，未降低容差。失败日志和修复前代码留档；首次失败未完成任何切片结果。
Python解释器为/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python，未新增依赖。执行bash code/run_analysis.sh，随后python code/validate.py与code/make_report.py；主要阶段和每片均经50GiB门禁，swap同步记录。无新图、无模型重训、旧RDS和既有报告未改；3G尚未启动。
交付：完整内部报告、独立小节草稿、主统计表、完整统计与复现ZIP。四份权威进度记录将更新于服务器，本地仅交付。

## References
{references}
"""
(P/'REPORT_CN.md').write_text(text,encoding='utf-8')
# Combined primary table keeps model and evidence type explicit.
a=primarysp.copy();a['analysis']='spatial_patient_direction';b=primarypair.copy();b['analysis']='same_library_patient_association'
pd.concat([a,b],ignore_index=True).to_csv(P/'tables/PRIMARY_EVIDENCE.tsv',sep='\t',index=False)
dest=P/'delivery'
for source,target in [('REPORT_CN.md','Phase257_full_report.md'),('SECTION_DRAFT_CN.md','Phase257_section_draft.md'),('tables/PRIMARY_EVIDENCE.tsv','Phase257_primary_statistics.tsv')]:shutil.copyfile(P/source,dest/target)
zpath=dest/'Phase257_statistics_reproducibility.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for folder in ['code','tables','inputs','qc']:
  for f in sorted((P/folder).rglob('*')):
   if f.is_file() and '__pycache__' not in str(f):z.write(f,f.relative_to(P))
 z.write(P/'scores/PAIRED_PATIENT_VALUES.tsv','scores/PAIRED_PATIENT_VALUES.tsv')
 z.write(P/'PLAN_CN.md','PLAN_CN.md')
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
manifest=[dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted(dest.iterdir()) if f.is_file()]
jwrite('DELIVERY_MANIFEST.json',manifest)
jwrite('STATUS.json',dict(status='COMPLETE_VALIDATED',time=stamp(),new_figures=0,old_report_modified=False,hotspot_started=False))
print(json.dumps(manifest,ensure_ascii=False,indent=2))
