from pathlib import Path
import json,re,copy,hashlib,shutil
ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
OLD=ROOT/'phase322_integrated_immune_report_20261005_v1'
SRC=OLD/'reader_revision_v2'
P=ROOT/'phase323_immune_manuscript_adaptation_20261005_v1'
P.mkdir(exist_ok=False)
for d in ['code','inputs','delivery','qc','logs','render','root_record_backups']: (P/d).mkdir()
for n in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:shutil.copy2(ROOT/n,P/'root_record_backups'/n)
old=json.loads((SRC/'inputs/DOCUMENT_CONTENT.json').read_text())
blocks=[]
def p(text):blocks.append({'type':'p','text':text})
def h(text,level=2):blocks.append({'type':'h','text':text,'level':level})
blocks.append({'type':'title','text':'TYMP与SPP1相关巨噬细胞的免疫特征'})
h('研究背景',1)
p('前述分析显示，肾癌恶性上皮细胞的TYMP相关转录状态伴随核苷、蛋氨酸及多胺代谢程序变化；配对组织代谢组提供了实测5′-甲硫腺苷（MTA）与SPP1表达之间的关联线索。髓系细胞分析进一步识别出SPP1富集、具有胆固醇处理和溶酶体特征的C4巨噬细胞状态。这些发现将研究问题从肿瘤代谢延伸至巨噬细胞的功能适应，需要进一步明确TYMP与SPP1各自对应何种免疫转录特征，以及这些特征与T细胞状态和组织分布的关系。')
p('SPP1相关巨噬细胞与成纤维细胞、细胞外基质及T细胞之间的联系，为上述问题提供了研究依据。结直肠癌单细胞与空间研究报道了SPP1阳性巨噬细胞与FAP阳性成纤维细胞的相互联系；实验研究进一步显示，TGF-β诱导的SPP1相关巨噬细胞可参与基质沉积和T细胞扩增调节。肾癌研究也报道了SPP1相关巨噬特征与免疫治疗结局的联系。[2–4] 因此，有必要在肾癌患者中同时考察巨噬细胞自身的免疫程序、同患者T细胞读数和局部基质环境。')
p('本部分首先检验TYMP与SPP1及巨噬细胞输出相关基因的联系，再在不预设免疫激活或抑制方向的条件下，比较两基因对应的干扰素反应、抗原呈递、吞噬及代谢程序。随后结合跨细胞关联、组织和空间读数、公开干预研究与临床队列，评价这些细胞内特征在组织和临床层面的对应。肿瘤细胞TYMP与髓系细胞内TYMP分别表征不同来源；C4用于定位SPP1富集状态，整体巨噬细胞用于评价更广泛的免疫功能联系。')
h('研究方法',1)
h('数据来源与独立生物学样本')
p('单细胞分析沿用前述5项研究、75个文库的整合集及冻结的细胞身份，保留C4原有成员。巨噬细胞内关联按每个基因对的表达资格确定患者：TYMP与SPP1整体巨噬分析为29名患者，C4内部分析为15名患者，其他基因对最多31名患者。原CD8接收端分析包含16名患者、19个同肿瘤文库和4项研究；Treg与NK端点各包含13名患者。')
p('外部单细胞数据采用Obradovic未治疗透明细胞肾癌队列的11名患者、13个CD45阳性文库，其中8名患者另有CD45阴性配套数据。该队列用于原单细胞集合之外的关联验证，并承载免疫程序探索及后续条件模型。HCRN GU16-260的13名透明细胞肾癌患者按基线8份与进展期6份文库分别分析，其中1名患者在两个时点重复取材。各分析按冻结资格使用有效患者，不合并治疗时点。')
p('组织转录组使用PRJNA1223824中经固定SNP面板核对供者身份、去除低信息观察样本并合并重复后的44个匿名供者。主要肾癌空间数据包括TJ-RCC的12名患者、GSE175540的23名患者24张切片，以及GLMF的2名患者；GLMF对应GSE210041中的两名患者。临床数据沿用全文的CheckMate和JAVELIN Renal 101队列。细胞和空间点用于形成样本读数，统计推断以独立患者、动物或独立混样池为单位。各端点报告其实际有效样本量。')
data_table=copy.deepcopy(next(b for b in old if b['type']=='table'));blocks.append(data_table)
h('巨噬细胞成员与免疫指标定义')
p('Obradovic分析使用两套已冻结的巨噬细胞成员。宽松集合纳入经复核的CD45阳性成熟巨噬簇2、3、5、6、8、9、10、15和17，共15,056个细胞，每名患者至少50个细胞；严格集合沿用逐细胞谱系复核后的巨噬名单，共13,065个细胞，每名患者至少100个细胞。原整合集的C4、CD8、常规CD4、Treg和NK身份沿用对应冻结版本，混合T/NK群单独保留。')
p('功能程序采用固定的50个Hallmark和6个GO基因集。6个GO程序分别为MHC-I抗原加工与呈递、MHC-II抗原加工与呈递、吞噬、淋巴细胞趋化、T细胞活化负调控和T细胞活化正调控。所有程序均移除TYMP和SPP1后计分。原始计数按细胞全基因总计数归一化至10,000，再作log1p变换；程序分数为可测成员基因归一化表达的等权均值。[1]')
p('纳入基因要求在各贡献文库中可测且在患者内至少10个细胞中非零。程序要求至少80%的成员可测，且至少20%的完整程序成员在不少于max〔3，向上取整（患者细胞数×1%）〕个细胞中检出。模型另要求残差维度至少8，目标基因及程序残差均具有变异。上述资格分别应用于宽松和严格集合。')
h('患者内关联与患者层统计汇总')
p('对TYMP或SPP1表达与目标程序分数分别取秩，并在每位患者内调整总计数秩、检出基因数秩、文库和冻结髓系亚群分类变量。两组调整后残差的Pearson相关系数作为该患者的秩残差相关ρ。跨患者报告ρ的中位数、正负方向及同向人数。免疫程序作为结局时采用上述核心模型；干扰素程序在这一分析中作为结局保留。')
p('TYMP—SPP1、TYMP—CTSB、TYMP—TGFB1、TYMP—RAB37、SPP1—CTSB和SPP1—RAB37六个基因对另设环境调整模型，在核心模型中加入干扰素α、干扰素γ和缺氧程序；这些调整程序先移除六个基因对涉及的全部目标基因。基因对与免疫程序分别检验。')
p('Obradovic患者内关联在至少10名具有非零有效效应的独立患者中，采用双侧精确二项符号检验，原假设为正负方向等概率；同时保留患者效应中位数。小于10名患者的该类分析仅描述方向和效应。原整合集及其他数据层沿用各自已实施的患者汇总与检验方法，完整模型和检验名称列于统计资料。')
h('相互条件调整与程序基因敏感性分析')
p('在既有10条关系中增加另一基因的表达秩作为协变量。具体包括TYMP和SPP1各自与糖酵解、干扰素α、干扰素γ的6条关系，TYMP与MHC-I及吞噬的2条关系，SPP1与MHC-II的1条关系，以及TYMP与总体炎症的1条关系。宽松和严格集合各计算一次，共20项条件模型。条件核查在取秩前将基因和程序数值统一至12位小数。')
p('程序敏感性分析分别移除糖酵解与所选免疫程序之间的共享成员，以及免疫程序彼此间的共享成员；删除后至少保留5个基因，并重新应用可测与检出资格。另对每个原程序成员逐一删除，记录患者内效应、汇总中位数和方向变化。依据宽松集合中单基因删除对原方向效应的衰减程度选定5个基因，合并删除后计算效应；严格集合使用同一名单。逐基因及联合5基因删除均作为描述性敏感性分析，不另计算基因筛选P值。')
h('跨细胞与组织免疫关联')
p('原单细胞集合在同一批合格肿瘤文库内配对巨噬来源指标和CD8、Treg或NK接收端读数，再汇总到患者。采用研究调整的偏Spearman相关、99,999次研究内置换和5,000次患者分层bootstrap。Obradovic的跨患者模型调整作者处理批次和pT分层；HCRN按治疗时点分别计算。肿瘤来源TYMP另按配套患者连接巨噬读数，来源和资格单列。')
p('CD8细胞毒程序由NKG7、GNLY、GZMB、PRF1和CTSW组成，功能障碍相关程序由PDCD1、CTLA4、LAG3、HAVCR2、TIGIT、TOX、CXCL13和ENTPD1组成。Treg比例以常规CD4与Treg之和为分母。探索性CD8程序采用全部50个Hallmark，CD8与合格巨噬总数中的CD8捕获比例另作组成读数。定向通信沿用冻结CellChatDB的9个SPP1受体映射，受体复合体要求在同一接收细胞中联合检出；相关定向程序和通信结果列入补充资料。')
p('组织表达采用log2（TPM或FPKM＋1），文献五基因评分为SPP1、FN1、CLEC5A、OLR1和IL1RN的表达均值。[4] 在44个供者中计算SPP1及五基因评分与组织免疫程序的Spearman相关，使用100,000次置换及5,000次供者bootstrap。SPP1的两个主要端点和五基因评分的两个对应端点分别校正。')
h('空间免疫分布与外部干预资料')
p('空间数据沿用全文的切片、细胞状态读数和组织坐标。界面分析将恶性细胞映射比例不低于0.20的区域定义为计算肿瘤表达区，采用400 μm为主要比较宽度，200及800 μm为敏感性宽度，分别计算内外总T及CD8读数，并评价其与局部SPP1相关状态的联系。深度、总巨噬及细胞组成调整按既定模型分别报告；先计算切片效应，同一患者的多切片再按原方案汇总。')
p('SPP1相关巨噬与成纤维、细胞外基质及胶原程序的联系采用预先固定的公开标记集合。在去除两状态共享基因后计算深度调整相关，并从基质程序中移除状态定义基因后再评价功能联系。SPP1单基因、旧空间SPP1相关状态和本项目C4保持各自定义，按队列分别汇总患者效应。高分辨率案例的来源定位与几何分析列入补充资料。')
p('外部功能资料按原论文的物种、肿瘤模型、干预对象及测量终点整理。可获得的源表按动物或独立混样池汇总，多个组织视野归属于相应动物。CRC宿主Spp1敲除的CAF计数采用原文线性模型复算，CD8中Ki67比例采用动物随机截距混合模型；CD8数量报告组均值和相对变化。肾癌离体阻断中IFNγ和GZMB结果引用原论文实验。源表复算、原始计数分析与文献引用分别记录。[3,4]')
h('临床转录评分与生存分析')
p('临床部分与全文的临床队列分析衔接。CheckMate包含311名RNA匹配患者，其中nivolumab组181名、everolimus组130名；JAVELIN Renal 101包含726名患者，其中avelumab联合axitinib组354名、sunitinib组372名。CheckMate文献五基因评分直接取发布标准化表达的均值，JAVELIN取log2（TPM＋1）均值。连续评分按治疗臂内每1个标准差缩放，文献分组采用完整队列固定中位数。')
p('Cox模型分析OS与PFS。CheckMate调整年龄、统一后的性别、试验及MSKCC分层；JAVELIN调整年龄和性别，扩展模型按原方案加入其他评分或总巨噬估计。使用各模型完整病例，治疗交互在CM025或JAVELIN各自随机两臂内检验。CheckMate使用原PFS和OS事件编码，JAVELIN的PFS事件为1−PFS_P_CNSR。')
p('C4转录评分沿用24名患者、5项研究的C4与其他成熟巨噬配对伪bulk发现结果，在临床数据中使用共同覆盖的固定19个基因；基因选择排除TYMP、SPP1及文献五基因，留一研究验证在各训练折内选基因。针对nivolumab PFS的非比例风险，采用β(t)＝β0＋β1log(1＋t)的时间变化系数模型，t以月计，以两自由度似然比检验C4增量关联，并报告预设时点HR及区间。')
h('多重检验与结果报告')
p('各分析报告独立生物学样本量、效应量、方向一致性、原始P值及对应FDR。原始P＜0.05作为名义显著标准，BH校正FDR＜0.05称为校正后显著。巨噬免疫发现的56个程序×2个目标基因×2套成员共224项属于同一家族；6个基因对×2种模型×2套成员共24项另为一家族；相互条件模型20项与去共享基因36项分别校正。上述家族各自保留，不合并为一个校正池。')
p('跨细胞、空间、组织和临床分析按其原先定义的同类假设家族分别校正。多因素模型使用目标基因、评分或评分与治疗交互项的P值。样本不足的描述性结果、不可估计端点及数据内选择的敏感性检查均保留相应记录。完整基因集、细胞和患者纳入名单、模型协变量、检验家族及全部P值与FDR见配套统计资料。')
# Preserve all results and discussion paragraphs and tables exactly, except small formal prose edits.
start=next(i for i,b in enumerate(old) if b['type']=='h' and b['text']=='2 结果')
end=next(i for i,b in enumerate(old) if b['type']=='h' and b['text']=='4 结论')
body=copy.deepcopy(old[start:end])
for b in body:
 if b['type']=='h':
  b['text']=re.sub(r'^\d+(?:\.\d+)?\s+','',b['text'])
  if b['text']=='结果':b['text']='分析结果'
 if b['type']=='p':
  b['text']=b['text'].replace('这项研究最值得保留的发现，是TYMP和SPP1虽然在巨噬细胞中正相关，却对应不同的免疫特征。','患者层分析显示，TYMP和SPP1在巨噬细胞中正相关，但对应不同的免疫转录特征。')
blocks.extend(body)
h('参考文献',1)
for b in old:
 if b['type']=='p' and re.match(r'^\[\d+\]',b['text']):blocks.append(copy.deepcopy(b))
for b in blocks:
 if b['type']=='table':b['title']=re.sub(r'^表(\d+)',r'表免疫\1',b['title'])
# Remove production-language from discussion without altering scientific content.
for b in blocks:
 if b['type']=='p':
  b['text']=b['text'].replace('报告可以明确讨论','这些结果支持讨论').replace('总体上，报告的核心结论','总体上，本部分的核心结论')
(P/'inputs/DOCUMENT_CONTENT.json').write_text(json.dumps(blocks,ensure_ascii=False,indent=2))
builder=(SRC/'code/build_docx.py').read_text()
builder=builder.replace("str(P.parent/'code/vendor')",repr(str(OLD/'code/vendor')))
builder=builder.replace("runpy.run_path(str(P/'code/report_content.py'))","")
builder=builder.replace("header.text='26210  TYMP and SPP1 immune analysis'","header.text=''")
# no forced section page breaks: keep chapter-ready flow
builder=builder.replace("z.paragraph_format.page_break_before=True","z.paragraph_format.page_break_before=False")
builder=re.sub(r"out=P/'delivery/[^']+';doc.save\(out\)","out=P/'delivery/26210_免疫章节_全文适配_v3.docx';doc.save(out)",builder)
(P/'code/build_docx.py').write_text(builder)
guide='''26210免疫章节全文适配说明

适配基准：Phase273正式全文与Phase322免疫报告修订版v2。
本次交付为免疫章节稿，尚未合并或覆盖45页完整报告。

1. 背景
新增背景接续恶性上皮TYMP相关代谢、实测MTA和SPP1富集巨噬状态，提出免疫程序及跨细胞问题。合并时将这三段置于全文研究设计的引入部分，与1.1科学问题衔接；不重复独立摘要。

2. 方法
“数据来源与独立生物学样本”补入全文1.2；巨噬成员、患者内关联、相互控制和删基因方法接1.10；空间补充并入1.12与1.15；临床段与原1.11合并去重；最后一节并入原1.19统计方法，保留原文其他检验方法。
本次小节采用Word标题层级，不预占全文编号，避免与1.1至1.19冲突。正式合并时统一更新编号、目录及交叉引用。

3. 结果和讨论
四项细胞内及组织免疫结果适合接原2.12之后；SPP1与基质/空间免疫段与原空间2.17至2.18合并，临床段与原2.13合并。均不可简单重复粘贴。
讨论分别融入原3.4、3.7、3.8和3.9；不另立第二份全文结论。
六项结果、七张表及影响核心结论的阴性和敏感性结果均保留。表采用“表免疫1至7”临时标识；合并时按最终位置统一编号。

4. 来源与统计
原结果的样本、效应、P值、FDR和检验家族不变；没有新增生物学计算。
参考文献保留本章节1至9的局部标识；正式并入时与全文已有文献去重，并同步正文引用。
完整统计继续使用Phase322“26210_TYMP_SPP1_完整统计与来源.zip”。

5. 版本
Phase273全文、Phase322原版及v2全部保留。当前背景和方法适配稿为本v3。
'''
(P/'delivery/全文衔接说明.txt').write_text(guide)
source_files=[SRC/'inputs/DOCUMENT_CONTENT.json', SRC/'delivery/26210_TYMP_SPP1_v2_免疫综合分析报告.docx',ROOT/'phase273_clinical_suppfig_wording_20260928_v1/report/report_materialized.json',ROOT/'phase314_TYMP_participation_immune_exploration_20261004_v1/resources/ANALYSIS_PROTOCOL.json',ROOT/'phase321_TYMP_SPP1_conditional_gene_robustness_20261005_v1/resources/ANALYSIS_PROTOCOL.json']
(P/'inputs/SOURCE_MANIFEST.json').write_text(json.dumps([{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in source_files],ensure_ascii=False,indent=2))
old_tables=[b for b in old if b['type']=='table'];new_tables=[b for b in blocks if b['type']=='table']
checks={'seven_tables_exact_values':len(new_tables)==7 and all(a['rows']==b['rows'] and a['headers']==b['headers'] for a,b in zip(old_tables,new_tables)),
'result_paragraphs_unchanged':all(b in blocks for b in old[start:next(i for i,b in enumerate(old) if b['type']=='h' and b['text']=='3 讨论')] if b['type']=='p'),
'no_abstract_or_separate_conclusion':not any(b.get('text') in ['摘要','4 结论'] for b in blocks),
'background_present':any(b.get('text')=='研究背景' for b in blocks)}
assert all(checks.values()),checks
(P/'qc/CONTENT_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':checks,'no_new_models':True,'not_merged_into_full_manuscript':True},ensure_ascii=False,indent=2))
print('PREPARED',P,'blocks',len(blocks))

