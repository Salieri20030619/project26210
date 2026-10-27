from pathlib import Path
import json,hashlib,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase164_SpaMTP_method_review_20260915_v1';O.mkdir(exist_ok=True)
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
review=dict(user_instruction='First assess feasibility; list a plan if feasible; do not implement now',status='FEASIBILITY_AND_PLAN_ONLY',client_pdf='C:/Users/31937/Downloads/SpaMTP.pdf',client_pdf_sha256='0890df980b57c90acec0fc5f217833e261e5da84fb72f7e577611ac0a7659fe2',pdf_review='6/6 pages visually inspected; Chinese introduction printout, many original figure placeholders not loaded, page6 blank; primary paper and official documentation used to verify methods',primary_paper='https://www.nature.com/articles/s41592-026-03140-8',official_repository='https://github.com/GenomicsMachineLearning/SpaMTP',official_docs='https://genomicsmachinelearning.github.io/SpaMTP/reference/index.html',docs_version_observed='1.1.0.9000; not installed or pinned for execution',available=dict(single_cell='Frozen Part1 malignant and Part2 macrophage/C4 reference',spatial_RNA='Existing GSE175540 24 sections (12 FFPE/12 frozen) and legacy Cell2location mappings; donor identity audit pending',bulk_metabolomics='Four public paired bulk RNA-metabolite cohorts,255 tumors; not same-coordinate spatial data',client_experiment='786-O figure only, no spatial raw data',paired_spatial_metabolomics='Not verified among the reviewed project inputs; no usable paired SM intensity+coordinate input currently established'),implemented=False,new_models_run=False,new_spatial_statistics=False,new_figures=False,installed_packages=False)
(O/'SOURCE_AND_FEASIBILITY_REVIEW.json').write_text(json.dumps(review,ensure_ascii=False,indent=2))
plan="""# SpaMTP方法可行性评估与待确认计划

**当前状态：只完成评估和计划，未实施。**

## 评估结论

建议优先采用方案A：借鉴SpaMTP的区域通路分析与空间数据组织思想，在现有空间RNA上开展代谢相关细胞状态和邻域联系的适配分析。它能补充“这些细胞状态出现在哪里、是否相邻”，不能提供空间实测MTA或确认MTA转运。

方案B为完整SpaMTP空间代谢组+转录组联合流程。它需要同一切片或可配准相邻切片的SM强度矩阵、空间坐标以及ST计数。当前审阅的项目输入尚未建立这一配对条件，因此B仅列为条件方案，不启动。

用户明确要求先评估、列计划、不立刻实施。本阶段没有安装SpaMTP、运行新模型、重算空间统计、重新投射C4或生成新分析图；现有分析结果保持原版本。

## 资料核查

客户提供的SpaMTP.pdf为6页图片型中文介绍稿。已逐页查看，部分原图在打印稿中为未加载占位框，第6页为空白。对应原始论文为Causer等，Nature Methods 2026，DOI:10.1038/s41592-026-03140-8，发表于2026-07-31。论文和官方文档支持的核心功能包括空间代谢物注释、空间对齐、区域/通路分析及多模态整合。[原始论文](https://www.nature.com/articles/s41592-026-03140-8)

官方包建立在Seurat/Cardinal框架上，分别存放空间RNA与空间代谢测量。现有RNA程序分数、scFEA预测或bulk MTA数值不能自动成为空间代谢实测层。[官方教程](https://genomicsmachinelearning.github.io/SpaMTP/articles/SpaMTP.html)

官方函数索引列有MapSpatialOmics、ApplySpatialAlignment、RunSpatialGraphPCA、CreatePathwayAssay、FindAllDEMs等。后续若实施，应锁定代码提交和环境，不把实时文档版本当作已经验证的安装版本。[官方函数索引](https://genomicsmachinelearning.github.io/SpaMTP/reference/index.html)

## 现有数据是否适配

| 输入 | 已知状态 | 可用于什么 | 当前不能替代什么 |
|---|---|---|---|
| 第一部分恶性细胞及第二部分巨噬/C4 | 已冻结、可追溯 | 候选基因与参考状态定义 | 空间代谢物测量 |
| GSE175540空间RNA | 已有24切片、分FFPE/冷冻各12 | 区域/邻域分析的基础 | 24独立患者；需核对患者映射 |
| 旧空间Cell2location参考 | 旧SPP1-state及其他状态 | 旧状态与总巨噬定位的参考 | 新C4；不能直接改名 |
| 公共255患者bulk RNA–代谢组 | 已有实测MTA等 | 独立组织层关联对照 | 同一空间切片的MTA分布 |
| 客户786-O湿实验图 | 仅作实验背景 | 解释假说与对应方向 | SM原始矩阵、空间坐标或精确代谢物注释 |
| 配对空间代谢组 | 在已审阅输入中尚未确认可用 | 若后续获得，可进入方案B | 不能由RNA或bulk值伪造 |

此次只核对已有文档和已知输入类型，没有全盘扫描服务器并据此宣称所有地方都没有SM数据。B的缺口是尚未建立可复用的配对输入，不是算法本身不可用。

## 借鉴范围及必须调整的环节

| SpaMTP思想/功能 | 本项目处理 |
|---|---|
| 统一空间坐标与区域组织 | A使用现有ST坐标；B在有SM后核验配准 |
| 区域通路分析 | A按RNA解释；合成、SAM供给、MTAP清除分别评分，排除TYMP/SPP1机械重叠 |
| GraphPCA空间结构表示 | 可列为敏感性，主要关联使用未平滑表达或预定区域汇总，防止平滑制造相关 |
| WNN联合整合 | 仅B有真实两模态输入时考虑；不把RNA预测分数包装成第二组学 |
| 代谢物注释 | 仅有m/z/强度时开展；保留多重候选、加合物和正交证据，不将质量匹配等同确证MTA |
| 随机像素池差异分析 | 不将像素池当独立生物重复；改用患者×区域层级的摘要和推断 |

原文的差异分析方案包含随机像素分池、默认3池。这可作算法内部探索，但本项目不能据此得到3名独立患者的统计证据，必须遵守既定患者单位要求。[原文Methods](https://www.nature.com/articles/s41592-026-03140-8)

此外，若用RNA通路一致性帮助选择代谢物注释，后续不能再把同一RNA关联当成该注释的独立验证。注释与结局检验应分开记录。

## 方案A：现有RNA层空间定位适配计划（推荐，待确认）

### A1 数据与身份审计

核对24切片的患者对应、制备方式、原始counts、坐标、组织图和Cell2location版本；同一患者多切片不增加独立n。确认旧SPP1-state与新C4的关系，保留旧标签。若无法建立新C4可靠投射，本轮空间主分析仍使用总巨噬和明确标注的旧状态，新C4只保留现有单细胞定位证据。

输出：患者/切片清单、输入哈希、参考版本对照、可分析人群及缺失项。此阶段只以身份和可复核性验收，不按关联强弱选择切片。

### A2 冻结空间问题与候选

主要问题：TYMP/SRM/SMS相关恶性区域，是否邻近SPP1相关巨噬生态位？

候选分别为：
- 合成端：ODC1、SRM、SMS；
- SAM供给：MAT2A、AMD1；
- 清除背景：MTAP；
- 并行候选：SLC16A3；
- 接收区域：整体巨噬SPP1相关状态；成纤维区域作为预列中转候选，C4作为条件可投射的次要定位。

保留单基因和既定模块，不按空间P值修改基因集。不能把ARG1/VEGFA/IL10/TIMP1组合直接命名为已验证的人肾癌MTA响应；Phase163的接收端阴性证据保留。

输出：预列假设、终点、协变量、检验家族及方向表。

### A3 少量切片技术试点

批准后按数据完整性和组织质量选择少量代表切片，兼顾制备类型；不看结局P值选样。核验坐标、区域边界、候选检出与空间尺度。该试点用于技术检查，不用少量切片宣称临床或患者级显著。

若计划投射新C4，必须先证明参考可区分并记录不确定性，不能只把旧SPP1-state标签替换为C4。若不满足，继续总巨噬主分析并保留C4单细胞结果。

输出：试点QC、覆盖和方向描述，不生成新分析图，除非后续另行确认。

### A4 预定区域及邻域比较

在通过技术检查的范围内，分别计算同区域与非自身邻域关系；明确来源区域和邻域的细胞组成协变量。主要尺度沿已有平台几何及预定邻域定义，较近/较远尺度只作预列敏感性，不挑最显著半径。

SPP1也在恶性细胞中表达，因此不能把spot总SPP1直接称为巨噬来源。采用巨噬富集区域与恶性/巨噬/成纤维丰度调整，并将其标为区域层代理。结果仅支持空间相容性，不自动证明细胞通信、MTA释放或串联中介。

输出：每患者/切片的效应、来源解释、预列敏感性和完整正负结果。

### A5 患者层汇总及旧结果对照

先在同患者内汇总多切片，再进行患者级统计；FFPE与冷冻分别呈现，不能把两种制备自动称作独立外部队列。有效n<10按探索性呈现；n≥10依既定规则报告效应、P与同类家族FDR。无有效重复时仅报告方向与效应。

控制空间自相关的置换/重抽样设计须在执行前固定，不能仅随机打散spot来获得更小P。区域和邻域终点、基因家族及外部验证分别校正。

与Phase75/78旧结果并列比较：旧分析未得到双制备支持。本方案的新增内容是合成/清除分离、来源解释及可核验的新状态定位，不是仅换算法再挑阳性结果。

输出：患者级统计表、敏感性、旧新口径对照、可保留正文结论。

### A6 验收后决定是否扩展

验收关注数据身份、统计单位、空间尺度、来源解释和重现性，不以P/FDR通过作为能否交付或保留的硬门槛。只有方案执行获用户确认后才进入A1–A5；绘图另按已约定的“先审文字统计”规则处理。

## 方案B：完整SpaMTP联合流程（条件计划，暂不启动）

1. 查找或确认配对SM–ST：同切片最佳，相邻切片需有组织及配准信息；同时记录独立患者数、测量平台及原始m/z矩阵。
2. 审核MTA是否确实可测：校验注释、极性、加合物、质量误差及可用MS/MS/标准品证据。未检测或注释不确定时，不用邻近峰或其他代谢物代替MTA。
3. 在隔离环境锁定SpaMTP及依赖后，进行导入、QC、归一化和空间配准；不修改共享环境。
4. 在真实两模态上评估联合区域/通路与单模态结果的一致性；联合聚类作为表征，不自动改写冻结C4。
5. 以患者×区域进行差异、邻域及多组学对应，保留模型不确定性与无效配对。
6. 如只找到脑/胶质瘤等其他癌种数据，可作方法试验或机制背景，不称为ccRCC独立验证。

条件B需要的新输入目前未建立，不能把公共bulk MTA复制到每个spot，也不能把scFEA预测量当空间实测强度填入SPM层。

## 执行顺序与当前停点

建议先审定A。批准后先做A1/A2和技术试点，确认数据与参考可用，再扩展患者层分析；B保留为有合适实测输入时的升级路线。若形成长任务，沿用服务器内存门槛和15分钟shell级完成等待规则。

当前停在**方案评估完成、等待用户决定**。无SpaMTP安装，无新模型，无新空间检验，无C4重映射，无新分析图。
"""
(O/'SpaMTP_FEASIBILITY_AND_PLAN_CN.md').write_text(plan,encoding='utf-8')
summary="""# SpaMTP方案评估（未实施）

建议方案A：借鉴区域通路与空间组织思想，在现有空间RNA中检查TYMP/SRM/SMS恶性区域与SPP1相关巨噬生态位的邻近关系；MTAP清除与成纤维中转单列，C4保持定位用途。该方案属于RNA层适配，不能称空间实测MTA分析。

完整方案B需要配对空间代谢组强度+坐标+空间RNA。现有已登记24张GSE175540切片为RNA空间输入，公共255患者MTA为bulk，客户为细胞系实验图；当前尚未建立可用的配对SM–ST输入，因此B只列条件计划。

批准后拟按：患者/切片与旧新参考审计 → 冻结候选和终点 → 少量切片技术试点 → 预定区域/邻域分析 → 患者级汇总与旧结果对照 → 验收后决定扩展。

关键改造：24切片不等于24患者；旧SPP1-state不能直接改名C4；spot总SPP1不直接归为巨噬；随机3个像素池不增加生物学n；RNA代理/预测MTA不能充当第二实测组学。FDR不作是否保留或继续的硬门槛，正负及旧分析未支持结果均保留。

用户明确要求“先评估列计划，不立刻实施”。当前仅完成文献和输入说明核查，未安装、未运行新模型、未重新投射C4、未生成新分析图。
"""
(O/'SUMMARY_CN.md').write_text(summary,encoding='utf-8')
status=dict(status='COMPLETE_FEASIBILITY_AND_PLAN_ONLY',user_decision_required_before_execution=True,implementation_started=False,model_runs=0,new_spatial_tests=0,new_analysis_figures=0,source_pdf_pages_reviewed=6,plan=str(O/'SpaMTP_FEASIBILITY_AND_PLAN_CN.md'),time=now)
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
note="""## 2026-09-15 Phase164 SpaMTP仅方案评估完成，未实施

用户提供Downloads/SpaMTP.pdf并明确“先评估方案，如果可行列出计划，不要立刻实施”。已逐页查看6页图片型介绍稿并追溯原始Nat Methods论文s41592-026-03140-8及官方文档。方案phase164_SpaMTP_method_review_20260915_v1/SpaMTP_FEASIBILITY_AND_PLAN_CN.md。
结论A可作为RNA层空间定位适配：合成SRM/SMS、MTAP清除、整体巨噬SPP1邻域与成纤维中转，C4仍为定位。完整B需配对SM强度/坐标/ST；现有已审阅输入尚未建立，不把bulk MTA或RNA预测当空间实测。
已知空间GSE175540共24切片12FFPE/12frozen；患者对应仍需执行前审计，不能当24独立患者。旧Phase17/75/78SPP1-state不能直接改名新C4；spot总SPP1不能直接归巨噬。随机像素3池不增加生物学n。借鉴区域通路/对齐思想，不照搬伪重复或用换算法掩盖旧阴性。
当前完成的是评估及计划文件：未安装SpaMTP、未跑新模型/空间统计、未投射C4、未生成新分析图。等待用户确认后才实施；此前Phase163结果已完成并保留。
"""
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=P/name;shutil.copy2(p,O/(name+'.before'));p.write_text(note+p.read_text(encoding='utf-8'),encoding='utf-8')
with (P/'RESULTS_MANIFEST.tsv').open('a',encoding='utf-8') as f:f.write('\t'.join(['phase164_SpaMTP_method_review_20260915_v1','Phase164','plan',str(O/'SpaMTP_FEASIBILITY_AND_PLAN_CN.md'),now,'REVIEW_COMPLETE_NOT_EXECUTED','6 source PDF pages reviewed; primary paper and official docs; existing spatial preregistrations','User explicitly requires plan review before implementation'])+'\n')
# Preserve completed analysis stage; add a pending plan pointer only.
p=P/'part2_rebuild_v1/STATUS.json';shutil.copy2(p,O/'PART2_STATUS.before.json');b=json.loads(p.read_text());b.update(latest_method_review='Phase164',proposed_next_plan=str(O/'SpaMTP_FEASIBILITY_AND_PLAN_CN.md'),proposed_plan_status='evaluated_not_executed_user_review_pending');p.write_text(json.dumps(b,ensure_ascii=False,indent=2))
manifest={x:hashlib.sha256((O/x).read_bytes()).hexdigest() for x in ['SUMMARY_CN.md','SpaMTP_FEASIBILITY_AND_PLAN_CN.md']};(O/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest));print('PLAN_ONLY_COMPLETE')

