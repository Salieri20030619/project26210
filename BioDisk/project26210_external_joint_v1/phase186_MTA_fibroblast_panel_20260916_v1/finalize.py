from pathlib import Path
import pandas as pd,json,hashlib,shutil,datetime,xml.etree.ElementTree as E
root=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');r=root/'phase186_MTA_fibroblast_panel_20260916_v1'
a=pd.read_csv(r/'PANEL_EFFECTS.tsv',sep='\t');v=json.loads((r/'ANALYSIS_VALIDATION.json').read_text())
p=a.pivot(index='gene',columns='method',values='delta_log2_median')
s=a[a.method=='median_ratio'].set_index('gene')
lines=['| 基因 | 主log2差 | 上四分位数法log2差 | 主分析组间分布 |','|---|---:|---:|---|']
for g in json.loads((r/'PRESET_PANEL.json').read_text())['genes']:
 label={'all_MTA_above':'所有处理文库高于所有对照','all_MTA_below':'所有处理文库低于所有对照','overlap':'范围重叠'}[s.loc[g,'complete_separation']]
 lines.append(f"| {g} | {p.loc[g,'median_ratio']:+.4f} | {p.loc[g,'upper_quartile']:+.4f} | {label} |")
table='\n'.join(lines)
report=f"""# Phase186：MTA成纤维扰动面板与SPP1中转依据

## 结论

**预设面板显示MTA处理后成纤维细胞IL6、HGF、WNT5A、PTGES表达升高；IL6具备继续衔接髓系SPP1的直接外部干预依据。** 本轮将IL6列为优先核查候选，不将其登记为已确认的MTA中介。旧MTA、A2B/STAT3及全部正负结果保留。

本轮取得完整计数、完成11基因描述性分析，并在原始文献中查到IL6中和降低髓系SPP1的直接实验。此前仅凭一般IL10响应猜测接收机制的范围已收窄为有明确两端依据的有限候选。

## Methods

数据来自[GSE249082](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249082)，HBF乳腺成纤维细胞系，MTA 10 µM处理24小时。原始featureCounts矩阵61544基因、8文库。C1_1/C2_1/C3_1/C4_2为对照，M1_1/M2_1/M3_1/M4_1为处理。保留原始列名；不将rep编号配成供体对子。逐GSM实验条件继承Phase185官方元数据。

主归一化：在全部8文库均为正计数的19116基因上计算几何均值参考，各文库取计数/参考中位数，再将因子几何均值标准化为1。敏感性：各文库正计数的75%分位数归一化，因子亦标准化为几何均值1。效应为处理组与对照组log2(归一化计数+1)中位数之差。

面板在计算前固定：HGF、TGFB1/2/3、IL6、CSF1、WNT5A、PTGS2、PTGES、SPP1，加MTAP清除背景。HGNC官方接口映射11个Ensembl ID，逐条归档。报告组间分布及两种归一化结果，未做P/FDR、差异基因筛选、患者签名训练或新图。每组删除一个文库的16种组合描述性敏感性另表保存。

## Results：IL6等输出相关转录本升高，两种归一化方向一致

{table}

IL6主效应+0.2424 log2单位，归一化计数中位数126.27→149.59，幅度较小。HGF和WNT5A同样上升，PTGES幅度最大。所有11基因的组间中位变化方向在两种归一化中一致。

关键反向和弱响应保留：TGFB3降低；TGFB1及MTAP变化小；PTGS2未升高。PTGES上升与PTGS2未升高并存，因此不能从RNA推断PGE2分泌增加。SPP1在该成纤维体系中计数低且组间范围重叠，不能作为巨噬SPP1上升证据。

文库计数总和约2505万–2976万，检出基因24485–25394；归一化全正参考基因的文库间相关最低0.9878。上述指标未见明显异常文库，但不能替代存活率、批次或独立培养验证。

## 文献核查：IL6有直接接收端依据，但MTA中转尚未闭合

### 1. 基质来源IL6中和可降低髓系SPP1

[van Baarle等，Nature Communications 2024](https://doi.org/10.1038/s41467-024-50438-2)，Fig.3g–i及Methods：肿瘤环境处理的肠神经胶质细胞上清含有较高IL6；上清作用于小鼠骨髓单核48小时，促进SPP1相关巨噬状态。以anti-IL6及磁珠清除上清IL6后，SPP1相关分化与Spp1/Arg1读出降低。该研究有蛋白测量和中和干预，超过单纯通讯预测。

这是**肠神经胶质—小鼠单核/巨噬**体系，不是乳腺成纤维或ccRCC，也没有MTA处理。可支撑“基质细胞IL6能够调节髓系SPP1”的生物学依据，不替代本项目的来源证明。原文GSE231804/GSE231709包含其他单细胞/批量实验，未将其直接登记为IL6中和转录组。

### 2. 直接IL6刺激与OPN反馈提供时序线索

[Yang等，Oncogenesis 2026](https://doi.org/10.1038/s41389-026-00596-3) 报告巨噬IL6刺激后24–48小时SPP1逐渐增加。THP-1条件培养液实验还显示，OPN中和降低接收细胞pSTAT3及SPP1，重组OPN补充可恢复相应读出。它提供“诱导与反馈维持可分开”的机制线索。

该文未测试MTA。pSTAT3与SPP1响应不是本项目STAT3转录代理的重新验证，也不凭此认定IL6—STAT3是本项目必经路径。其GSE311076为Rab37背景单细胞数据，不冒充直接IL6处理矩阵。此次精读正文，未独立重算其补图效应或重复数。

### 3. 其他候选保持各自证据等级

- WNT5A：本轮RNA方向支持，但[Qi等CRC研究](https://doi.org/10.1038/s41467-022-29366-6)中的WNT5A—FZD2主要是通讯推断；不能作为WNT5A直接诱导SPP1的实验证明。
- HGF：检索到的HGF—OPN直接诱导主要在上皮/肿瘤或成骨细胞；没有据此补成巨噬证据。
- PTGES/PGE2：本轮只测合成相关RNA。已知OPN促进巨噬COX2/PGE2的方向不能倒置成PGE2诱导SPP1。暂不升格。
- CSF1、TGFB2：上升方向保留但组间重叠；本轮不另行扩大搜索或筛选。

## 作者原始qPCR包核查

通过ZIP末尾范围请求解析170个中央目录条目，按成员位置定向提取原始文件。成功取得并通过CRC/长度校验的6份工作簿：Fig4C、Fig4E、Fig5D、S20、S21、S24。扫描其XML字符串、端点和条件标签，包含NOS2、ARG1、TNFa、CD163及内参，**未找到SPP1/OPN端点**。

Fig4F两次返回HTTP504，仍未取得；整包未完整下载，补充PDF未完整检查。不将“这6份文件没有SPP1”扩大为所有作者数据均无SPP1。工作簿只读，未改写、未重算客户湿实验。目录与范围请求、CRC及字符串清单保留在服务器。

## Discussion：当前可以走到哪一步

当前有两端分段依据：
**MTA处理 → 成纤维IL6 RNA升高；基质来源IL6 → 髓系SPP1响应。**

仍缺少相同体系中的分泌蛋白、MTA处理上清的IL6必要性，以及本项目ccRCC中对应来源/接收端。8个文库不等于8名供体；独立培养批次和供体身份未确认，所以本轮仅作探索性描述，不报告“4/4患者复现”或统计显著性。组内一致也不等于跨供体复现。不能仅凭两篇研究将跨体系证据写成已证实中介。

本轮缩小了后续核查范围，尚未确认TYMP—MTA—IL6—STAT3—SPP1完整路线。第一部分结论、整体巨噬STAT3—SPP1关联、C4定位及成熟巨噬IL10/STAT3干预的关键反向结果全部保留。

## 下一步固定任务

1. 在冻结75文库中先确认IL6的细胞来源；对成纤维与其他来源分别汇总，核查是否有足够独立患者，而不默认IL6由肿瘤或成纤维单独产生。
2. 在整体巨噬中核查IL6R/IL6ST与STAT3代理、SPP1的测量资格，复用Phase182已有受体结果；低RNA不解释为无蛋白，也不自动引入可溶性受体转信号。
3. 只有匹配患者/文库足够时，检验成纤维IL6—巨噬SPP1与STAT3的预设关联，报告效应、样本一致性及相应统计；不跨患者拼中介，不训练新评分追阳性。免疫与空间终点继续暂缓并保持分开。
"""
(r/'MTA_FIBROBLAST_PANEL_RESULTS_CN.md').write_text(report)
now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
validation={'time':now,'input':v,'independent_recalculation':json.loads((r/'INDEPENDENT_RECALCULATION.json').read_text()),'qpcr_members_validated':6,'Fig4F':'HTTP504','active_jobs':[],'inference':'descriptive only; donor independence unresolved'}
(r/'FINAL_VALIDATION.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2))
backup=r/'project_record_backups';backup.mkdir(exist_ok=True)
note=f"""## Phase186 MTA成纤维预设面板完成（{now}）
GSE249082完整61544基因8库；22/22范围块及gzip校验，HGNC11基因映射。主中位比/上四分位数11/11同向。IL6 +0.2424、HGF +0.2219、WNT5A +0.3264、PTGES +0.7529 log2；TGFB3 -0.3491、PTGS2 -0.0326，成纤维SPP1低计数重叠。无独立供体证明，无P/FDR；22效应独立复算通过，未删旧结果。
新直接依据：PMC11271635基质（肠神经胶质）IL6中和降低小鼠髓系SPP1；PMC12823564 IL6时序/OPN中和补充—pSTAT3/SPP1。两篇非MTA/ccRCC，IL6仅优先候选，不确认串联中介。WNT5A通讯推断不冒充干预。
Zenodo170目录条目、6份qPCR工作簿成员CRC通过，未见SPP1端点；Fig4F仍504，补充PDF未全查，整体ZIP未验证。
报告：{r}/MTA_FIBROBLAST_PANEL_RESULTS_CN.md
下一断点：冻结75库IL6来源与整体巨噬IL6R/IL6ST测量资格，复用Phase182；匹配独立样本够才作预设跨区室检验。无新图/评分，免疫空间暂缓且分开；无活跃任务。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']:
 pp=root/n;shutil.copy2(pp,backup/n)
 if n!='RESULTS_MANIFEST.tsv':pp.write_text(note+'\n'+pp.read_text())
 else:
  with open(pp,'a') as f:f.write('\n'+'\t'.join(['phase186_MTA_fibroblast_panel','phase186','report',str(r/'MTA_FIBROBLAST_PANEL_RESULTS_CN.md'),now,'PASS_WITH_SOURCE_BOUNDARIES','22 effect recalculations; 22 range chunks; 6 qPCR member CRC','descriptive only; Fig4F504; no complete causal chain'])+'\n')
d=r/'delivery';d.mkdir(exist_ok=True)
for n in ['MTA_FIBROBLAST_PANEL_RESULTS_CN.md','PANEL_EFFECTS.tsv','PANEL_SAMPLE_VALUES.tsv','LIBRARY_QC.tsv','LEAVE_ONE_LIBRARY_SENSITIVITY.tsv','FINAL_VALIDATION.json','QPCR_ENDPOINT_AUDIT.json']:
 shutil.copy2(r/n,d/n)
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in d.iterdir()}
(d/'SHA256.json').write_text(json.dumps(hashes,indent=2))
print('DELIVERY',len(hashes),'REPORT_CHARS',len(report));print(json.dumps(hashes,indent=2))
