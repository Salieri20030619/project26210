from pathlib import Path
import json,copy,hashlib
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');old=r/'phase189_integrated_TYMP_MTA_macrophage_report_20260917_v1';d=r/'phase190_mechanism_focused_report_20260917_v1'
for sub in ['report','scripts','qc','figures','supplement']:(d/sub).mkdir(parents=True,exist_ok=True)
a=json.loads((old/'report/report_materialized.json').read_text());b=copy.deepcopy(a)
def p(t):return {'type':'p','text':t}
def h(t):return {'type':'h2','text':t}
edits={
142:'组织TYMP与MTA的未调整关联未获支持；SMS及多胺短路径的间接效应区间跨零（表7）。与直接TYMP关系相比，公共数据对MTA的支持主要来自合成与清除背景的分段联系。',
147:'组织层联系、直接髓系扰动与本项目细胞状态发现分别列示，实际实验对象和读出见表中对应列。',
148:'组织层MTA-SPP1关联与直接髓系作用依据共同指向后续问题：本项目中哪些巨噬细胞状态集中表达SPP1，其代谢与细胞内处理程序具有何种特征？下述分析通过独立分群和患者比较进行定位。',
154:'C4内TYROBP检出100%，CD68约98.9%，C1QC约98.4%，GPNMB约97.6%，SPP1约47.8%。SPP1相对富集与广泛存在的巨噬骨架共同构成该群身份；C1QC仍广泛表达。谱系审核中识别的上皮/髓系混合候选已排除出成熟巨噬比较。',
155:'C4核心在不同种子及重抽样中的Jaccard分别为0.8574和0.7589，均达到预设候选核心门槛，与多基因身份及患者一致性相互支持。',
157:'排除TYMP和SPP1后，C4的KEGG胆固醇代谢和溶酶体程序仍在26/26名患者升高，中位评分差分别0.62468和0.34919，P均8.80e-6、原FDR均2.00e-5。抗原呈递程序在25/26名患者降低；炎症和干扰素相关程序整体较低。吞噬体程序为14/26正向、P=0.54216。',
159:'功能程序以独立患者配对比较；SPP1贡献采用原始counts分解。同文库16患者与C4配对26患者分别使用其实际分母。',
160:'同文库16名患者中，C4占成熟巨噬的比例中位数为7.4%，贡献捕获SPP1原始RNA的中位54.6%。C4因而构成SPP1信号集中的细胞状态，而整体成熟巨噬保留了该信号在其他状态中的分布。两类读出分别用于描述状态定位和总体表达联系。',
164:'主模型具有有效拟时序的C4为442/2674细胞（16.53%），敏感性模型为1743/2674（65.18%）。各模型使用其实际有效患者集合。',
165:'可估计轨迹中的C4位于单核之后，与其他成熟巨噬的相对排序没有统一方向。图7的C1-C4为恶性细胞基因曲线显示簇，本节C4为冻结巨噬候选，二者采用不同命名体系。',
166:'上述结果将SPP1相关信号定位于具有脂质处理及溶酶体表达特征的巨噬状态，并将其与TYMP相关代谢发现及MTA候选联系纳入同一研究框架。',
}
for i,t in edits.items():b[i]['text']=t
discussion=[
{'type':'h1','text':'3 讨论'},
h('3.1 TYMP相关代谢状态将核苷与蛋氨酸多胺环节联系起来'),
p('TYMP的直接生化功能位于核苷代谢，可催化胸苷与无机磷酸生成胸腺嘧啶和2-脱氧核糖-1-磷酸的可逆反应[7]。本项目发现的TYMP相关恶性上皮状态还伴随蛋氨酸及部分多胺程序改变，说明后续研究需要从单一酶反应扩展到该细胞状态的代谢组织方式。F1、通路评分以及ODC1/SRM表达提供了不同层面的对应，使“核苷相关状态与氨基酸/多胺背景共同变化”成为有数据支撑的解释。'),
p('这一联系中，蛋氨酸/SAM供给与多胺合成承担不同环节。SAM相关供给为氨丙基转移提供代谢背景，SRM/SMS参与的合成反应则直接伴随MTA生成。因而MTA并非由某一个相关基因临时引出的标志物，而是位于本项目所涉及代谢环节交汇处的具体产物。讨论重点应放在这些环节如何构成候选代谢背景，而不是将全部代谢酶视为必须同时升高的一组指标。'),
h('3.2 MTA的生成清除与局部暴露共同影响髓系响应'),
p('MTA在本研究中的机制意义包含两个层次：一是多胺合成的产物，二是具有直接髓系作用依据的候选代谢介质。SMS与实测MTA的正向联系及MTAP与MTA的负向背景分别对应生成和清除。稳态水平取决于这些过程的共同作用，因此解释MTA时需要同时考虑生成、摄取与分解，而不是仅依据一项酶RNA的升降。'),
p('细胞内、胞外和组织平均水平还反映不同的代谢区室。Barekatain等在胶质母细胞瘤模型中观察到MTAP缺失细胞释放MTA，而MTAP表达基质细胞可代谢胞外MTA[6]。该结果提供了一个与本研究直接相关的解释：周围细胞的清除能力会改变局部暴露，使培养体系、细胞内测量与组织平均读数出现差异。这里引用的是区室与代谢交换原理。'),
p('MTA转向髓系研究的依据来自直接扰动。Hansen等通过MTA处理及受体/STAT3干预，显示髓系表型响应具有信号依赖[4]；Keyel等显示MTA与TLR背景共同改变炎症和共刺激响应[5]。两类研究说明，MTA的生物学意义不止于代谢产物标记，其效应还取决于接收细胞和刺激环境。由此，本项目进一步考察具有明确分子特征的巨噬状态，而不是仅比较全组织中的单个标志物。'),
h('3.3 SPP1富集与脂质处理溶酶体程序共同定位巨噬状态'),
p('C4的解释价值来自多基因身份与程序的一致组合。GPNMB、CCL18、CYP27A1、LPL和APOE富集，与去除TYMP/SPP1后仍保留的胆固醇及溶酶体程序升高相互对应。单一SPP1升高可以出现在不同细胞背景，而这一组合将本项目的SPP1信号进一步定位于脂质处理和细胞内降解相关的巨噬状态。'),
p('胆固醇处理和溶酶体程序共同提供了一种代谢适应解释：巨噬细胞状态不仅涉及炎症信号强弱，也涉及对局部物质的处理和细胞内加工。在本项目中，这两类程序的患者一致性，为讨论SPP1富集状态的功能取向提供了比单个M1/M2标志更具体的依据。抗原呈递、炎症/干扰素程序及吞噬体的实际方向也应与这一组合同时解释，形成该状态的完整表达特征。'),
p('数量与信号贡献进一步区分了状态定位和总体关系。C4仅占成熟巨噬的一部分，却集中贡献捕获的SPP1 RNA；与此同时，SPP1也存在于其他成熟巨噬状态。因此，本研究以整体巨噬承载连续表达关系，以C4定位信号集中的状态。这种分工既保留细胞异质性，也使后续免疫和空间问题具有明确的观察对象。'),
h('3.4 工作模型连接代谢背景与巨噬状态'),
p('基于上述结果，本研究提出的工作模型是：TYMP相关恶性细胞状态伴随核苷及蛋氨酸/多胺代谢组织变化，MTA作为其中具有明确生成和清除途径、并可作用于髓系的候选代谢物，将肿瘤代谢环境与SPP1相关巨噬状态联系起来。公共MTA-SPP1的名义正向关联提供组织层线索，单细胞分群则把SPP1信号定位于以C4为代表的脂质处理和溶酶体相关状态。'),
p('这个模型的研究价值在于把代谢物候选与细胞状态落到可对应的层面：上游关注代谢环境，中间关注髓系响应，接收端关注SPP1相关状态的身份及其程序。它为后续分别研究免疫细胞功能和空间组织提供了统一起点。'),
h('3.5 证据适用范围与主线待验证环节'),
p('上述工作模型属于分段证据支撑的机制解释。TYMP经典反应不直接生成MTA；公共TYMP-MTA直接关联不稳定，同文库肿瘤TYMP与整体巨噬SPP1亦未形成稳定联系，已有串联模型的间接效应区间跨零。主线中尚待解决的是相同体系内TYMP相关代谢改变、胞外MTA暴露与髓系SPP1响应的连续对应。该缺口不改变已观察到的代谢程序、SMS-MTA分段联系及SPP1巨噬状态差异。'),
p('组织MTA与组织SPP1分别混合了多种来源，不直接等同肿瘤细胞分泌或巨噬专属读出。外部髓系扰动的ARG1/VEGFA等读出也不是本项目的SPP1诱导验证；受体或STAT3功能依赖不能直接写成普遍的磷酸化升高。不同癌种的MTA交换研究在此用于机制解释，不作为ccRCC局部暴露的测量。'),
p('C4发现与marker比较使用同一数据，属于分群后的条件性发现；其核心稳定性得到支持，完整17簇体系的全局稳定性未通过。患者和研究构成不均衡、轨迹有效覆盖随模型改变，也限制了普遍亚型或终末分化解释。横断面排序不证明细胞转化；C4与旧SPP1规则状态不同，不能直接继承旧标签的下游结论。'),
p('通路、调控网络、蛋白活性和代谢代理分别反映其模型测量层次，不能替代酶促速率、代谢物浓度、蛋白分泌或实测免疫功能。CNV树不等同经突变验证的祖先关系，F1的选择使用过TYMP相关信息；恶性细胞研究构成、主克隆敏感性及各端点有效患者数均按原结果保留。统计推断以独立患者为单位，名义显著与原家族FDR分别解释，主要反向结果不因本次文字重组而删除。'),
p('现有客户图仅对应TYMP敲低后的代谢变化，没有展示胞外MTA、MTA补回或巨噬SPP1干预终点。实验对应作为后续验证层面讨论，不参与公共数据候选筛选或效应估计。'),
]
wet=copy.deepcopy(b[187:193])
wet[-1]['text']='客户敲低结果中Thymidine升高与Thymine下降，对应TYMP相关核苷反应；Methionine、SAM和Spermidine下降则与本项目蛋氨酸/多胺背景相对应。MTA下降为候选代谢环节提供实验方向上的联系。由生信发现、具体代谢关系到细胞状态定位，再与已有干预方向对应，构成当前报告的论证顺序。'
conclusion=[h('3.7 阶段结论'),p('本研究显示TYMP相关恶性上皮状态具有核苷、蛋氨酸及部分多胺代谢特征，并据生化网络和实测分段证据引出MTA候选。公共组织联系与直接髓系作用依据支撑进一步研究SPP1相关状态；单细胞结果将该信号定位于具有脂质处理和溶酶体特征的C4巨噬状态。该框架构成后续免疫功能与空间组织研究的共同基础。')]
out=b[:167]+discussion+wet+conclusion+b[196:]
out[1]['text']='项目 26210　机制讨论整合修订版　2026-09-17'
for x in out:
 if x.get('text','').startswith('本报告保留第一部分2026'):
  x['text']='本报告以第一部分既有结果为基础展开至巨噬状态定位。2026-09-17修订版重写机制讨论，将主线边界集中呈现，并将中间机制探索移至独立补充材料。原图、所有表格数值及P/FDR保持不变。'
# New verified mechanism references.
out += [p('[6] Barekatain Y et al. Homozygous MTAP deletion in primary human glioblastoma is not associated with elevation of methylthioadenosine. Nature Communications. 2021;12:4228. https://doi.org/10.1038/s41467-021-24240-3'),p('[7] Reactome. TYMP catalysed thymidine phosphorolysis. https://reactome.org/content/detail/R-HSA-112265')]
assert [x for x in a if x['type']=='table']==[x for x in out if x['type']=='table']
assert sum(x['type']=='figure' for x in out)==7
main='\n'.join(x.get('text','') for x in out[:next(i for i,x in enumerate(out) if x.get('text','').startswith('4 统计'))])
assert all(w not in main for w in ['IL10','LGALS3','IL6扩展','PRMT5'])
(d/'report/report_materialized.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
supp="""# 中间机制探索补充材料

本材料保存主线之外的机制调查，不构成正文的新增论证步骤。所有结果来自已完成分析，本次未重算。

## S1 STAT3相关状态
整体巨噬STAT3代理-SPP1主31患者29正，中位rho=0.08087、P=4.32e-6；公共MTA-STAT3主rho=-0.12332、P=0.10520、FDR=0.15780，短路径区间跨零。TYMP-STAT3肿瘤内及公共组织联系保留，同文库跨区室主结果未获支持。以上为辅助状态证据，不登记为已确认中介。完整来源Phase179、180。

## S2 IL10与LGALS3
IL10扰动方向随体系变化，成熟人巨噬验证未支持通用IL10-STAT3-SPP1正向链。LGALS3在C4富集，但其与SPP1关系受组成影响，STAT3-LGALS3未支持正向必经联系。原阳性、反向和未支持结果均保留。完整来源Phase181-184。

## S3 成纤维中转与IL6
GSE249082的MTA处理后成纤维IL6 RNA小幅上升，两种归一化同向；独立供体关系未确认，仅描述。本项目肿瘤Fibroblast标签39细胞、1患者，来源关联不可估计。整体巨噬IL6-STAT3在19/19患者正向，但IL6-SPP1主rho=0.000944、P=0.7323；公共255例MTA-IL6主rho=-0.09149、P=0.1596、FDR=0.1596。不能拼成已确认MTA中转。完整来源Phase185-187。

## S4 上游短路径及其他调查
SMS-MTA分段结果与TYMP-SMS/多胺短路径未获稳定支持的结果同时保留。PRMT5及其他候选调查不进入正文机制模型。不得因本次正文取舍删除服务器原始统计、失败记录或候选清单。

## 完整材料入口
"""
for pattern in ['phase174*','phase175*','phase179*','phase180*','phase181*','phase182*','phase183*','phase184*','phase185*','phase186*','phase187*']:
 dd=next(r.glob(pattern));supp+='\n- '+str(dd)
(d/'supplement/MECHANISM_EXPLORATION_SUPPLEMENT_CN.md').write_text(supp)
(d/'qc/CONTENT_VALIDATION.json').write_text(json.dumps({'tables_unchanged':31,'figures_unchanged':7,'new_statistics':False,'exploration_removed_from_main':True,'main_boundary_section':'3.5','parent_source_sha256':hashlib.sha256((old/'report/report_materialized.json').read_bytes()).hexdigest()},indent=2))
renderer=(old/'scripts/render_report_v3.py').read_text().replace("phase189_integrated_TYMP_MTA_macrophage_report_20260917_v1","phase190_mechanism_focused_report_20260917_v1").replace("macrophage_integrated_report_v1","macrophage_integrated_report_v2")
(d/'scripts/render_report.py').write_text(renderer)
print('CONTENT_READY',len(out))
