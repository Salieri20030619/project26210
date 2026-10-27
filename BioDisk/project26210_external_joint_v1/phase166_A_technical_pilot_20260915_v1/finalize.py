
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil,datetime
from scipy.sparse import load_npz
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');O=P/'phase166_A_technical_pilot_20260915_v1';V=O/'strict_region_v2'
g=pd.read_csv(V/'GRAPH_QC.tsv',sep='\t')
g=g.rename(columns={'malignant_threshold':'malignant_fraction_median','macrophage_threshold':'macrophage_fraction_median'})
g['source_min_fraction']=.2;g['source_dominant_required']=True;g['receiver_min_fraction']=np.maximum(.1,g.macrophage_fraction_median);g.to_csv(V/'GRAPH_QC.tsv',sep='\t',index=False)
checks=[]
for gsm in g.section:
 A=load_npz(V/(gsm+'_adjacency.npz'));s=pd.read_csv(V/(gsm+'_SPOT_TECHNICAL_QC.tsv.gz'),sep='\t')
 checks.extend([dict(section=gsm,check='symmetric_no_self_nonnegative',passed=bool((A!=A.T).nnz==0 and A.diagonal().sum()==0 and (A.data>0).all())),dict(section=gsm,check='spot_degree_alignment',passed=bool(A.shape[0]==len(s) and np.array_equal(np.asarray(A.sum(axis=1)).ravel(),s.degree))),dict(section=gsm,check='minimum_three_receivers',passed=bool((s.loc[s.eligible_receiver_link,'eligible_macrophage_neighbors']>=3).all()))])
pd.DataFrame(checks).to_csv(O/'VALIDATION.tsv',sep='\t',index=False);assert all(x['passed'] for x in checks)
report="""# A方案技术试点：冷冻切片通过，FFPE邻域终点不适用
本轮完成两张预选切片的空间邻接、区域覆盖和设计矩阵核验。未计算生物学关联、P值、中介模型或绘图。

| 项目 | FFPE GSM5924034 | 冷冻 GSM5924045 |
|---|---:|---:|
| 总spot | 4734 | 1692 |
| 原六邻接中移除的过长连接 | 3.27% | 10.62% |
| 连通区数 | 1 | 2 |
| 严格源区域spot | 104 | 1076 |
| 有至少3个合格巨噬邻居的源spot | 2 | 282 |
| 上游设计可识别 | 通过 | 通过 |
| 接收端设计可识别 | 未通过 | 通过 |

上述spot数只用于技术QC，不是独立患者n。冷冻片最大连通区占99.82%，独立3-spot小块未被强行连接。距离规则与Visium阵列第一圈拓扑完全一致，坐标仿射最大残差均低于0.71像素。6项图结构核验通过。未做病理图像分区，不能称组织学边界已验收。

## 区域规则在查看关联前修正并冻结
初版按恶性丰度中位数分区，FFPE门槛仅4.53%，不能据此称肿瘤富集。
严格版源区：Malignant归一化q05丰度至少20%，并高于每一其他类别与合并成熟巨噬。
接收区：C1QC-state、Other、SPP1-state三类成熟巨噬合并归一化q05丰度至少10%，同时不低于本片中位数；单核单列。
邻域最多六个非自身最近邻，距离至多本片最近邻距离中位数的1.3倍，且符合阵列第一圈拓扑；至少3个合格接收邻居才纳入。
源和接收区域允许重叠，邻接不含自身。归一化q05丰度不是测定的细胞纯度，以上为操作阈值。
首版作为诊断保留；严格版同时修正上游设计预测变量为TYMP。没有根据关联方向调参。严格版FFPE接收端失败未触发阈值放宽。

## 模型口径
上游拟以源区TYMP预测SRM、SMS、SLC16A3，分别检验。
接收端拟以源区候选表达对应邻居SPP1平均表达，调整源区TYMP。邻居平均是明确空间终点，不对全片表达先平滑。
基础调整含源区深度、恶性和合并巨噬归一化丰度；接收端增加邻域对应变量。扩展调整增加基质COL1A1/COL1A2/DCN/LUM和缺氧CA9/EGLN3/VEGFA/BNIP3/PDK1评分，不含目标基因。
严格冷冻矩阵满秩，条件数最高13.27；FFPE上游满秩，接收端仅2源spot，不能识别9–13列设计。
本轮仅核验未秩变换线性设计；正式模型前仍要复核秩变换设计、终点方差和邻域重叠。条件数不代表空间有效样本量。表达检出详见strict_region_v2/REGION_GENE_COVERAGE.tsv。

## 参考及身份
TYMP、SPP1、SRM、SMS、SLC16A3均在旧13状态参考表中。需后续检查目标基因剔除敏感性；本轮没有重拟合解卷积。旧参考没有新C4、也未单列成纤维细胞。基质评分不写成成纤维细胞比例；spot表达不直接赋予细胞来源。
GEO未提供患者ID；原文页面访问失败，患者映射仍未解决。本轮未联系作者。
来源：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175540
原文：https://www.cell.com/immunity/fulltext/S1074-7613(22)00081-4

## 后续决策
下一步对24切片统一应用固定严格规则，先形成端点可用性清单，不以关联正负或P选择样本。
源区或有效接收端源spot少于50、设计不满秩或条件数至少30，则不拟合对应端点。所有技术合格切片均纳入，不仅保留本次冷冻片。FFPE/frozen分层。
患者身份未解决时只可做逐片效应描述，不能把24切片当24患者汇总P/FDR。
MTA多胺轴保留，MTAP单列清除背景，SLC16A3并行，C4延续单细胞状态定位。本轮没有新增支持或反对TYMP—X—SPP1关系的生物学证据。
"""
(O/'TECHNICAL_PILOT_CN.md').write_text(report)
s=json.loads((O/'STATUS.json').read_text());s.update(status='PILOT_COMPLETE_CONDITIONAL_EXPANSION',all_designs_pass=False,strict_source_pass_sections=2,strict_receiver_pass_sections=1,validation_checks=6,validation_passed=6,primary_output='strict_region_v2',next_step='24 section endpoint eligibility QC')
(O/'STATUS.json').write_text(json.dumps(s,ensure_ascii=False,indent=2))
pd.DataFrame([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in list(O.glob('*.py'))+list(V.glob('*.tsv'))+[O/'TECHNICAL_PILOT_CN.md']]).to_csv(O/'ARTIFACT_HASHES.tsv',sep='\t',index=False)
back=O/'authority_backups';back.mkdir(exist_ok=True)
note='\n## Phase166 A技术试点完成 2026-09-15\n两片图结构通过；严格冷冻接收端282源spot可识别，FFPE仅2个不可拟合。无关联/P值/图。下一步24片统一端点可用性；患者映射未解决。报告：'+str(O/'TECHNICAL_PILOT_CN.md')+'\n'
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=P/n;shutil.copy2(f,back/n);f.write_text(note+'\n'+f.read_text())
f=P/'RESULTS_MANIFEST.tsv';shutil.copy2(f,back/f.name)
with f.open('a') as h:
 for n in ['TECHNICAL_PILOT_CN.md','STATUS.json','VALIDATION.tsv','strict_region_v2/GRAPH_QC.tsv','strict_region_v2/DESIGN_IDENTIFIABILITY.tsv','strict_region_v2/REGION_GENE_COVERAGE.tsv']:
  h.write('\t'.join(['phase166_'+n.replace('/','_'),'Phase166','technical_pilot',str(O/n),datetime.datetime.now().astimezone().isoformat(),'CONDITIONAL','6 structure checks; 1 of 2 receiver designs eligible','No association tests'])+'\n')
print(hashlib.sha256((O/'TECHNICAL_PILOT_CN.md').read_bytes()).hexdigest())

