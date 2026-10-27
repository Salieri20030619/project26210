from pathlib import Path
import json,hashlib,datetime
import pandas as pd
P=Path(__file__).resolve().parents[1];R=P.parent
s=P/'code/draw_figure.py';t=s.read_text().replace('[113,124,131]','[113,124,126]').replace("'SPP1 score and survival'","'SPP1 TAM score and survival'").replace("'SPP1 score: treatment interaction across clinical endpoints'","'SPP1 TAM score: treatment interaction across clinical endpoints'");s.write_text(t)
f=pd.read_csv(P/'data/FORMAL_SUMMARIES.tsv',sep='\t');b=pd.read_csv(P/'data/PANEL_B_STATISTICS.tsv',sep='\t')
cap='''26210｜TYMP–SPP1巨噬细胞免疫程序与临床关联综合图（独立预览 v1）

图件规格：纵向A4，210 × 297 mm；PDF为矢量图，PNG为300 dpi。A–E为独立综合图面板编号，尚未插入正式文章。

A｜共同代谢背景下的免疫程序差异。以巨噬细胞示意图概括Obradovic未经治疗RCC队列内的表达关联：TYMP与IFN-γ响应、MHC-I抗原呈递及吞噬程序正相关；SPP1与IFN-γ响应及MHC-II抗原呈递程序负相关；两者均与糖酵解程序正相关。实线及“+ / −”表示关联方向，不表示基因干预或因果传递。细胞、囊泡及膜表面结构是功能解释用的示意形态，并非显微成像或实测亚细胞定位。所有程序评分均排除TYMP和SPP1。

B｜患者层面的完整免疫程序分布。宽松巨噬定义共15,056个细胞，独立统计单位为11名患者；横轴1–8及A–C对应冻结患者标识，两基因区块使用完全相同的患者顺序。纵轴保留预设15项免疫展示程序，并补充既有糖酵解结果，未按本次显著性重新挑选通路。色块为每名患者的秩残差相关系数ρ：在患者内控制测序深度、检测基因数、文库及髓系簇，红色为正、蓝色为负，统一色标−0.4至0.4，不按行标准化。右侧为11名患者ρ的中位数，粗体表示原224项免疫检验家族的BH FDR < 0.05；常规字体包括名义信号及非显著结果，精确原始P/FDR均保留在服务器配套表。T细胞激活正/负调节两行是巨噬细胞中相应GO基因程序的表达关联，不是直接测得的T细胞功能。

A所示五项免疫方向在相互控制另一基因及去除免疫程序重叠基因后仍保留。IFN-α及炎症程序的信号对基因集重叠处理更敏感，因此B保留原始结果，但A不将其列为稳健的功能主线。糖酵解的相关性不等同于代谢通量；删去贡献最大的五个基因后TYMP关联明显减弱，不能概括为整个代谢网络均匀改变。

C｜SPP1相关巨噬五基因评分与总生存。评分为SPP1、FN1、CLEC5A、OLR1及IL1RN表达均值，以完整311例CheckMate队列中位数分组；显示181例nivolumab患者，高组84例、低组97例，共123个死亡事件。阶梯曲线及色带为Kaplan–Meier估计和95%置信区间，“+”为删失，下方为在险人数。高评分组OS较短，log-rank P = 0.014836，FDR = 0.059343，为名义显著。本面板评分不同于单基因SPP1表达，也不同于C4状态评分。

D｜基线C4状态评分与PFS的时间依赖关联。使用冻结19基因评分，排除TYMP及上述全部五基因；CheckMate nivolumab完整模型纳入164例、145个PFS事件，控制年龄、统一性别、试验、MSKCC、五基因评分及总巨噬估计。纵轴为基线C4评分每增加1 SD的条件HR（对数刻度），实线和色带为HR及点态95%置信区间；固定时间函数β(t)=β0+β1log(1+t)，未搜索时间截点。两自由度增量检验P = 0.0007173，FDR = 0.001435。12个月后灰底标示在险人数减少，属于显示提示；此曲线描述风险关联随随访时间变化，不是C4细胞功能发生逆转。当前未建立C4独立OS关联或相对治疗获益预测。

E｜五基因评分与治疗的交互效应。点和横线为交互HR比值及95%置信区间，虚线为1，横轴对数刻度。CM025比较nivolumab与everolimus并调整年龄、性别、MSKCC；JAVELIN比较avelumab联合axitinib与sunitinib并调整年龄、性别。HR比值大于1表示评分升高时免疫治疗臂相对对照臂的风险比增加。保留CM025 OS阳性交互、CM025 PFS及JAVELIN PFS非显著结果；后两项当前不支持相应PFS治疗交互。FDR沿用各队列/终点内单基因和五基因指标的原校正家族，不将不同队列终点混池。

A–B来自单细胞发现队列，C–E来自独立临床转录组队列，并非相同患者的一条连续因果路径。临床评分结果用于联系组织转录状态与结局，不能直接换算为某一巨噬细胞亚群丰度。

主要免疫结果的精确统计（宽松定义；患者内相关，跨患者符号检验；未重算）：
'''
keys=[('TYMP','Interferon Gamma Response'),('TYMP','GO:0002474__antigen processing and presentation of peptide antigen via MHC class I'),('TYMP','GO:0006909__phagocytosis'),('SPP1','Interferon Gamma Response'),('SPP1','GO:0002495__antigen processing and presentation of peptide antigen via MHC class II'),('TYMP','Glycolysis'),('SPP1','Glycolysis')]
for gene,target in keys:
 z=b[(b.source==gene)&(b.target==target)].iloc[0];c=f[(f.condition=='L2_broad')&(f.source==gene)&(f.target==target)&(f.variant=='conditional')].iloc[0]
 cap+=f'{gene} / {target.split("__")[-1]}：原始中位ρ={z.median_rho:.6f}，正/负患者={z.positive}/{z.negative}，P={z.P:.9f}，FDR={z.FDR:.9f}（224项家族）；相互控制后中位ρ={c.median_rho:.6f}，P={c.P:.9f}，FDR={c.FDR:.9f}（20项家族）。\n'
cap+='\n来源：Phase314冻结患者效应及224项汇总；Phase321条件/去重稳健性汇总；Phase273临床综合图原始绘图表。此次仅重绘与综合排版，无新增检验、模型、分群或样本筛选。完整输入哈希、作图脚本、统计表与校验在服务器Phase328目录保留。\n'
(P/'delivery/26210_immune_clinical_A4_v1_图注.txt').write_text(cap,encoding='utf-8-sig')
print('CAPTION_READY')
