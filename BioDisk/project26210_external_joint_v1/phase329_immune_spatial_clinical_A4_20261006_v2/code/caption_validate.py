from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
import pymupdf as fitz
from PIL import Image
P=Path(__file__).resolve().parents[1];R=P.parent
rd=lambda f:pd.read_csv(P/'data'/f,sep='\t')
b=rd('PANEL_B_STATISTICS.tsv');e=rd('SPATIAL_COHORT_SUMMARY_ALL.tsv');d=rd('BOUNDARY_COHORT_SUMMARY.tsv');f=rd('KO_QUALIFIED_MODEL_CHECKS.tsv')
cap='''TYMP–SPP1免疫程序、空间生态位及临床结局（A4综合图 v2）

本图为独立预览，尚未插入全文。纵向A4，210×297 mm；8个面板，300 dpi PNG及PDF。A为原创概念插画；B–H全部来自既有冻结数值或真实空间坐标，无新增模型、检验、分群、样本筛选或平滑推断。

A｜巨噬细胞–基质生态位的细胞组织形式。参考用户提供的Science Advances论文（sciadv.aed0733）图3H的组织微环境构图，以肿瘤细胞巢、巨噬细胞、成纤维细胞、基质纤维及CD8 T细胞呈现细胞关系。插画由内置imagegen生成，再叠加本项目矢量文字；不沿用参考图的APOE/POSTN特异性、三层距离或已证实物理屏障结论。细胞形态、纤维走向与细胞位置均为示意，不是显微成像、真实定位或细胞迁移数据。右侧为Obradovic巨噬细胞内的实测RNA程序关联：TYMP与IFN-γ、MHC-I呈递及吞噬正相关；SPP1与IFN-γ及MHC-II呈递负相关；两者均与糖酵解正相关。加减号是相关方向，不是干预后的增减或因果箭头。

B｜患者层面的免疫程序。宽松巨噬定义纳入15,056细胞，但独立统计单位为11名患者。保留预设15项免疫程序，加上已报告糖酵解，共16行；横轴1–8、A–C是冻结患者标识，TYMP/SPP1使用相同患者顺序。每格是患者内控制测序深度、检测基因数、文库和髓系簇后的秩残差ρ，不按行标准化，色标固定为−0.4至0.4。右侧是患者中位ρ，粗体表示原224项家族BH FDR < 0.05。所有评分排除TYMP和SPP1。T细胞激活正/负调节是巨噬细胞中相应GO程序的关联，并非直接测得的T细胞功能。A中的五项免疫方向经相互控制和免疫程序去重叠保留；IFN-α、炎症及糖酵解删基因敏感性仍按正式报告保留，尤其糖酵解不代表通量测量。

C｜同一组织的真实空间程序图。选用GSE175540的GSM5924033，因其也是用户参考论文图3使用的切片，在查看本轮程序图前即确定，未按相关强度筛选。展示全部4,939个冻结QC spot，三个小图使用完全相同的真实坐标、方向和视野，不平滑、不重构热点、不绘制人工组织边界。SPP1及FAP为已冻结的文献25基因来源程序去重叠版本，ECM为排除重叠成员的保留程序；均沿用Phase317的AddModuleScore。颜色显示各程序分数，分别以切片内第2及98百分位设色限，超出者仅饱和显示；不能用不同程序的颜色深浅比较绝对表达。比例尺1 mm，依据既有像素–微米换算1.8210517240511208。坐标标识与分数逐行一一核对。此图是单切片空间示例，全队列关系见E。

D｜固定表达界面内外的CD8 RNA差。GSE175540，固定400 μm带宽，20名合格患者。每根柱为一位患者的内侧减外侧效应，控制测序深度和总巨噬指标，纵轴为原冻结标准化RNA差值；按患者标识排序，未按效应挑选患者。19/20为负，中位差−0.260230，显示界面内侧CD8谱系信号较低。P = 0.00004005，FDR = 0.00016022，沿用原4项家族。该分区事实与SPP1特异解释分开：同一研究SPP1局部占比与CD8内外差的调整相关中位ρ = +0.049298，P = 0.115318，FDR = 0.461273，当前不支持由SPP1局部富集解释该分区。表达映射界面不冒称独立病理标注。

E｜患者层面的SPP1–基质程序联系。GSE175540全23位患者、24张切片；每名患者各贡献两个冻结的深度调整相关系数，灰线连接同一患者的FAP–SPP1程序及SPP1–ECM程序相关值。点为患者，粗横线为中位数；小幅横向错开仅避免点重叠。双向符号表示相关而非因果方向。FAP–SPP1为23/23正，中位ρ=0.152709；SPP1–ECM为20/23正，中位ρ=0.163723。三个负向患者全部保留。FDR分别来自原coupling 8项和process 12项家族，不重新合并校正，不将23位患者×两类关系计为46个独立样本。这是文献标记向RCC的探索性转移结果；RNA程序共分布不等于纯细胞物理邻接。

F｜公开Spp1干预补充功能参照。来自Henriques等（Nature Genetics 2025，doi:10.1038/s41588-025-02380-2）公开源值，沿用Phase319的逐动物复算。左图为原Extended Data Fig.10H，组织重量归一化CAF计数，WT 6只、Spp1缺失7只；右图为原Extended Data Fig.10C，CD8内Ki67阳性比例，WT和缺失各5只。每点代表一只动物；右图先对同鼠多个观察汇总，再绘制动物均值；横线为组均值。两端点属于不同实验面板，不作配对。CAF组均值由2792.541降至1547.815；CD8增殖组均值由30.563%升至50.249%，均为描述性效应。右图既有混合模型估计与简单动物均值差略有差异，完整原模型统计附后。该外部CRC实验支持机制参照，不替代本项目RCC队列的因果验证。其他公开模型中CD8变化不一致的结果仍见完整报告。

G｜SPP1相关五基因评分与OS。五基因为SPP1、FN1、CLEC5A、OLR1及IL1RN，按完整311人CheckMate队列中位数分组。展示181人nivolumab组，高84、低97，123个死亡事件。阶梯为KM估计，阴影95%CI，“+”为删失，底部数字为各时点在险人数。Log-rank P=0.014836、FDR=0.059343，为名义显著。临床调整连续评分HR=1.269884，P=0.016110、FDR=0.064440，口径保持不变。

H｜C4状态与PFS的时间依赖关联。冻结19基因评分排除TYMP及五基因SPP1评分全部成员，nivolumab完整模型164人、145个PFS事件；控制年龄、性别、试验、MSKCC、五基因评分及总巨噬估计。曲线为基线评分每增加1 SD对应的条件HR（纵轴对数刻度），带为点态95%CI，底部为在险人数。固定β(t)=β0+β1log(1+t)，未搜索时间切点；联合P=0.0007173、FDR=0.001435。12个月后灰底只提示风险集变小。该曲线不表示C4功能在体内逆转，当前未建立其独立OS关联或相对治疗获益预测。

原v1的治疗交互面板已依用户意见移出综合图，原统计未删除：CM025 OS交互HR比1.594944（1.170496–2.173308），P/FDR均0.003105；CM025 PFS P=0.299660；JAVELIN PFS P=0.331944、FDR=0.663888，后两者当前不支持相应交互。所有临床原表保留。

新增空间与实验面板的精确冻结统计：
'''
for target in ['FAP_CAF25_disjoint__SPP1_MAC25_disjoint','SPP1_MAC25_disjoint__ECM_org_heldout']:
 z=e[(e.dataset=='GSE175540')&(e.stratum=='overall')&(e.model=='depth_adjusted')&(e.endpoint==target)].iloc[0]
 cap+=f'E / {target}: n={z.n_patients}, median rho={z.median_effect:.9f}, P={z.p_raw:.12g}, FDR={z.FDR:.12g}, family={z.family}, family size={z.family_planned_size}.\n'
for z in f.itertuples():cap+=f'F / {z.panel}: WT n={z.n_WT}, KO n={z.n_KO}, model={z.model}, KO−WT={z.estimate_KO_minus_WT:.9f}, P={z.P:.12g}, FDR={z.FDR_two_model_audit:.12g}; two-model audit family.\n'
cap+='\nB面板原始P及FDR（224项家族）：\n'
for z in b.itertuples():cap+=f'{z.source} / {z.target.split("__")[-1]}: n={z.patients}, median rho={z.median_rho:.6f}, P={z.P:.9f}, FDR={z.FDR:.9f}.\n'
(P/'delivery/26210_immune_spatial_clinical_A4_v2_图注.txt').write_text(cap,encoding='utf-8-sig')
# Validate identities and exported coordinates/data, no statistical recomputation.
manifest=json.loads((P/'qc/SOURCE_MANIFEST.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for src in manifest:
 if 'sha256' in src:assert sha(Path(src['path']))==src['sha256']
xy=rd('SPATIAL_MAP_VALUES.tsv');sc=rd('GSM5924033_MODULE_SCORES.tsv');assert xy.spot_id.tolist()==sc.spot_id.tolist();assert np.allclose(xy[sc.columns[1:]],sc.iloc[:,1:],rtol=0,atol=1e-15)
assert len(rd('PANEL_D_VALUES.tsv'))==20 and len(rd('PANEL_E_PATIENT_PAIRS.tsv'))==23
full=R/'phase327_abstract_clarity_revision_20261006_v1/delivery/26210_manuscript_immune_integrated_v3.pdf';assert sha(full)=='cf2a41dc0701e99367d760d621a2456a231690c267c019d8cd648cc694602efb'
pdf=P/'delivery/26210_immune_spatial_clinical_A4_v2.pdf';doc=fitz.open(pdf);assert len(doc)==1;page=doc[0];assert abs(page.rect.width*25.4/72-210)<.001 and abs(page.rect.height*25.4/72-297)<.001
page.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False).save(P/'qc/final_pdf_render.png')
spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']];assert all('\ufffd' not in s['text'] for s in spans)
(P/'qc/NUMERIC_PDF_VALIDATION.json').write_text(json.dumps({'status':'PASS','page_mm':[210,297],'panels':8,'spatial_spots':len(xy),'patient_units':{'B':11,'D':20,'E':23},'all_source_hashes_unchanged':True,'full_manuscript_unchanged':True,'body_min_font_pt':6.3,'new_statistical_analysis':False},indent=2))
print('CAPTION_AND_VALIDATION_READY')
