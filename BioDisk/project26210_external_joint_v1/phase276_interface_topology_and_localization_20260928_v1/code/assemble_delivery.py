from pathlib import Path
import pandas as pd,numpy as np,json,datetime,shutil,hashlib,zipfile
P=Path(__file__).resolve().parents[1]
s=pd.read_csv(P/'tables/COHORT_SUMMARY.tsv',sep='\t');primary=s[s['quantile']==.75].copy()
def fmt(x):
 if pd.isna(x):return '—'
 return x if isinstance(x,str) else f'{x:.6g}'
def table(d,cols,labels):
 out=['| '+' | '.join(labels)+' |','| '+' | '.join(['---']*len(cols))+' |']
 for _,r in d.iterrows():out.append('| '+' | '.join(fmt(r[k]) for k in cols)+' |')
 return '\n'.join(out)
cols=['dataset','endpoint','model','n','positive','negative','median','P','FDR']
labels=['队列','端点','模型','患者n','正向','负向','中位效应','主P','FDR']
report='''# 26210：界面结构、缺口与独立定位尝试（Phase276）

完成日期：2026-09-28。依据用户“1、2先进行尝试”开展空间结构与特异关联计算，并核查独立病理/多重定位资产。本轮不出图、不修改正式稿；既有Phase275结果保留。

## 结论与实际完成范围

**结构层面的补充检验仍未建立稳定的SPP1特异屏障证据。** 已完成38张切片、37位患者的边界重建、连续段/覆盖/低覆盖缺口描述，以及患者层面的相关和同轮廓匹配检验。SPP1覆盖率与CD8内外差的预期负关联未在主要队列获得支持。

主阈值下，GSE175540的严格缺口匹配保留一条描述性线索：9位可评估患者中7位缺口相对覆盖段的CD8内外差更高，中位差+0.05497评分单位。但n<10不作推断检验，而且预设50%/80%阈值下不稳定。不能把这一筛选后可评估小子集写成总体屏障成立。

**独立病理边界/多重定位未能实际运行。** 本地已查资产有38张配套H&E、GLMF两份旧AI边界草稿（另有相同副本），但没有审核病理边界或可用的SPP1/巨噬/CD3/CD8原始多通道图像/细胞坐标表。已经完成核查、定位测量方案与数据交接模板；没有把这些准备工作写成染色实验或蛋白定位验证完成。

## 新计算具体做了什么

沿用Phase234冻结的恶性映射比例≥0.20表达分区。只在三个原始100 µm邻接边均存在的三角形内连接肿瘤/外侧跨区边中点，形成表达交界网络；不跨缺失spot、组织空洞或连通分量，不人为闭合轮廓。网络最大度数为2。由几何坐标确定起点和方向，按400 µm弧长分段，末段不足200 µm保留在QC但不进入相关检验。

38片合计39611个跨区边中点、38424条边界连接、5985个边界分段；3930段同时满足几何与内外各至少10个spot的覆盖要求。所有段都不是独立患者，未用段数或spot数放大统计样本量。

SPP1高状态以每片巨噬比例≥0.05的spot内SPP1/总巨噬份额第75百分位为主阈值，50/80百分位事先固定为敏感性；边界中点用两端spot的SPP1 q05之和/巨噬q05之和，要求合并巨噬比例≥0.05。C1QC按同一规则作状态对照。两者是相对映射富集，不是蛋白阳性细胞分类。

主要关联为：边界段高状态覆盖率与该段400 µm虚拟图邻域“内侧−外侧”的总T log2映射丰度/CD8谱系RNA差异。负相关才符合覆盖增强、内部免疫相对降低的预期。
- 未调整/核心调整各需至少15个合格段；
- 核心调整包括平均深度、深度内外差、总巨噬丰度、平均恶性比例、恶性比例内外差及内外spot数量比；
- 扩展调整再加入内皮比例和成纤维细胞RNA，至少30段；
- 秩变量完整性、变化范围、设计矩阵与条件数标准冻结，不因结果放宽。

同一患者多切片先合并，相关用Fisher-z平均、差值用算术平均；队列分开。n≥10采用双侧精确符号检验，Wilcoxon并列保留；n<10仅描述。主覆盖关联家族为12项（两状态×两读出×三模型），敏感性24项；缺口家族主2项、敏感性4项。全126队列端点、12个FDR家族完整保留。

## SPP1覆盖与免疫内外差：主要结果

'''+table(primary[primary.endpoint.str.startswith('SPP1__')],cols,labels)+'''

GSE175540扩展调整后CD8 n19，中位rho=+0.04939，10正9负，主P=1、FDR=1；总T n22，中位rho=+0.08907，15正7负，P=0.1338、FDR=0.7975。预期的稳定负向并未出现。TJ-RCC CD8有效n8，仅描述；扩展调整6正2负，中位+0.06326。

GLMF两例CD8关联均为负向：GLMF1未调整−0.73099、核心调整−0.32312、扩展调整−0.15686；GLMF2分别−0.09332、−0.04035、−0.13399。这是可保留的病例方向，但只有两位患者，并且C1QC扩展调整也均负向，不能据此宣称SPP1特异性或替代外部队列结果。

## C1QC对照与主要反向结果

'''+table(primary[primary.endpoint.str.startswith('C1QC__')],cols,labels)+'''

C1QC也可出现负向空间关联。第80百分位敏感性中，GSE175540核心调整C1QC-CD8 n17，13负4正，rho=−0.06535，P=0.04904、FDR=0.5885；C1QC-总T n20，15负5正，rho=−0.05695，P=0.04139、FDR=0.5885。两项属于名义显著探索信号，不因FDR未过而删除，也不能移用为SPP1证据。

SPP1第80百分位扩展调整CD8为正向（rho=+0.08389，12正7负）；主P=0.3593，伴随Wilcoxon P=0.04937/FDR=0.3373。这是反预期、检验方法依赖的名义信号，完整保留。

C1QC/SPP1同属总巨噬内组成份额，具有组成依赖；此对照用于检查状态特异性，不能称完全独立阴性对照。成纤维RNA也不等于独立CAF细胞计数。

## 同一轮廓的覆盖段—缺口匹配

以SPP1覆盖率≥75%的段为覆盖段，≤25%的段为低覆盖“缺口”。在同一条表达轮廓内，用核心6项协变量标准化后的距离进行无放回最小成本匹配，固定卡钳2；每切片至少3对才合格。匹配不使用免疫结局数值来选择方向。

结局为“缺口处内外差−覆盖段内外差”，正值才符合缺口对应更高内部免疫信号。这里的“缺口”是SPP1低覆盖区，尚不是显微镜确认的物理通道。

'''+table(s[s.family.str.startswith('matched')],['dataset','quantile','endpoint','n','positive','negative','median','P','FDR'],['队列','富集分位阈值','读出','患者n','正向','负向','中位差','P','FDR'])+'''

主阈值GSE175540 CD8 7/9位方向符合预期，但50百分位12位中6正6负、80百分位6位中3正3负。总T主阈值n10，6正4负，中位差+0.05676，P=0.7539/FDR=1。TJ-RCC主阈值只有3位；GLMF两例均不满足至少3对的资格，未为了主案例展示而放宽匹配。

因此，“缺口允许更多CD8进入”目前只有小样本、阈值不稳定的描述性线索，尚未补足SPP1屏障证据。全部匹配对、未匹配原因、共同患者阈值对照和留一患者结果随表交付。

## 覆盖、连续性与横向范围

按主阈值、同患者切片等权平均后，SPP1覆盖率的队列描述性中位数约为GSE175540 18.67%、TJ-RCC 19.20%、GLMF 24.71%；最长连续高状态段的对应中位数为750、550、1825 µm。它们描述了SPP1相对高值带在表达界面的组织方式，不是检验“比随机更多”的P值，也不意味着完整围墙形成。

闭合轮廓数量是表达分区轮廓的几何属性，不是SPP1完整包围环数。高状态阈值依赖切片内分布，覆盖率不能直接当屏障面积比例。横向范围用分配给该段的高状态spot至界面的90百分位图距离描述，限定400 µm内，内外各至少3个高状态spot；不能称病理层厚度。原始逐段描述表完整交付。

## 第2步独立定位：查到了什么、缺什么

在四个明确项目/公共空间资料目录中扫描3368个文件，候选核对如下：
- GLMF1/2两个边界JSON及各自副本：标注not_pathologist_reviewed，来源为AI辅助且曾参考计算边界，不能作为独立审核轮廓；
- 38张H&E配套记录：图像路径/内嵌RGB维度存在，支持组织定位，不含所需独立标记通道；
- TLS资料KC1/KC2/KC3三个TIFF：单页RGB（YXS、3样本通道），没有可验证的SPP1/巨噬/CD3/CD8原始通道及配套细胞表，未纳入多重定位验证。

本结论仅限明确扫描范围，不等同于全互联网不存在可用数据。未重启以前暂停的大文件下载，未联系作者、未提交受控申请、未重新分析客户湿实验。实际独立定位仍待合格图像/坐标及病理审核记录。

已交付“独立定位测量与交接方案.md”、图像/细胞/ROI空模板及38切片病理审核任务表。它们用于把缺失输入具体化，不代表实验完成。审核过程应对SPP1/CD8结果盲法，保留不确定区域、组织坏死和切片边缘；不得用AI自动标注补成病理确认。

## 证据定位与验收

本轮增加的是“连续结构与缺口”的直接问题，未重复搜索既有局部相关阈值。主要结果和固定敏感性仍不支持稳定的SPP1特异屏障关联；此前Phase275的CD8内外差异成立，不能被本轮阴性反向抹去。可以继续保留免疫空间分区主线；SPP1屏障属于待独立定位与功能验证的假说。

1596条切片效应、1554条患者效应、126条队列端点均已计算。独立核对中点几何、长度守恒、无人工闭合、覆盖率、DFS连续段、患者合并、精确二项P和BH；376项独立偏秩相关最大误差8.89e−16，597条匹配记录的同轮廓、卡钳、无复用及差值均通过。原输入SHA未变，进程退出0。完整验证见FINAL_VALIDATION.json。

本轮图件0，正式稿修改0。计算结果已完成；独立病理/多重定位验证因输入缺失未完成。
'''
(P/'delivery/界面结构与定位尝试报告.md').write_text(report,encoding='utf-8')
protocol='''# 独立定位测量与交接方案（准备稿，尚无合格数据）

目标：验证SPP1来源及其与巨噬、CD8、病理界面的空间位置。区别组织内RNA区域与真实细胞定位；本文件是数据与测量方案，不是已完成实验。

## 优先资料

1. 同一切片的原始多通道图像及标记—通道表、像素/µm比例、扫描/染色批次、单染和背景控制。优先可回溯细胞分割与坐标，不以文章截图或合成RGB代替。
2. 核心身份：核标记、巨噬身份标记、SPP1、CD3/CD8，以及病理确认肿瘤区。具体抗体/平台和阳性阈值由实验方依据验证与控制确定，不从预期结论反推。
3. 如果归因“巨噬来源SPP1”，需确认细胞内SPP1 RNA/蛋白与巨噬身份对应；细胞外SPP1沉积与细胞内表达分开，不能把邻近荧光直接归到巨噬来源。
4. 基质分支可加FAP、αSMA、胶原读出；血管及周细胞结构单独标注，避免把血管壁平滑肌信号当CAF。
5. 连续切片可以做区域对应，但未完成可靠配准时不能声称同一个细胞共表达。供体/患者、组织块、切片、ROI的层级必须有明确ID。

## 病理审核

以独立病理图像确定肿瘤/间质/坏死/不确定区，审核者不查看SPP1/CD8统计方向。已有GLMF草稿仅供工具/坐标参考，不是已审核真值；最好重新盲法审核并记录更改。提供reviewer、review_date、审核状态、polygon坐标系及源图像SHA。未知区域不自动归为正常/间质。

病理审核任务表列出38片原始图像位置；未预填任何审核完成记录。不得把空模板作为已取得批准或已完成人工审核。

## 预定义测量问题

- 身份：SPP1信号在巨噬与其他细胞/细胞外区域的分布；
- 位置：SPP1阳性巨噬在病理边界内外的实际密度（细胞/mm²）、覆盖与连续性；
- 免疫：同一段边界两侧CD3+CD8+细胞密度、至肿瘤边界的有符号距离和内部浸润；
- 特异：同患者、类似总巨噬/基质/血管背景下SPP1阳性与阴性巨噬相关边界段的差异；
- 缺口：预设定义的低覆盖边界段与高覆盖段比较，未发现低覆盖段时记不可评估；
- 不把静态相邻或内外分布差异直接称为物理阻挡，不用Ki67/存活变化代替迁移测量。

正式图像分析前，依据实际分辨率和组织尺度另行冻结细胞阳性规则、边界带宽、覆盖连接距离与主要端点，不能机械沿用本轮Visium第75百分位作为蛋白阳性阈值。先做分割/染色可行性与测量重复性检查，检验和功效按独立患者而非细胞/ROI数量规划。推断继续遵守患者n、原P/FDR与小样本规则。

## 交接表

IMAGE_MANIFEST_TEMPLATE.tsv：每图一行，保留患者/切片/通道/标尺/批次/审核/许可来源。
CELL_TABLE_TEMPLATE.tsv：每细胞一行，保留坐标、分割ID、分类与不确定性；SPP1巨噬和CD3+CD8+身份必须可追溯。
ROI_TABLE_TEMPLATE.tsv：每ROI一行，关联组织面积、病理轮廓和审核记录。
这些模板只有表头、没有虚构患者或细胞数据。

## 方法背景

相关研究用于界定所需证据层次，不替代26210自身实验：
- 肝癌SPP1/CAF屏障研究：https://pubmed.ncbi.nlm.nih.gov/36708811/
- 肾透明细胞癌SPP1与CD8功能研究：https://pubmed.ncbi.nlm.nih.gov/40680256/
本轮未执行功能干预，不新增实验结论。
'''
(P/'delivery/独立定位测量与交接方案.md').write_text(protocol,encoding='utf-8')
schemas={
 'IMAGE_MANIFEST_TEMPLATE.tsv':['patient_id','sample_id','section_id','image_id','assay','source_image_path','sha256','channel_marker_map','microns_per_pixel_x','microns_per_pixel_y','batch_id','same_section_or_serial','registration_record','pathology_review_status','reviewer','review_date','source_permission'],
 'CELL_TABLE_TEMPLATE.tsv':['patient_id','section_id','image_id','roi_id','cell_id','x_um','y_um','segmentation_id','cell_class','classification_confidence','macrophage_positive','SPP1_intracellular_positive','SPP1_extracellular_flag','CD3_positive','CD8_positive','quality_flag'],
 'ROI_TABLE_TEMPLATE.tsv':['patient_id','section_id','image_id','roi_id','polygon_path','coordinate_space','tissue_area_mm2','tumor_area_mm2','stroma_area_mm2','necrosis_area_mm2','uncertain_area_mm2','reviewer','review_date','review_status']
}
for f,columns in schemas.items():pd.DataFrame(columns=columns).to_csv(P/'delivery'/f,sep='\t',index=False)
reg=pd.read_csv(P/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t');im=pd.read_csv(P/'tables/IMAGE_AND_LOCALIZATION_READINESS.tsv',sep='\t')
todo=reg[['sample','dataset','patient']].merge(im[['sample','image_path','width','height']],on='sample',validate='one_to_one')
todo['review_status']='NOT_REVIEWED';todo['reviewer']='';todo['review_date']='';todo['approved_polygon_path']=''
todo.to_csv(P/'delivery/PATHOLOGY_REVIEW_TASKS.tsv',sep='\t',index=False)
for f in (P/'tables').glob('*.tsv'):shutil.copy2(f,P/'delivery'/f.name)
for f in ['PROTOCOL_BEFORE_TESTING.md','FREEZE.json']:shutil.copy2(P/f,P/'delivery'/f)
for f in ['FINAL_VALIDATION.json','LOCALIZATION_AUDIT.json','SOFTWARE.json']:shutil.copy2(P/'qc'/f,P/'delivery'/f)
for f in ['run_tests.py','validate.py','audit_localization.py']:shutil.copy2(P/'code'/f,P/'delivery'/f)
primary.to_csv(P/'delivery/PRIMARY_Q75_SUMMARY.tsv',sep='\t',index=False)
validation=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());audit=json.loads((P/'qc/LOCALIZATION_AUDIT.json').read_text())
assert validation['status']=='PASS' and audit['usable_SPP1_macrophage_CD8_multiplex_datasets']==0
for f in schemas:assert len(pd.read_csv(P/'delivery'/f,sep='\t'))==0
assert len(todo)==38 and todo['sample'].is_unique
files=sorted(f for f in (P/'delivery').iterdir() if f.is_file() and f.name!='DELIVERY_SHA256.tsv')
pd.DataFrame([dict(file=f.name,bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]).to_csv(P/'delivery/DELIVERY_SHA256.tsv',sep='\t',index=False)
with zipfile.ZipFile(P/'phase276_delivery.zip','w',compression=zipfile.ZIP_DEFLATED) as zz:
 for f in sorted((P/'delivery').iterdir()):zz.write(f,arcname=f.name)
dump=dict(at=datetime.datetime.now().astimezone().isoformat(),status='COMPUTATION_VALIDATED_LOCALIZATION_INPUT_MISSING',delivery_files=len(files),figures=0,background_jobs=0)
(P/'STATUS.json').write_text(json.dumps(dump,indent=2))
print(json.dumps(dump,indent=2));print('ZIP_BYTES',(P/'phase276_delivery.zip').stat().st_size)
