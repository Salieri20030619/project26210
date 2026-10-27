from common import *
import shutil
def fmt(x):
 if pd.isna(x):return '—'
 if isinstance(x,(int,np.integer)):return str(x)
 if isinstance(x,(float,np.floating)):return f'{x:.4g}'
 return str(x).replace('|','/').replace('\n',' ')
def mdtable(d):
 out=['| '+' | '.join(d.columns)+' |','| '+' | '.join(['---']*len(d.columns))+' |']
 out += ['| '+' | '.join(fmt(v) for v in row)+' |' for row in d.itertuples(index=False,name=None)]
 return '\n'.join(out)
def main():
 b=pd.read_csv(P/'tables/A1_COHORT_SUMMARY.tsv',sep='\t',float_precision='round_trip')
 n=pd.read_csv(P/'tables/A2_A4_COHORT_SUMMARY.tsv',sep='\t',float_precision='round_trip')
 raw=pd.read_csv(P/'tables/A1_SECTION_EFFECTS.tsv',sep='\t',float_precision='round_trip')
 regs=pd.read_csv(P/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t')
 cohorts=['GSE175540','TJRCC']
 qcols=['dataset','exposure','target','scale','model','n','positive','negative','median','CI_low','CI_high','P','P_wilcoxon','FDR']
 names={'dataset':'来源','exposure':'来源变量','target':'目标','scale':'尺度','model':'调整','n':'患者n','positive':'正向','negative':'负向','median':'效应中位数','CI_low':'CI下限','CI_high':'CI上限','P':'主P（符号）','P_wilcoxon':'伴随P（Wilcoxon）','FDR':'主FDR'}
 def table(d,cols=qcols):return mdtable(d[cols].rename(columns=names))
 w=b[(b.layer=='broad')&(b.scale=='W0')&(b.dataset.isin(cohorts))]
 chunks=['# 第四部分重订方案执行结果：广域共变化与局部空间结构分层\n\nPhase234｜2026-09-21｜基于Phase233锁定方案｜服务器权威版本\n',
 '''## 执行结论

多胺/MTA转录背景与SPP1相关巨噬绝对丰度呈跨来源一致的区域共变化，且在固定400μm邻域仍保留；巨噬内SPP1状态比例未同步升高。SPP1与血管的关系主要表现为绝对丰度的共同分布，在构成比例或完整条件模型中并未获得一致支持。当前数据不支持将主线写为普遍的SPP1—CAF/ECM邻近增强。

本轮完成38片、114640个QC spot的输入与物理标定，70个宽层端点、固定800/400/200μm与第一圈诊断、400μm资格及逐层调整、表达交界距离分析、9个SPP1及2个POSTN配受体表达兼容性分析，并分来源整理结果。GLMF两例保持主体病例展示地位；GSE175540按23名登记患者、TJ-RCC按12名登记患者单独汇总。不新增图，不修改旧RDS、既有C4成员或Phase208正文。

## Methods：数据、变量与统计

输入复用Phase229已核对的38片计数层、条码、图像范围和邻接图。GLMF1/2与GSE210041为同一两例；GSE175540的c_2双制备按患者合并。旧13状态映射的绝对q05从原后验表读取，未由比例反推。V1=SPP1-state绝对q05，V2=SPP1-state/13状态总和，V3=SPP1-state/三类成熟巨噬总和，V4=去SPP1的冻结表达程序。血管目标同时保留q05和13状态归一化比例；总巨噬、C1QC及APOE程序为比较对象。

程序采用原始counts按每spot总量归一到10000后log1p，在每片内逐基因计算总体标准差Z分数，再取等权均值。原冻结基因集保留，四条KEGG及custom代谢程序均移除TYMP。KEGG与custom分开校正，旧GSVA结果继续按旧方法独立保留；本轮均明确标为Phase234平均Z程序。程序资格为至少半数且至少3个基因可用（两基因集至少2个），至少20个spot检出至少2个程序基因、评分方差非零；完整实际基因和缺失项记录于资格表。

W0全组织、W1坐标秩分位4×4/5×5块中位数分别计算相关。至少50个有效spot；粗块至少8个占用块；相关两变量各至少3个不同取值。W2沿原组织邻接图取800/400/200μm可达邻域，W3为原第一圈，均排除自身，至少3邻居与50来源spot。Q0为宽层有限分母，Q1巨噬构成≥0.05，Q2≥0.10。尺度诊断按Q1固定；代谢的恶性富集组织层另以原恶性构成≥0.20定义。

所有38片均通过原始array_row/array_col六角阵列边及图像仿射配套核验。物理标定采用标准Visium阵列水平间距100μm、行间距50√3μm，第一圈边长100μm；物理图路径不能跨断开的组织。图像仿射残差最大0.527μm，轴比例约1.0046。原spot_diameter_fullres元数据换算导致约85μm邻点间距的差异已记录，未将其作为本轮半径单位。依据为[10x Visium官方技术说明](https://cdn.10xgenomics.com/image/upload/v1758134784/support-documents/CG000686_VisiumHDSpatialApplicationsPerformance_TechnicalNote_RevB.pdf)。完整资格见audit/PHYSICAL_CALIBRATION.tsv及tables/SCALE_ELIGIBILITY.tsv。

400μm调整按A0无协变量→A1来源/邻域深度→A2加入组成→A3加入来源位置同一目标程序/变量；POSTNCAF/ECM另有A4坐标二阶项。V1组成使用排除SPP1的其他巨噬构成，V3使用总巨噬构成，均加入邻域对应构成及恶性构成。采用秩残差相关；协变量满秩、条件数<30。阈值和模型差异分别报告同患者配对效应变化与bootstrap区间，不用P跨过0.05代替效应差异检验。

每片一个效应；同患者相关系数Fisher-z等权合并，差值/斜率算术等权合并。患者≥10时以双侧精确符号检验为主、Wilcoxon为固定伴随，10000次患者bootstrap给出中位效应区间；不足10只作方向和效应描述。FDR在来源×同类端点家族×分析层内计算，平行变量和该层全部尺度/模型纳入；不可估计行在校正向量中占P=1但显示P/FDR为空。配对变化亦按原端点家族分开校正。全部P、伴随P、FDR、有效n及不可估计原因均保留。

表达交界沿冻结恶性构成分区确定，分别比较内/外侧≤800/400/200μm带与同侧>800μm固定参照，每区至少20spot；每200μm分箱连续曲线另存。配受体名单继承冻结CellChatDB：9个SPP1配对及2个POSTN对照；所有受体亚基均需存在且在至少20spot共同检出，配体至少20spot检出，以最小亚基归一表达代表复合体，检验来源配体与固定邻域受体的表达兼容性。
''',
 '## Results 1：TYMP相关区域保留嘧啶与甘油磷脂程序联系\n\nGSE175540四条KEGG程序均呈正向全组织关联；TJ-RCC对嘧啶和甘油磷脂的方向一致性得到校正后支持，另两条未获得同等支持。效应总体较小，体现组织内连续共变化。\n\n'+table(w[w.family=='KEGG_four_axes']),
 '## Results 2：多胺/MTA背景对应SPP1绝对丰度，而非一致的巨噬内状态偏移\n\nGSE175540全组织23/23患者、TJ-RCC12/12患者的MTA-context—SPP1绝对q05效应为正，中位rho分别0.2256与0.1634；主P分别2.384e-7与0.0004883，FDR分别7.153e-7与0.005859。粗区域结果同方向，TJ-RCC的5×5为名义显著。固定400μm邻域仍为23/23和11/12正向。V4程序在全组织两来源亦有支持。\n\nV3巨噬内比例在全组织GSE175540仅14/23正、TJ-RCC3/12正，主P均未达0.05。GSE175540该比例关联在FFPE中10/12正、frozen中5/12正，显示制备/来源差异；c_2跨制备仍是一名患者。\n\n'+table(b[(b.layer=='broad')&(b.family=='MTA_context_SPP1')&b.dataset.isin(cohorts)])+'\n\n固定400μm：\n\n'+table(n[(n.layer=='scale')&(n.family=='MTA_context_SPP1')&(n.scale=='400')&n.dataset.isin(cohorts)]),
 '## Results 3：SPP1与血管的绝对数量关系不同于构成比例关系\n\nGSE175540的SPP1绝对q05与内皮、周细胞绝对q05均为19/23正，中位rho0.3539/0.3868，FDR均0.01248。TJ-RCC相应10/12正，为名义显著，FDR0.2057。总巨噬绝对丰度也与血管绝对丰度共同变化，说明该特征并非已经证明SPP1状态特异。\n\n同一SPP1绝对丰度换成内皮构成比例目标时，TJ-RCC12/12负、中位rho−0.1227，FDR0.007812；这种反向结果必须与数量层面的正相关同时呈现。400μm完整条件模型未建立统一正向联系。\n\n'+table(w[(w.family=='SPP1_vascular')&w.exposure.isin(['V1_SPP1_q05','V3_SPP1_macro_share'])])+'\n\n400μm无调整与A3：\n\n'+table(n[(n.layer=='adjustment')&(n.family=='SPP1_vascular')&n.model.isin(['A0','A3'])&n.dataset.isin(cohorts)]),
 '## Results 4：CAF/ECM伙伴未显示一致的SPP1丰度或状态正向邻近\n\nGSE175540与TJ-RCC的SPP1绝对丰度、巨噬内比例对POSTNCAF/ECM的全组织和固定400μm结果均未获得一致正向支持。GSE175540的V4—POSTNCAF全组织为19/23正、中位rho0.0880，P0.002599、FDR0.06239，属于程序层名义信号，不能替代丰度或空间细胞身份结果。\n\nGLMF两病例的V3—POSTNCAF/ECM全组织均为负，其主体展示地位不预设正向伙伴。完整A0→A4及同患者效应变化在附件，不据某一层P变化宣称中介或混杂已确认。\n\n'+table(w[w.family=='SPP1_stroma_program'])+'\n\n400μm V3逐层调整：\n\n'+table(n[(n.layer=='adjustment')&(n.family=='SPP1_stroma_program')&(n.exposure=='V3_SPP1_macro_share')&n.dataset.isin(cohorts)]),
 '## Results 5：SPP1程序与免疫、侵袭背景共变化，巨噬内比例呈不同功能方向\n\nV4去SPP1程序与髓系免疫抑制、T细胞功能障碍和侵袭/EMT程序在两个来源的全组织层均呈校正后支持的正向关联。绝对丰度也与部分功能程序同向，但细胞毒性和功能障碍程序可同时增强，说明这里反映组织生态背景，不能直接等同T细胞被抑制。\n\nV3与细胞毒性程序在GSE17554019/23负，rho−0.0984、P0.002599/FDR0.01040；TJ-RCC10/12负，rho−0.0499、P0.03857/FDR0.08050，为名义信号。GSE175540的V3与T细胞功能障碍程序也为18/23负，而非功能障碍升高。\n\n'+table(w[(w.family=='SPP1_function')&w.exposure.isin(['V3_SPP1_macro_share','V4_SPP1_program'])]),
 '## Results 6：表达交界附近出现基质增强与SPP1巨噬内比例降低的探索性方向\n\n固定肿瘤侧400μm带相对>800μm核心，GSE175540合格9患者中V3为8负1正，中位差−0.1483（−14.83百分点）；ECM8正1负，中位差0.3625评分单位。TJ-RCC对应7患者V3为4负3正，中位差−0.0570；ECM6正1负，中位差0.1937。GLMF两例同向表现为V3较低、POSTNCAF/ECM较高。各端点患者不足10，本轮只报告探索性方向和效应，不给新推断P。\n\n外侧同侧>800μm参照更少；缺少固定核心的切片未通过改距离补入。全部200/400/800μm及200μm分箱记录保留。\n\n'+table(n[(n.layer=='expression_boundary')&(n.scale=='400')&n.exposure.isin(['V3_SPP1_macro_share','POSTNCAF','ECM4'])]),
 '## Results 7：SPP1配受体兼容性呈轴与来源差异，CD44未获统一支持\n\n400μm SPP1—CD44在GSE175540与TJ-RCC均未获得患者层支持。SPP1—ITGAV/ITGB1在TJ-RCC11/12正，中位rho0.2082，P0.006348、FDR0.05713；GSE17554016/23正、P0.09314，未达到名义显著。SPP1—ITGAV/ITGB3在GSE175540可估计16患者中13正，P0.02127、FDR0.1784；TJ-RCC因完整复合体检出不足不可估计，不能当成零效应。相反方向的ITGA8/ITGB1结果完整保留。\n\nPOSTN对照不是SPP1链的替代；固定400μm其ITGAV/ITGB5在GSE175540得到支持，TJ-RCC为名义显著。全部9+2候选和4尺度均登记，包括未通过检出资格者。\n\n'+table(n[(n.layer=='LR_compatibility')&(n.scale=='400')&n.dataset.isin(cohorts)]),
 '''## Results 8：旧证据保持原定义，与新问题并列

- Phase222/223旧SPP1-state比例的第一圈血管主终点阴性保留；Phase230/231第一圈POSTNCAF/ECM及调整诊断原结果保留。本轮全组织和固定物理半径没有改写这些端点。
- Phase232确认旧GLMF绝对丰度、TJ-RCC复合niche及TYMP-MTA core的显著/名义结果保留。复合niche已包含血管，T-cell balance的其他免疫项含SPP1，不能从复合关联推导组成成分两两正向。
- Phase219旧方法去共享LGALS3后，SPP1—血管程序条件关联衰减、不再获支持；本轮统一平均Z评分且无调整的广域程序关系另列，不能反向宣称旧条件模型恢复。新去共享基因广域结果、原含共享基因及原条件模型均保留，具体算法与调整层明确分开。
- TLS严格注释补充继续引用12患者10负、中位−7.27百分点、P0.03857；主要来自BIONIKK，不改为CAF屏障、CD8排斥或TLS形成抑制证据。
- 旧TYMP-high/low参考、COMMOT运输量及MST1/ICAM5等原阴性/不可估计记录留在原目录，不冒充本轮新13状态或新SPP1候选结果。
''',
 '## Discussion：结论边界与叙事定位\n\n本轮支持将第四部分主线写为“代谢转录背景与SPP1相关巨噬的区域共分布，并伴随免疫/侵袭程序及特定尺度结构差异”。绝对数量、组织构成和巨噬内状态是不同事实；数量关联不能自动推为状态诱导。广域与局部关系也不能互相否定。多胺/MTA-context为MAT2A/MAT2B/AMD1/SRM/SMS转录背景，不是实测MTA、代谢净通量或空间中介。客户湿实验仅作实验背景，不进入本轮公共效应或方向筛选。\n\n所有相关来自混合spot、既有映射与转录程序；未验证CAF/CD8独立空间身份，不能声称细胞屏障、物理排斥或确定发送/接收细胞。边界为表达定义，不是已审核病理侵袭前沿。GLMF每治疗层仅一例，仅作病例背景；未知治疗不当作未治疗。登记来源去重不等于跨来源分子身份查重，不能称本轮建立全新独立验证队列。\n\n绝对q05仅在同片估效应，跨平台不比较绝对丰度均值。邻域模型均为关联；加入来源程序可移除共享组织结构，也可能调整真实组织过程，不能一概视为唯一真实效应。新方案属于已知历史结果基础上的锁定探索性执行。有限层级结束后不因P未显著扩展半径、换标签或重选患者。\n',
 '''## 后续图件内容方案（尚未绘图）

图组1：GLMF两例完整组织上的TYMP、四代谢程序、MTA-context及SPP1 V1/V3/V4；补充GSE175540/TJ-RCC患者效应与全部方向。GLMF以完整视野为主，不挑阳性ROI。

图组2：同一目标下绝对丰度/构成/巨噬内比例的并列结果，W0、4×4、5×5和固定800/400/200/第一圈尺度变化；CAF/ECM阴性及血管反向保留。400μm A0→A3用同患者效应变化呈现。

图组3：功能程序、表达交界200μm距离剖面及配受体候选资格/兼容性。病理边界、CAF/CD8身份不合格者不画成已验证细胞结构。TLS独立补充。每张图按A4实际单页尺寸构思，长图注另存；内容审核后再绘图。

## 交付目录与审计入口

- A1_SECTION/PATIENT_EFFECTS、A1_COHORT_SUMMARY：全部宽层与恶性富集层。
- A2_A4_SECTION/PATIENT_EFFECTS、A2_A4_COHORT_SUMMARY：全部尺度、资格、调整、表达边界及LR。
- PAIRED_LAYER_CHANGES、PAIRED_DENOMINATOR_CHANGES：共同患者的逐层/变量定义变化。
- INPUT_SECTION_REGISTER、FROZEN_ENDPOINTS、PROGRAM_LR_ELIGIBILITY、SCALE_ELIGIBILITY、FROZEN_LR_CANDIDATES：输入、基因/复合体和技术资格。
- GSE175540_STRATIFIED_DESCRIPTION及两类PATIENT_EFFECTS：制备和原来源分层，主推断仍为去重患者。
- EXPRESSION_BOUNDARY_BIN_PROFILES、LR_SPATIAL_COVERAGE：距离连续描述和表达兼容性覆盖。
- TECHNICAL_SCOPE_AND_UNESTIMABLE：不可估计及独立补充边界；正式病理、CAF/CD8特异细胞定位和细胞类型LR均未获资格。MISTy为计划中可选补充，未新增拟合；旧COMMOT结果保留，不将本轮兼容性改名为COMMOT。
- legacy_tables及LEGACY_EVIDENCE_REGISTER：原定义统计与原报告索引。
- FINAL_VALIDATION.json、audit/INDEPENDENT_PATIENT_SIGN_CHECKS.tsv、audit/NUMERICAL_EFFECT_CHECKS.tsv、audit/GRAPH_INTEGRITY.tsv：验收证据；scripts含全部本轮代码。
''']
 (P/'PART4_RESULTS_CN.md').write_text('\n\n'.join(chunks),encoding='utf-8')
 brief='''# 26210第四部分：本轮分析摘要
Phase234｜2026-09-21

已按重订方案完成可执行的A0–A5分析与结果整理；无新图。38片、37个登记患者单位，GLMF2例/GSE17554023患者/TJ-RCC12患者分开。

1. **多胺/MTA转录背景与SPP1绝对丰度的区域关联最一致。** 全组织GSE17554023/23正，rho中位0.2256、P2.384e-7/FDR7.153e-7；TJ-RCC12/12正，rho0.1634、P0.0004883/FDR0.005859。400μm亦同向。这里是转录背景，不是实测MTA。
2. **巨噬内状态比例未同步支持。** 同一MTA背景对SPP1巨噬内比例，两来源主P均未显著；不能把绝对丰度关联写成SPP1状态诱导。
3. **血管关联取决于所测变量。** SPP1与血管绝对丰度可共同增加，换构成比例或完整条件模型后减弱/反向。CAF/ECM未显示一致的SPP1丰度或状态正向邻近。
4. **功能与边界有可解释的分层。** SPP1程序与免疫抑制/侵袭背景共变化；巨噬内比例与细胞毒性程序偏负。表达交界附近探索性显示基质程序较高、SPP1巨噬内比例较低；相关有效n<10，只作方向和效应描述。
5. **LR保持候选层级。** CD44未获统一支持；部分整合素有来源相关名义线索。完整受体复合体检出不足明确写不可估计，不当作阴性。当前只测混合spot兼容性，不确认发送/接收细胞。

物理尺度已通过原始阵列坐标和图像配套核对；旧spot直径元数据约85μm间距差异留档，本轮固定标准阵列100μm间距。原阴性、旧复合阳性、共享LGALS3衰减和TLS补充均保留。

完整Methods/Results/Discussion及全部主P、伴随Wilcoxon P、FDR、效应、有效n见PART4_RESULTS_CN.md和统计附件。病理审核、独立CAF/CD8空间身份等未获资格项目单列。MISTy可选补充未新增；图件须内容审核后再制作。
'''
 (P/'SUMMARY_CN.md').write_text(brief,encoding='utf-8')
 # Standalone unified section4 prose for later manuscript integration; old manuscript remains unchanged.
 (P/'MANUSCRIPT_SECTION4_DRAFT_CN.md').write_text('\n\n'.join(chunks[:11]),encoding='utf-8')
 # Result register includes every estimable and unestimable endpoint, no significance filtering.
 allsum=pd.concat([b,n],ignore_index=True);allsum['evidence_level']=np.where(allsum.n<10,'exploratory_descriptive',np.where(allsum.P<.05,np.where(allsum.FDR<.05,'FDR_supported_association','nominal_association'),'current_data_not_supported'))
 tab('tables/ALL_ENDPOINT_EVIDENCE_REGISTER.tsv',allsum)
 print('REPORT_WRITTEN',len(allsum),'summaries',flush=True)
if __name__=='__main__':main()
