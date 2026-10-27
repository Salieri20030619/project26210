from pathlib import Path
import pandas as pd,numpy as np,json,datetime,shutil,hashlib,zipfile
P=Path('phase275_immune_interface_tests_20260928_v1')
s=pd.read_csv(P/'tables/COHORT_SUMMARY.tsv',sep='\t')
def fmt(x):
 if pd.isna(x):return '—'
 if isinstance(x,str):return x
 return f'{x:.6g}'
def table(d,cols,labels):
 out=['| '+' | '.join(labels)+' |','| '+' | '.join(['---']*len(cols))+' |']
 for _,r in d.iterrows():out.append('| '+' | '.join(fmt(r[k]) for k in cols)+' |')
 return '\n'.join(out)
primary=s[s.width_um.eq(400)].copy()
names={'band_primary_readouts__T_abundance_log2':'总T映射丰度（log2）','band_primary_readouts__CD8_RNA':'CD8谱系RNA','band_auxiliary_readouts__T_RNA':'总T谱系RNA','band_auxiliary_readouts__CD8_codetection':'CD8B/T标记共检出','SPP1_coupling__T_abundance_log2':'SPP1与总T内外差','SPP1_coupling__CD8_RNA':'SPP1与CD8内外差'}
models={'unadjusted':'未调整','depth_macro_adjusted':'深度/总巨噬调整','depth_macro_composition_adjusted':'深度/巨噬/组成/几何调整'}
for d in [s,primary]:
 d['读出']=d.endpoint.map(names);d['模型']=d.model.map(models)
 d['内低人数或负相关人数']=d.negative.astype(int).astype(str)+'/'+d.n.astype(int).astype(str)
cols=['dataset','读出','模型','n','内低人数或负相关人数','median','P','FDR']
labels=['队列','读出','模型','有效患者数','负向/总数','中位效应','原始P','FDR']
text='''# 26210：免疫界面检验结果（Phase275）

完成日期：2026-09-28。任务范围：先完成统计检验，视结果决定图件。本轮未出图、未修改正文；正式稿仍为Phase273。

## 核心结论

**现有数据支持肿瘤表达区内侧CD8/T谱系信号较低的空间分区，但当前数据不支持“界面SPP1巨噬越富集，T/CD8进入越少”的特异关联。** 可以形成免疫空间分布结果，尚不能把主图的机制结论定为“SPP1巨噬构建屏障并阻挡T细胞”。

最稳定的结果是GSE175540的CD8谱系RNA内外差异。400 µm界面带内，未经调整16/20位患者内侧较低，中位差−0.2821，P=0.011818、FDR=0.023636；调整深度和总巨噬丰度后19/20位内侧较低，中位差−0.2602，P=0.00004005、FDR=0.0001602。200/800 µm同样19/20位方向一致。这里的单位是固定基因各自标准化后均值的评分单位，不是CD8细胞比例，也不是完整评分再标准化后的一个标准差。

SPP1特异性检验未支持预期方向：GSE175540调整后SPP1局部巨噬占比与CD8内外差的中位相关为+0.0493（P=0.1153，FDR=0.4613），与总T内外差为+0.0420（P=0.4049，FDR=0.8097）。预期“SPP1高、进入少”需要负相关；目前既缺乏显著证据，也没有稳定负向。TJ-RCC亦未支持该关联。

## 数据、资格与预设方案

共38张切片、37位患者：GSE175540 24张/23位（c_2的两种制备按同一患者合并）；TJ-RCC 12张/12位；GLMF 2张/2位。各队列单独总结。复用Phase234冻结的恶性细胞映射比例≥0.20表达分区、原物理图和连通分量；主距离400 µm，敏感性200/800 µm。未移动分区阈值、未挑选最佳距离或最佳切片。

400 µm内外差与局部耦合的总T有效患者数分别均为23/11/2（GSE175540/TJ-RCC/GLMF）；CD8为20/9/2。TJ-RCC的R51_T外侧带仅15个spot，低于每侧至少20个的标准；CD8另有全切片标记覆盖不合格，均保留在资格表中。局部中心每侧至少5个spot，每张切片至少50个合格中心。不能把这些spot/中心数量当独立生物学样本量。

独立患者数≥10时，事先指定双侧精确符号检验为主，Wilcoxon为伴随检验；不足10只作方向与效应描述。患者中位数CI来自10000次患者重抽样。主宽度按队列×分析模块分别校正4个假设（2读出×2模型）；敏感性把200和800 µm共8个假设放在对应家族。不可估计假设以P=1保留于校正分母，最终P/FDR显示为空。共108个队列端点、18个校正家族，未混入既有Phase234/246的检验池。

详细冻结方案见PROTOCOL_BEFORE_TESTING.md；新结果已知前保存时间与SHA见FREEZE.json。本研究是已知前期结果基础上的探索性扩展，不能称前瞻性预注册研究。

## 界面内外的主要读出

负差值表示肿瘤表达区内侧低于外侧。总T为映射q05丰度取log2后比较，未使用组成比例作主要结局；CD8为CD3D/CD3E/TRAC/CD8A/CD8B固定RNA评分。内外带每侧至少20个spot、合计至少50个；调整模型为每切片结局对内外侧指示变量、log1p深度、log1p总巨噬q05的线性模型，取侧别系数后再按患者汇总，不使用spot层面的回归P值。

'''+table(primary[primary.endpoint.str.startswith('band_primary')],cols,labels)+'''

GSE175540总T未经调整的方向为14/23位内侧低，P=0.4049；调整后17/23位内侧低，中位log2系数−0.4758，对应条件性几何丰度比约0.719，P=0.03469、FDR=0.04625。这是依赖协变量调整的证据，不能写成未经调整的总T已普遍减少。TJ-RCC调整后9/11位同向，中位log2系数−0.2380，主检验P=0.06543、FDR=0.2617；伴随Wilcoxon P=0.01367、FDR=0.05469，为方法依赖的名义信号，不能替代未通过的预设主检验。

TJ-RCC调整后CD8 9/9位内侧更低，中位差−0.1385；按n<10规则仅作描述。GLMF两例均显示总T/CD8内侧较低，但不足以独立推断，也没有处理组间检验。

## 辅助读出

CD8共检出要求同一spot检出CD8B及至少一个CD3D/CD3E/TRAC；沿用全切片既有标记资格，不要求内外两侧各有最低阳性数，以免把真实低信号区域排除。此指标是混合spot的转录本共检出，仍不是确认同一个CD8细胞。调整后的二元线性模型系数表示条件性检出比例差，不能外推为细胞占比。

'''+table(primary[primary.endpoint.str.startswith('band_auxiliary')],cols,labels)+'''

GSE175540调整后CD8共检出20/20位内侧较低，中位差−0.03888（约−3.89个百分点），P=1.907e−6、FDR=7.629e−6；总T RNA为22/23位同向。两者与CD8程序共同支持内外RNA空间分区，但来自同一表达数据，不能称独立实验验证。

## SPP1是否解释较低的内侧免疫信号

在冻结的表达界面、总巨噬比例≥0.05的中心上，用200 µm源邻域内SPP1 q05之和/总巨噬q05之和定义SPP1巨噬占比；结局是中心周围400 µm内侧减外侧的总T log2丰度或CD8评分。SPP1与内外差为负相关才符合预期排斥方向。调整固定包括邻域平均深度、内外深度差、源邻域总巨噬丰度、恶性比例均值及内外差、内外spot数量比；每张切片得到秩相关或偏秩相关，同患者Fisher-z等权合并。

'''+table(primary[primary.endpoint.str.startswith('SPP1')],cols,labels)+'''

TJ-RCC未调整CD8局部耦合9/9位是正向，中位rho=+0.0833；调整后5正4负，中位+0.0810。这是必须保留的反向描述，不可包装为SPP1排斥。GLMF1调整后CD8耦合−0.0602，而GLMF2为+0.2725；两例不一致。不能选GLMF1单例替代整体证据。

## 固定尺度与患者稳健性

'''+table(s[(s.dataset=='GSE175540')&(s.model!='unadjusted')&(~s.endpoint.str.contains('auxiliary'))],['读出','width_um','n','内低人数或负相关人数','median','P','FDR'],['读出','宽度µm','患者数','负向/总数','中位效应','P','FDR'])+'''

CD8调整后的200/400/800 µm中位差依次−0.2396/−0.2602/−0.2845，均19/20位同向；三个宽度的共同患者集相同。逐个剔除患者后，400 µm CD8中位效应始终为负（−0.2644至−0.2560），主P最大7.629e−5。总T调整后方向也在各尺度保持，但400 µm留一患者主P可达0.05248，其证据强度弱于CD8。

SPP1耦合在固定宽度间不稳定：GSE175540调整后CD8 rho在200/400/800 µm分别+0.0051/+0.0493/−0.0607，均无主检验支持。限制到三宽度共同患者集（CD8 n19、总T n22）仍发生方向变化，并非仅由宽度资格变化造成。200 µm总T调整后伴随Wilcoxon P=0.02292、FDR=0.1661，但为正向且主P=0.05248；800 µm未调整总T伴随P=0.04152、FDR=0.1661，但主P=0.09314。完整保留这些名义信号，不改换主检验来构造SPP1屏障结论。

完整主/伴随P与FDR、CI在COHORT_SUMMARY.tsv；全部患者、切片效应、留一法与共同患者尺度分析分别见对应表。

## 与之前SPP1/CD8结果的关系

Phase246的GSE175540结果曾显示：SPP1高区域到CD8富集区域距离增加，中位差+300 µm，n16，P=0.03516/FDR=0.05273；400 µm CD8邻域负相关调整后减弱至rho约−0.0916，P=0.2632/FDR=0.3948。旧结果衡量的是SPP1区域与CD8区域的相对位置，本轮衡量的是表达分区内外差及这种差异是否随局部SPP1变化，二者不是同一个命题。

新增检验补上了“内外侧”的证据，但没有补上“SPP1特异驱动”的证据。旧距离结果不能替代当前局部耦合阴性；既有CAF/ECM与总T绝对丰度的非一致结果同样保留，不拼接成已证实的物理屏障。原P/FDR原样保存在EXISTING_*表，不重新缩小校正家族。

## 推断边界与后续定位

边界来自同一转录组的恶性细胞映射，尚无审核病理轮廓；恶性比例、免疫映射和谱系RNA来自同一混合spot数据。组成/深度调整降低部分混杂，不能消除分区定义的依赖或证明真实T细胞被挡住。CD8B/T共检出亦受捕获深度影响，局部邻域相互重叠；本轮不对spot或中心作独立样本检验。

当前最可成立的结果定位是“肿瘤表达区内侧CD8/T谱系信号降低，形成免疫空间分区”。若后续做机制主图，需要以审核病理轮廓和独立CD8/CD3、SPP1/巨噬标记定位检验界面结构；物理阻挡或SPP1因果作用需要相应干预证据。本轮不追加这些尚未授权或尚无数据的实验，不以C4名称替代旧空间SPP1-state。

## 计算与交付验收

38/38切片正常完成，退出状态0；1368条切片效应、1332条患者效应、108条队列端点。独立数学精确二项式P、手工BH、同患者合并、statsmodels侧别系数、QR偏秩相关、直接邻域几何及标记重建通过。406项数值交叉核验最大误差3.78e−15，标记重建最大误差8.89e−16。原输入派生文件SHA与原H5AD元数据检查通过。

复核过程中发现CSV默认浮点读取会改变个别并列秩，导致约4e−7的复算差异；验证器改用round_trip精度读取后与原分析一致。未修改原分析、端点或统计阈值，记录见VALIDATION_READER_FIX.txt。FINAL_VALIDATION.json为完整验收记录。

所有新结果保存在服务器Phase275目录；本地交付包含报告、冻结方案、完整统计/资格/QC表、分析和独立复核脚本、原有证据对照。大体积逐中心中间量及原始数据仍保留服务器。图件0，主稿修改0。
'''
(P/'delivery/免疫界面检验报告.md').write_text(text,encoding='utf-8')
for f in (P/'tables').glob('*.tsv'):shutil.copy2(f,P/'delivery'/f.name)
for f in ['PROTOCOL_BEFORE_TESTING.md','FREEZE.json']:shutil.copy2(P/f,P/'delivery'/f)
for f in ['FINAL_VALIDATION.json','VALIDATION_READER_FIX.txt']:shutil.copy2(P/'qc'/f,P/'delivery'/f)
for f in ['run_tests.py','validate.py']:shutil.copy2(P/'code'/f,P/'delivery'/f)
for f in ['EXISTING_CD8_STROMA_EVIDENCE.tsv','EXISTING_TCELL_BALANCE_EVIDENCE.tsv']:
 shutil.copy2(Path('phase274_client_immune_exclusion_priority_20260928_v1/delivery')/f,P/'delivery'/f)
primary.to_csv(P/'delivery/PRIMARY_400UM_SUMMARY.tsv',sep='\t',index=False)
files=sorted(f for f in (P/'delivery').iterdir() if f.is_file() and f.name!='DELIVERY_SHA256.tsv')
pd.DataFrame([dict(file=f.name,bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]).to_csv(P/'delivery/DELIVERY_SHA256.tsv',sep='\t',index=False)
with zipfile.ZipFile(P/'phase275_delivery.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted((P/'delivery').iterdir()):z.write(f,arcname=f.name)
(P/'STATUS.json').write_text(json.dumps(dict(status='VALIDATED_AWAITING_LOCAL_DELIVERY',at=datetime.datetime.now().astimezone().isoformat(),figures=0,files=len(files)),indent=2))
print('DELIVERABLES',len(files),'ZIP_BYTES',(P/'phase275_delivery.zip').stat().st_size)
