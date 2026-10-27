from analyze import *
def table(d):
 cols=list(d.columns);lines=['| '+' | '.join(cols)+' |','|'+'|'.join(['---']*len(cols))+'|']
 for row in d.itertuples(index=False,name=None):
  vals=[]
  for v in row:
   if isinstance(v,(float,np.floating)):vals.append('—' if not np.isfinite(v) else f'{v:.6g}')
   else:vals.append(str(v))
  lines.append('| '+' | '.join(vals)+' |')
 return '\n'.join(lines)
def main():
 assert (P/'logs/exact.exit').read_text().strip()=='0'
 for f in ['qc/INDEPENDENT_VALIDATION.json','qc/EXACT_NULL_VALIDATION.json']:assert json.loads((P/f).read_text())['status']=='PASS'
 su=read(P/'tables/SPECTRAL_COHORT.tsv');ex=read(P/'tables/EXACT_NULL_SECTION.tsv');old=read(P/'tables/SPECTRAL_SECTION.tsv');ma=read(P/'tables/MATCH_COHORT.tsv')
 cmp=ex.merge(old[['sample','radius','model','endpoint','spatial_p','spatial_FDR']],on=['sample','radius','model','endpoint'],suffixes=('_exact','_approx'))
 tab(cmp,'NULL_MODEL_COMPARISON.tsv')
 rows=[]
 for keys,g in ex.groupby(['dataset','radius','model','endpoint']):
  v=g[g.effect.notna()]
  rows.append(dict(zip(['dataset','radius','model','endpoint'],keys),sections=len(v),patients=v.patient.nunique(),positive=int((v.effect>0).sum()),negative=int((v.effect<0).sum()),positive_spatial_nominal=int(((v.effect>0)&(v.spatial_p<.05)).sum()),negative_spatial_nominal=int(((v.effect<0)&(v.spatial_p<.05)).sum()),positive_spatial_FDR=int(((v.effect>0)&(v.spatial_FDR<.05)).sum()),negative_spatial_FDR=int(((v.effect<0)&(v.spatial_FDR<.05)).sum())))
 diag=pd.DataFrame(rows);tab(diag,'SPATIAL_NULL_COHORT_DIAGNOSTIC.tsv')
 pa=read(P/'tables/SPECTRAL_PATIENT.tsv');chain=[]
 for keys,g in pa.groupby(['dataset','radius','model']):
  wide=g.pivot(index='patient',columns='endpoint',values='effect')
  for end in ['joint_to_T_distribution','joint_to_CD8_distribution']:
   z=wide[['tumorTYMP_to_joint',end]].dropna()
   for patient,r in z.iterrows():chain.append(dict(zip(['dataset','radius','model'],keys),patient=patient,endpoint=end,upstream=float(r.iloc[0]),downstream=float(r.iloc[1]),directions_compatible=bool(r.iloc[0]>0 and r.iloc[1]<0)))
 ch=pd.DataFrame(chain);tab(ch,'SAME_PATIENT_DIRECTION_DETAILS.tsv')
 cs=ch.groupby(['dataset','radius','model','endpoint']).agg(patients=('patient','size'),compatible=('directions_compatible','sum')).reset_index();tab(cs,'SAME_PATIENT_DIRECTION_SUMMARY.tsv')
 q=read(P/'tables/MATCH_QUALIFICATION.tsv');qc=read(P/'tables/EXACT_NULL_QC.tsv')
 primary=su[(su.radius==200)&(su.model=='core')]
 primarydiag=diag[(diag.radius==200)&(diag.model=='core')]
 nominal=su[su.p_patient<.05]
 comp=cs[(cs.radius==200)&(cs.model=='core')]
 tj=ex[(ex.dataset=='TJRCC')&(ex.radius==200)&(ex.model=='core')&(ex.endpoint=='tumorTYMP_to_joint')&(ex.reason=='OK')]
 glm=ex[(ex.dataset=='GLMF_GSE210041')&(ex.radius==200)&(ex.model=='core')]
 mq=q.groupby(['radius','model']).agg(sections=('sample','size'),qualified=('eligible','sum'),median_pairs=('pairs','median'),maximum_pairs=('pairs','max')).reset_index()
 report=f"""# TYMP相关生态位：空间结构零模型与匹配区域补充检验
2026-09-30；Phase291。完整审计材料，全部方向和不可评估项保留。

## 本轮结论
两种补充方法没有把现有TYMP—SPP1相关生态位线索提升为稳定的免疫屏障参与证据。
TJRCC在200 μm主分析仍有9例中8例肿瘤侧TYMP—联合生态位正向，患者效应中位ρ=0.123034；9例只作探索性描述，不能因较前轮减少一例就迁移原n=10的P值。本轮{int((tj.spatial_p<.05).sum())}/9例通过严格协变量约束空间零模型的原始P<0.05，空间模型最小P={tj.spatial_p.min():.6g}。这说明方向线索仍在，但尚未明确超出本轮空间背景模型。
SPP1状态—间质活化在GSE175540的400 μm血管调整模型有名义正向关联（n=15，12例正向，中位ρ=0.224690，患者P=0.0353394，FDR=0.106018）。这是尺度和模型限定的协作线索，未同时获得TYMP上游与T细胞进入减少的对应支持。
严格非重叠匹配仅200 μm核心7片、血管4片通过平衡资格，400 μm均不足；各队列仅1–4例可评价主要匹配读出。匹配后的SPP1状态/联合指标有效对不足，不能把不可评估写成阴性。

## 与前轮的区别
1. 将近邻聚集结构纳入零模型。谱替身保留预设图上的Moran二次型和表达残差功率；严格版本还保留协变量正交约束，而非逐spot随机打乱。
2. 使用几何稀疏化减少重叠边界重复占权；匹配区进一步严格无共享spot。匹配暴露为肿瘤侧TYMP高低，而非重复Phase276的SPP1覆盖/缺口对照。
3. 总T采用映射后验下分位丰度代理的内外log差，不只依赖细胞比例；原CD8表达内外差独立保留。内外分布差是静态空间读出，不是直接迁移测量。
边界仍为Phase289计算边界；SPP1映射仍是旧SPP1相关巨噬状态，不改称当前C4。CAF活化为混合间质程序，不等同纯CAF数量。

## 患者层面的200 μm主要结果
effect为患者内偏秩相关，再在患者间取中位数；各队列独立报告。
{table(primary[['dataset','endpoint','patients','positive','negative','median_effect','p_patient','FDR_patient']])}

主要解释：
- TJRCC总T内外分布n=9，5负4正，中位ρ=-0.016630；CD8 n=4，3负1正，中位ρ=-0.073390。没有稳定的进入减少信号。
- GSE175540总T n=15，9负6正，ρ=-0.067562，P=0.330261；CD8 n=13，7负6正，ρ=-0.024279，P=0.946045。
- TYMP与CAF活化本身：TJRCC n=10、6正4负，ρ=0.009983、P=0.695312；GSE175540 n=20、9正11负，ρ=-0.006986、P=0.409098。不能把联合指标的变化自动解释为TYMP增加两个组成部分。
- GLMF两病例总T均负、TYMP—联合均正，但只有2病例；CD8一正一负。作为案例描述保留。

## 全部尺度/模型中患者原始P<0.05的结果
{table(nominal[['dataset','radius','model','endpoint','patients','positive','negative','median_effect','p_patient','FDR_patient']])}
完整{len(su)}行谱方法队列汇总见SPECTRAL_COHORT.tsv。患者统计只计算一次：两种零模型共享完全相同的观察效应和患者P，不算两个验证。

## 严格空间零模型结果
以下计数为切片内空间模型的名义/FDR信号，不是额外患者样本量；不据空间P选择患者，再做患者P。
{table(primarydiag[['dataset','endpoint','sections','patients','positive_spatial_nominal','negative_spatial_nominal','positive_spatial_FDR','negative_spatial_FDR']])}

GLMF两个案例全部主要结果：
{table(glm[['sample','endpoint','effect','spatial_p','spatial_FDR','regions']])}
全部正向和反向空间结果见EXACT_NULL_SECTION.tsv；原近似版本在SPECTRAL_SECTION.tsv，二者对照在NULL_MODEL_COMPARISON.tsv。任何单片空间P均不得称为生物学处理效应P。

## 同患者两段方向
这里仅检查同一患者是否同时TYMP—联合为正、联合—T/CD8内外差为负；没有做因果中介或合并P。
{table(comp)}

## 非重叠区域匹配
高低TYMP阈值固定为每切片几何合格区的上下四分位。先选中心距离>=2倍半径的区域，同组织连通分量做一对一无放回匹配，每协变量差<=1 SD、RMS差<=0.75；>=3对且最大配对后绝对SMD<=0.25才纳入。完整资格失败不放宽阈值。
{table(mq)}

200 μm核心匹配患者结果（差值=高TYMP减低TYMP）：
{table(ma[(ma.radius==200)&(ma.model=='core')][['dataset','endpoint','patients','positive','negative','median_effect']])}
TJRCC 4例匹配后CAF活化3例降低、1例升高；T内外分布差3例升高、1例降低，不支持一致的CAF增强/T减少。GLMF仅1例，GSE175540仅2例，总T方向相反。SPP1状态/联合端点因匹配对两端的髓系富集区可评估性不足，各队列无>=3有效对的合格患者，属于测量/匹配可行性失败，不作为反向证据。

## Methods
输入33个已合格切片，跨GLMF、TJRCC、GSE175540；原38片资格表在Phase289，不重新选择边界。200/400 μm的区域成员、六基因CAF面板、旧SPP1映射保持。仅从Phase234原Tcell_q05重聚合T读出，未读取/修改源H5AD计数。
谱分析采用坐标确定性贪心稀疏化（中心距离>=radius）；以有效组织ring图物理最短路径、4r截点、2r高斯带宽生成对称W。至少20区域、至少8残差维。核心模型控制两侧深度、巨噬和恶性比例、侧spot数比、x/y线性趋势、图连通分量；血管模型另加两侧内皮/周细胞比例。总T/CD8模型增加外侧相应读出，CD8另控制外侧总T丰度代理。
所有连续暴露/结局和组成协变量取秩；x/y为标准化坐标；D包含常数及图分量固定效应。严格版本以Q=null(D')构造V=Q eigenvectors(Q'WQ)，对a=V'xr独立符号翻转999次。以残差相关检验空间对齐，双侧空间P=(1+|null|>=|observed|-1e-12的次数)/1000。此条件谱实现保存协变量约束、范数和Moran二次型。
患者层面：同患者多片相关以Fisher-z等权汇总；匹配差值直接等权平均；各队列分别双侧Wilcoxon，n<10不作患者推断。组织家族3项、免疫家族2项按方法/队列/半径/模型分别BH；不可估计P不补造，登记家族规模与有效检验数均保留。空间P单独在切片/半径/模型/家族内BH，不与患者P混池。
空间谱思路参考[Wagner & Dray 2015](https://doi.org/10.1111/2041-210X.12407)及[adespatial官方文档](https://adeverse.github.io/adespatial/reference/msr.html)。本轮约束谱零模型为明确记录的Python实现，不声称直接运行原R函数。

## 质量核验与讨论
完整输入SHA保持；原分析与校验退出码均0。独立核验503个观察相关、每模型999空间替身、114040个区域T聚合及399个匹配对，患者聚合/P/FDR、匹配无重复spot及协变量平衡均通过。
严格零模型全部503模型：最大Moran误差{qc.max_moran_error.max():.3g}、最大协变量正交误差{qc.max_cov_error.max():.3g}、最大范数相对误差{qc.max_norm_error.max():.3g}。观察相关及患者P与原版本一致。
初版最短路径截点出现浮点非对称，数值闸门阻止产出；已按1e-6 μm容差修正并保留日志。原近似谱替身重投影后Moran有偏移（最大0.46245），因此补登记全模型严格约束版；两条记录均保留。极少数等效符号翻转与观察量相等，统一1e-12容差保守计入尾部，数值复核通过。
本轮为同批数据上已知前轮结果后的方法补充；方向可用于探索性组织定位，不能作为独立队列验证。空间零模型依赖图、谱对称及协变量设置，不证明生物因果。静态T内外差不能区分进入、迁出、存活、增殖和技术检出。原边界不是病理金标准，肿瘤侧TYMP和间质程序均为混合区域信号；未检出可重复负向关系也不等于否定所有TYMP参与机制。
当前可写：TYMP与SPP1相关生态位存在队列限定的局部方向线索；SPP1状态—间质活化有模型/尺度限定的协作信号。当前数据不支持进一步写成TYMP促进SPP1—CAF屏障、导致T细胞进入减少。后续若换证据，应优先改变数据或测量层级，而非继续调整同批相关模型的阈值。

## 交付
完整代码、分析前方案、质量驱动补登记、全患者/切片/匹配统计及校验记录随附；大体积区域与谱缓存保留权威服务器inputs。无新图、未修改Phase273正式稿、原C4成员及源数据不变。
"""
 (P/'补充方法检验与完整结果.md').write_text(report,encoding='utf-8')
 dump(P/'qc/COMPLETION.json',dict(status='COMPLETE_VERIFIED',analysis_exit=0,validation_exit=0,exact_null_exit=0,figures=0,manuscript_unchanged='Phase273',source_data_unchanged=True))
 print('REPORT_READY')
 print(primarydiag.to_string(index=False))
 print(comp.to_string(index=False))
if __name__=='__main__':main()
