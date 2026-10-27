from pathlib import Path
import pandas as pd,numpy as np,json,datetime,hashlib,shutil,zipfile
P=Path(__file__).resolve().parents[1];E=P/'output_function_extension';D=P/'delivery'
s=pd.read_csv(P/'tables/WITHIN_PATIENT_SUMMARY.tsv',sep='\t');src=pd.read_csv(P/'tables/INHERITED_OUTPUT_ENRICHMENT.tsv',sep='\t').set_index('gene')
rec=pd.read_csv(P/'tables/RECEIVER_COMPLEX_LIBRARY_QUALIFICATION.tsv',sep='\t');db=pd.read_csv(P/'tables/DATABASE_COVERAGE.tsv',sep='\t').set_index('gene');fn=pd.read_csv(E/'tables/OUTPUT_CD8_FUNCTION_SUMMARY.tsv',sep='\t')
rows=[]
for r in s[(s.group=='C4')&(s.model=='depth_complexity_adjusted')].itertuples():
 old=src.loc[r.gene];positive=pd.notna(r.P) and r.P<.05 and r.median>0
 row=dict(program=r.program,gene=r.gene,source_n=int(old.n),source_delta=old.median_delta,source_P=old.P,source_FDR=old.FDR,source_nominal_enriched=bool(old.source_nominal_enriched),C4_coupling_n=r.n,C4_coupling_positive=r.positive,C4_coupling_rho=r.median,C4_coupling_P=r.P,C4_coupling_FDR=r.FDR,within_C4_nominal_positive=positive,source_side_candidate=bool(positive and old.source_nominal_enriched),database_status=db.loc[r.gene,'status'])
 for lineage in ['CD8_T','Fibroblast']:
  rr=rec[(rec.ligand==r.gene)&(rec.receiver==lineage)&rec.receptor_eligible&rec.same_library_C4_ge20&(rec.C4_ligand_detected_cells>0)]
  row[lineage+'_compatible_patient_n']=rr.patient.nunique() if db.loc[r.gene,'database_interactions']>0 else np.nan
 for ep in ['CD8_cytotoxicity','CD8_dysfunction']:
  rr=fn[(fn.gene==r.gene)&(fn.endpoint==ep)&(fn.model=='study_adjusted')].iloc[0];row[ep+'_rho']=rr.rho;row[ep+'_P']=rr.P;row[ep+'_FDR']=rr.FDR
 rows.append(row)
priority=pd.DataFrame(rows);priority.to_csv(P/'tables/CANDIDATE_EVIDENCE_REGISTER.tsv',sep='\t',index=False)
nominated=priority[priority.source_side_candidate];assert set(nominated.gene)=={'APOE','CTSB','CTSD','LGALS3','SPP1'} and len(nominated)==6
leave=pd.read_csv(P/'tables/NOMINAL_C4_LEAVE_ONE_STUDY.tsv',sep='\t')
for i,row in leave.iterrows():
 n=int(s[(s.group=='C4')&(s.model=='depth_complexity_adjusted')&(s.program==row.program)&(s.gene==row.gene)].iloc[0]['n']);leave.at[i,'omitted_effective_patients']=n-int(row.n)
leave.to_csv(P/'tables/LEAVE_STUDY_EFFECTIVE_SCOPE.tsv',sep='\t',index=False)
def fmt(x):
 if pd.isna(x):return '—'
 return x if isinstance(x,str) else f'{x:.6g}'
def table(d,cols,labels):
 out=['| '+' | '.join(labels)+' |','| '+' | '.join(['---']*len(cols))+' |']
 for _,r in d.iterrows():out.append('| '+' | '.join(fmt(r[c]) for c in cols)+' |')
 return '\n'.join(out)
main=s[(s.group=='C4')&(s.model=='depth_complexity_adjusted')&(s.program=='lysosome_disjoint')&s.gene.isin(['APOE','CTSB','CTSD','LGALS3','SPP1'])]
report='''# 26210：代谢状态—候选输出—接收细胞深化结果（Phase278）

2026-09-28。继续围绕“代谢相关巨噬状态参与并调节免疫屏障”深化。本轮从C4内部的代谢变化寻找候选输出，再核查同文库CD8功能和接收端受体资格；没有重复扩大空间阈值搜索，没有重分C4，没有绘图或改正式稿。

## 新结果与研究重点

**本轮形成了一组有患者一致性的来源端候选：C4溶酶体相关程序与LGALS3、CTSB、CTSD、APOE正向关联，SPP1为名义显著对应。** 它们在原C4对其他巨噬的配对比较中亦富集，因此比单纯重复边界内外分布更具体地定位了代谢状态相关输出。

优先机制方向可收敛为：C4的溶酶体相关状态与LGALS3及基质处理候选分子的联系，SPP1作为协同线索保留。这里“输出”指候选分子的RNA表达，尚非已测定的胞外分泌/蛋白活性。来源端关联成立，并不意味着已经证明这些分子调节T细胞进入。

**接收端尚未形成完整支持。** 同文库11患者的CD8细胞毒/功能障碍读出未支持核心候选的稳定免疫抑制关联；独立肿瘤成纤维细胞仅1位患者39个，不能建立CAF响应的人群推断。这些是当前要补的环节，不能反向抹去来源端已确立的关联。

## 预设候选与非循环评分

预先固定15个输出：SPP1、APOE、LGALS3、TGFB1、VEGFA、IL10、IL1B、CCL2、CSF1、PDGFB、MMP9、MMP14、CTSB、CTSD、CD274。保留全部方向和资格失败，不按结果改变列表。

两个代谢暴露为Phase277冻结的去共享溶酶体程序、多胺/MTA转录上下文；再剔除全部15个输出。溶酶体最终123基因（进一步移除CTSB、CTSD），多胺仍为MAT2A/MAT2B/AMD1/SRM/SMS。**LGALS3、APOE、SPP1及CTSB/CTSD都不参与溶酶体评分，避免基因自重叠造成机械相关。**

使用原22444个成熟巨噬、固定2674个C4成员，肿瘤部分覆盖36患者。每患者在全体肿瘤巨噬内标准化log1pCP10k基因表达并平均得到程序分数，再在C4内部或整体巨噬中计算关联：
- C4至少50细胞，整体巨噬至少100细胞；
- 输出基因至少10阳性细胞，代谢程序至少20细胞检出两个以上程序基因；
- 主要调整为log深度和检测基因数；整体巨噬再控制原亚群标签；
- 每患者得到一个相关系数，细胞不作为独立统计样本；
- n≥10以精确双侧符号检验为主、Wilcoxon伴随；每个群体×代谢程序家族包含15输出×2模型=30项BH；n<10仅描述。

完整120项输出关联、4320条患者资格/效应记录均交付。样本不合格或缺少有效重复的条目不补造P值。

## 核心：C4内部的溶酶体相关输出

'''+table(main,['gene','n','positive','negative','median','P','FDR'],['输出候选','患者n','正向','负向','中位偏rho','原始P','FDR'])+'''

LGALS3为15/16位正向，中位rho=0.21755，P=0.0005188、FDR=0.002223；CTSB/CTSD均16/16位正向，中位rho=0.3865/0.4104，P均3.052e−5、FDR均0.0002289；APOE为14/16位正向，rho=0.1562，P=0.004181/FDR=0.01568。

SPP1为12/15位正向，rho=0.08214，P=0.03516、FDR=0.07233，按用户原P筛选标准保留为名义显著候选，不能称FDR显著。效应弱于前述LGALS3/蛋白酶关联；不因它是客户关注分子而夸大。

这些有效C4关联主要来自3项研究；不是全部5研究都有足够C4细胞。剔除最大研究GSE207493后，LGALS3仍6/7位正向、CTSB/CTSD为7/7正向，SPP1为5/6正向；这些留一结果只作方向稳健性，不以剩余小样本重新宣称显著。LEAVE_STUDY_EFFECTIVE_SCOPE.tsv明确标出实际移除人数，避免把本来没有合格患者的研究当独立重复验证。

整体巨噬控制亚群标签后，溶酶体与LGALS3为30/31正向（rho0.1712）、CTSB30/31（0.2825）、CTSD31/31（0.3775）、APOE31/31（0.3130），均通过FDR。这是同一资料内的组成敏感性核验，不能称新增独立队列。整体巨噬SPP1 rho仅0.0080、P0.06143，未通过主门槛；不能用不同P值直接证明C4关联显著强于其他亚群。

## 原C4富集证据与候选分级

沿用Phase200肿瘤C4对其他巨噬的原全基因P/FDR，不把15个候选另缩成更有利的差异表达校正池。LGALS3、CTSB、CTSD、APOE、SPP1的原配对中位log2CPM差依次约1.985、1.106、1.861、1.285、1.447，均具有原家族统计支持。

按事前规则“原C4名义富集＋C4内部调整后的正向代谢关联原P<0.05”，共有5个不同基因、6个程序—输出组合通过：溶酶体对应上述5个基因，另有多胺程序—CTSD（12/15正向，rho0.04971，P0.03516/FDR0.23145）。这只是来源端候选标准，不是完整跨细胞机制标准。

多胺—IL10在C4中10/12位正向，rho0.05349，P0.03857/FDR0.23145，也是需要保留的名义信号；但原C4的IL10水平相对其他巨噬降低，因此它不进入“C4富集输出”清单，仍作为状态内部变化线索单列。不能把来源富集门槛当作所有生物学作用的必要条件。多胺—SPP1未获支持，转录上下文也不等于MTA浓度。

MMP14与溶酶体的主符号P=0.07681，伴随Wilcoxon P=0.001678/FDR=0.004578；完整保留检验方法依赖性，按预设主门槛不升级为核心候选。

## 关键反向：不能将C4写成普遍促纤维化输出增强

C4相对其他巨噬的TGFB1、VEGFA、PDGFB、IL1B和IL10均降低。C4内部溶酶体程序与IL1B、PDGFB、VEGFA还出现显著负向关联；TGFB1/IL10的溶酶体关联未得到支持。因此不能直接套用文献中“所有TGFβ/VEGF促纤维化输出增强”的机制。

本轮优先指向溶酶体相关的LGALS3与蛋白酶候选。基质处理可能包含降解与重塑，不能直接等同胶原沉积增加或屏障增厚；需由实际胞外活性、组织结构及迁移读出确定方向。

## 同文库CD8功能：分别检验代谢程序与输出分子

采用Phase213原冻结的CD8细胞毒/功能障碍程序和C4资格，精确按同患者、同文库匹配，C4至少20细胞；同患者多文库先平均。有效11患者，其中9位来自GSE207493，另两项研究各1位，不能称均衡多队列验证。

代谢程序：2暴露×2功能×2模型共8项独立家族；研究内置换99999次、5000次患者bootstrap。研究调整后溶酶体—CD8细胞毒rho0.5325、P0.1469，溶酶体—功能障碍rho0.2955、P0.4457；多胺对应rho0.3923/−0.0622，P0.3029/0.8884。均无原P<0.05支持，也不能把细胞毒正相关写成抑制。

为补足“候选输出—接收功能”的直接环节，另在运行前冻结扩展方案，保留全部15个原输出，不只挑来源端阳性：15×2功能×2模型=60项同一BH家族。原8项不合并重算、不被替换。核心5候选的研究调整后结果如下：

'''+table(fn[fn.gene.isin(['LGALS3','CTSB','CTSD','APOE','SPP1'])&(fn.model=='study_adjusted')],['gene','endpoint','n','rho','P','FDR'],['C4输出RNA','CD8功能','患者n','研究调整rho','P','FDR'])+'''

扩展中只有CSF1—CD8功能障碍（研究调整rho−0.8088，P0.01098/FDR0.6588）及未调整VEGFA—功能障碍（rho−0.6818，P0.02504/FDR0.7512）达到名义显著，方向是较低功能障碍，必须保留，不能作为免疫抑制支持。核心来源端候选的接收功能关联尚未建立；这也没有直接检验空间迁移，因此不能据此否定或确认T细胞进入机制。

## 接收端受体资格

从锁定的CellChatDB导出原15输出的全部27条已收录相互作用，不按当前表达选择受体。原始计数提取29177个肿瘤CD8细胞（36患者、50文库）和39个肿瘤成纤维细胞（1患者、1文库）。此前参考总体有2位成纤维患者，限制肿瘤后只剩1位；不能把总体n2当肿瘤CAF独立n2，也不能把普通成纤维标签直接当POSTN⁺CAF。

每文库至少20个接收细胞，完整受体所有亚基须在同一细胞中共同检出，至少5细胞且≥10%。仅有单个整合素亚基不算合格。

SPP1-CD44在30位CD8患者满足接收资格，其中20位同时有同文库≥20个C4且检出SPP1；SPP1-ITGA4/ITGB1对应16位接收合格、15位同文库C4匹配。其余SPP1复合体未满足本轮资格。该结果是RNA相容性，不能证明配体结合或受体激活，也不能替代功能阴性。

APOE-TREM2/TYROBP在所检查CD8/成纤维细胞中均不合格。LGALS3、CTSB、CTSD等在当前锁定库中无条目：这表示数据库未覆盖，不能解释成它们没有生物学作用。完整1377条文库—受体资格记录随附。

## 为什么优先LGALS3及基质处理方向

已有直接研究显示，胞外Galectin-3可结合基质和糖基化IFNγ，限制IFNγ扩散及趋化因子梯度；拮抗Galectin后T细胞募集改善。该机制提供“基质中的信号可及性”这一可检验方向，不要求仅靠CAF数量或完整物理围墙解释所有排斥。[Gordon-Alonso等，Nature Communications2017](https://www.nature.com/articles/s41467-017-00925-6)

该论文主要采用肿瘤来源Galectin-3等实验体系，并未证明本项目C4向胞外分泌Galectin-3。本项目新增的是去自重叠、患者一致的C4内转录关联；胞外定位、来源与作用仍待验证。

CTSB/CTSD的巨噬与纤维化联系有实验文献依据，但其胞内溶酶体功能和胞外作用需要分开，不把蛋白酶RNA升高直接当成基质屏障增强。[Jiang等，Experimental & Molecular Medicine2023](https://www.nature.com/articles/s12276-023-01027-8)

因此后续问题可以具体化为：代谢条件改变是否影响C4相关巨噬的胞外LGALS3及基质处理活性，并改变基质中的趋化信号分布和T细胞跨界面迁移。优先测量实际蛋白定位/分泌与功能，维持SPP1作为协同候选；本轮没有执行或宣称完成这些实验。

## 验收及与前期结果的关系

本轮不改变Phase275/276/277的空间或协同结果，也未把来源端阳性拼成已经确认的代谢—屏障因果链。保留原C4成员、原差异表达全基因FDR、所有候选、原始P和本轮校正家族。

来源与受体计算均exit0。2400项独立偏秩相关最大误差1.96e−15，1377条受体完整复合体资格独立复核通过；8项代谢功能与60项输出功能的置换P、BH及文库/患者聚合复核通过。输出扩展SPP1伪bulk还与此前冻结文库数值逐一一致，60项独立相关误差<3.4e−16。原输入SHA/元数据未变。

完整报告、候选分级、120项来源端统计、8＋60项接收功能统计、资格、程序、代码和验收文件随附。细胞级中间量和置换分布保留服务器。新图0、正式稿修改0。
'''
(D/'代谢状态候选输出与接收端深化报告.md').write_text(report,encoding='utf-8')
# Full all-candidate register, not only shortlisted hits.
for srcdir,label in [(P,'source_and_receiver'),(E,'output_function_extension')]:
 target=D/label
 for sub in ['tables','code','qc']:(target/sub).mkdir(parents=True,exist_ok=True)
 for f in (srcdir/'tables').glob('*.tsv'):shutil.copy2(f,target/'tables'/f.name)
 for f in (srcdir/'code').iterdir():
  if f.is_file() and f.suffix in ['.py','.R'] and f.name!='assemble_delivery.py':shutil.copy2(f,target/'code'/f.name)
 for f in (srcdir/'qc').glob('*.json'):shutil.copy2(f,target/'qc'/f.name)
for f in ['PROTOCOL_BEFORE_TESTING.md','PROGRAMS_FROZEN.json','FREEZE.json']:shutil.copy2(P/f,D/f)
shutil.copy2(E/'ADDENDUM_BEFORE_OUTPUT_FUNCTION.json',D/'output_function_extension/ADDENDUM_BEFORE_OUTPUT_FUNCTION.json')
v=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());ev=json.loads((E/'qc/FINAL_VALIDATION.json').read_text());assert v['status']==ev['status']=='PASS'
combined=dict(status='PASS',at=datetime.datetime.now().astimezone().isoformat(),source_receiver=v,output_function=ev,source_side_genes=sorted(nominated.gene.unique()),candidate_program_gene_pairs=len(nominated),figures=0,manuscript_changes=0)
(P/'qc/COMBINED_VALIDATION.json').write_text(json.dumps(combined,indent=2));(D/'COMBINED_VALIDATION.json').write_text(json.dumps(combined,indent=2))
files=sorted(f for f in D.rglob('*') if f.is_file() and f.name!='DELIVERY_SHA256.tsv')
pd.DataFrame([dict(file=f.relative_to(D).as_posix(),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]).to_csv(D/'DELIVERY_SHA256.tsv',sep='\t',index=False)
with zipfile.ZipFile(P/'phase278_delivery.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(D.rglob('*')):
  if f.is_file():z.write(f,arcname=f.relative_to(D).as_posix())
(P/'STATUS.json').write_text(json.dumps(dict(status='VALIDATED_AWAITING_DELIVERY',at=datetime.datetime.now().astimezone().isoformat(),files=len(files),figures=0,background_jobs=0),indent=2))
print('DELIVERY_FILES',len(files),'ZIP_BYTES',(P/'phase278_delivery.zip').stat().st_size)
