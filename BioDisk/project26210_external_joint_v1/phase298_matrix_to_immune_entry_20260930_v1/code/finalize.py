from pathlib import Path
import pandas as pd,numpy as np,json,csv,io,hashlib,datetime,os,shutil,subprocess
from scipy import stats
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
with (P/'logs/memory_gate.log').open('a') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,check=True)
s=pd.read_csv(P/'tables/SECTION_RESULTS.tsv',sep='\t');p=pd.read_csv(P/'tables/PATIENT_RESULTS.tsv',sep='\t');c=pd.read_csv(P/'tables/COHORT_RESULTS.tsv',sep='\t')
checks=0
for r in c.itertuples():
 z=p[(p.dataset==r.dataset)&(p.radius==r.radius)&(p.model==r.model)&(p.program==r.program)&(p.endpoint==r.endpoint)].effect.dropna().values
 assert len(z)==r.n and int((z<0).sum())==r.negative and int((z>0).sum())==r.positive
 assert np.isclose(np.median(z),r.median_rho,equal_nan=True)
 if len(z)>=10:assert np.isclose(stats.wilcoxon(z,alternative='two-sided',method='auto').pvalue,r.P)
 else:assert pd.isna(r.P) and pd.isna(r.FDR)
 checks+=1
for _,g in c.groupby(['dataset','radius','model']):
 q=stats.false_discovery_control(g.P.fillna(1));assert np.allclose(q[g.P.notna()],g.FDR[g.P.notna()]);checks+=1
for _,g in s.groupby(['sample','radius','model']):
 q=stats.false_discovery_control(g.spatial_P.fillna(1));assert np.allclose(q[g.spatial_P.notna()],g.spatial_FDR[g.spatial_P.notna()]);checks+=1
assert (P/'logs/run.exit').read_text().strip()=='0'
report="""# Phase298：ECM组织程序与相对肿瘤内T细胞减少的方向线索
按用户最新取舍：前半段以已有前两部分结果描述，本轮只补“基质重塑—免疫进入受限”；没有重新要求肿瘤TYMP—SPP1患者串联显著。

## TJ-RCC提供方向一致的探索性空间线索
固定200μm主尺度，控制外侧T细胞供给、内外测序深度、恶性/巨噬比例、区域数量比及空间坐标/连通分量后，外侧ECM组织程序与T细胞内外log信号差为负向：9患者7负2正，中位rho−0.192010。
进一步控制内外血管/壁细胞比例，9/9患者均为负向，中位rho−0.144059。全部合格区域外侧T q05估计>0，避免将完全没有外侧供给误称阻挡。此为映射估计非实测T细胞计数，不能据>0赋予病理“免疫热”标签。
400μm core与vascular各9例7负，中位rho−0.161769/−0.174626，方向仍负但不是全部患者同向。n9按既定规则不作患者推断检验。
这条线索适合补充“基质组织状态与T细胞相对肿瘤内不足相联系”的描述；不能将控制外侧信号等同排除了所有供给和组成混杂。

## 完整对照显示证据尚未稳定跨队列
TJ-RCC 200μm血管调整的9个切片严格空间零模型P均≥0.05（最小0.057），因此未显示超出该空间背景模型的单片名义支持；不能只凭9/9方向称显著。
TJ-RCC CD8仅4患者，200μm ECM组织—CD8为3/4负向；不推断。ECM_remodeling程序与总T主要偏正，与ECM_organization不可混称一个普遍排斥程序。
GSE175540主200μm模型：ECM组织—总T n15、rho−0.017816、P0.890381/FDR0.890381；—CD8 n13、rho0.107930、P0.414307/FDR0.851562，均无支持。
该队列400μm core ECM组织—CD8出现反向于排斥的名义信号：n12、9正3负、rho0.169469、P0.042480/FDR0.169922；血管调整后不保留（P0.733398）。完整反向证据保留，不用TJ-RCC方向替代跨队列结论。
GLMF两病例200μm ECM组织—CD8均负，而总T均正；ECM重塑—总T均负。作为案例描述，不能称全免疫一致排斥。

## 文献可补机制依据，不能充当本项目验证
Salmon等在人肺肿瘤活切片中直接观察基质结构与T细胞运动，并通过胶原酶处理增加T细胞接近/接触癌细胞，为“基质结构限制进入”的解释提供实验依据。原始研究：https://www.jci.org/articles/view/45817 ，DOI10.1172/JCI45817。
本项目为ccRCC的ECM转录程序与静态空间分布对应，组织、技术和证据层级不同。没有直接纤维取向、孔隙、T细胞轨迹或本项目基质干预；不把他人实验转述成本项目因果结果。

## 可纳入主线的文字
“进一步的空间分析在TJ-RCC中观察到，外侧ECM组织程序与相对肿瘤内T细胞信号呈负向联系，且在控制外侧免疫供给和血管相关组成后保持一致方向。这一探索性结果与基质组织参与限制免疫细胞进入的机制相符，为前述SPP1相关髓系—间质重塑联系提供了下游空间线索。”
随后在讨论集中说明：9例描述性结果、单片严格零模型未通过、其他队列未一致重复以及RNA/静态分布的测量边界。不得删去这些会改变整体解释的证据。

## 方法与验收
完整复用Phase289计算边界和Phase291预设几何/协变量；外侧暴露为既定ECM_organization和ECM_remodeling。结局总T为q05丰度代理区域均值的内外log1p差，CD8为既定表达内外差。
主要200μm core，200μm vascular与400μm两个模型均为敏感性；没有结果后选尺度/程序。每队列×尺度×模型4个假设分别BH，原P<0.05为名义门槛。患者少于10只方向/效应，不用spot补n。
使用既有D/W/yr/V缓存逐edge精确对齐，新暴露重新投影，999次协变量约束谱翻转；保存全部528模型，340可评估。独立核查结局残差、rho、3个代理的协变量正交/范数/空间二次型不变量。患者合并、P/FDR及完整方向另验收。无新空间P筛患者。
属于既往分析后的定向探索，并非独立验证；新结果仅对最后一环提供有限补充。没有重算前两部分、修改旧对象、标签、边界或文章，没有新图。

输出：tables/COHORT_RESULTS.tsv、PATIENT_RESULTS.tsv、SECTION_RESULTS.tsv；qc/MODEL_VALIDATION.tsv、FINAL_VALIDATION.json；PROTOCOL_FROZEN.md；code/与logs/。
"""
reportpath=P/'delivery/基质重塑与免疫进入_补充证据.md';reportpath.write_text(report)
now=datetime.datetime.now().astimezone().isoformat()
entry=f"""## Phase298（2026-09-30）：用户收窄为仅补基质重塑—免疫进入
前半段沿用前两部分描述，不再强求肿瘤TYMP—SPP1患者串联。报告：{P.name}/delivery/基质重塑与免疫进入_补充证据.md。
固定现有边界直接检验外侧ECM到T/CD8内外差。TJ-RCC ECM组织—总T：200μm core9例7负rho−0.192010；vascular9/9负rho−0.144059；400μm均7/9负。n9只描述，严格单片P全部未过0.05。为可保留的方向线索，非功能屏障证实。
GSE175540主模型不支持，400μm core ECM组织—CD8反向名义正关联P0.042480/FDR0.169922完整保留。ECM重塑与组织程序不混同。528模型/340可评估及患者/P/FDR检查完成。
前半段不重算，TGFB1不另起主线；旧数据/C4/边界/文章不改，无新图。本轮完成，无后台。

"""
names=['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv'];old={n:(R/n).read_bytes() for n in names}
for n,b in old.items():
 with (P/'root_record_backups'/n).open('xb') as f:f.write(b)
new={}
for n,b in old.items():
 if n.endswith('.md'):new[n]=entry.encode()+b
 else:
  fields=b.decode('utf-8-sig').splitlines()[0].split('\t');o=io.StringIO();w=csv.DictWriter(o,fieldnames=fields,delimiter='\t',lineterminator='\n')
  for rid,typ,path in [('phase298_report','report',reportpath),('phase298_cohort','table',P/'tables/COHORT_RESULTS.tsv'),('phase298_section','table',P/'tables/SECTION_RESULTS.tsv')]:
   w.writerow(dict(result_id=rid,stage='Phase298',artifact_type=typ,absolute_path=str(path),generated_at=now,validation_status='VALIDATED_EXPLORATORY',validation_evidence=str(P/'qc/FINAL_VALIDATION.json'),notes='SHA256='+sha(path)))
  new[n]=b+(b'' if b.endswith(b'\n') else b'\n')+o.getvalue().encode()
for n in names:assert (R/n).read_bytes()==old[n]
for n,b in new.items():
 tmp=R/(n+'.phase298.tmp')
 with tmp.open('xb') as f:f.write(b)
 shutil.copymode(R/n,tmp);os.replace(tmp,R/n)
assert all((R/n).read_bytes()==new[n] for n in names)
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(dict(status='PASS',at=now,aggregate_checks=checks,spatial_models_checked=340,section_models=528,new_figures=0,report_sha256=sha(reportpath),root_records_updated=names),indent=2))
print('PHASE298_COMPLETE',checks,'aggregate checks')
