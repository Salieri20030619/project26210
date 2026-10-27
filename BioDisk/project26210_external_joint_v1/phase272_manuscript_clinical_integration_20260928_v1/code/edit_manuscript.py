from pathlib import Path
import json,re,hashlib,shutil,copy
import pandas as pd
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
base=r/'phase267_manuscript_abstract_logic_20260927_v1';p=r/'phase272_manuscript_clinical_integration_20260928_v1'
for d in ['code','report','figures','qc','logs','delivery','delivery/clinical_statistics']:(p/d).mkdir(parents=True,exist_ok=True)
z=json.loads((base/'report/report_materialized.json').read_text());b=copy.deepcopy(z)
c=r/'phase271_C4_state_clinical_extension_20260928_v1';a=r/'phase270_SPP1_literature_clinical_revalidation_20260927_v1'
def rd(root,name):return pd.read_csv(root/'tables'/name,sep='\t')
cm=rd(c,'C4_SURVIVAL_MODELS.tsv');tv=rd(c,'C4_PFS_TIME_VARYING_MODELS.tsv');th=rd(c,'C4_PFS_TIME_SPECIFIC_HR.tsv');sens=rd(c,'C4_TIME_VARIATION_TRIAL_SENSITIVITY.tsv');lit=rd(a,'ALL_SURVIVAL_MODELS.tsv');lc=rd(c,'PHASE270_ADJUSTED_SEX_CORRECTION_AUDIT.tsv');li=rd(a,'TREATMENT_INTERACTIONS.tsv')
def row(t,**kw):
 for k,v in kw.items():t=t[t[k]==v]
 assert len(t)==1
 return t.iloc[0]
def n(x):return f'{x:.4g}'
def hr(x,name='HR'):return f'{x[name]:.3f} ({x.CI_low:.3f}-{x.CI_high:.3f})'
q=row(tv,model_id=12);q6=row(th,model_id=12,month=6);q12=row(th,model_id=12,month=12)
ls=row(lit,dataset='CheckMate',endpoint='OS',arm='NIVOLUMAB',score='Jiang5',model='continuous')
la=row(lc,endpoint='OS',arm='NIVOLUMAB',score='Jiang5')
lm=row(lit,dataset='CheckMate',endpoint='OS',arm='NIVOLUMAB',score='Jiang5',model='median_cohort')
ix=row(li,dataset='CheckMate',endpoint='OS',score='Jiang5')
sc=row(sens,parent_model_id=12,sensitivity='CM025_only');st=row(sens,parent_model_id=12,sensitivity='stratify_trial')
su=row(cm,model_id=24)
def P(t):return {'type':'p','text':t,'_new_clinical':True}
def H(t):return {'type':'h2','text':t,'_new_clinical':True}
def N(t):return {'type':'note','text':t,'_new_clinical':True}
methods=[H('1.11 临床队列与C4状态转录评分'),
P('临床分析纳入CheckMate公开RNA与随访资料的311名患者，其中nivolumab组181名、everolimus组130名；JAVELIN Renal 101纳入具有RNA数据的726名患者，其中avelumab联合axitinib组354名、sunitinib组372名[21,22]。按患者身份连接表达与临床记录，性别统一为女性或男性。CheckMate的PFS_CNSR和OS_CNSR按原始事件编码使用，JAVELIN PFS事件为1-PFS_P_CNSR。'),
P('文献SPP1巨噬评分由SPP1、FN1、CLEC5A、OLR1及IL1RN组成[19]。CheckMate使用已标准化发布表达的算术均值；JAVELIN使用log2(TPM+1)的均值。文献分组比较采用完整队列中位数，高组为评分不低于中位数；连续评分按治疗臂内1 SD缩放。SPP1单基因比较按既有文献定义单独计算[20]。'),
P('C4状态评分基于肿瘤组织24名患者、5项研究的C4与其他成熟巨噬配对伪bulk；每组至少20个细胞。筛选要求配对Wilcoxon P<0.05、中位log2CPM差异≥0.5、至少70%患者同向、至少80%研究中位差异为正且C4中位检测比例≥0.1。排除SPP1、TYMP、文献五基因及线粒体/核糖体基因，按研究平衡效应及患者效应排序取20个。每个留一研究折只在训练研究选择基因。临床队列共同覆盖其中19个基因，未覆盖的CCL18不计分；两队列按相同成员及上述表达均值规则计算C4评分。'),
P('Cox模型分析OS与PFS。CheckMate调整年龄、性别、试验及MSKCC；JAVELIN调整年龄和性别，再分别加入SPP1、文献五基因及官方总巨噬估计。按各模型完整病例进行拟合，并在同一批患者中比较嵌套模型。用Schoenfeld残差检验比例风险；针对CheckMate PFS的非比例风险，登记并拟合β(t)=β0+β1log(1+t)，t以月计，采用两自由度似然比检验C4的增量关联。治疗交互在CM025随机两臂及JAVELIN两臂分别检验。按试验分层基线风险、CM025子集及逐患者遗漏进行敏感性分析。')]
results=[H('2.13 SPP1相关转录特征与较短OS相关，C4状态显示时间依赖性PFS关联'),
P(f'公开临床队列的回顾性再分析显示，文献SPP1巨噬评分与nivolumab治疗后的较短OS相关。在181名患者、123个死亡事件中，高评分组84名、低评分组97名，log-rank P={n(lm.logrank_P)}、FDR={n(lm.logrank_FDR)}，与原报告方向及P值相符[19]。连续评分每升高1 SD的HR为{hr(ls)}，P={n(ls.P)}、FDR={n(ls.FDR)}；调整年龄、性别、试验和MSKCC后，164名患者中HR为{hr(la)}，P={n(la.P)}、FDR={n(la.FDR)}，保留名义显著关联（附表S19）。'),
P(f'在CM025随机两臂的248名完整病例、189个死亡事件中，文献评分与治疗的OS交互HR比值为{hr(ix,"HR_ratio")}，P={n(ix.P)}、FDR={n(ix.FDR)}。这一方向提示，评分升高与nivolumab相对everolimus的OS优势减弱相关，为SPP1相关转录特征的治疗情境差异提供探索性证据。'),
P('本项目C4状态评分进一步提供了不同于文献指标的状态读出。留一研究内部验证中，24/24名患者的C4评分高于同患者其他成熟巨噬，五项研究方向一致；移除临床未覆盖的CCL18后仍保持该方向。19基因C4评分与SPP1单基因的相关性较弱，CheckMate和JAVELIN的rho分别为0.123和0.135；与文献五基因评分的rho分别为0.426和0.476，显示两类转录特征存在区分。'),
P(f'在nivolumab组，C4与PFS的关联随随访时间变化。调整临床因素、文献五基因及总巨噬估计后，164名患者、145个事件的时间模型提供增量关联信息，两自由度似然比为{q.joint_LR:.2f}，P={n(q.joint_P)}、FDR={n(q.joint_FDR)}。每升高1 SD评分的条件进展风险HR在6个月为{hr(q6)}，12个月为{hr(q12)}，对应在险人数分别为62和27。按试验分层基线风险时P={n(st.joint_P)}、FDR={n(st.joint_FDR)}；仅保留119名CM025 nivolumab患者时P={n(sc.joint_P)}、FDR={n(sc.joint_FDR)}，逐患者遗漏后时间系数方向保持一致（附表S20）。'),
P(f'不同治疗背景下，C4评分还呈现不同的条件关联。在JAVELIN sunitinib组，控制文献评分和总巨噬估计后，较高C4评分与较低PFS风险相关，371名患者中HR为{hr(su)}，P={n(su.P)}、FDR={n(su.FDR)}。上述结果将C4的脂质处理状态与临床病程联系起来，并显示其临床读出有别于SPP1总体转录特征。')]
neg=P('临床终点与治疗差异。文献五基因评分在CheckMate的OS关联得到复现，但其PFS并未同步显著。C4在nivolumab组的独立OS关联未获支持，完整调整后P=0.2660（以附表S20为准）；CM025的C4恒定OS交互P=0.7479，允许时间变化后的PFS治疗交互P=0.4180、FDR=0.5886，当前数据不支持C4预测相对免疫治疗获益。JAVELIN联合组的C4恒定PFS关联未获支持；sunitinib较低风险信号在加入文献评分后出现，属于条件关联。SPP1单基因在JAVELIN连续调整模型中为名义显著（HR=1.218，P=0.02745，FDR=0.1098），但中位数分组未精确复现原文HR=1.61、P=0.004[20]，具体阈值或随访实现的差异尚未确定。')
cn=row(cm,model_id=4);neg['text']=neg['text'].replace('P=0.2660（以附表S20为准）',f'HR={cn.HR:.3f}，P={n(cn.P)}、FDR={n(cn.FDR)}')
lim=P('临床分析基于公开队列的回顾性再分析。文献复现和不同模型检查不等同新增独立临床队列验证；C4评分的留一研究评估基于既有冻结标签。时间模型由比例风险诊断触发，属于探索性扩展。随时间变化的是基线评分的条件风险系数，不是直接观察到C4功能逆转；12及24个月在险人数分别为27及17，晚期HR受风险集构成影响。组织RNA评分不能直接等同C4比例，亦未建立预测准确性或临床决策效用。')
disc=P('临床再分析将SPP1相关巨噬特征与生存结局相联系，同时揭示状态内部的差异。文献五基因评分对应较短OS及CM025的OS治疗交互线索；由本项目患者配对差异确定的C4评分则呈现时间依赖性PFS关联，并在调整文献评分和总巨噬估计后保留增量信息。结合C4的脂质处理和溶酶体表型，这些结果支持从具体功能状态解释SPP1富集巨噬的临床异质性。')
stats=P('临床模型使用目标评分或评分×治疗项的P值，不以其他协变量P值筛选。文献比较按队列、终点、模型和阈值定义同类家族；C4恒定模型按队列×终点×调整方式的两治疗臂进行BH校正。时变模型按两治疗臂、时变治疗交互按两调整模型，试验敏感性按两种固定检查分别校正。时点HR和区间描述同一时间函数，完整原P、FDR、模型公式及有效患者数见附表和配套统计文件。')
out=[];mapping={}
for i,x in enumerate(b):
 if x['type']=='h2':
  m=re.match(r'([12])\.(\d+) ',x['text'])
  if m and ((m[1]=='1' and int(m[2])>=11) or (m[1]=='2' and int(m[2])>=13)):
   old=x['text'];x['text']=re.sub(r'^[12]\.\d+',f'{m[1]}.{int(m[2])+1}',old);mapping[old]=x['text']
 if i==48:out.extend(methods)
 if i==166:out.extend(results)
 if i==6:
  out.append(P(f'临床转录组再分析显示，文献SPP1巨噬评分与nivolumab治疗后的较短OS相关（181名患者，每1 SD HR={ls.HR:.3f}，P={n(ls.P)}，FDR={n(ls.FDR)}）。排除该文献基因后的C4状态评分在留一研究评估中24/24名患者同向；在164名nivolumab患者中，控制文献评分及总巨噬估计后，C4与PFS呈时间依赖性关联（增量P={n(q.joint_P)}，FDR={n(q.joint_FDR)}）。'))
 if i==3:x['text']+=' 临床转录组分析另纳入CheckMate的311名及JAVELIN Renal 101的726名患者。'
 if i==1:x['text']=x['text'].replace('2026-09-27','2026-09-28')
 if i==10:x['text']=x['text'].replace('巨噬细胞状态及其空间分布','巨噬细胞状态、临床结局及其空间分布');x['text']+=' 临床生存资料采用公开治疗队列进行回顾性再分析。'
 if i==214:x['text']='3.4 SPP1富集巨噬的脂质处理特征与临床异质性'
 if i==218:out.append(disc)
 if i==236:out.append(neg)
 if i==246:
  x['text']='空间区域特征与临床结局未在同一批患者中联合测定；临床分析检验的是独立来源的组织转录评分，不能将空间邻近关系直接解释为临床效应。'
  out.append(lim)
 if i==240:x['text']=x['text'].replace('本研究采用横断面多组学设计','本研究整合横断面多组学及回顾性临床资料')
 if i==251:x['text']=x['text'].replace('配对代谢组和实验方向将MTA保留为候选代谢联系。','配对代谢组和实验方向将MTA保留为候选代谢联系。临床再分析显示SPP1相关巨噬转录特征与生存结局相联系，C4状态提供了不同于文献评分的时间依赖性进展风险信息。')
 if i==75:out.append(stats)
 if i==345:
  out.append(H('4.15 临床转录评分与生存模型'))
  rows=[]
  for label,v in [('文献五基因 OS 中位数分组',lm),('文献五基因 OS 连续评分',ls),('文献五基因 OS 临床调整',la)]:
   rows.append([label,f'{int(v.n)}/{int(v.events)}',hr(v),n(v.P),n(v.FDR)])
  rows.append(['文献五基因 CM025 OS治疗交互','248/189',hr(ix,'HR_ratio'),n(ix.P),n(ix.FDR)])
  js=row(lit,dataset='JAVELIN',endpoint='PFS',arm='Avelumab+Axitinib',score='SPP1',model='continuous_adjusted')
  rows.append(['JAVELIN联合组 SPP1连续PFS','354/157',hr(js),n(js.P),n(js.FDR)])
  out.append({'type':'table','caption':'附表S19 文献SPP1转录特征的临床再分析','headers':['模型','n/事件','HR或HR比值 (95%CI)','P','FDR'],'rows':rows,'_new_clinical':True})
  out.append(N('前三行均为CheckMate nivolumab组OS；中位数行列Cox P，log-rank P=0.01484、FDR=0.05934。连续HR按每1 SD，调整后值使用统一性别标签。治疗交互HR比值表示两治疗臂评分斜率之比，不是治疗HR；文献交互FDR保留SPP1与五基因两个指标家族。所有结果为公开队列再分析。'))
  rows=[]
  for id,label in [(4,'Nivolumab OS 恒定HR'),(12,'Nivolumab PFS 恒定HR'),(20,'JAVELIN联合 PFS 恒定HR'),(24,'Sunitinib PFS 恒定HR')]:
   v=row(cm,model_id=id);rows.append([label,f'{int(v.n)}/{int(v.events)}',hr(v),n(v.P),n(v.FDR)])
  rows += [['Nivolumab PFS 时变C4增量','164/145',f'LR={q.joint_LR:.3f}; df=2',n(q.joint_P),n(q.joint_FDR)],['时变模型 按试验分层','164/145','联合2 df',n(st.joint_P),n(st.joint_FDR)],['时变模型 仅CM025 nivo','119/104','联合2 df',n(sc.joint_P),n(sc.joint_FDR)]]
  ct=rd(c,'C4_TREATMENT_INTERACTIONS.tsv');cti=rd(c,'C4_TIME_VARYING_TREATMENT_INTERACTIONS.tsv')
  v=row(ct,interaction_id=2);rows.append(['C4 CM025 OS治疗交互','248/189',hr(v,'HR_ratio'),n(v.P),n(v.FDR)])
  v=row(cti,interaction_id=4);rows.append(['C4 CM025 时变PFS交互','248/220',f'LR={v.joint_LR:.3f}; df=2',n(v.joint_interaction_P),n(v.FDR)])
  out.append({'type':'table','caption':'附表S20 C4状态评分的生存关联及敏感性','headers':['模型','n/事件','效应量或统计量','P','FDR'],'rows':rows,'_new_clinical':True})
  out.append(N('均为临床因素、文献五基因及官方总巨噬估计调整模型；交互模型同时加入文献评分×治疗项。恒定HR按每1 SD；时变行采用两自由度似然比检验，不能解读为单一HR。6/12个月HR分别0.785 (0.642-0.958)/0.621 (0.473-0.816)，在险人数62/27；时点区间为同一模型的描述性区间。'))
  genes=json.loads((c/'qc/CLINICAL_SCORE_FREEZE.json').read_text())['genes_used']
  out.append(N('临床C4评分19个成员：'+', '.join(genes)+'。发现端原20基因另含CCL18；两队列按共同覆盖统一评分。24项恒定模型、6项治疗交互、4项时变模型、2项时变交互及完整敏感性保留于配套clinical_statistics文件夹，包含各模型原始P、FDR和家族。'))
 out.append(x)
 if i==363:
  for text in [
  '[19] Jiang W et al. SPP1+ Tumor-Associated Macrophages Drive Immunotherapy Resistance via CD8+ T-cell Dysfunction in Clear-Cell Renal Cell Carcinoma. Cancer Immunology Research. 2025;13:1533-1546. https://doi.org/10.1158/2326-6066.CIR-24-1146',
  '[20] Xu X et al. SPP1 expression indicates outcome of immunotherapy plus tyrosine kinase inhibition in advanced renal cell carcinoma. Human Vaccines & Immunotherapeutics. 2024;20:2350101. https://doi.org/10.1080/21645515.2024.2350101',
  '[21] Braun DA et al. Interplay of somatic alterations and immune infiltration modulates response to PD-1 blockade in advanced clear cell renal cell carcinoma. Nature Medicine. 2020;26:909-918. https://doi.org/10.1038/s41591-020-0839-y',
  '[22] Motzer RJ et al. Avelumab plus axitinib versus sunitinib in advanced renal cell carcinoma: biomarker analysis of the phase 3 JAVELIN Renal 101 trial. Nature Medicine. 2020;26:1733-1741. https://doi.org/10.1038/s41591-020-1044-8']:out.append(P(text))
assert [t for t in out if t['type']=='table' and not t.get('_new_clinical')]==[t for t in z if t['type']=='table']
for section,last in [(1,19),(2,21)]:
 nums=[int(re.match(rf'{section}\.(\d+)',x['text'])[1]) for x in out if x['type']=='h2' and re.match(rf'{section}\.\d+',x['text'])]
 assert nums==list(range(1,last+1)),nums
(p/'report/report_materialized.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
(p/'qc/SECTION_NUMBER_MAP.json').write_text(json.dumps(mapping,ensure_ascii=False,indent=2))
for f in (base/'figures').glob('*.pdf'):shutil.copy2(f,p/'figures'/f.name)
script=(base/'code/render_report.py').read_text().replace(base.name,p.name).replace('26210_manuscript_abstract_logic_revised_v1','26210_manuscript_with_clinical_section_v1').replace('摘要补全 全文逻辑修订 阴性结果集中讨论 含补充图S1','SPP1与C4临床小节整合 含补充图S1')
script=script.replace("if x['caption'].startswith('附表S18 '):weights=[2.0,.45,.65,1.2,1.6]","if x['caption'].startswith('附表S18 '):weights=[2.0,.45,.65,1.2,1.6]\n  if x['caption'].startswith(('附表S19 ','附表S20 ')):weights=[2.7,.8,2.3,1.0,1.0]")
(p/'code/render_report.py').write_text(script)
v=(base/'code/validate_render.py').read_text().replace('26210_manuscript_abstract_logic_revised_v1','26210_manuscript_with_clinical_section_v1').replace('tables=32','tables=34').replace('len(b)-32-6-1','len(b)-34-6-1')
(p/'code/validate_render.py').write_text(v)
files=[base/'report/report_materialized.json',base/'report/26210_manuscript_abstract_logic_revised_v1.pdf',*list((base/'figures').glob('*.pdf')),a/'tables/ALL_SURVIVAL_MODELS.tsv',a/'tables/TREATMENT_INTERACTIONS.tsv']
for f in ['C4_SURVIVAL_MODELS.tsv','C4_TREATMENT_INTERACTIONS.tsv','C4_PFS_TIME_VARYING_MODELS.tsv','C4_PFS_TIME_SPECIFIC_HR.tsv','C4_TIME_VARYING_TREATMENT_INTERACTIONS.tsv','C4_TIME_VARIATION_TRIAL_SENSITIVITY.tsv','C4_TIME_VARIATION_TRIAL_HR.tsv','C4_TIME_VARIATION_LEAVE_ONE_PATIENT.tsv','PHASE270_ADJUSTED_SEX_CORRECTION_AUDIT.tsv','CLINICAL_GENE_COVERAGE.tsv','FROZEN_C4_SIGNATURE.tsv','C4_CONSTRUCT_CORRELATIONS.tsv']:
 files.append(c/'tables'/f);shutil.copy2(c/'tables'/f,p/'delivery/clinical_statistics'/f)
shutil.copy2(a/'tables/TREATMENT_INTERACTIONS.tsv',p/'delivery/clinical_statistics/LITERATURE_TREATMENT_INTERACTIONS.tsv')
lit[lit.model!='continuous_adjusted'].to_csv(p/'delivery/clinical_statistics/LITERATURE_UNADJUSTED_AND_MEDIAN_MODELS.tsv',sep='\t',index=False)
lit[(lit.dataset=='JAVELIN')&(lit.model=='continuous_adjusted')].to_csv(p/'delivery/clinical_statistics/JAVELIN_LITERATURE_ADJUSTED_MODELS.tsv',sep='\t',index=False)
(p/'qc/INPUTS.json').write_text(json.dumps([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files],indent=2))
(p/'qc/EDIT_VALIDATION.json').write_text(json.dumps({'old_tables_unchanged':32,'new_tables':2,'methods_sections':19,'results_sections':21,'new_result_position':'after prior2.12 and Figure4 captions, before former spatial2.13 (now2.14)','new_analyses':0},indent=2))
print('AUTHORING_DONE',len(out),'blocks')
