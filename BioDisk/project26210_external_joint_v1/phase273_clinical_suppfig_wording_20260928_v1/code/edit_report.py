from pathlib import Path
import json,copy,re,hashlib,shutil
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
base=r/'phase272_manuscript_clinical_integration_20260928_v1';p=r/'phase273_clinical_suppfig_wording_20260928_v1'
b=json.loads((base/'report/report_materialized.json').read_text());before=copy.deepcopy(b);changes=[]
for i,x in enumerate(b):
 for key in ['text','caption']:
  if key in x:
   old=x[key];new=old.replace('公开临床队列的回顾性再分析显示','临床队列分析显示')
   new=new.replace('临床分析基于公开队列的回顾性再分析。','临床分析使用已发表公开队列的回顾性资料。')
   new=new.replace('临床再分析','临床队列分析').replace('回顾性再分析','回顾性分析').replace('再分析','分析')
   if key=='caption' and new.startswith('附表S19 '):new='附表S19 SPP1相关转录特征的临床关联'
   if new!=old:changes.append({'index':i,'field':key,'before':old,'after':new});x[key]=new
 for key in ['text']:
  if key not in x:continue
  t=x[key]
  if t.startswith('临床队列分析显示，文献SPP1巨噬评分'):t=t.replace('（附表S19）','（补充图S2A；附表S19）')
  if t.startswith('在CM025随机两臂的248名'):t=t.replace('探索性证据。','探索性证据（补充图S2B）。')
  if t.startswith('在nivolumab组，C4与PFS'):t=t.replace('（附表S20）','（补充图S2C-D；附表S20）')
  x[key]=t
caption=[
 {'type':'caption','text':'补充图 S2｜SPP1相关转录特征与C4状态的临床关联'},
 {'type':'p','text':'（A）CheckMate nivolumab组的总生存（OS）Kaplan-Meier曲线。SPP1相关巨噬评分为SPP1、FN1、CLEC5A、OLR1及IL1RN的表达均值，按完整311名患者队列的中位数分组。图示181名nivolumab患者，其中高组84名、低组97名，共123个死亡事件。横轴为治疗后月数，纵轴为总生存概率；阶梯线为估计值，浅色带为95%置信区间，“+”标记删失，底部列各时点在险人数。高评分组OS较短，log-rank P=0.01484、FDR=0.05934，为名义显著；连续评分模型的效应见正文及附表S19。'},
 {'type':'p','text':'（B）连续SPP1相关巨噬评分与治疗的交互效应。点为交互HR比值，横线为95%置信区间，虚线为1；横轴使用对数刻度。CM025比较nivolumab与everolimus，调整年龄、性别及MSKCC；JAVELIN比较avelumab联合axitinib与sunitinib，调整年龄和性别。CM025 OS交互HR比值1.595（1.170-2.173），P=0.003105、FDR=0.003105；CM025 PFS与JAVELIN PFS的交互P分别为0.2997和0.3319。比值大于1表示评分升高时，免疫治疗臂相对对照臂的风险比增加。所列FDR保留各队列及终点内SPP1单基因和五基因两个指标的校正家族。'},
 {'type':'p','text':'（C）C4状态评分与PFS的时间依赖性关联。采用不含TYMP及上述五基因的固定19基因评分，CheckMate nivolumab完整模型纳入164名患者、145个事件，控制年龄、统一性别、试验、MSKCC、文献五基因评分及总巨噬估计。横轴为治疗后月数，纵轴为基线C4评分每升高1 SD对应的条件风险比，采用对数刻度；实线和浅色带分别为HR及点态95%置信区间，虚线为HR=1。固定时间函数为β(t)=β0+β1log(1+t)，未搜索时间截点。两自由度增量检验P=0.0007173、FDR=0.001435；6及12个月HR分别为0.785（0.642-0.958）和0.621（0.473-0.816）。底部列在险人数；12个月后的灰色背景提示晚期风险集减少，仅为显示提示，不是分组或建模截点。'},
 {'type':'p','text':'（D）C4时间模型的稳健性。显示完整调整、按试验分层基线风险及仅CM025患者三个模型在固定12个月的HR与95%置信区间；分数缩放保持相同，横轴为对数刻度。右列P/FDR检验各模型完整的两自由度时间效应，不是仅对12个月HR的检验。三个模型的联合P分别为0.0007173、0.0008573和0.002245。后两项是同一队列的内部敏感性分析。'},
 {'type':'note','text':'临床资料来自已发表公开队列。C4的时间系数描述组织评分与条件风险的关系，不代表观察到C4功能随时间逆转；晚期HR取决于仍在险患者。当前未建立C4独立OS关联或相对治疗获益预测。完整模型、检验家族及来源见附表S19、S20和配套统计文件。'}
]
b.append({'type':'supplementary_figure','figure':'S2','source':str(p/'figures/Supplementary_Figure_S2.pdf')});b.extend(caption)
assert not any('再分析' in x.get('text','') or '再分析' in x.get('caption','') for x in b)
# Numeric table cells remain exactly identical.
ot=[x for x in before if x['type']=='table'];nt=[x for x in b if x['type']=='table']
assert len(ot)==len(nt)==34
for x,y in zip(ot,nt):assert x['headers']==y['headers'] and x['rows']==y['rows']
(p/'report/report_materialized.json').write_text(json.dumps(b,ensure_ascii=False,indent=2))
(p/'qc/WORDING_CHANGES.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2))
(p/'report/Supplementary_Figure_S2_caption_CN.md').write_text('\n\n'.join(x['text'] for x in caption),encoding='utf-8')
for f in (base/'figures').glob('*.pdf'):shutil.copy2(f,p/'figures'/f.name)
shutil.copytree(base/'delivery/clinical_statistics',p/'delivery/clinical_statistics',dirs_exist_ok=True)
s=(base/'code/render_report.py').read_text().replace(base.name,p.name).replace('26210_manuscript_with_clinical_section_v1','26210_manuscript_clinical_suppfig_v1')
s=s.replace("figpaths['S1']=d/'figures/Supplementary_Figure_S1.pdf'","figpaths['S1']=d/'figures/Supplementary_Figure_S1.pdf';figpaths['S2']=d/'figures/Supplementary_Figure_S2.pdf'")
s=s.replace("'S1':set('ABC')}","'S1':set('ABC'),'S2':set('ABCD')}")
s=s.replace(r'补充图\s*S1',r'补充图\s*S[12]').replace(r'补充图\s*(S1)',r'补充图\s*(S[12])')
s=s.replace("n='S1' if m[1] else int(m[2])","n=m[1] if m[1] else int(m[2])")
s=s.replace("'补充图 S1' if n=='S1' else f'图{n}'","f'补充图 {n}' if isinstance(n,str) else f'图{n}'")
s=s.replace("level=0 if n=='S1' else 1","level=0 if isinstance(n,str) else 1")
s=s.replace("if n=='S1':\n   anchors=","if n=='S2':\n   anchors=[('A',10,29),('B',110,29),('C',10,145),('D',110,145)]\n  elif n=='S1':\n   anchors=")
s=s.replace('figures=6','figures=7').replace('SPP1与C4临床小节整合 含补充图S1','SPP1与C4临床关联 含补充图S1-S2')
(p/'code/render_report.py').write_text(s)
v=(base/'code/validate_render.py').read_text().replace('26210_manuscript_with_clinical_section_v1','26210_manuscript_clinical_suppfig_v1').replace('len(b)-34-6-1','len(b)-34-7-1')
(p/'code/validate_render.py').write_text(v)
inputs=[base/'report/report_materialized.json',base/'report/26210_manuscript_with_clinical_section_v1.pdf',*list((base/'figures').glob('*.pdf'))]
for ph,name in [
 ('phase270_SPP1_literature_clinical_revalidation_20260927_v1','CHECKMATE_ANALYSIS_INPUT.tsv'),
 ('phase270_SPP1_literature_clinical_revalidation_20260927_v1','ALL_SURVIVAL_MODELS.tsv'),
 ('phase270_SPP1_literature_clinical_revalidation_20260927_v1','TREATMENT_INTERACTIONS.tsv'),
 ('phase271_C4_state_clinical_extension_20260928_v1','C4_MODEL_MEMBERSHIP.tsv'),
 ('phase271_C4_state_clinical_extension_20260928_v1','C4_PFS_TIME_VARYING_MODELS.tsv'),
 ('phase271_C4_state_clinical_extension_20260928_v1','C4_PFS_TIME_SPECIFIC_HR.tsv'),
 ('phase271_C4_state_clinical_extension_20260928_v1','C4_TIME_VARIATION_TRIAL_HR.tsv'),
 ('phase271_C4_state_clinical_extension_20260928_v1','C4_TIME_VARIATION_TRIAL_SENSITIVITY.tsv')]:inputs.append(r/ph/'tables'/name)
(p/'qc/INPUTS.json').write_text(json.dumps([{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in inputs],indent=2))
print('EDIT_COMPLETE',len(changes),'wording replacements',len(b),'blocks')
