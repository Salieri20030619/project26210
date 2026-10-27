from pathlib import Path
import copy,json,re,shutil
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
old=r/'phase191_discovery_order_report_20260917_v1'
d=r/'phase201_three_candidates_C4_main_report_20260918_v4'
a=json.loads((old/'report/report_materialized.json').read_text());b=copy.deepcopy(a)
# Keep all original figures and tables; add one bounded candidate-comparison section.
for p in (old/'figures').glob('*.pdf'):shutil.copy2(p,d/'figures'/p.name)
i=next(i for i,x in enumerate(b) if x.get('type')=='h2' and x.get('text','').startswith('2.9 '))
for x in b:
 if x.get('type')=='h2':
  for n in [9,10,11,12]:
   if x.get('text','').startswith(f'2.{n} '):
    x['text']=x['text'].replace(f'2.{n} ',f'2.{n+1} ',1);break
new=[
 {'type':'h2','text':'2.9 三个SPP1高表达候选群共同拓宽信号来源覆盖'},
 {'type':'p','text':'在12个成熟巨噬候选群的并列鉴定中，C4、C10和C15在患者等权表达汇总中呈现较高的SPP1表达。三群同时具有不同的标志组合：C4以GPNMB、CCL18和CYP27A1为代表，C10呈USP53/ALOX5AP特征，C15呈MKI67/TOP2A增殖特征（图8A）。因此，SPP1分布于多个巨噬状态，三群均保留在候选谱系和来源比较中。'},
 {'type':'p','text':'为补充评估SPP1来源覆盖，将原C4、C10和C15成员合计为联合候选群，共4,502个细胞。在原先固定的16名同文库配对患者中，联合群占成熟巨噬细胞的比例中位数为9.76%，高于C4的7.41%；其贡献的捕获SPP1原始RNA计数份额中位数为59.63%，高于C4的54.59%。联合描述因此提供了更完整的SPP1来源视角。'},
 {'type':'p','text':'三群的跨患者支持范围不同（图8B）。C4出现在37名患者和五项研究，26名患者满足群内与其余成熟巨噬均至少20个细胞的配对条件；C10仅2名患者满足该条件，且97.64%的细胞来自同一患者；C15有9名有效配对患者，并具有明确的增殖状态特征。依据候选鉴定阶段的患者覆盖、多基因身份及研究内方向一致性，后续功能、调控和状态排序继续以C4为主分析对象，C10和C15作为补充候选保留。'},
 {'type':'figure','text':'Figure 8: figures/Figure8_three_SPP1_candidates.pdf; candidate comparison before C4 main analysis'},
 {'type':'caption','text':'图8 三个SPP1高表达候选群的标志特征与患者支持','keep_with_next':True},
 {'type':'note','text':'A：横轴为C4、C10、C15，纵轴为代表marker；颜色为原12个成熟巨噬候选群患者等权表达的行Z分数。SPP1在三群中表达，GPNMB/CCL18/CYP27A1、USP53/ALOX5AP及MKI67/TOP2A分别区分其状态特征。B：横轴为各群细胞在五项研究中的构成比例，纵轴为候选群；右侧n为两侧均至少20个细胞的有效配对患者数。C4具有更广泛的患者支持，C10高度集中于单一研究及患者，C15为增殖相关候选。图中研究构成按细胞计数描述，统计样本量仍以患者为单位。'}
]
b[i:i]=new
j=next(i for i,x in enumerate(b) if x.get('type')=='h2' and x.get('text','').startswith('2.10 '))
b[j]['text']='2.10 C4的多基因身份与患者一致性支持后续主分析'
b[j+1]['text']='基于上述候选比较，C4作为后续主分析状态。该群包含2,674个细胞，出现在37名患者及五项研究。GPNMB、CCL18、CYP27A1、LPL、APOE与SPP1的共同富集构成其多基因身份；SPP1由分群后的标志鉴定识别。'
# Add the descriptive union aggregation method without introducing a new inferential endpoint.
mi=next(i for i,x in enumerate(b) if x.get('type')=='h2' and x.get('text','').startswith('1.11 '))
b.insert(mi,{'type':'p','text':'候选谱系比较并列描述C4、C10和C15的代表marker、研究构成及有效患者数。补充来源分析在固定16名同文库配对患者内，分别计算C4及三群联合的细胞占比和SPP1原始RNA计数份额，再按患者取中位数；分母均为该文库的成熟巨噬细胞或其SPP1原始计数总量。'})
# Align the summary and discussion with the same primary-analysis choice.
b[1]['text']='项目 26210　候选群比较与C4主分析增补版　2026-09-18'
for x in b:
 t=x.get('text','')
 if t.startswith('在肿瘤代谢程序及具体生化反应的基础上'):
  x['text']=t.replace('独立的无监督髓系细分识别出2,674个细胞的SPP1富集C4；','无监督髓系鉴定及候选比较将2,674个细胞的SPP1富集C4确定为主分析状态；')
 if t.startswith('数量与信号贡献进一步区分了状态定位和总体关系'):
  x['text']='数量与信号贡献进一步区分了状态定位和总体关系。C4集中贡献捕获的SPP1 RNA，C10和C15等状态也具有较高SPP1表达；三群联合可补充呈现信号来源覆盖。考虑C4的患者支持、多基因身份与研究内方向一致性，本研究以整体巨噬承载连续表达关系，以C4开展主要功能、调控和状态排序分析，其余候选用于补充解释细胞异质性。'
# Original numerical results remain intact.
oldtables=[x for x in a if x['type']=='table'];newtables=[x for x in b if x['type']=='table']
assert len(oldtables)==len(newtables)==31
assert all(x==y for x,y in zip(oldtables,newtables))
assert sum(x['type']=='figure' for x in b)==8
(d/'report/report_materialized.json').write_text(json.dumps(b,ensure_ascii=False,indent=2))
(d/'qc/CONTENT_VALIDATION.json').write_text(json.dumps({'original_tables_identical':31,'original_figures_retained':7,'new_figures':1,'added_section':'2.9','C4_remains_primary':True,'no_new_inferential_results_added_to_main':True,'unchanged_union_full_analysis_archived':True},indent=2))
shutil.copy2(r/'phase200_union_C4_C10_C15_reassessment_20260918_v1/FULL_REASSESSMENT_CN.md',d/'supplement/三候选联合复评_完整补充记录.md')
shutil.copy2(old/'supplement/MECHANISM_EXPLORATION_SUPPLEMENT_CN.md',d/'supplement/原机制探索补充材料.md')
# Reuse the original typesetting engine and append Figure8.
s=(old/'scripts/render_report.py').read_text().replace('phase191_discovery_order_report_20260917_v1','phase201_three_candidates_C4_main_report_20260918_v4').replace('macrophage_integrated_report_v3','macrophage_C4_primary_report_v4')
s=s.replace("'Figure7_trajectory_candidate_composite_v1.pdf']","'Figure7_trajectory_candidate_composite_v1.pdf','Figure8_three_SPP1_candidates.pdf']")
s=s.replace("figi+=1;story +=[Spacer(1,8),Fig(figi),Spacer(1,5)]","figi+=1; fobj=Fig(figi); gap=Spacer(1,5)\n  if figi==8:fobj.keepWithNext=True;gap.keepWithNext=True\n  story +=[Spacer(1,8),fobj,gap]")
s=s.replace("story.append(para(x['text'],typ));md+=","pobj=para(x['text'],typ)\n  if x.get('keep_with_next'):pobj.keepWithNext=True\n  story.append(pobj);md+=")
(d/'scripts/render_report.py').write_text(s)
print('REPORT_CONTENT_READY',len(b))
