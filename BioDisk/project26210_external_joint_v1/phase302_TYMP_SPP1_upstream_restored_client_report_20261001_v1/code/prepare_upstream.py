from pathlib import Path
import ast
W=Path(__file__).resolve().parent
s=(W/'base_builder.py').read_text(encoding='utf-8')
start=s.index("para('SPP1相关巨噬与基质免疫屏障的空间证据'")
data=r'''
Q=R/'phase280_TYMP_source_niche_function_20260928_v1'
K=R/'phase282_C4_specificity_detection_regulation_20260929_v1'
up_paths=[Q/'tables/WITHIN_MACRO_SUMMARY.tsv',Q/'tables/CONTEXT_SENSITIVITY_SUMMARY.tsv',K/'tables/TYMP_INTENSITY_SPP1_DETECTION_SUMMARY.tsv',K/'tables/CONDITIONAL_REMAINING_SUMMARY.tsv',Q/'PROTOCOL_BEFORE_TESTING.md']
for path in up_paths:src[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
(P/'inputs/SOURCE_MANIFEST.json').write_text(json.dumps(src,indent=2))
within,context,detection,conditional=[pd.read_csv(path,sep='\t',float_precision='round_trip') for path in up_paths[:4]]
up=[]
def add_up(scope,model,row,family,source):
 up.append(dict(scope=scope,model=model,n=int(row['n']),positive=int(row['positive']),rho=float(row['median']),P=float(row['P']),FDR=float(row['FDR']),family=family,source=str(source)))
for scope,code in [('C4','C4'),('整体巨噬','WholeMacro')]:
 for model,label in [('unadjusted','未调整'),('adjusted','测序调整')]:
  row=within[within.group.eq(code)&within.endpoint.eq('SPP1')&within.model.eq(model)].iloc[0]
  add_up(scope,label,row,'F1',up_paths[0])
 row=context[context.group.eq(code)&context.endpoint.eq('SPP1')].iloc[0]
 add_up(scope,'环境调整',row,'F2',up_paths[1])
for scope,code in [('C4','C4'),('其他巨噬','OtherMacro')]:
 for model,label in [('core','强度与检出'),('context','强度与检出加环境')]:
  row=detection[detection.scope.eq(code)&detection.model.eq(model)].iloc[0]
  add_up(scope,label,row,'F3',up_paths[2])
for rg in ['ILF2(+)','SPI1(+)']:
 row=conditional[conditional.scope.eq('C4')&conditional.regulon.eq(rg)].iloc[0]
 add_up('C4','再调'+rg[:-3],row,'F4',up_paths[3])
(P/'inputs/UPSTREAM_MATERIALIZED.json').write_text(json.dumps(up,ensure_ascii=False,indent=2))
def vf(x):return f'{x:.3e}' if abs(x)<.001 else f'{x:.6f}'
'''
s=s[:start]+data+s[start:]
s=s.replace('SPP1相关巨噬与基质免疫屏障的空间证据','TYMP与SPP1相关免疫屏障的分段证据')
s=s.replace('肾癌组织中基质活动与T细胞分布的探索性分析','巨噬细胞内表达关联与基质及T细胞空间联系')
s=s.replace('这三项结果为何与屏障有关','这些结果为何与屏障有关')
s=s.replace('证据来自静态空间关联，尚未直接证明T细胞受到物理阻挡或迁移能力下降。','屏障的组织层证据来自静态空间关联，尚未直接证明T细胞受到物理阻挡或迁移能力下降。')
s=s.replace("h('摘要')", "h('摘要')\npara('单细胞数据表明，TYMP表达与SPP1表达在巨噬细胞内存在患者一致的正向联系。C4中，15例患者有13例正向，P=0.007385、FDR=0.029541；进一步调整干扰素和缺氧状态后15/15例正向，P=6.104×10⁻⁵、FDR=3.662×10⁻⁴。整体巨噬也保留正关联。这些结果将TYMP与SPP1相关巨噬状态联系起来，为其下游基质免疫屏障研究提供了具体的上游表达证据。')",1)
s=s.replace("['哪些细胞与基质活动相联系'", "['TYMP与SPP1巨噬有何联系','C4及整体巨噬内TYMP与SPP1正关联通过FDR','提供TYMP与SPP1状态相联系的上游表达证据'],\n['哪些细胞与基质活动相联系'",1)
s=s.replace("h2('1.1 围绕肿瘤边界比较外侧基质与内外T细胞分布')", "h2('1.1 单细胞和空间数据分别提供表达及组织证据')\npara('巨噬细胞内表达分析沿用已冻结C4及整体成熟巨噬，主要TYMP—SPP1模型分别有15例和29例有效患者。逐患者计算细胞内偏秩相关，再以患者为独立单位汇总并作双侧符号检验，P及FDR继承原家族；完整资格、调整和模型对照见附表S7。以下空间分析独立沿用9例患者，不将三种人数相加。')")
# Insert the upstream findings as the first results page, preserving downstream paragraphs and tables.
for old,new in [('2.4 逐例','2.5 逐例'),('2.3 巨噬','2.4 巨噬'),('2.2 基质','2.3 基质'),('2.1 七位','2.2 七位')]:s=s.replace(old,new)
insert=r'''
newpage();h('2 结果')
h2('2.1 TYMP与SPP1在巨噬细胞内呈患者一致的正关联')
para('在C4这一已识别的SPP1富集巨噬状态中，同一患者内TYMP表达较高的细胞，SPP1表达也倾向于较高。调整测序深度和复杂度后，15例中13例为正向，中位ρ=0.134872，P=0.007385、FDR=0.029541，达到校正后显著。')
para('进一步纳入干扰素α、干扰素γ及缺氧程序后，15/15例为正向，中位ρ=0.156478，P=6.104×10⁻⁵、FDR=3.662×10⁻⁴。所测共同炎症和缺氧背景未解释掉这项表达联系。')
para('整体成熟巨噬中，在控制测序指标及冻结亚群后，29例中24例为正向，中位ρ=0.049133，P=0.000546、FDR=0.002184；加入环境程序后仍为24/29例正向。整体巨噬的5项研究中位方向均为正；去掉最大研究后仍13/17例正向，C4中仍5/6例正向。')
para('表1 TYMP与SPP1表达联系的主要患者级结果','Caption')
main_up=[x for x in up if x['family'] in ['F1','F2'] and x['model']!='未调整']
table(['细胞范围','模型','患者n','正向人数','中位ρ','原始P','FDR'],[[x['scope'],x['model'],str(x['n']),f"{x['positive']}/{x['n']}",f"{x['rho']:.6f}",vf(x['P']),vf(x['FDR'])] for x in main_up],[1.15,1.25,.7,1,1.25,1.45,1.45])
note('测序调整：C4控制深度和复杂度，整体巨噬同时控制亚群；环境调整在此基础上加入干扰素α 干扰素γ和缺氧。P来自患者效应的双侧符号检验。原FDR家族完整保留，见附表S7。')
para('这项正面结果说明，TYMP与SPP1相关巨噬状态已有明确的表达层面联系。后续空间结果进一步检验SPP1相关巨噬、外侧基质和T细胞分布是否对应，从而形成上游细胞状态与下游组织位置的分段证据。')
newpage()
'''
s=s.replace("newpage();h('2 结果')",insert,1)
s=s.replace("para('表1 血管调整后SPP1", "para('表2 血管调整后SPP1")
s=s.replace('表1逐例列出','表2逐例列出')
old="整体研究以TYMP相关代谢背景及SPP1相关巨噬发现为起点。本部分补充其下游的空间问题：SPP1相关巨噬状态是否与基质环境相联系，以及同一基质环境是否对应T细胞内外分布变化。这样可以将TYMP研究背景与免疫屏障研究问题衔接，但这9例空间分析没有建立TYMP干预经SPP1和基质改变而影响T细胞进入的因果链。髓系TYMP相关结果也不能替代肿瘤来源TYMP的作用。"
new="TYMP与SPP1之间已有新的正面表达结果：C4及整体巨噬的患者内正关联通过FDR，且在所测环境程序调整后保留。因此，TYMP在本部分有具体的上游证据；SPP1相关巨噬与基质、基质与T细胞位置的空间联系则提供下游证据。不同分析具有各自的患者集合和读出，不能把这些分段结果当作同病例的因果中介。明确的髓系细胞内联系继续保留，也不替代肿瘤来源TYMP的作用。"
assert old in s;s=s.replace(old,new)
# Retain the relevant contrasts while explaining why they do not cancel the positive result.
needle="para('现有读出分别是映射状态"
pos=s.index(needle)
s=s[:pos]+"para('TYMP—SPP1正关联的效应较小，表明表达状态相联系，尚不确定因果方向。进一步的TYMP表达强度与SPP1检出分析，C4环境模型13/14例正向并通过FDR，基础模型未达名义标准；仅在双阳性细胞中比较表达强度时，C4两模型均未达名义标准。C4与其他巨噬的关联差异也未获支持，这些对照不否定已确立的患者内正关联。跨患者肿瘤TYMP—巨噬SPP1模型未获支持，外源TYMP蛋白刺激的公开培养体系SPP1呈下降方向；它们检验的是不同来源和暴露，不能迁移为髓系细胞内联系不存在。')\n"+s[pos:]
# A focused, complete endpoint appendix including positive and non-significant models.
append=r'''
newpage();h('附表S7 TYMP与SPP1表达联系及模型对照')
para('C4沿用冻结成员，不以本轮相关性重新分群。原分析要求C4每患者至少50个细胞，整体巨噬至少100个细胞，两基因各至少10个检出细胞。每位患者内先计算表达秩的条件相关，再以独立患者作双侧符号检验；细胞数不代替患者n。整体巨噬额外控制冻结亚群。')
para('环境程序为固定干扰素α 干扰素γ和缺氧，剔除TYMP及各结局或代谢程序基因。强度与检出模型比较连续TYMP强度秩与SPP1是否检出，并非概率或比值比。条件模型在原环境模型上再纳入去结局基因的ILF2或SPI1靶表达程序。')
table(['范围','模型','n','正向\n人数','中位ρ','原始P','FDR','家族'],[[x['scope'],x['model'],str(x['n']),f"{x['positive']}/{x['n']}",f"{x['rho']:.6f}",vf(x['P']),vf(x['FDR']),x['family']] for x in up],[1.1,2,.5,.8,1.2,1.4,1.4,.5],compact=True)
note('F1为原患者内表达分析，每个细胞范围的6端点×2模型共12项BH；F2为另行登记的环境模型，6端点×2范围共12项BH；F3为TYMP强度与SPP1检出分析，2范围×2模型共4项BH；F4为原靶程序条件模型的28项BH。本表原P及FDR直接继承相应完整家族，不在本表重新校正。')
note('F4列出用于讨论的ILF2及SPI1条件结果；两者靶基因高度重叠，不作为两条独立机制。细胞内正关联不等于TYMP上调SPP1的干预效应。新C4单细胞身份与旧SPP1空间状态分别沿用原定义，不互换，也不将有重叠的患者加成新的样本量。')
'''
s=s.replace("(P/'inputs/DOCUMENT_CONTENT.json').write_text",append+"\n(P/'inputs/DOCUMENT_CONTENT.json').write_text",1)
s=s.replace('26210_SPP1相关基质免疫屏障_客户易读修订版.docx','26210_TYMP与SPP1免疫屏障_补全上游证据版.docx')
s=s.replace('intended_pages=11','intended_pages=13')
ast.parse(s)
(W/'build_upstream_complete.py').write_text(s,encoding='utf-8')
print('Prepared upstream result section, abstract integration, discussion correction and S7')
