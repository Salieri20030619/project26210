from pathlib import Path
import json,re,shutil,hashlib
P=Path(__file__).resolve().parents[1];R=P.parent;O=R/'phase327_abstract_clarity_revision_20261006_v1';F=R/'phase330_figure_F_tissue_immune_20261006_v3'
b=json.loads((O/'report/report_materialized.json').read_text());original=json.loads(json.dumps(b));changes=[]
for i,x in enumerate(b):
 x['_phase327_index']=i
 for k in ['text','caption','source']:
  if isinstance(x.get(k),str):
   old=x[k];x[k]=re.sub(r'图\s*5(?!\d)',lambda m:m[0].replace('5','6'),old)
   if x[k]!=old:changes.append({'index':i,'change':'old figure 5 -> 6','field':k})
 if x.get('type')=='figure' and x.get('figure')==5:x['figure']=6;x['text']='Figure6.pdf'
# Section numbers and all numerical tables remain intact.
adds={202:'细胞生态位及免疫程序概览见图5A。',208:'患者层面的完整免疫程序分布见图5B。',218:'两基因的患者层面糖酵解关联见图5B。',224:'两类程序的供者散点分布见图5F。',254:'同一切片的程序空间分布及全队列患者关联分别见图5C和图5E。',263:'各患者的界面内外差见图5D。'}
for i,s in adds.items():b[i]['text']+=s;changes.append({'index':i,'change':'panel cross-reference','added':s})
b[230]['text']=b[230]['text'].replace('（补充图S2A；','（图5G；补充图S2A；')
b[233]['text']=b[233]['text'].replace('（补充图S2C-D；','（图5H；补充图S2C-D；')
i=next(i for i,x in enumerate(b) if x.get('text','').startswith('主图1-5分别展示'))
b[i]['text']='主图1-6依次展示整合细胞图谱、恶性上皮调控状态、代谢及巨噬信号来源、C4功能和髓系排序、免疫与空间及临床综合结果，以及代谢相关组织空间对应。图5A-H分别对应细胞生态位示意、巨噬免疫程序、同切片空间分布、CD8界面差异、基质程序联系、组织免疫散点、生存曲线和C4时变关联。补充图S1展示两例空间代谢通路及外部ATAC足迹活性投射，补充图S2保留临床交互与敏感性结果。完整统计仍见正文及补充表；图注说明坐标、颜色、统计单位和对应证据层级。'
notes=[
'A：肿瘤细胞、巨噬细胞、成纤维细胞、基质纤维及CD8 T细胞的概念示意，采用AI辅助原创插画；形态、位置及纤维走向不是实测组织影像。右侧概括本项目巨噬细胞内RNA关联，加减号表示正负关联，不代表干预或因果箭头。',
'B：Obradovic队列11名患者的16项程序，宽松巨噬定义共15,056细胞；独立统计单位为患者。每格为患者内控制测序深度、检出基因数、文库及髓系簇后的秩残差相关ρ；两基因区块患者顺序一致，色标-0.4至0.4，不按行标准化。程序均排除TYMP/SPP1，右侧为患者中位ρ，粗体为原224项家族FDR<0.05。T细胞调节行表示巨噬细胞中的对应GO程序。相互控制和去重叠结果见表8；糖酵解删基因敏感性见2.15节。',
'C：GSE175540的GSM5924033整片4,939个QC spot，依次显示冻结的SPP1、FAP去重叠程序及保留ECM程序。三图使用完全相同的真实坐标、方向和范围，不平滑、不重构热点。颜色按各程序片内第2-98百分位显示，不能跨程序比较绝对量；比例尺1 mm，像素换算1.82105 μm。切片用于空间示例，队列汇总见E；程序共分布不等于纯细胞物理邻接。',
'D：GSE175540的20名合格患者，固定400 μm表达映射界面，纵轴为控制深度及总巨噬指标后的CD8谱系RNA内侧减外侧标准化差；每柱一位患者，按标识排序。19/20内侧较低，中位差-0.260230，P=0.00004005、FDR=0.00016022（原4项家族）。局部SPP1富集与该内外差未获支持，P=0.115318、FDR=0.461273；界面不等同独立病理标注。',
'E：GSE175540共23名患者、24切片，灰线连接同一患者的FAP-SPP1及SPP1-ECM深度调整ρ，粗横线为中位数；每患者各贡献两个不同关系。前者23/23正，中位ρ=0.152709，P=2.38419e-7、FDR=4.76837e-7；后者20/23正，中位ρ=0.163723，P=FDR=0.000488281。分别沿用原8项耦合及12项过程家族。三个负向患者保留；双向连接表示相关。',
'F：公开组织bulk的44名合格独立供者，每点一位供者，全部保留。横轴SPP1 log2(TPM+1)，纵轴为既定细胞毒或功能障碍程序成员的log2(TPM+1)等权均值。Spearman ρ分别为0.670754和0.491896，原双侧置换P分别为9.99990e-6和0.000789992，两个主要端点BH FDR分别为1.99998e-5和0.000789992。组织层两类程序共同升高不表示细胞毒功能下降，也不直接定位到CD8或巨噬细胞。',
'G：CheckMate nivolumab组181人、123个死亡事件。SPP1相关五基因评分为SPP1、FN1、CLEC5A、OLR1、IL1RN表达均值，以完整311人队列中位数分组，高84人、低97人。阶梯和阴影为KM及95%CI，加号为删失，下方为在险人数。Log-rank P=0.014836、FDR=0.059343，为名义显著。连续模型及交互保留在正文和补充图S2。',
'H：冻结19基因C4评分排除TYMP及五基因SPP1评分成员；CheckMate nivolumab完整模型164人、145个PFS事件，控制临床变量、五基因评分及总巨噬估计。纵轴为每升高1 SD的条件HR（对数刻度），曲线及带为估计和点态95%CI，下方为在险人数；固定时间函数β(t)=β0+β1log(1+t)，不搜索截点。联合P=0.0007173、FDR=0.001435；12个月后灰底提示风险集减少，不表示C4功能逆转。',
'各面板分别对应单细胞、空间、组织bulk和临床队列，人数不相加，不将分段关联拼为同患者因果链。绘图沿用既有冻结数值、程序和坐标，未新增检验或筛选。'
]
blocks=[{'type':'figure','figure':5,'text':'Figure5.pdf','_new_figure':'Phase330 approved v3'},{'type':'caption','text':'图5 TYMP–SPP1免疫程序、空间生态位与临床关联','keep_with_next':True}]+[{'type':'note','text':t} for t in notes]
pos=next(i for i,x in enumerate(b) if x.get('text','').startswith('2.18 '));b[pos:pos]=blocks
(P/'report/report_materialized.json').write_text(json.dumps(b,ensure_ascii=False,indent=2))
for f in (O/'figures').iterdir():shutil.copy2(f,P/'figures'/('Figure6.pdf' if f.name=='Figure5.pdf' else f.name))
shutil.copy2(F/'delivery/26210_immune_spatial_clinical_A4_v3.pdf',P/'figures/Figure5.pdf')
s=(O/'code/render_report.py').read_text().replace('phase327_abstract_clarity_revision_20261006_v1',P.name).replace('26210_manuscript_immune_integrated_v3','26210_manuscript_immune_figure_integrated_v4')
s=s.replace('for n in [1,2,3,4,5]','for n in [1,2,3,4,5,6]').replace("5:set('ABCDEFGHI')","5:set('ABCDEFGH'),6:set('ABCDEFGHI')").replace('[1-5]','[1-6]')
s=s.replace('elif n==5:',"elif n==5:\n   anchors=[('A',9,17),('B',9,73),('C',9,145),('D',116,145),('E',9,197),('F',116,197),('G',9,244),('H',116,244)]\n  elif n==6:",1)
s=s.replace('if n==5 else','if n==6 else').replace('figures=7','figures=8')
s=s.replace("else 125))", "else 80 if x['text'].startswith('2.26 ') else 125))")
s=s.replace("  tt=Table(data,", "  assert len(weights)==n and all(v>0 for v in weights), (x['caption'],n,weights)\n  tt=Table(data,")
(P/'code/render_report.py').write_text(s)
assert [x for x in b if x['type']=='table']==[dict(x,_phase327_index=i) for i,x in enumerate(original) if x['type']=='table']
assert [x['text'] for x in b if x['type'] in ['h1','h2']]==[x['text'] for x in original if x['type'] in ['h1','h2']]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(P/'qc/INTEGRATION.json').write_text(json.dumps({'source_manuscript':str(O/'report/26210_manuscript_immune_integrated_v3.pdf'),'source_sha256':sha(O/'report/26210_manuscript_immune_integrated_v3.pdf'),'new_figure_sha256':sha(P/'figures/Figure5.pdf'),'insert_before':'2.18','figure_number_map':{'old 5':6,'new':5},'existing_tables_unchanged':41,'headings_unchanged':True,'changes':changes},ensure_ascii=False,indent=2));print('INTEGRATED',len(b),'blocks')
