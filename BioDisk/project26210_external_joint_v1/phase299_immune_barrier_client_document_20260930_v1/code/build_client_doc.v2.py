from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from xml.sax.saxutils import escape
import pandas as pd,numpy as np,json,hashlib,subprocess
P=Path(__file__).resolve().parents[1];R=P.parent
with (P/'logs/memory_gate.log').open('a') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,check=True)
sp=R/'phase298_matrix_to_immune_entry_20260930_v1/tables'
paths=[sp/'COHORT_RESULTS.tsv',sp/'SECTION_RESULTS.tsv',R/'phase213_immune_restart_execution_20260919_v1/I2_main_statistics/PRIMARY_12_HYPOTHESES.tsv']
src={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
co=pd.read_csv(paths[0],sep='\t');sec=pd.read_csv(paths[1],sep='\t');im=pd.read_csv(paths[2],sep='\t')
(P/'inputs/SOURCE_MANIFEST.json').write_text(json.dumps(src,indent=2))
W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'
ns=f'xmlns:w="{W}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
parts=[];audit=[];pending_page=False
def run(t,bold=False,size=None):
 pr=('<w:b/>' if bold else '')+(f'<w:sz w:val="{size}"/>' if size else '')
 ss=[]
 for j,line in enumerate(str(t).split('\n')):
  if j:ss.append('<w:br/>')
  ss.append('<w:t xml:space="preserve">'+escape(line)+'</w:t>')
 return '<w:r><w:rPr>'+pr+'</w:rPr>'+''.join(ss)+'</w:r>'
def para(t,style='Normal',align=None,keep=False):
 global pending_page
 pr=f'<w:pStyle w:val="{style}"/>'+('<w:jc w:val="'+align+'"/>' if align else '')+('<w:keepNext/>' if keep else '')
 if pending_page:pr+='<w:pageBreakBefore/>';pending_page=False
 parts.append('<w:p><w:pPr>'+pr+'</w:pPr>'+run(t)+'</w:p>');audit.append(dict(type=style,text=str(t)))
def h(t):para(t,'Heading1')
def h2(t):para(t,'Heading2')
def note(t):para(t,'Note')
def newpage():
 global pending_page
 pending_page=True
def table(headers,rows,weights,size=20):
 total=10080;widths=[int(total*x/sum(weights)) for x in weights];widths[-1]+=total-sum(widths)
 borders=''.join(f'<w:{side} w:val="single" w:sz="4" w:color="D9D9D9"/>' for side in ['top','left','bottom','right','insideH','insideV'])
 tbl='<w:tbl><w:tblPr><w:tblW w:w="10080" w:type="dxa"/><w:tblLayout w:type="fixed"/><w:tblBorders>'+borders+'</w:tblBorders><w:tblCellMar><w:top w:w="75" w:type="dxa"/><w:left w:w="100" w:type="dxa"/><w:bottom w:w="75" w:type="dxa"/><w:right w:w="100" w:type="dxa"/></w:tblCellMar></w:tblPr><w:tblGrid>'+''.join(f'<w:gridCol w:w="{x}"/>' for x in widths)+'</w:tblGrid>'
 for ri,row in enumerate([headers]+rows):
  tbl+='<w:tr><w:trPr><w:cantSplit/>'+('<w:tblHeader/>' if ri==0 else '')+'</w:trPr>'
  for ci,t in enumerate(row):
   fill='DDEBF7' if ri==0 else 'FFFFFF'
   alg='left' if ci==0 else 'center'
   tbl+=f'<w:tc><w:tcPr><w:tcW w:w="{widths[ci]}" w:type="dxa"/><w:shd w:fill="{fill}"/><w:vAlign w:val="center"/></w:tcPr><w:p><w:pPr><w:jc w:val="{alg}"/><w:spacing w:before="0" w:after="0" w:line="260" w:lineRule="auto"/></w:pPr>'+run(t,ri==0,size)+'</w:p></w:tc>'
  tbl+='</w:tr>'
 tbl+='</w:tbl>';parts.append(tbl);parts.append('<w:p><w:pPr><w:spacing w:after="40" w:line="100" w:lineRule="auto"/></w:pPr></w:p>')
 audit.append(dict(type='table',headers=headers,rows=rows))
def val(x,d=5):return '—' if pd.isna(x) else f'{x:.{d}f}'
# 1
para('SPP1相关免疫微环境与基质屏障的探索性分析','Title')
para('免疫接收端分析及T细胞空间分布补充报告','Subtitle')
para('项目26210    2026年9月30日','Note')
h('摘要')
para('本研究围绕SPP1相关髓系状态、间质基质重塑与免疫细胞肿瘤内分布开展分层分析。空间结果探索性提示，在TJ-RCC队列中，外侧细胞外基质组织程序越高，T细胞信号相对于肿瘤外侧越低；该负向联系在控制外侧免疫供给及血管相关组成后，见于全部9位可评价患者。')
para('200 μm主要模型中，9位患者有7位呈负向关联，中位偏相关系数为−0.192；血管调整模型中9位患者全部负向，中位系数为−0.144。400 μm分析亦保持总体负向方向。此前SPP1相关巨噬状态与间质活化的空间联系，以及髓系TYMP与间质胶原成熟的关联，为这一末端空间发现提供了髓系与基质过程背景。')
para('综合证据支持将“SPP1相关髓系状态—间质基质组织—T细胞相对肿瘤内不足”作为探索性机制框架。当前证据以患者方向一致性和静态空间关联为主；单片严格空间检验及跨队列重复的限制集中见讨论。')
h('主要发现')
table(['分析层级','主要结果','证据定位'],[
['髓系与间质过程','髓系TYMP—胶原成熟\n13例 ρ=0.810','患者关联\n原P=0.00895 FDR=0.02685'],
['空间状态协作','SPP1状态—间质活化\n15例中12例正向','400 μm血管调整\n原P=0.03534 FDR=0.10602'],
['基质与T细胞分布','TJ-RCC血管调整\n9例中9例负向','200 μm 中位ρ=−0.144\n描述性方向一致']], [2,3.5,3.7])
note('ρ为相应分析的相关效应；不同分析的患者和测量层级分别统计。“全部同向”专指9位患者的关联方向，不表示所有免疫细胞类型或全部区域的绝对信号均下降。')
#2
newpage();h('1 方法')
h2('数据与独立统计单位')
para('免疫分析采用冻结的成熟巨噬、C4及淋巴细胞身份。CD8功能纳入16位患者，Treg组成及NK功能纳入13位患者；C4群内表达按端点纳入10至11位。患者为独立单位，同患者合格文库等权汇总。')
para('空间分析采用TJ-RCC、GSE175540和GLMF（GSE210041）数据及计算肿瘤边界。33张候选切片沿用既定几何、区域和自由度资格。各端点分别计n；TJ-RCC总T为9位、CD8为4位，GLMF作病例描述。')
h2('免疫功能与基质程序')
para('巨噬SPP1采用原始计数伪汇总的log2(CPM+1)表达；C4比例以成熟巨噬为分母。CD8及NK细胞毒基因为NKG7、GNLY、GZMB、PRF1、CTSW；CD8功能障碍基因为PDCD1、CTLA4、LAG3、HAVCR2、TIGIT、TOX、CXCL13、ENTPD1。Treg组成以常规CD4 T细胞和Treg之和为分母。')
para('空间暴露采用冻结的外侧ECM组织及ECM重塑程序，基因集合分别为FN1、COL1A1、COL3A1、LOX、LOXL2、PLOD2，以及MMP2、MMP14、TIMP1、TIMP2、ADAM12。沿用既有归一化表达程序，按固定边界区域汇总。')
h2('肿瘤内外免疫分布读出')
para('总T读出采用细胞映射丰度后验第5百分位估计值（q05），分别计算肿瘤内外区域均值，经log1p转换后取内侧减外侧的差值。CD8采用既定表达程序的内外差。负向关联表示外侧基质程序较高时，相对内侧免疫信号较低。')
para('200 μm为主要尺度，400 μm为敏感性分析。基础模型控制内外深度、巨噬与恶性比例、两侧区域数量比、坐标及连通分量；T模型另控制外侧T供给，CD8模型另控制外侧CD8和T供给。血管模型再控制内外血管及壁细胞比例。')
h2('统计分析')
para('免疫功能关联采用研究调整的秩偏相关及99,999次研究内置换。空间分析采用秩残差相关和999次受协变量约束的谱随机化；同患者多切片以Fisher变换后均值合并。独立患者达到10例时进行双侧Wilcoxon检验；少于10例报告方向、效应量和患者一致性。')
para('原始P<0.05用于名义信号识别；多重比较采用同类假设家族内的Benjamini–Hochberg校正。基质—免疫比较在每个队列、尺度和模型内包含两种ECM程序与两个免疫终点，共4项假设。空间模型P与患者层面P分别保留。')
#3
newpage();h('2 结果')
h2('SPP1相关髓系状态具有免疫接收背景')
para('在具备C4与其他成熟巨噬配对资格的患者中，C4的SPP1表达总体较高：CD8接收端配对范围内，11位患者中9位为正向，中位log2(CPM+1)差为1.105。部分免疫接收细胞具有CD44或ITGA4/ITGB1完整受体的联合检出资格，提供了SPP1接收的表达背景。')
para('免疫功能端点未呈现统一的研究调整后关联。整体成熟巨噬SPP1与CD8细胞毒、CD8功能障碍相关程序、Treg/CD4及NK细胞毒的原始P均大于0.05。C4比例与CD8功能障碍的未经调整关联为正（ρ=0.597，P=0.01461），研究调整后为ρ=0.325、P=0.35629、FDR=0.84644；留一患者与留一研究时保留正方向。完整12项结果见附表S1。')
h2('髓系状态与间质胶原成熟具有过程联系')
para('13位同患者同文库对象中，髓系TYMP与间质胶原成熟程序呈正向关联（ρ=0.810，P=0.00895，FDR=0.02685）。独立GSE264586基线4位患者中，巨噬TYMP与严格成纤维候选胶原成熟的描述性ρ为1.0，SPP1对应ρ为0.8。该外部对应保留了髓系与基质过程的同向关系。')
para('13位患者中，巨噬SPP1与胶原成熟的关联为正（ρ=0.667，P=0.13583，FDR=0.33958），13次留一患者均为正向。空间层面，GSE175540在400 μm血管调整模型中，SPP1相关巨噬状态与间质活化呈名义正向关联：15位患者中12位正向，中位ρ=0.225，P=0.03534，FDR=0.10602。')
para('表1 髓系与间质过程的主要关联','Caption')
table(['关系','独立患者n','效应ρ','原始P','FDR'],[
['髓系TYMP与胶原成熟','13','0.8096','0.00895','0.02685'],
['SPP1与胶原成熟','13','0.6666','0.13583','0.33958'],
['独立病例TYMP与胶原成熟','4','1.0000','—','—'],
['SPP1状态与空间间质活化','15','0.2247','0.03534','0.10602']], [3.4,1.1,1.1,1.2,1.2])
note('第一项保留原胶原成熟分析的P/FDR；第二项来自统一13患者的分段关联，BH家族不同，二者不构成独立重复。空间状态协作项限于400 μm血管调整模型。4例结果为描述性。')
h2('间质表达变化与既有ECM特征方向对应')
para('接收端状态比较的7位患者中，原始文献[4]的27个ECM差异或趋势基因有21个中位方向对应；排除用于定义已有状态或程序的基因后，剩余20个中有15个对应。进一步限定原文P<0.05的非定义基因，13/18个方向对应。该结果为间质变化的基质过程解释提供补充。')
#4
newpage();h2('TJ RCC中外侧ECM组织与相对内侧T细胞信号呈负向联系')
para('200 μm主要模型中，TJ-RCC 9位患者有7位呈负向关联，中位ρ为−0.192。进一步控制血管及壁细胞组成后，9位患者全部呈负向关联，中位ρ为−0.144。全部合格区域均检出正的外侧T细胞映射丰度估计。')
para('这一结果显示，在外侧免疫供给及主要组成因素得到调整后，外侧ECM组织程序与T细胞相对肿瘤内不足具有患者间一致的方向。各患者效应见表2。')
para('表2 TJ RCC患者的血管调整结果','Caption')
z=sec[(sec.dataset=='TJRCC')&(sec.radius==200)&(sec.model=='vascular')&(sec.program=='ECM_organization')&(sec.endpoint=='T')&sec.reason.eq('OK')].copy()
assert len(z)==9 and z.effect.lt(0).all()
rows=[[str(r.patient).replace('TJRCC:',''),val(r.effect,4),val(r.spatial_P,3),val(r.spatial_FDR,3)] for r in z.itertuples()]
table(['患者','偏相关ρ','空间模型P','空间模型FDR'],rows,[1.6,1.6,2,2])
note('独立生物学样本量为9，故不进行患者层面推断检验。表内P/FDR为每张切片的严格空间零模型结果；FDR在该切片、尺度和模型的4项基质—免疫假设内计算。')
para('表3 TJ RCC总T结果的尺度与模型对照','Caption')
zs=co[(co.dataset=='TJRCC')&(co.program=='ECM_organization')&(co.endpoint=='T')].sort_values(['radius','model'])
table(['尺度','模型','患者n','负向患者','中位ρ'],[[str(r.radius)+' μm','基础' if r.model=='core' else '血管调整',str(r.n),f'{r.negative}/{r.n}',val(r.median_rho,4)] for r in zs.itertuples()],[1.4,1.8,1,1.5,1.5])
note('200 μm基础模型为主分析，其余为敏感性分析。负向表示ECM组织程序越高，T细胞内侧减外侧的信号差越低。')
#5
newpage();h('3 讨论')
h2('患者方向一致性支持基质参与免疫分布的探索性解释')
para('本研究将免疫接收背景、髓系与间质过程联系，以及末端空间分布结果衔接起来。最直接的末端线索来自TJ-RCC：外侧ECM组织程序与相对肿瘤内T细胞信号在主要模型中总体负向，血管调整后9例方向完全一致。该观察为SPP1相关髓系—间质重塑参与免疫空间组织提供了探索性依据。')
para('原始研究已在人肺肿瘤活切片中观察到，肿瘤巢周围基质纤维的密度与方向影响T细胞运动，胶原酶处理可增加T细胞与癌细胞接触[1]；DDR1研究进一步将胶原纤维排列与免疫排斥联系起来[2]。SPP1相关巨噬与成纤维细胞协作也有其他癌种的单细胞及空间依据[3]。这些文献提供机制背景，本项目提供ccRCC中的表达和空间方向对应。')
h2('信号的测量层级与适用范围')
para('“9例全部同向”是9个患者内偏相关系数均为负，不表示9例T细胞都发生纵向下降，也不表示每个区域或所有免疫亚型均减少。总T使用映射丰度代理，CD8使用表达程序；二者均为静态分布读出。当前尚无本项目纤维取向、孔隙、实时迁移或基质干预测量，因此结论定位为探索性空间关联。')
para('TJ-RCC的9例未作患者推断检验；严格切片空间零模型的原始P均未小于0.05，最小为0.057。计算边界未经独立病理金标准确认，spot混合成分和映射误差仍可影响结果。内外差及外侧供给调整增强了“相对内侧不足”的解释，但不能完全区分进入、滞留、局部增殖和存活。')
h2('跨队列差异和关键反向结果')
para('GSE175540未在200 μm主模型中重复这一关系。其400 μm基础模型中，ECM组织与CD8内外差呈反向于排斥假设的名义正关联（n=12，ρ=0.169，P=0.04248，FDR=0.16992）；血管调整后不保留。TJ-RCC的ECM重塑程序与总T多数为正向，说明基质组织与基质重塑程序不能合并解释为普遍排斥。GLMF仅2例，总T与CD8方向也不完全一致。')
para('免疫功能分析的12项研究调整检验均未获名义支持，受体检出亦不等于功能激活。13例间质主关联的接收池主要为广义壁细胞；独立4例为转移部位的严格成纤维候选。SPP1状态—间质空间协作限于特定尺度与模型，7例ECM文献对应并非独立人群验证。上述边界不否定已观察到的方向一致性，但限定了其外推范围。')
h2('综合结论')
para('现有结果支持“SPP1相关髓系状态与基质过程相联系，基质组织在部分ccRCC病例中对应T细胞相对肿瘤内不足”的探索性框架。其中TJ-RCC血管调整后9例同向，是可保留的下游空间线索；完整免疫屏障及因果链尚未由当前数据确立。')
#6
newpage();h('参考文献')
for t in [
'[1] Salmon H, et al. Matrix architecture defines the preferential localization and migration of T cells into the stroma of human lung tumors. J Clin Invest. 2012;122:899–910. DOI 10.1172/JCI45817. https://www.jci.org/articles/view/45817',
'[2] Sun X, et al. Tumour DDR1 promotes collagen fibre alignment to instigate immune exclusion. Nature. 2021;599:673–678. DOI 10.1038/s41586-021-04057-2.',
'[3] Qi J, et al. Single-cell and spatial analysis reveal interaction of FAP+ fibroblasts and SPP1+ macrophages in colorectal cancer. Nat Commun. 2022;13:1742. DOI 10.1038/s41467-022-29366-6.',
'[4] Tomko LA, et al. Targeted matrisome analysis identifies thrombospondin-2 and tenascin-C in aligned collagen stroma from invasive breast carcinoma. Sci Rep. 2018;8:12941. DOI 10.1038/s41598-018-31126-w.']:
 para(t,'Reference')
h('附表S1 免疫功能关联的完整主要结果')
predictors={'mature_SPP1_pseudobulk':'整体巨噬SPP1','C4_fraction':'C4比例','C4_SPP1_pseudobulk':'C4内SPP1'}
endpoints={'CD8_cytotoxicity':'CD8细胞毒','CD8_dysfunction':'CD8功能障碍','Treg_fraction_within_CD4':'Treg/CD4','NK_cytotoxicity':'NK细胞毒','Treg_fraction':'Treg/CD4','Treg_CD4':'Treg/CD4'}
print('IMMUNE_LABELS',im.predictor.unique(),im.endpoint.unique(),flush=True)
ir=[]
for r in im.itertuples():
 pred=predictors.get(r.predictor,r.predictor);ep=endpoints.get(r.endpoint,r.endpoint)
 ir.append([pred,ep,str(r.n_patients),val(r.rho,4),val(r.p_value),val(r.FDR)])
table(['预测量','免疫端点','n','调整ρ','原始P','FDR'],ir,[2.0,2.0,.5,1.0,1.2,1.2],20)
note('采用研究调整的患者级关联，双侧研究内置换。BH分别在预先定义的整体SPP1 CD8两项、C4 CD8四项、C4 Treg两项及C4 NK两项家族内计算。整体SPP1的Treg和NK为单独预设终点，不作额外多重校正；FDR以原结果记录为准。')
#7-9
for idx,(ds,name) in enumerate([('TJRCC','TJ RCC'),('GSE175540','GSE175540'),('GLMF_GSE210041','GLMF GSE210041')],start=2):
 newpage();h(f'附表S{idx} {name}基质与免疫分布的完整结果')
 para('每行汇总独立患者内的空间偏相关。负向表示外侧基质程序较高时，免疫信号的相对内侧分布较低。基础模型已控制外侧免疫供给；血管模型进一步调整血管及壁细胞组成。','Note')
 df=co[co.dataset.eq(ds)].sort_values(['radius','model','program','endpoint'])
 rows=[]
 for r in df.itertuples():
  rows.append([str(r.radius),'基础' if r.model=='core' else '血管','组织' if r.program=='ECM_organization' else '重塑',r.endpoint,str(r.n),f'{r.negative}/{r.n}',val(r.median_rho,4),val(r.P),val(r.FDR)])
 table(['尺度\nμm','模型','ECM\n程序','免疫\n读出','n','负向\n人数','中位ρ','原始P','FDR'],rows,[.9,.9,.9,.8,.5,.8,1.3,1.3,1.3],20)
 note('200 μm基础模型为主分析。患者n<10时不计算推断P或FDR，表中以“—”表示。n≥10时采用双侧Wilcoxon检验；BH家族为同队列、同尺度、同模型下的4项比较。全零或不可评价与无统计支持分别记录。')
 if ds=='TJRCC':
  para('总T端点9位患者，CD8端点4位患者；两种读出不可混用。ECM组织程序在总T中的负向一致性，未在ECM重塑程序中呈现相同模式。','Note')
 elif ds=='GSE175540':
  para('400 μm基础模型中ECM组织—CD8为正向名义显著（P=0.04248、FDR=0.16992），方向与排斥假设相反；该结果及其血管调整后的减弱均保留。','Note')
 else:
  para('GLMF两个病例分别对应GSE210041中的Patient1及Patient7。400 μm仅1位患者满足模型资格，不作为新增患者或独立队列。','Note')
# archive narrative and tabular representation for reproducibility
(P/'inputs/DOCUMENT_CONTENT.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
# Define native Word styles; no decorative rules.
def st(sid,name,size,bold=False,after=100,before=0,keep=False):
 return f'<w:style w:type="paragraph" w:styleId="{sid}"><w:name w:val="{name}"/><w:basedOn w:val="Normal"/><w:pPr><w:spacing w:before="{before}" w:after="{after}" w:line="300" w:lineRule="auto"/>'+('<w:keepNext/>' if keep else '')+f'</w:pPr><w:rPr><w:color w:val="000000"/><w:sz w:val="{size}"/>'+('<w:b/>' if bold else '')+'</w:rPr></w:style>'
styles=f'<w:styles {ns}><w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="Liberation Serif" w:hAnsi="Liberation Serif" w:eastAsia="Noto Sans CJK SC"/><w:color w:val="000000"/><w:sz w:val="22"/><w:lang w:val="en-US" w:eastAsia="zh-CN"/></w:rPr></w:rPrDefault><w:pPrDefault><w:pPr><w:spacing w:after="100" w:line="300" w:lineRule="auto"/><w:widowControl/></w:pPr></w:pPrDefault></w:docDefaults><w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/></w:style>'
styles+=st('Title','Title',38,True,150,0,True)+st('Subtitle','Subtitle',24,False,120,0,True)+st('Heading1','heading 1',28,True,120,180,True)+st('Heading2','heading 2',24,True,100,120,True)+st('Note','Note',19,False,90)+st('Reference','Reference',19,False,100)+st('Caption','Caption',21,True,80,100,True)+'</w:styles>'
# Letter portrait with 0.75in margins matches skill baseline.
sect='<w:sectPr><w:headerReference w:type="default" r:id="rId2"/><w:footerReference w:type="default" r:id="rId3"/><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1000" w:right="1080" w:bottom="1000" w:left="1080" w:header="420" w:footer="420"/><w:docGrid w:linePitch="312"/></w:sectPr>'
document=f'<w:document {ns}><w:body>'+''.join(parts)+sect+'</w:body></w:document>'
header=f'<w:hdr {ns}><w:p><w:pPr><w:jc w:val="right"/></w:pPr>'+run('项目26210  免疫与基质屏障补充分析',size=17)+'</w:p></w:hdr>'
footer=f'<w:ftr {ns}><w:p><w:pPr><w:jc w:val="center"/></w:pPr>'+run('第 ',size=18)+'<w:fldSimple w:instr="PAGE"/>'+run(' 页',size=18)+'</w:p></w:ftr>'
rels='<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/header" Target="header1.xml"/><Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" Target="footer1.xml"/></Relationships>'
ct='<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/><Override PartName="/word/header1.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.header+xml"/><Override PartName="/word/footer1.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/></Types>'
rootrels='<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/></Relationships>'
dest=P/'delivery/26210_免疫微环境与基质屏障探索性分析_客户交付版.docx'
with ZipFile(dest,'w',ZIP_DEFLATED) as zip:
 for fn,content in {'[Content_Types].xml':ct,'_rels/.rels':rootrels,'word/document.xml':document,'word/styles.xml':styles,'word/header1.xml':header,'word/footer1.xml':footer,'word/_rels/document.xml.rels':rels}.items():zip.writestr(fn,'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'+content)
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in src.items())
(P/'qc/BUILD.json').write_text(json.dumps(dict(status='BUILT_AWAITING_RENDER',document=str(dest),source_unchanged=True,intended_pages=9,body_elements=len(audit),tables=sum(x['type']=='table' for x in audit)),indent=2))
print('BUILT',dest,flush=True)
