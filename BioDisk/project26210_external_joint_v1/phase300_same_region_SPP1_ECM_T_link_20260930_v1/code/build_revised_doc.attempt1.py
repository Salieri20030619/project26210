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

edge=pd.read_csv(P/'tables/SAME_REGION_EDGES.tsv',sep='\t')
align=pd.read_csv(P/'tables/PATIENT_LINK_ALIGNMENT.tsv',sep='\t')
summ=pd.read_csv(P/'tables/ALIGNMENT_SUMMARY.tsv',sep='\t')
for path in [P/'tables/SAME_REGION_EDGES.tsv',P/'tables/PATIENT_LINK_ALIGNMENT.tsv',P/'tables/ALIGNMENT_SUMMARY.tsv']:
 src[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
(P/'inputs/SOURCE_MANIFEST.json').write_text(json.dumps(src,indent=2))
v=align[align.model.eq('vascular')].sort_values('patient')
assert len(v)==9 and int(v.two_links_compatible.sum())==7 and int(v.all_three_compatible.sum())==5
para('SPP1相关巨噬状态与基质组织及T细胞分布','Title')
para('同患者同区域的探索性关联分析','Subtitle')
para('项目26210    2026年9月30日    修订版','Note')
h('摘要')
para('本报告承接TYMP与SPP1相关髓系状态的研究背景，直接检验SPP1相关巨噬状态、同一个ECM组织程序和T细胞相对肿瘤内分布，是否在同一批患者及同一空间区域中呈现相接的方向。分析固定TJ-RCC原有9位患者，使用200 μm区域及既定基础、血管调整模型。')
para('血管调整后，SPP1相关巨噬状态与外侧ECM组织程序在7/9位患者中呈正向关联，中位ρ为0.271；ECM组织程序与T细胞内侧减外侧信号差在9/9位患者中呈负向关联，中位ρ为−0.144。两段方向同时符合“SPP1状态较高—ECM组织较高—T细胞相对内侧较低”的患者为7/9位。')
para('在这7位两段方向相接的患者中，5位同时表现为SPP1状态与T细胞相对内侧信号负向关联。基础模型中，两段方向相接为5/9位，三项方向同时符合为3/9位。因此，现有数据提供同患者、同区域的中段连接线索，其证据强度仍为探索性。')
h('连接对象与对应结果')
table(['连接环节','血管调整结果','同患者方向'],[
['SPP1相关巨噬状态与ECM组织','中位ρ=0.2714','7/9正向'],
['同一ECM组织与T细胞内外差','中位ρ=−0.1441','9/9负向'],
['SPP1状态与T细胞内外差','中位ρ=−0.0818','6/9负向'],
['前两段同时符合预期方向','逐患者核对两条关系','7/9符合'],
['两段及端到端关系均符合','同时核对三条关系','5/9符合']], [3.8,3.2,2.0])
note('9/9只指ECM组织与T细胞分布这一条末端关系；整条两段关联为7/9方向相接。方向相接不等于因果中介，9例不进行患者层面推断检验。')
h('同一对象的定义')
para('M代表巨噬细胞中SPP1相关空间状态的映射占比；E代表外侧ECM组织程序；T代表总T细胞映射信号的肿瘤内外差。本报告的核心判断始终使用相同的M、E和T，不以其他队列的间质活化或胶原成熟结果替代E。','Note')
newpage();h('1 方法')
h2('固定患者与空间区域')
para('固定此前ECM组织—T细胞分析中全部9位可评价TJ-RCC患者，不新增或排除患者。主要尺度为200 μm，同时报告基础模型和血管调整模型。每位患者在每个模型内，三条关系采用完全相同的区域行、协变量矩阵和空间邻接结构。')
para('患者为R114、R29、R51、X49、X98、Y12、Y27、Y7和Z43，每人对应一张切片。原计算肿瘤边界、几何稀疏化及区域资格保持固定。基础与血管模型可具有不同的有效区域集合，但每个模型内部三条关系的区域严格一致。')
h2('三个读出的固定定义')
para('M沿用原SPP1相关巨噬空间状态占比，在冻结的外侧巨噬富集区域内汇总。它是映射状态读出，区别于SPP1单基因表达、蛋白量及单细胞C4身份。E采用同一外侧ECM组织程序，基因为FN1、COL1A1、COL3A1、LOX、LOXL2、PLOD2，不重新选基因。')
para('T采用总T细胞映射丰度后验第5百分位估计值（q05），先计算内外区域均值，再取log1p转换后的内侧减外侧差。M—E正向与E—T负向表示两条条件关联方向符合候选模型；另记录M—T是否负向。')
h2('共同调整与空间检验')
para('基础模型共同控制内外测序深度、巨噬与恶性细胞比例、两侧区域数量比、坐标、连通分量及外侧T细胞供给；血管模型再加入内外血管和壁细胞比例。三条关系均对相同协变量作秩残差化，再计算残差相关。')
para('每条关系进行999次受协变量约束的谱随机化。每位患者、每个模型的M—E、E—T及M—T三项构成一个BH校正家族；原始P和FDR均保留。既有末端专项分析的4项家族结果不覆盖，两套结果的效应量逐项核对，P值的随机种子和家族分别记录。')
para('9位患者仅报告中位效应、方向人数及逐患者对照，不以区域数替代患者数，不计算患者总体P或中介效应。两段一致定义为同患者M—E>0且E—T<0；三关系一致再要求M—T<0。全部9位患者均保留。')
h2('背景结果与原专项结果')
para('免疫功能分析沿用已冻结的12项患者关联，完整列于附表S1。此前三个空间队列的ECM组织及重塑与T/CD8分布结果完整保留于附表S4至S6。这些结果为背景及范围对照，不与当前9例合并为新增样本或独立验证。')
newpage();h('2 结果')
h2('SPP1相关巨噬状态与同一区域ECM组织呈总体正向联系')
para('在固定9位患者中，M—E关系的基础模型和血管调整模型均为7例正向。中位ρ分别为0.2451和0.2714。血管调整后，R114、R29、R51、X49、X98、Y27和Y7为正向，Y12和Z43为负向。')
para('这一步直接建立了SPP1相关空间状态与末端分析中同一个ECM组织读出的对应，未使用其他队列或其他基质程序替代。')
h2('同一ECM组织与T细胞相对内侧信号呈负向联系')
para('E—T基础模型为7/9例负向，中位ρ=−0.1920；血管调整后为9/9例负向，中位ρ=−0.1441。复算效应与原末端专项分析一致。')
h2('七位患者同时具备两段方向对应')
para('血管调整下，7/9位患者同时具有M—E正向及E—T负向；基础模型为5/9位。在血管调整的7位两段一致者中，R114、R29、R51、X49、X98这5位的M—T也为负向。Y27及Y7的M—T为正向，因此没有将其写成三项关系均一致。')
para('表1 同患者同区域的三条关联','Caption')
table(['患者','M—E\nρ','E—T\nρ','M—T\nρ','两段\n一致','三关系\n一致'],[[r.patient.replace('TJRCC:',''),val(r.M_E,4),val(r.E_T,4),val(r.M_T,4),'是' if r.two_links_compatible else '否','是' if r.all_three_compatible else '否'] for r in v.itertuples()],[1.1,1.4,1.4,1.4,1.2,1.3])
note('表1为200 μm血管调整模型。M—E为SPP1相关巨噬状态与ECM组织；E—T为ECM组织与T细胞内外差；M—T为SPP1状态与T细胞内外差。两段一致：前者正向、后者负向。三关系一致：在两段一致基础上，M—T也负向。')
para('同患者方向核对将末端的9/9负向定位到完整两段模型：其中7/9具有前后两段对应，5/9进一步具有端到端对应。这构成局部状态与基质及免疫分布相联系的探索性证据。')
newpage();h('3 讨论')
h2('同一基质读出提供了中间连接')
para('本轮将SPP1相关巨噬状态和T细胞分布共同定位到外侧ECM组织程序。血管调整后，7位患者同时具有SPP1状态较高对应ECM组织较高，以及同一ECM组织较高对应T细胞相对内侧较低的方向。与分别列出髓系、基质和免疫结果相比，该对照明确了哪些患者在同一空间框架下呈现相接的关联。')
para('基质影响T细胞定位具有原始实验依据：人肺肿瘤活切片研究显示纤维密度和取向影响T细胞运动，胶原酶处理可增加T细胞与癌细胞接触[1]；DDR1研究将纤维排列与免疫排斥联系起来[2]。其他癌种也报道了SPP1相关巨噬与成纤维细胞的空间协作[3]。这些机制背景帮助解释当前同区域方向对应，不能替代本项目干预证据。')
h2('方向相接与统计证据强度')
para('本轮是在观察既有ECM—T结果并发现连接缺口后开展的定向探索，所用患者不是独立验证集。两段一致由两个条件相关的符号共同定义，并非因果中介；三条关系共同调整外侧T供给也可能引入条件依赖。基础模型为5/9两段一致，血管调整为7/9，说明对应程度依赖模型。')
para('新检验中仅个别切片的M—E或M—T达到原始P<0.05。基础模型M—E有1例通过三项家族FDR；血管调整M—E的最低P为0.039、对应FDR为0.0585，属名义信号。两模型E—T均无切片原始P<0.05。因此主要价值在患者方向与同区域对应，不表述为完整链条已获统计或因果确认。')
h2('原免疫及跨队列结果限定外推范围')
para('此前12项免疫功能研究调整关联均未达到名义标准，SPP1受体检出不能直接证明免疫抑制。GSE175540也未在200 μm主模型中重复ECM组织—T细胞相对内侧降低；400 μm基础模型中ECM组织—CD8为反向名义正关联（P=0.04248，FDR=0.16992）。这些结果及GLMF案例差异完整保留。')
para('M为空间映射状态，E为转录程序，T为静态丰度代理。当前没有实测纤维结构、实时迁移或同病例干预。SPP1状态不能直接改称C4或SPP1蛋白，ECM组织也不能与间质活化、胶原成熟及ECM重塑互换。计算边界、spot混合成分和空间自相关仍限制解释。')
h2('综合结论')
para('现有数据探索性提示，在TJ-RCC部分患者中，SPP1相关巨噬状态与外侧ECM组织呈正向联系，同一ECM组织又与T细胞相对肿瘤内不足相联系。血管调整后7/9位患者具备两段方向对应，其中5/9同时具备端到端对应。该结果为SPP1相关髓系—基质—免疫空间分布提供了明确的同对象连接线索。')
newpage();h('参考文献')
for t in [
'[1] Salmon H, et al. Matrix architecture defines the preferential localization and migration of T cells into the stroma of human lung tumors. J Clin Invest. 2012;122:899–910. DOI 10.1172/JCI45817. https://www.jci.org/articles/view/45817',
'[2] Sun X, et al. Tumour DDR1 promotes collagen fibre alignment to instigate immune exclusion. Nature. 2021;599:673–678. DOI 10.1038/s41586-021-04057-2.',
'[3] Qi J, et al. Single-cell and spatial analysis reveal interaction of FAP+ fibroblasts and SPP1+ macrophages in colorectal cancer. Nat Commun. 2022;13:1742. DOI 10.1038/s41467-022-29366-6.']:
 para(t,'Reference')
h('附表S1 免疫功能分析背景')
pred={'mature_SPP1_pseudobulk':'整体巨噬SPP1','C4_fraction':'C4比例','C4_SPP1_pseudobulk':'C4内SPP1'}
ends={'CD8_cytotoxicity':'CD8细胞毒','CD8_dysfunction':'CD8功能障碍','Treg_CD4':'Treg/CD4','NK_cytotoxicity':'NK细胞毒'}
table(['预测量','免疫端点','n','调整ρ','原始P','FDR'],[[pred[r.predictor],ends[r.endpoint],str(r.n_patients),val(r.rho,4),val(r.p_value),val(r.FDR)] for r in im.itertuples()],[2,2,.5,1,1.2,1.2])
note('研究调整的患者级关联；各假设继承原检验家族。整体SPP1的CD8两项、C4的CD8四项、Treg两项及NK两项分别BH；整体SPP1的Treg及NK各为单独预设终点，不作额外FDR校正。此表不是同区域连接的替代证据。')
for idx,model in enumerate(['core','vascular'],start=2):
 newpage();h(f'附表S{idx} 同区域连接的'+('基础' if model=='core' else '血管调整')+'模型')
 note('固定200 μm、同9位患者。M—E为SPP1状态与ECM组织，E—T为ECM组织与T细胞内外差，M—T为端到端关系。')
 dd=edge[edge.model.eq(model)].sort_values(['patient','edge'])
 lab={'M_E':'M—E','E_T':'E—T','M_T':'M—T'}
 table(['患者','关系','ρ','空间模型P','空间模型FDR'],[[r.patient.replace('TJRCC:',''),lab[r.edge],val(r.rho,4),val(r.spatial_P,3),val(r.spatial_FDR,4)] for r in dd.itertuples()],[1.2,.9,1.4,1.7,1.9])
 note('每位患者、每模型的三项关系为一个BH家族。表内P/FDR是切片空间模型统计，不能替代患者总体P。末端专项分析采用不同随机种子及4项家族，其结果原样保留于附表S4至S6；效应量与本轮逐项一致。')
for idx,(ds,name) in enumerate([('TJRCC','TJ RCC'),('GSE175540','GSE175540'),('GLMF_GSE210041','GLMF GSE210041')],start=4):
 newpage();h(f'附表S{idx} {name}原末端专项结果')
 note('原ECM组织及重塑与T/CD8分布的完整比较。每行以独立患者汇总，保留原始P/FDR和检验家族，不用同区域连接结果覆盖原检验。')
 dd=co[co.dataset.eq(ds)].sort_values(['radius','model','program','endpoint'])
 table(['尺度\nμm','模型','ECM\n程序','免疫\n读出','n','负向\n人数','中位ρ','原始P','FDR'],[[str(r.radius),'基础' if r.model=='core' else '血管','组织' if r.program=='ECM_organization' else '重塑',r.endpoint,str(r.n),f'{r.negative}/{r.n}',val(r.median_rho,4),val(r.P),val(r.FDR)] for r in dd.itertuples()],[.9,.9,.9,.8,.5,.8,1.3,1.3,1.3])
 note('200 μm基础模型为主要分析，其他为敏感性分析。患者n<10不计算推断P或FDR；n≥10用双侧Wilcoxon检验，BH在同队列、尺度和模型的4项比较内计算。不同尺度和模型不作为独立验证。')
 if ds=='GSE175540':note('400 μm基础模型的ECM组织—CD8为正向名义信号，与排斥假设方向相反；血管调整后不保留。')
 if ds=='GLMF_GSE210041':note('GLMF两例对应GSE210041的Patient1及Patient7；400 μm仅1位满足模型资格，不重复计为新患者或队列。')
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
dest=P/'delivery/26210_SPP1基质与T细胞空间关联_客户修订版.docx'
with ZipFile(dest,'w',ZIP_DEFLATED) as zip:
 for fn,content in {'[Content_Types].xml':ct,'_rels/.rels':rootrels,'word/document.xml':document,'word/styles.xml':styles,'word/header1.xml':header,'word/footer1.xml':footer,'word/_rels/document.xml.rels':rels}.items():zip.writestr(fn,'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'+content)
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in src.items())
(P/'qc/BUILD.json').write_text(json.dumps(dict(status='BUILT_AWAITING_RENDER',document=str(dest),source_unchanged=True,intended_pages=10,body_elements=len(audit),tables=sum(x['type']=='table' for x in audit)),indent=2))
print('BUILT',dest,flush=True)
