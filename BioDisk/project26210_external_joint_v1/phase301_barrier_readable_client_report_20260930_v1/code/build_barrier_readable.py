from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from xml.sax.saxutils import escape
import pandas as pd,numpy as np,json,hashlib,subprocess
P=Path(__file__).resolve().parents[1];R=P.parent
O=R/'phase300_same_region_SPP1_ECM_T_link_20260930_v1'
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
def table(headers,rows,weights,size=20,compact=False):
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
 
 if compact:tbl=tbl.replace('w:w="75"', 'w:w="55"')
 tbl+='</w:tbl>';parts.append(tbl);parts.append('<w:p><w:pPr><w:spacing w:after="40" w:line="100" w:lineRule="auto"/></w:pPr></w:p>')
 audit.append(dict(type='table',headers=headers,rows=rows))
def val(x,d=5):return '—' if pd.isna(x) else f'{x:.{d}f}'

edge=pd.read_csv(O/'tables/SAME_REGION_EDGES.tsv',sep='\t')
align=pd.read_csv(O/'tables/PATIENT_LINK_ALIGNMENT.tsv',sep='\t')
summ=pd.read_csv(O/'tables/ALIGNMENT_SUMMARY.tsv',sep='\t')
for path in [O/'tables/SAME_REGION_EDGES.tsv',O/'tables/PATIENT_LINK_ALIGNMENT.tsv',O/'tables/ALIGNMENT_SUMMARY.tsv']:
 src[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
(P/'inputs/SOURCE_MANIFEST.json').write_text(json.dumps(src,indent=2))
v=align[align.model.eq('vascular')].sort_values('patient')
assert len(v)==9 and int(v.two_links_compatible.sum())==7 and int(v.all_three_compatible.sum())==5

para('SPP1相关巨噬与基质免疫屏障的空间证据','Title')
para('肾癌组织中基质活动与T细胞分布的探索性分析','Subtitle')
para('项目26210    客户修订版','Note')
h('摘要')
para('现有空间数据探索性提示，部分肾癌患者呈现与基质免疫屏障相符的空间联系：肿瘤外侧SPP1相关巨噬状态占比较高的区域，基质组织程序也较高；同一基质程序较高的区域，肿瘤内侧相对外侧的T细胞信号差较低。血管调整后，9位患者中7位同时符合这两段方向，其中5位还表现为SPP1相关巨噬状态与T细胞内外分布差负向关联。')
h('免疫屏障在组织中对应什么问题')
para('这里关注的免疫屏障，是肿瘤周围的细胞和基质环境限制T细胞接近肿瘤细胞的过程。它需要回答一个空间问题：T细胞即使出现在肿瘤外侧，是否仍难以进入肿瘤细胞所在区域。因而，基质与肿瘤内外T细胞分布的联系，比单独观察T细胞总量更贴近本研究问题。')
h('这三项结果为何与屏障有关')
table(['屏障问题','本研究观察','对应意义'],[
['哪些细胞与基质活动相联系','SPP1相关巨噬状态与外侧基质组织程序在7/9例中正向关联','定位与候选屏障环境相联系的巨噬状态'],
['这种基质活动是否伴随T细胞分布变化','基质组织程序与T细胞内外差在血管调整后9/9例中负向关联','基质活动较高对应较低的内侧相对信号'],
['上述联系是否出现在同一患者','7/9例同时符合前两段方向  5/9例同时符合三条关系','将巨噬状态 基质和T细胞分布连接到同一对象']
],[2.4,4.4,3.2],size=20)
para('这些结果支持“SPP1相关巨噬与基质活动相联系，而基质活动与T细胞向肿瘤内侧的分布不足相联系”的探索性解释。证据来自静态空间关联，尚未直接证明T细胞受到物理阻挡或迁移能力下降。','Note')

newpage();h('1 方法')
h2('1.1 围绕肿瘤边界比较外侧基质与内外T细胞分布')
para('分析沿用TJ-RCC全部9位可评价患者，每位对应一张切片，编号为R114、R29、R51、X49、X98、Y12、Y27、Y7和Z43。主要尺度为200 μm。沿用原计算肿瘤边界、区域资格和几何稀疏化，在同一空间框架中汇总肿瘤外侧巨噬状态、外侧基质组织程序，以及边界两侧的T细胞信号。')
h2('1.2 用三个固定读出描述候选屏障环境')
para('SPP1相关巨噬状态：用既有空间映射得到的SPP1相关状态占比，在已冻结的外侧巨噬富集区域汇总。该指标描述巨噬细胞的状态构成。')
para('基质组织程序：用FN1、COL1A1、COL3A1、LOX、LOXL2和PLOD2组成的原固定程序，描述外侧区域中与基质形成、组织及加工有关的转录活动。ECM是细胞外基质的英文缩写。')
para('T细胞内外分布差：对总T细胞空间映射丰度后验第5百分位估计值q05，分别计算内外区域均值，取log1p转换后的内侧减外侧差。这个差值越低，表示内侧相对外侧的T细胞信号越低。')
h2('1.3 在同一区域检验三条关系')
para('分别检验巨噬状态与基质、基质与T细胞内外差、巨噬状态与T细胞内外差。每位患者在同一模型下，三条关系采用完全相同的区域行、协变量和空间邻接结构。基础模型控制两侧测序深度、巨噬及恶性细胞比例、两侧区域数量比、坐标、连通分量及外侧T细胞供给；血管调整模型再加入两侧血管和壁细胞比例。基础与血管模型分别保留。')
h2('1.4 以独立患者汇总方向和效应')
para('三条关系经相同协变量秩残差化后计算相关系数ρ，并进行999次受协变量约束的谱随机化。每位患者每个模型的三项检验分别构成一个BH校正家族，保留原始P和FDR，共18个家族。9位患者汇总中位ρ、方向人数及逐例关系，不作患者总体推断检验。区域数不计为患者数。原12项免疫功能结果及三个空间队列的完整比较保留于附表S1至S6；原专项四项校正家族与本次三项家族分别记录。')

newpage();h('2 结果')
h2('2.1 七位患者的SPP1相关巨噬状态与外侧基质活动同向变化')
para('9位患者中，7位的外侧SPP1相关巨噬状态占比与同一区域的基质组织程序呈正向关联。基础模型中位ρ为0.2451，血管调整后为0.2714。这表明，在这些患者各自的组织切片内，SPP1相关巨噬状态较集中的外侧区域，也倾向于具有较高的基质组织相关转录活动。')
para('这一结果把SPP1相关巨噬状态与候选屏障的基质环境联系起来：研究对象不只是一类巨噬细胞，还包括与其共同变化的肿瘤周围基质活动。')
h2('2.2 基质活动较高对应较低的T细胞内侧相对信号')
para('外侧基质组织程序与T细胞内侧减外侧信号差，在基础模型中7/9例为负向，中位ρ为−0.1920；血管调整后9/9例均为负向，中位ρ为−0.1441。')
para('负向关联的具体含义是：同一患者中，基质组织程序较高的区域，T细胞在肿瘤内侧相对外侧的信号差较低。这个观察直接对应屏障问题中的“内外分布”，提示较高的外侧基质活动伴随T细胞向肿瘤内侧分布不足的方向。血管调整结果说明，这种方向在纳入血管和壁细胞比例后仍保留。')
h2('2.3 巨噬与基质的联系和基质与T细胞的联系出现在同一批病例')
para('血管调整后，R114、R29、R51、X49、X98、Y27和Y7这7位患者，同时表现为“SPP1相关巨噬状态与基质活动正向关联”以及“同一基质活动与T细胞内外差负向关联”。其中前5位还表现为SPP1相关巨噬状态与T细胞内外差负向关联。')
para('因此，现有数据将候选巨噬状态、候选基质环境与T细胞分布变化连接到了同一患者和同一区域。9/9描述的是基质与T细胞分布这一项关系；7/9描述前两项关系方向同时符合；5/9描述三项关系方向同时符合。')

newpage();h('2.4 逐例对照显示哪些患者符合屏障相关方向')
para('表1逐例列出同一切片中的三条关系。若前两项同时符合，表示该患者既有SPP1相关巨噬与基质的联系，也有同一基质与较低T细胞内外差的联系；第三项再为负向，则巨噬状态也直接对应较低的T细胞内外差。')
para('表1 血管调整后SPP1相关巨噬 基质与T细胞分布的逐例关系','Caption')
table(['患者','巨噬状态与\n基质活动\nρ','基质活动与\nT细胞内外差\nρ','巨噬状态与\nT细胞内外差\nρ','前两项方向\n同时符合','三项方向\n同时符合'],[
[r.patient.replace('TJRCC:',''),val(r.M_E,4),val(r.E_T,4),val(r.M_T,4),'是' if r.two_links_compatible else '否','是' if r.all_three_compatible else '否'] for r in v.itertuples()
],[.85,1.65,1.7,1.7,1.3,1.3])
note('巨噬状态指SPP1相关巨噬空间状态占比；基质活动指固定ECM组织程序；T细胞内外差为肿瘤内侧减外侧信号。前两项方向同时符合要求第一项ρ>0且第二项ρ<0；三项方向同时符合再要求第三项ρ<0。是或否仅表示方向，不表示统计显著。全部9位患者保留。')
h2('基础模型保留同样的联系方向  对应病例数有所不同')
para('基础模型中，5/9位患者符合前两项方向，3/9位符合全部三项。两模型的巨噬状态与基质活动均为7/9例正向；基质与T细胞内外差则由基础模型7/9例负向，变为血管调整后9/9例负向。两模型均指向巨噬相关基质活动与T细胞分布的联系，其病例一致程度受调整条件影响。')

newpage();h('3 讨论')
h2('3.1 这些空间联系为何支持基质免疫屏障这一研究方向')
para('本研究得到的具体线索是：SPP1相关巨噬状态与肿瘤外侧基质活动相联系，而同一基质活动又与较低的T细胞内侧相对信号相联系。它为“巨噬相关的基质环境参与T细胞空间排斥”提供了探索性支持，并将研究问题从免疫细胞功能评分推进到细胞与基质所在的组织位置。')
para('这一解释具有已发表机制背景。人肺肿瘤活切片实验观察到，肿瘤周围基质结构影响T细胞运动，胶原酶处理可增加T细胞与肿瘤细胞接触[1]；DDR1研究将胶原纤维排列与免疫排斥相联系[2]；结直肠癌研究报道了SPP1巨噬与FAP成纤维细胞的空间协作[3]。这些文献解释为什么基质与T细胞位置值得共同研究，项目中的联系仍按本身的测量和统计证据判断。')
h2('3.2 TYMP与这部分屏障结果怎样衔接')
para('整体研究以TYMP相关代谢背景及SPP1相关巨噬发现为起点。本部分补充其下游的空间问题：SPP1相关巨噬状态是否与基质环境相联系，以及同一基质环境是否对应T细胞内外分布变化。这样可以将TYMP研究背景与免疫屏障研究问题衔接，但这9例空间分析没有建立TYMP干预经SPP1和基质改变而影响T细胞进入的因果链。髓系TYMP相关结果也不能替代肿瘤来源TYMP的作用。')
h2('3.3 当前证据支持到哪一步')
para('当前支持的是屏障相关的空间关联线索。9例只作方向与效应描述；本次为已观察末端结果后的定向探索。严格空间检验中，两模型的基质与T细胞内外差均没有切片达到P<0.05。巨噬状态与基质关系仅个别切片有统计支持，血管模型R29为P=0.039、FDR=0.0585，属于名义信号；基础模型Y27为P=0.010、FDR=0.030。完整54项空间模型见附表S2和S3。')
para('现有读出分别是映射状态、基质转录程序及静态T细胞信号，未实测胶原厚度、纤维排列、SPP1蛋白或T细胞迁移。较低的内外差不等于所有内侧区域都缺乏T细胞，也不等于T细胞总数下降。原计算边界、混合spot、空间自相关及共同调整外侧T供给形成的条件依赖，限制了功能解释。')
para('此前12项免疫功能患者关联未达到研究调整后的名义标准。GSE175540未在200 μm主模型中重复末端方向，400 μm基础模型的基质组织与CD8呈反向名义正关联，P=0.04248、FDR=0.16992；GLMF病例差异也完整保留。因此，现有结果尚不能确认为稳定跨队列的功能屏障。')
newpage();h('参考文献')
for t in [
'[1] Salmon H, et al. Matrix architecture defines the preferential localization and migration of T cells into the stroma of human lung tumors. J Clin Invest. 2012;122:899–910. DOI 10.1172/JCI45817. https://www.jci.org/articles/view/45817',
'[2] Sun X, et al. Tumour DDR1 promotes collagen fibre alignment to instigate immune exclusion. Nature. 2021;599:673–678. DOI 10.1038/s41586-021-04057-2.',
'[3] Qi J, et al. Single-cell and spatial analysis reveal interaction of FAP+ fibroblasts and SPP1+ macrophages in colorectal cancer. Nat Commun. 2022;13:1742. DOI 10.1038/s41467-022-29366-6.']:
 para(t,'Reference')
h('附表S1 SPP1相关状态与免疫细胞功能的患者级结果')
pred={'mature_SPP1_pseudobulk':'整体巨噬SPP1','C4_fraction':'C4比例','C4_SPP1_pseudobulk':'C4内SPP1'}
ends={'CD8_cytotoxicity':'CD8细胞毒','CD8_dysfunction':'CD8功能障碍','Treg_CD4':'Treg/CD4','NK_cytotoxicity':'NK细胞毒'}
table(['预测量','免疫端点','n','调整ρ','原始P','FDR'],[[pred[r.predictor],ends[r.endpoint],str(r.n_patients),val(r.rho,4),val(r.p_value),val(r.FDR)] for r in im.itertuples()],[2,2,.5,1,1.2,1.2])
note('研究调整的患者级关联；各假设继承原检验家族。整体SPP1的CD8两项、C4的CD8四项、Treg两项及NK两项分别BH；整体SPP1的Treg及NK各为单独预设终点，不作额外FDR校正。此表不是同区域连接的替代证据。')
for idx,model in enumerate(['core','vascular'],start=2):
 newpage();h(f'附表S{idx} 同区域连接的'+('基础' if model=='core' else '血管调整')+'模型')
 note('固定200 μm和同9位患者。M表示SPP1相关巨噬状态，E表示外侧基质组织程序，T表示肿瘤内侧减外侧的T细胞信号差。M—E检查巨噬状态与基质联系，E—T检查基质与T细胞分布联系，M—T检查巨噬状态与T细胞分布联系。ρ为正表示同向变化，为负表示反向变化。')
 dd=edge[edge.model.eq(model)].sort_values(['patient','edge'])
 lab={'M_E':'M—E','E_T':'E—T','M_T':'M—T'}
 table(['患者','关系','ρ','空间模型P','空间模型FDR'],[[r.patient.replace('TJRCC:',''),lab[r.edge],val(r.rho,4),val(r.spatial_P,3),val(r.spatial_FDR,4)] for r in dd.itertuples()],[1.2,.9,1.4,1.7,1.9],compact=True)
 note('每位患者、每模型的三项关系为一个BH家族。表内P/FDR是切片空间模型统计，不能替代患者总体P。末端专项分析采用不同随机种子及4项家族，其结果原样保留于附表S4至S6；效应量与本轮逐项一致。')
for idx,(ds,name) in enumerate([('TJRCC','TJ RCC'),('GSE175540','GSE175540'),('GLMF_GSE210041','GLMF GSE210041')],start=4):
 newpage();h(f'附表S{idx} {name}原末端专项结果')
 note('原ECM组织及重塑与T/CD8分布的完整比较。每行汇总一个队列中的独立患者，保留原始P/FDR和检验家族，不用同区域连接结果覆盖原检验。')
 dd=co[co.dataset.eq(ds)].sort_values(['radius','model','program','endpoint'])
 table(['尺度\nμm','模型','ECM\n程序','免疫\n读出','n','负向\n人数','中位ρ','原始P','FDR'],[[str(r.radius),'基础' if r.model=='core' else '血管','组织' if r.program=='ECM_organization' else '重塑',r.endpoint,str(r.n),f'{r.negative}/{r.n}',val(r.median_rho,4),val(r.P),val(r.FDR)] for r in dd.itertuples()],[.9,.9,.9,.8,.5,.8,1.3,1.3,1.3])
 note('200 μm基础模型为主要分析，同时保留400 μm尺度及血管调整的比较。患者n<10不计算推断P或FDR；n≥10用双侧Wilcoxon检验，BH在同队列、尺度和模型的4项比较内计算。不同尺度和模型不作为独立验证。')
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
dest=P/'delivery/26210_SPP1相关基质免疫屏障_客户易读修订版.docx'
with ZipFile(dest,'w',ZIP_DEFLATED) as zip:
 for fn,content in {'[Content_Types].xml':ct,'_rels/.rels':rootrels,'word/document.xml':document,'word/styles.xml':styles,'word/header1.xml':header,'word/footer1.xml':footer,'word/_rels/document.xml.rels':rels}.items():zip.writestr(fn,'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'+content)
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in src.items())
(P/'qc/BUILD.json').write_text(json.dumps(dict(status='BUILT_AWAITING_RENDER',document=str(dest),source_unchanged=True,intended_pages=11,body_elements=len(audit),tables=sum(x['type']=='table' for x in audit)),indent=2))
print('BUILT',dest,flush=True)
