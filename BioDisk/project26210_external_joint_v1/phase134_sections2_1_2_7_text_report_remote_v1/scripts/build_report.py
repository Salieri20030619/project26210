from pathlib import Path
import json,csv,hashlib,html,datetime,sys,re
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,KeepTogether,CondPageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER,TA_JUSTIFY,TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

w=Path(__file__).resolve().parents[1]
E=json.loads((w/'qc/EVIDENCE.json').read_text())
blocks=json.loads((w/'report/report_content.json').read_text())
pdfmetrics.registerFont(TTFont('Song',str(w/'fonts/simsun.ttc'),subfontIndex=0))
pdfmetrics.registerFont(TTFont('Hei',str(w/'fonts/simhei.ttf')))
pdfmetrics.registerFontFamily('Song',normal='Song',bold='Hei',italic='Song',boldItalic='Hei')
S={}
S['title']=ParagraphStyle('title',fontName='Hei',fontSize=16,leading=23,alignment=TA_CENTER,spaceAfter=9,keepWithNext=True,textColor=colors.black)
S['date']=ParagraphStyle('date',fontName='Song',fontSize=10,leading=15,alignment=TA_CENTER,spaceAfter=15,keepWithNext=True)
S['h1']=ParagraphStyle('h1',fontName='Hei',fontSize=13,leading=19,spaceBefore=12,spaceAfter=6,keepWithNext=True)
S['h2']=ParagraphStyle('h2',fontName='Hei',fontSize=11.2,leading=16.5,spaceBefore=9,spaceAfter=4,keepWithNext=True)
S['p']=ParagraphStyle('body',fontName='Song',fontSize=10.4,leading=15.8,alignment=TA_JUSTIFY,firstLineIndent=20.8,spaceAfter=5,wordWrap='CJK',allowWidows=0,allowOrphans=0)
S['caption']=ParagraphStyle('caption',fontName='Hei',fontSize=9.7,leading=14,spaceBefore=7,spaceAfter=5,keepWithNext=False)
S['note']=ParagraphStyle('note',fontName='Song',fontSize=8.6,leading=12.3,spaceAfter=7,wordWrap='CJK')
S['cell']=ParagraphStyle('cell',fontName='Song',fontSize=8.3,leading=11.8,wordWrap='CJK')
S['cellc']=ParagraphStyle('cellc',parent=S['cell'],alignment=TA_CENTER)
S['cellh']=ParagraphStyle('cellh',parent=S['cell'],fontName='Hei')
S['cellhc']=ParagraphStyle('cellhc',parent=S['cellh'],alignment=TA_CENTER)
story=[];material=[];table_registry=[]
WIDTH=A4[0]-90

def clean(t):
    t=str(t).replace('−','-').replace('–','-').replace('—','-').replace('‑','-')
    supers=str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹','-0123456789')
    return re.sub(r'×10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)',lambda m:'e'+m.group(1).translate(supers),t)
def para(t,style='p'):
    t=clean(t)
    latin='Times-Bold' if style in ['title','h1','h2','caption','cellh','cellhc'] else 'Times-Roman'
    parts=re.split(r'([\x20-\x7e]+)',t)
    rich=''.join('<font name="'+latin+'">'+html.escape(x)+'</font>' if x and all(32<=ord(c)<=126 for c in x) else html.escape(x).replace('\n','<br/>') for x in parts)
    return Paragraph(rich,S[style])
def add(t,style='p'):
    story.append(para(t,style));material.append({'type':style,'text':clean(t)})
def num(v,digits=5):
    if v is None or str(v) in ['','NA','nan','None']:return '—'
    v=float(v)
    if v==0:return '0'
    if abs(v)<0.0001:return format(v,'.3g')
    return format(v,'.'+str(digits)+'f').rstrip('0').rstrip('.')
def sig(p,q=None):
    if q not in [None,'','NA'] and float(q)<.05:return 'FDR显著'
    if p not in [None,'','NA'] and float(p)<.05:return '名义显著'
    return '未显著'
def direction(r,pk='n_positive',nk='n_negative',zk='n_zero'):
    return '/'.join(str(r.get(k,'—'))for k in [pk,nk,zk])
def T(caption,headers,rows,widths=None,note='',leftcols=(0,),small=False):
    story.append(CondPageBreak(75))
    add(caption,'caption')
    if widths is None:widths=[WIDTH/len(headers)]*len(headers)
    else:widths=[WIDTH*x/sum(widths)for x in widths]
    data=[[para(x,'cellh'if j in leftcols else'cellhc')for j,x in enumerate(headers)]]
    data += [[para(x,'cell'if j in leftcols else'cellc')for j,x in enumerate(row)]for row in rows]
    table=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
    style=[('VALIGN',(0,0),(-1,-1),'MIDDLE'),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dfeaf0')),
           ('LINEABOVE',(0,0),(-1,0),.6,colors.black),('LINEBELOW',(0,0),(-1,0),.45,colors.black),
           ('LINEBELOW',(0,-1),(-1,-1),.6,colors.black),('LEFTPADDING',(0,0),(-1,-1),5),
           ('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]
    for i in range(1,len(data)):
        if i%2==0:style.append(('BACKGROUND',(0,i),(-1,i),colors.HexColor('#f2f5f6')))
    table.setStyle(TableStyle(style));story.append(table);story.append(Spacer(1,6))
    material.append({'type':'table','caption':caption,'headers':headers,'rows':rows})
    table_registry.append({'caption':caption,'rows':len(rows),'columns':len(headers)})
    if note:add(note,'note')
def h(t):
    story.append(CondPageBreak(200 if t.startswith('4.2 ') else 105))
    add(t,'h2')

pathnames={'AUC_Pyrimidine':'嘧啶','AUC_CysMet':'Cys/Met','AUC_ArgPro':'Arg/Pro','AUC_Glycerophospholipid':'甘油磷脂','Pyrimidine_flux':'嘧啶通量','CysMet_flux':'Cys/Met通量','ArgPro_flux':'Arg/Pro通量','Choline_proxy':'胆碱局部代理','F1':'F1','TYMP_log1p_CP10K_v1':'TYMP RNA','infercnv_burden_over_threshold_v1':'CNV负荷','progression_anchor_v1':'起点锚定指标'}
methodnames={'Monocle3_pseudotime_patient_v1':'Monocle3','Slingshot_pseudotime_patient_v1':'Slingshot','scTour_pseudotime_patient_v1':'scTour','CytoTRACE2_Score':'CytoTRACE2'}
endpointnames={'continuous_F1_spearman':'连续F1','F1_Q4_minus_Q1_median_difference':'F1 Q4−Q1','F1_dominant_median_difference':'F1主导组','F1_continuous':'连续F1','F1_dominant_group':'F1高主克隆','F1_patient_quartile':'F1 Q4−Q1','trajectory_terminal_shift':'末端−早期'}
axnames={'TYMP_axis_primary':'TYMP严格','TYMP_axis_broad_sensitivity':'TYMP广义','RRM2B_axis_primary':'RRM2B严格','RRM2B_axis_broad_sensitivity':'RRM2B广义'}
def main_table(key,cap):
    if key=='cohorts':
        rows=[[x['study'],str(x['libraries']),str(x['patients']),str(x['malignant_patients']),format(x['malignant_cells'],',')]for x in E['cohorts']]
        rows.append(['合计','75','40','26','11,950'])
        T(cap,['队列','纳入文库','纳入个体','恶性细胞患者','恶性上皮细胞'],rows,[1.6,1,1,1.3,1.4],
          '文库数和个体数对应原冻结纳入范围；恶性细胞患者数对应保守 CNV 判定后的正式分析集。')
    elif key=='cnv':
        mapping={'malignant_consensus':'双方法恶性一致','malignant_infercnv_only':'仅inferCNV支持且另一方法不可用','non_malignant_consensus':'双方法非恶性一致','method_discordant':'方法冲突','indeterminate':'未定'}
        order=['malignant_consensus','malignant_infercnv_only','non_malignant_consensus','method_discordant','indeterminate']
        rows=[]
        for name in order:
            x=next(a for a in E['cnv_counts']if a['malignancy_status_v1']==name)
            rows.append([mapping[name],format(int(x['n_cells']),','),'是'if x['malignant_conservative_v1']=='True'else'否'])
        T(cap,['判定','细胞数','纳入恶性主分析'],rows,[3.8,1,1.8])
    elif key=='reg':
        T(cap,['算法','调控子','n','rho','正/负/零','P','FDR'],
          [[x['method'],x['tf'],x['n_eligible_patients'],num(x['median_patient_effect']),direction(x),num(x['wilcoxon_p_two_sided']),num(x['BH_FDR'])]for x in E['reg_cont']],
          [1.35,1,0.45,.85,1,.95,.95],leftcols=(0,1),
          note='P 为双侧 Wilcoxon；FDR 在 GENIE3 或 GRNBoost2 的完整连续调控子家族内计算。')
    elif key=='tymp_activity':
        rows=[]
        for target in ['CASP9','CHMP1A','PRSS8','PXN']:
            for model in ['unadjusted','technical_adjusted']:
                x=next(z for z in E['tymp_activity']if z['target']==target and z['model']==model)
                rows.append([target,'未调整'if model=='unadjusted'else'技术调整',x['n'],num(x['median_effect']),direction(x,'positive','negative','zero'),num(x['P']),num(x['FDR'])])
        T(cap,['目标','模型','n','rho','正/负/零','P','FDR'],rows,[1,1.1,.4,.9,.9,.95,.95],leftcols=(0,1),
          note='技术调整控制原始文库量和检出基因数的秩；两个模型分别在四蛋白家族内校正。')
    elif key=='auc_pair':
        rows=[]
        for x in E['auc_pair']:
            rows.append(['TYMP阳性'if x['contrast']=='TYMP-positive_vs_negative'else'F1高主克隆',pathnames[x['pathway']],x['n_patients'],num(x['median_patient_delta']),str(round(float(x['fraction_patients_positive'])*int(x['n_patients'])))+'/'+x['n_patients'],num(x['paired_wilcox_PValue']),num(x['FDR'])])
        T(cap,['比较','通路','n','Δ AUCell','正向/n','P','FDR'],rows,[1.25,1.1,.4,.95,.7,.85,.85],leftcols=(0,1),
          note='分别比较 TYMP 阳性−阴性、F1 高主克隆−其他主克隆；每个比较内四通路独立 BH。Cys/Met 为半胱氨酸和蛋氨酸，Arg/Pro 为精氨酸和脯氨酸。')
    elif key=='terminal_key':
        rows=[[x['method'],x['gene'],x['n'],num(x['median_effect']),direction(x,'positive','negative','zero'),num(x['P']),num(x['FDR'])]for x in E['terminal_enzyme']if float(x['P'])<.05]
        T(cap,['方法','基因','n','Δ RNA','正/负/零','P','FDR'],rows,[1.2,.8,.4,.9,.9,.9,.9],leftcols=(0,1),
          note='Δ RNA 为患者末端−早期 log2(CPM+1) 的中位数；四项为全部 39 端点中 P<0.05 的结果，均未通过 39 项 BH。')
    elif key=='wetlab':
        T(cap,['客户报告的敲低变化','生信对应','证据定位'],[
          ['Thymidine升高；Thymine下降','TYMP与F1、嘧啶AUCell及嘧啶相关功能代理同向','核苷相关状态对应'],
          ['Methionine和SAM下降','TYMP阳性及F1高主克隆的Cys/Met通路升高','通路层对应；未直接测SAM'],
          ['Spermidine下降','ODC1和SRM与F1正相关；Slingshot的SRM末端名义升高','基因与局部轨迹支持'],
          ['5′-MTA下降','Cys/Met和Arg/Pro通路、多胺相关酶提供邻近代谢背景','MTA未被显式建模'],
          ['客户图未给出其他端点方向','scTour的SMS和Slingshot的PCYT1A末端反向结果保留','不扩展为未测定实验变化']
        ],[1.6,2.5,1.5],leftcols=(0,1,2),
        note='湿实验列仅重述客户提供的方向；不据本表推导代谢物变化倍数或因果效应。')
    elif key=='parameters':
        T(cap,['分析','关键参数和统计单位'],[
          ['原始计数质控','UMI≥200；基因≥200；线粒体≤25%；文库内MAD过滤；Scrublet预期双细胞率0.05。'],
          ['整合','5队列、75文库、40个体；scVI 8,000 HVG、30潜在维度、31 epochs；15近邻。'],
          ['恶性上皮与克隆','11,950细胞、26患者；阈值0.25；主克隆≥20细胞；14位主克隆配对患者。'],
          ['cNMF与pseudobulk','2,500基因；rank 7；20次运行；14患者×2组；edgeR检验15,308基因。'],
          ['差异排序基因ORA','各方向按P排序top250；背景15,308基因；每方向201通路独立BH。'],
          ['pySCENIC','GENIE3 299、GRNBoost2 314调控子；连续家族按算法，分组家族按算法×程序。'],
          ['metaVIPER','KIRC/KIRP；rank；minsize=25；pleiotropy=FALSE；≥30细胞/患者；四MR蛋白各端点独立BH。'],
          ['代谢评分','AUCell前5%秩；四重点通路；85通路扩展独立BH；scFEA M168、100 epochs、四代理轴。'],
          ['酶RNA与轨迹','12酶F1关联独立BH；13基因×3拟时序末端pseudobulk，共39项BH；早/末20%，每臂≥5细胞。'],
          ['统计','患者为独立重复；保留原精确/近似Wilcoxon实现；小于10例按探索性；名义P与FDR分开。']
        ],[1.15,4.8],leftcols=(0,1))
def appendices():
    h('4.1 F1 与 TYMP 的比较及上调方向富集')
    T('附表 S1A F1 主克隆 TYMP 比较与阈值敏感性',['阈值','n','Δ RNA','正向/n','P','FDR'],
      [[x['distance_threshold'],x['n_patients'],num(x['median_patient_tymp_delta']),str(round(float(x['positive_patient_fraction'])*int(x['n_patients'])))+'/'+x['n_patients'],num(x['wilcoxon_two_sided_p']),num(x['wilcoxon_BH_FDR_within_threshold'])]for x in E['clone_tymp']if x['factor']=='F1'],
      [.6,.4,1,1,1,1],note='效应为 TYMP log1p CP10K 差值；FDR 为每一距离阈值下 7 程序家族。主阈值为 0.25，不能以敏感性阈值替换主结果。')
    T('附表 S1B F1 高主克隆上调排序基因集的全部 FDR 显著通路',['通路','重叠数','P','FDR'],
      [[x['pathway'].replace('Reactome_metabolism::','Reactome: ').replace('KEGG_metabolism::','KEGG: ').replace('Hallmark_v5.2::HALLMARK_','Hallmark: ').replace('_',' '),x['overlap_size'],num(x['PValue']),num(x['FDR'])]for x in E['ora']if float(x['FDR'])<.05],
      [4.1,.55,.9,.9],note='各方向按edgeR P值排序取250基因，背景15,308基因；上调17/201显著，下调0/201显著。原始完整402项结果随统计源表保留。')
    h('4.2 候选调控子的探索性分组结果')
    T('附表 S2 F1 主导细胞与其他细胞的调控子比较',['算法','调控子','n','Δ AUCell','正/负/零','P','FDR'],
      [[x['method'],x['tf'],x['n_eligible_patients'],num(x['median_patient_effect_internal']),direction(x),num(x['wilcoxon_p_two_sided']),num(x['BH_FDR'])]for x in E['reg_group']],
      [1.3,.8,.4,.95,.9,.9,.9],leftcols=(0,1),note='9位独立患者，全部按探索性结果解释；效应为患者内F1主导细胞减去其他细胞的AUCell差值中位数。')
    h('4.3 候选蛋白的完整 F1 端点')
    T('附表 S3 metaVIPER 活性与 F1 的三个端点',['蛋白','端点','n','中位效应','正/负/零','P','FDR'],
      [[x['protein'],endpointnames[x['endpoint']],x['n_eligible_patients'],num(x['median_patient_effect']),direction(x),num(x['wilcoxon_p_two_sided']),num(x['BH_FDR'])]for x in E['activity_f1']],
      [.8,1.2,.4,.9,.95,.9,.9],leftcols=(0,1),
      note='连续效应为rho，分组/四分位效应为Δ NES。四MR蛋白按每个端点独立BH；HIF1A为单独背景目标。n=9的分组结果保持探索性定位。')
    h('4.4 代谢通路的连续关联与通量配对结果')
    rows=[]
    for key,label in [('auc_cor','AUCell'),('flux_cor','scFEA')]:
        for x in E[key]:
            rows.append([label,pathnames[x['pathway']],'F1'if x['predictor']=='F1'else'TYMP',x['n_patients'],num(x['median_patient_rho']),str(round(float(x['fraction_positive'])*int(x['n_patients'])))+'/'+x['n_patients'],num(x['signed_rank_PValue']),num(x['FDR'])])
    T('附表 S4A 四通路与四通量轴的完整连续关联',['类型','通路/轴','预测量','n','rho','正向/n','P','FDR'],rows,
      [.8,1.4,.7,.35,.8,.7,.85,.85],leftcols=(0,1),note='AUCell与scFEA分别计算；每个预测量对应四特征家族。')
    T('附表 S4B scFEA 四轴的患者配对结果',['比较','代谢轴','n','Δ评分','正向/n','P','FDR'],
      [['TYMP阳性'if x['contrast']=='TYMP-positive_vs_negative'else'F1高主克隆',pathnames[x['pathway']],x['n_patients'],num(x['median_patient_delta']),str(round(float(x['fraction_patients_positive'])*int(x['n_patients'])))+'/'+x['n_patients'],num(x['paired_wilcox_PValue']),num(x['FDR'])]for x in E['flux_pair']],
      [1.2,1.35,.35,.8,.65,.85,.85],leftcols=(0,1),note='每种比较的四轴独立BH；推断评分不是实测代谢物浓度。AUCell完整配对结果见正文表5。')
    h('4.5 预设酶及转运体的完整 F1 关联')
    T('附表 S5 十二基因患者内 RNA 关联',['基因','n','rho','正向/n','P','FDR'],
      [[x['gene'],x['n_patients'],num(x['median_rho']),x['positive_patients']+'/'+x['n_patients'],num(x['P']),num(x['FDR'])]for x in E['enzyme_f1']],
      [1,.4,.9,.8,.9,.9],note='12基因家族独立BH，采用原近似Wilcoxon；TYMP效应与正文2.2相同，但其检验实现和校正家族不同。')
    h('4.6 三种拟时序及分化潜能的完整患者关联')
    for method in ['Monocle3_pseudotime_patient_v1','Slingshot_pseudotime_patient_v1','scTour_pseudotime_patient_v1','CytoTRACE2_Score']:
        T('附表 S6 '+methodnames[method]+' 的十二端点',['指标','n','rho','正向/n','P','FDR'],
          [[pathnames[x['target']],x['n_patients'],num(x['median_patient_rho']),str(round(float(x['fraction_positive'])*int(x['n_patients'])))+'/'+x['n_patients'],num(x['signed_rank_PValue']),num(x['FDR'])]for x in E['trajectory']if x['method']==method],
          [1.9,.4,.9,.85,.95,.95],note='各方法原12端点独立BH；CytoTRACE2为分化潜能，不作为方向相同的拟时序替代。'if method=='CytoTRACE2_Score'else'各方法原12端点独立BH；CNV负荷为19位患者，其余20位。')
    h('4.7 候选活性及 RNA 的原末端比较')
    T('附表 S7 原细胞中位数末端端点',['方法','目标','n','Δ','正/负/零','P','FDR'],
      [[x['method'],x['feature'],x['n_eligible_patients'],num(x['median_patient_effect']),direction(x),num(x['wilcoxon_p_two_sided']),num(x['BH_FDR'])]for x in E['activity_terminal']],
      [1.1,1.1,.35,.9,.9,.9,.9],leftcols=(0,1),
      note='四MR蛋白每方法独立BH，HIF1A单独；TYMP/RRM2B RNA每方法为两端点家族。蛋白为Δ NES，RNA为Δ log1p CP10K。该端点保留，不以S8的pseudobulk替换。')
    h('4.8 末端 pseudobulk 的全部酶基因结果')
    for method in ['Monocle3','Slingshot','scTour']:
        T('附表 S8 '+method+' 的十三酶末端 RNA',['基因','n','Δ RNA','正/负/零','P','FDR'],
          [[x['gene'],x['n'],num(x['median_effect']),direction(x,'positive','negative','zero'),num(x['P']),num(x['FDR'])]for x in E['terminal_enzyme']if x['method']==method],
          [1,.4,.9,1,.95,.95],note='全部13基因×3方法统一39项BH；Δ为log2(CPM+1)的末端−早期差值。')
    h('4.9 TYMP 和 RRM2B 功能代理的完整汇总')
    for ax in ['TYMP','RRM2B']:
        T('附表 S9 '+ax+' 严格与广义代理',['代理','端点/方法','n','中位效应','正/负/零','P','FDR'],
          [[axnames[x['feature']],x['method']if x['endpoint']=='trajectory_terminal_shift'else endpointnames[x['endpoint']],x['n_eligible_patients'],num(x['median_patient_effect']),direction(x),num(x['wilcoxon_p_two_sided']),num(x['BH_FDR'])]for x in E['axis']if x['feature'].startswith(ax)],
          [1.1,1.2,.35,.85,.9,.9,.9],leftcols=(0,1),
          note='连续端点为rho，其余为标准化代理Δ。F1高主克隆比较为14位患者，区别于9位患者的F1主导细胞比较；原BH家族保持。')

for x in blocks:
    typ=x['type']
    if typ=='table':main_table(x['id'],x['caption'])
    elif typ=='appendix':appendices()
    else:
        if typ=='h1':story.append(CondPageBreak(95))
        add(x['text'],typ)
(w/'report/report_materialized.json').write_text(json.dumps(material,ensure_ascii=False,indent=2))
md=[]
for x in material:
    if x['type']=='table':
        md+=['\n'+x['caption'],'| '+' | '.join(x['headers'])+' |','| '+' | '.join(['---']*len(x['headers']))+' |']
        md+=['| '+' | '.join(str(v).replace('|','/')for v in row)+' |'for row in x['rows']]
    elif x['type']in['title','h1','h2']:md+=['\n'+{'title':'# ','h1':'## ','h2':'### '}[x['type']]+x['text']]
    else:md+=['\n'+x['text']]
(w/'report/26210_TYMP_report_no_figures_v1.md').write_text('\n'.join(md))
output=w/'report/26210_TYMP_report_no_figures_v1.pdf'
page_meta=[]
def footer(c,doc):
    c.saveState();c.setFont('Song',9);c.drawCentredString(A4[0]/2,23,str(doc.page));c.restoreState()
    page_meta.append(doc.page)
doc=SimpleDocTemplate(str(output),pagesize=A4,rightMargin=45,leftMargin=45,topMargin=40,bottomMargin=39,
 title='肾癌恶性上皮细胞 TYMP 相关程序与代谢特征',author='项目26210',subject='2.1 2.2 2.4 2.5 2.6 2.7 无图综合报告')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
(w/'qc/TABLE_REGISTRY.json').write_text(json.dumps(table_registry,ensure_ascii=False,indent=2))
(w/'qc/BUILD.json').write_text(json.dumps({'pages':len(page_meta),'tables':len(table_registry),'images_inserted':0,'sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'bytes':output.stat().st_size},indent=2))
print(json.dumps({'output':str(output),'pages':len(page_meta),'tables':len(table_registry),'bytes':output.stat().st_size}))
