from pathlib import Path
import json,re,hashlib,zipfile,csv,io,html,math
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pymupdf as fitz
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image,ImageDraw

w=Path(__file__).resolve().parent;q=w/'qc';out=w/'delivery';figs=w/'figures'
for p in [q,out,figs]:p.mkdir(exist_ok=True)
mods=json.loads((w/'MODULES.json').read_text(encoding='utf-8'))
excluded={'E12':'旧SPP1-state轨迹/通信与当前C4定义不同，且当前主稿已有新的轨迹结果；本轮不重复展开。','E29':'该独立病例主结果已在完整稿2.25呈现，新增信息有限，不再单列章节。','E30':'版本失败、资格不足与旧资产只进入盘点记录，不作为科学补充结果。'}
mods=[m for m in mods if m['id'] not in excluded]
manifest=json.loads((q/'STATISTICAL_SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
assert all(hashlib.sha256((w/'statistics'/x['file']).read_bytes()).hexdigest()==x['sha256'] for x in manifest)
def table(phase,end):
 a=[x for x in manifest if x['phase']==phase and x['file'].endswith(end)]
 assert a,(phase,end)
 x=a[0];return pd.read_csv(w/'statistics'/x['file'],sep='\t'),x
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
plot_sources=[]
# S8: actual eight cell-line baseline values; no fitted model or new P value.
nci,x=table(103,'renal_line_level_v1.tsv');plot_sources.append(x)
fig,axs=plt.subplots(1,2,figsize=(7.48,4.8));fig.subplots_adjust(left=.10,right=.96,bottom=.25,top=.90,wspace=.35)
for ax,col,title in zip(axs,['tymp_log2_fpkm_plus1','mtap_log2_fpkm_plus1'],['TYMP','MTAP']):
 ax.scatter(nci[col],nci.mta_mean_raw,s=28,c='#97516D');ax.axhline(.086,c='#AAA',ls='--',lw=.7)
 floor_count=0
 for num,row in enumerate(nci.itertuples(),start=1):
  if row.mta_at_global_floor:
   offset=(4+floor_count*5,8+floor_count*12);floor_count+=1
  else:offset=(4,5)
  ax.annotate(str(num),(getattr(row,col),row.mta_mean_raw),xytext=offset,textcoords='offset points',fontsize=7,arrowprops={'arrowstyle':'-','lw':.4,'color':'#777'} if row.mta_at_global_floor else None)
 ax.set(xlabel=title+' log2(FPKM+1)',ylabel='Measured MTA (relative abundance)');ax.set_title(title+' vs MTA | renal n=8',loc='left',weight='bold')
fig.text(.06,.095,'   |   '.join(f'{k+1}: {s}' for k,s in enumerate(nci.official_cell_line) if k<4),fontsize=7)
fig.text(.06,.045,'   |   '.join(f'{k+1}: {s}' for k,s in enumerate(nci.official_cell_line) if k>=4),fontsize=7)
fig.savefig(figs/'S8.pdf');plt.close(fig)
# S9: all eight original decomposition combinations, including descriptive n4.
dec,x=table(282,'DECOMPOSITION_SUMMARY.tsv');plot_sources.append(x)
fig,axs=plt.subplots(1,2,figsize=(7.48,4.9));fig.subplots_adjust(left=.22,right=.96,bottom=.14,top=.87,wspace=.50)
for ax,model in zip(axs,['core','context']):
 z=dec[dec.model==model].sort_values(['scope','part']);y=np.arange(len(z));ax.axvline(0,c='#999',lw=.8)
 ax.scatter(z['median'],y,c=['#A55A72' if s=='C4' else '#508798' for s in z.scope],s=32)
 for k,row in enumerate(z.itertuples()):ax.annotate(f'n={row.n}',(row.median,k),xytext=(5,6),textcoords='offset points',size=7)
 ax.set_yticks(y,[s+' / '+p.replace('co_detection','co-detection').replace('double_positive_intensity','positive intensity') for s,p in zip(z.scope,z.part)],fontsize=7)
 ax.set(xlabel='Median patient effect',xlim=(-.01,.22));ax.invert_yaxis();ax.set_title('Core model' if model=='core' else 'Environment-adjusted',weight='bold')
fig.savefig(figs/'S9.pdf');plt.close(fig)
# S10: all prespecified source x endpoint comparisons in the same 11 patients.
rc,x=table(306,'RECIPIENT_ASSOCIATIONS_ALL.tsv');plot_sources.append(x)
fig,axs=plt.subplots(1,2,figsize=(7.48,4.8));fig.subplots_adjust(left=.19,right=.96,bottom=.18,top=.84,wspace=.40)
names={'PDCD1_ENTPD1_double_detected_fraction':'PD-1 / CD39 detection','TOX_detected_fraction':'TOX detection','IFNG_GZMB_effector_mean_logCP10K':'IFNG / GZMB program'}
for ax,source in zip(axs,['SPP1','TYMP']):
 z=rc[(rc.dataset=='OBRADOVIC')&(rc.source==source)]
 for j,(model,label,col) in enumerate([('batch_stage_adjusted','Batch + stage','#467F9C'),('batch_stage_CD8depth_adjusted','+ CD8 depth','#AB6179')]):
  zz=z[z.model==model].set_index('endpoint').loc[list(names)];ax.scatter(zz.rho,np.arange(3)+(j-.5)*.15,label=label,color=col,s=30)
 ax.set_yticks(range(3),list(names.values()),fontsize=7);ax.invert_yaxis();ax.axvline(0,c='#888',ls='--',lw=.7);ax.set(xlim=(-.85,.3),xlabel='Partial Spearman rho');ax.set_title(source+' source | n=11',weight='bold')
axs[1].legend(frameon=False,fontsize=7,loc='lower left');fig.savefig(figs/'S10.pdf');plt.close(fig)
# S11: individual independent mixed pools, all values retained, no inferential P.
pool,x=table(308,'POOL_ENDPOINTS_ALL.tsv');plot_sources.append(x)
endpoints=[('Fibroblast','collagen_maturation','Collagen maturation'),('Fibroblast','focal_adhesion','Focal adhesion'),('Fibroblast','actomyosin','Actomyosin'),('Fibroblast','myCAF_markers','myCAF program'),('Fibroblast','steady_state_markers','Steady-state program'),('Macrophage','CXCL9_10','Macrophage CXCL9/10')]
fig,axs=plt.subplots(2,3,figsize=(7.48,6.1));fig.subplots_adjust(left=.09,right=.98,bottom=.12,top=.93,hspace=.55,wspace=.45)
for ax,(lineage,endpoint,title) in zip(axs.ravel(),endpoints):
 z=pool[(pool.lineage==lineage)&(pool.endpoint==endpoint)&pool.eligible]
 for j,group in enumerate(['WT','KO']):
  vals=z[z.group==group].sort_values('library').value.values;assert len(vals)==3
  ax.scatter(j+np.array([-.08,0,.08]),vals,c='#548799' if j==0 else '#A25C78',s=26);ax.plot([j-.18,j+.18],[np.median(vals)]*2,c='#333',lw=1)
 ax.set_xticks([0,1],['WT','Spp1 KO']);ax.set(xlim=(-.5,1.5),ylabel='Mean log1p(CP10K)');ax.set_title(title,fontsize=8,weight='bold')
fig.savefig(figs/'S11.pdf');plt.close(fig)
(q/'NEW_FIGURE_SOURCES.json').write_text(json.dumps(plot_sources,ensure_ascii=False,indent=2),encoding='utf-8')

# Disposition index: distinguish directory registration from result-level review.
entries=[x for x in json.loads((w/'PHASE_FILE_INVENTORY.json').read_text(encoding='utf-8')) if x['phase']<=336]
inmods={}
for m in mods:
 for n in m['phases']:inmods.setdefault(n,[]).append(m['id'])
report_phases={int(f.name[:3]) for f in (w/'evidence_reports').iterdir() if f.is_file()}
presentation={22,23,40,41,42,43,44,45,46,56,68,89,90,91,92,93,94,95,96,100,106,108,110,111,112,113,115,116,117,118,119,120,122,125,127,128,129,130,131,132,134,135,136,137,138,139,140,141,142,143,146,149,152,155,156,164,165,171,172,173,181,185,188,189,190,191,192,193,194,195,196,197,198,199,201,202,203,204,205,206,207,208,209,210,211,212,215,216,217,218,221,224,225,227,228,232,233,236,237,239,241,245,247,250,251,252,253,254,255,258,259,260,261,262,263,264,265,266,267,268,272,273,274,279,281,283,285,286,288,294,297,299,300,301,302,304,310,315,316,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336}
ledger=[]
for e in entries:
 n=e['phase'];name=Path(e['root']).name
 if n in inmods:decision='纳入当前补充';why='对应'+','.join(inmods[n])+'；正文摘要与完整统计分层保留';depth='结果摘要/统计源表核对'
 elif n in [83,84,85,69,77,79,114]:decision='历史或资格限制，不纳入';why='旧版本、状态定义或资格未满足；不以当前科学结论引用';depth='历史结果与版本规则核对'
 elif n in [87,88,109,248,257]:decision='不重复展开';why=excluded['E12'] if n in [87,88,109] else excluded['E29'];depth='结果摘要及正文对照'
 elif n in presentation:decision='实现/整理/主稿已有';why='该阶段以方案、身份整理、资产恢复或文图修订为主；不按阶段号重复新增分析';depth='阶段记录/文件索引对照'
 else:decision='登记，未纳入本轮证据';why='未作为本轮选定模块的有效新增结果；保留来源，不据目录名宣称分析未做或已通过';depth='目录/文件索引核查'
 ledger.append([n,name,e['root'],len(e['files']),decision,','.join(inmods.get(n,[])),depth,why])
for n in [1,59,146,147,148]:
 if n not in {x[0] for x in ledger}:ledger.append([n,'无同名阶段目录',('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1' if n in [146,147,148] else ''),'','非同名目录/阶段记录','', '交接记录核查','146–148对应part2_rebuild_v1，不代表未实施；1/59未见同名目录，不作完成判断'])
ledger.sort(key=lambda x:(x[0],x[1]));header=['阶段','目录名','服务器路径','浅层登记文件数','处置','补充模块','核查层次','依据与说明']
buf=io.StringIO();cw=csv.writer(buf,delimiter='\t',lineterminator='\n');cw.writerow(header);cw.writerows(ledger)
(out/'逐阶段盘点与取舍.tsv').write_text(buf.getvalue(),encoding='utf-8')

# Select deliverable statistics, keep cell/pixel caches and superseded attempts on server.
reject=[];kept=[]
for x in manifest:
 name=Path(x['file']).name.lower();reason=None
 if x['module'] in excluded:reason='该模块不进入当前补充正文'
 elif any(t in name for t in ['pixel_qc','cell_pseudotime','profiles_v2','state_membership','pseudobulk.tsv','all_model_inputs','spatial_readouts','cell_identity','gene_coverage','marker_coverage','cluster_markers','refined_markers','qc_by_','source_audit','source_access','filelist','identity_profile','download','source_column_units','all_numeric_source_cells','all_observations','panel_inventory','input_audit','coordinate_map','refined_identity']):reason='逐细胞/像素/表达缓存、身份及技术来源辅助资料保留服务器'
 elif x['phase']==307 and ('cd44_response_associations_all' in name or 'cd44_response_association_pairs' in name):reason='有效结果采用患者伪bulk V2，初版资格失败表不作当前模型'
 elif x['phase']==293 and name.endswith('perturbation_all_genes.tsv.gz'):reason='采用QUALIFIED版本，原未限定版本不重复'
 elif x['bytes']<5:reason='空结果文件，不当有效统计'
 if reason:reject.append(dict(x,reason=reason))
 else:kept.append(x)
(q/'FINAL_TABLE_SELECTION.json').write_text(json.dumps({'included':kept,'excluded':reject},ensure_ascii=False,indent=2),encoding='utf-8')

# Main reader document.
pdfmetrics.registerFont(TTFont('Song',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('Hei',r'C:\Windows\Fonts\simhei.ttf'))
styles={k:ParagraphStyle(k,fontName='Hei' if k in ['title','h1','h2'] else 'Song',fontSize=fs,leading=ld,spaceAfter=after,wordWrap='CJK',keepWithNext=k in ['h1','h2']) for k,fs,ld,after in [('title',18,26,15),('h1',13,19,9),('h2',10.5,16,5),('body',10.2,15.5,8),('note',8,11.5,7),('cell',8.5,12,0)]}
def P(s,k='body'):
 assert '敏感性' not in s
 return Paragraph(html.escape(s).replace('\n','<br/>'),styles[k])
story=[P('TYMP相关代谢、SPP1巨噬状态及基质免疫联系','title'),P('扩充补充材料 v2｜与完整稿v6配套','h1'),P('本册依据历史阶段盘点重新整理。范围限定于当前研究的代谢背景、髓系状态、候选输出与接收端、临床联系和空间结构。正文已充分呈现且无新增信息的独立病例重复、旧状态轨迹/通信，以及纯实现和排版记录不再单列。'),P(f'共整理{len(mods)}个分析模块，并保留原补充图S3–S7、新增S8–S11。新增图全部取自既有结果，没有重新拟合模型或产生新的推断P值。完整统计表按模块分组，原始数据、逐细胞缓存和代码仍留在服务器。'),P('统计及版本原则','h2'),P('独立生物学样本是推断单位；细胞、spot、区域和技术重复不增加n。既有原P、FDR、效应量与检验家族分别保留，不把不同目的家族中较低的FDR择优替换。小样本、混样池、单供者和模型扰动按各自证据类型解释；不可评估不填零或P=1。'),P('与完整稿的关系','h2'),P('完整稿v6和主图2保持不变；已有补充图S1/S2及附表S1–S20不重新编号。这里补充未展开的模型、方向、输出和限制，不把关键阴性从正文移走。采用最终版本；初版失败、错误编码和后续不再适用的身份定义只在盘点清单中登记。'),PageBreak(),P('分析目录与正文对应','h1')]
for m in mods:story.append(P(m['id']+'  '+m['title']+' ｜ '+m['main'],'note'))
story +=[PageBreak()]
toc=[]
for m in mods:
 h=P(m['id']+'  '+m['title'],'h1');h.bookmark=m['id'];story.append(h)
 story +=[P('分析对象与方法','h2'),P(m['method']),P('主要结果','h2'),P(m['result'])]
 rows=[[P(x,'cell') for x in ['比较/读出','实际单位或范围','原结果与统计']]]+[[P(x,'cell') for x in row] for row in m['rows']]
 tb=Table(rows,colWidths=[125,135,255],repeatRows=1,hAlign='LEFT');tb.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#DFEAF0')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F3F5F7')]),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story+=[tb,Spacer(1,8),P('结果解释与适用范围','h2'),P(m['boundary']),P('数据来源：Phase '+', '.join(map(str,m['phases']))+'；完整统计见同名模块目录。'+('既有D01–D03表保留在BASE文件夹。' if m['id']=='E01' else '既有D05表保留在BASE文件夹。' if m['id']=='E28' else ''),'note'),Spacer(1,12)]
story +=[PageBreak(),P('盘点范围、版本取舍及数据索引','h1'),P('历史目录索引覆盖403个目录、331个不同阶段编号（截至Phase336），并补充登记没有同名目录的阶段。Phase146–148的实施入口是part2_rebuild_v1。阶段号不等于独立研究或独立统计分析。'),P('本轮先对阶段记录和文件索引登记去向，再对纳入模块的结果报告及统计文件作内容核对。逐阶段清单明确区分“结果/源表核对”与“仅目录登记”，不把未深读的原始对象写成已审计通过。'),P('本次取舍','h2')]
for s in excluded.values():story.append(P(s))
story +=[P('此外，失败bulk分解、未合格MTA转移签名、历史客户实验方向指纹，以及与当前问题不相容的旧数据资产不进入科学结果正文。它们保留在来源目录及逐阶段登记表。'),P(f'统计交付包含{len(kept)}份选定文件，另保留上一版D01–D06数据作为BASE。文件按字节复制并核对SHA256；完整来源索引记录原路径、阶段、大小和校验值。逐细胞/像素缓存、表达矩阵及资格失败的替代版本另列排除原因，保留服务器原件。'),P('补充图索引','h2')]
figlist=[('S3','cNMF代表基因载荷'),('S4','cNMF运行稳定性与来源构成'),('S5','F1高主克隆的配对差异表达'),('S6','控制另一基因及删去重叠基因后的免疫程序关联'),('S7','去LGALS3及参考可识别性'),('S8','肾来源细胞系TYMP/MTAP与MTA'),('S9','SPP1检出与双阳性表达强度'),('S10','CD8接收端的原模型与深度调整模型'),('S11','宿主Spp1缺失后的逐池基质程序')]
for n,t in figlist:story.append(P(n+'  '+t,'note'))
class Doc(SimpleDocTemplate):
 def afterFlowable(self,f):
  if hasattr(f,'bookmark'):toc.append([1,f.getPlainText(),self.page])
def foot(c,doc):
 c.setFont('Song',8);c.setFillColor(colors.HexColor('#637582'));c.drawString(40,A4[1]-24,'26210｜扩充补充材料 v2');c.drawCentredString(A4[0]/2,22,str(doc.page))
base=w/'text_layout.pdf';Doc(str(base),pagesize=A4,leftMargin=40,rightMargin=40,topMargin=43,bottomMargin=40).build(story,onFirstPage=foot,onLaterPages=foot)
doc=fitz.open(base);textpages=len(doc)
old=fitz.open(r'C:\Users\31937\Downloads\26210\Phase336_补充材料_S3-S7\26210_Supplementary_Materials_S3-S7_v1.pdf')
for idx,(n,title) in enumerate(figlist[:5],start=1):
 doc.insert_pdf(old,from_page=idx,to_page=idx);p=doc[-1];p.add_redact_annot(fitz.Rect(35,43,570,72),fill=(1,1,1));p.apply_redactions(images=0,graphics=0)
 # A small ReportLab overlay preserves legible Chinese typography.
 from reportlab.pdfgen import canvas
 buf=io.BytesIO();cv=canvas.Canvas(buf,pagesize=A4);cv.setFont('Hei',12);cv.drawString(40,A4[1]-61,'补充图'+n+'  '+title);cv.save();overlay=fitz.open(stream=buf.getvalue(),filetype='pdf');p.show_pdf_page(p.rect,overlay,0)
 toc.append([1,'补充图'+n+' '+title,len(doc)])
captions={
'S8':'8个肾来源NCI60细胞系的基线RNA与实测MTA。横轴分别为TYMP或MTAP的log2(FPKM+1)，纵轴为原代谢相对丰度；虚线为观测下限0.086。每点为一个细胞系的技术重复均值，4/8位于下限。TYMP–MTA rho=0.3298，MTAP–MTA rho=-0.09570；不作小样本推断P、FDR或处理效应解释。完整全谱系与组织调整结果见E02。',
'S9':'原八个范围×读出×模型组合全部显示：共同检出与双阳性表达强度分别计算。点为原患者效应中位数，标示每项实际n；C4共同检出仅4患者，保留描述性效应而不补P。其他成熟巨噬的双阳性强度联系通过原家族FDR，C4对应结果未达原P<0.05。此图不按P值大小推算机制贡献比例，完整P/FDR见E15。',
'S10':'同11位Obradovic患者的全部TYMP/SPP1来源与三个CD8接收端比较。蓝点为批次/分期调整，红点为另加CD8测序深度；横轴原偏Spearman rho，无重新拟合。PD-1/CD39共同检出的部分负向名义结果随模型改变，TOX和IFNG/GZMB未提供支持性名义结果。两个模型来自同一患者，不视为独立复制；全部P/FDR及HCRN描述结果见E21。',
'S11':'宿主Spp1缺失的PDAC实验，每个点为一个独立混样池（每池2个肿瘤），每组3池；横线为组中位数。显示预设胶原成熟、黏着斑、收缩、myCAF、稳态及巨噬CXCL9/10全部六项，不新增推断检验。基质程序与CXCL9/10的方向分别解释；其余组成、SPP1及资格不足的CD8端点完整保留于E22，不能把基质变化等同于免疫功能恢复。'}
for n,title in figlist[5:]:
 p=doc.new_page(width=A4[0],height=A4[1]);src=fitz.open(figs/f'{n}.pdf');fh=(A4[0]-56)*src[0].rect.height/src[0].rect.width
 buf=io.BytesIO();cv=canvas.Canvas(buf,pagesize=A4);cv.setFont('Hei',13);cv.drawString(40,A4[1]-57,'补充图'+n+'  '+title)
 para=P(captions[n]);ww,hh=para.wrap(A4[0]-80,200);para.drawOn(cv,40,A4[1]-(90+fh+24)-hh);cv.save();over=fitz.open(stream=buf.getvalue(),filetype='pdf');p.show_pdf_page(p.rect,over,0)
 p.show_pdf_page(fitz.Rect(28,90,A4[0]-28,90+fh),src,0)
 toc.append([1,'补充图'+n+' '+title,len(doc)])
for k,p in enumerate(doc):
 p.draw_rect(fitz.Rect(220,A4[1]-37,375,A4[1]-10),fill=(1,1,1),color=None);p.insert_text((280,A4[1]-20),f'SM-{k+1}',fontsize=8)
doc.set_toc([[1,'概览与目录',1]]+toc)
pdf=out/'26210_Expanded_Supplement_v2.pdf';doc.save(pdf,garbage=4,deflate=True)
fulltext='\n'.join(p.get_text() for p in doc);assert '敏感性' not in fulltext
for k,p in enumerate(doc):p.get_pixmap(matrix=fitz.Matrix(1.15,1.15)).save(q/f'page_{k+1:02}.png')
for st in range(0,len(doc),6):
 sheet=Image.new('RGB',(1190,1680),'#e6e6e6')
 for j,k in enumerate(range(st,min(st+6,len(doc)))):
  im=Image.open(q/f'page_{k+1:02}.png').convert('RGB');im.thumbnail((390,535));x=(j%3)*396;y=(j//3)*840;sheet.paste(im,(x,y+28));ImageDraw.Draw(sheet).text((x+5,y+5),f'Page {k+1}',fill='black')
 sheet.save(q/f'contact_{st//6+1}.png')
readme=f'''扩充补充材料 v2\n包含{len(mods)}个筛选后的分析模块、S3–S11九张补充图。\n逐阶段登记以403个历史目录/331个编号为索引范围；146–148对应part2_rebuild_v1。清单明确区分索引核查和结果源表核对。\n当前直接相关的科学结果进入正文，旧定义、重复正文、失败和纯实现版本不按阶段数堆入。\n全文不使用用户禁用的中文术语；原始统计字段、来源路径和旧文件名为溯源信息，保持原值。\n完整稿仍Phase335 v6；旧补充v1保留。\n'''
(out/'阅读说明.txt').write_text(readme,encoding='utf-8')
with zipfile.ZipFile(out/'26210_Expanded_Statistics_v2.zip','w',zipfile.ZIP_DEFLATED) as z:
 for x in kept:z.write(w/'statistics'/x['file'],x['file'])
 z.writestr('SOURCE_MANIFEST.json',json.dumps(kept,ensure_ascii=False,indent=2).encode('utf-8'))
 z.writestr('MODULE_INDEX.json',json.dumps(mods,ensure_ascii=False,indent=2).encode('utf-8'))
 z.write(out/'逐阶段盘点与取舍.tsv','逐阶段盘点与取舍.tsv');z.writestr('阅读说明.txt',readme.encode('utf-8'))
 with zipfile.ZipFile(r'C:\Users\31937\Downloads\26210\Phase336_补充材料_S3-S7\26210_Supplementary_Data_D01-D06_v1.zip') as prior:
  for name in prior.namelist():z.writestr('BASE/'+name,prior.read(name))
summary={'modules':len(mods),'included_phase_numbers':sorted(inmods),'historical_directories':403,'distinct_historical_phases':331,'registered_rows':len(ledger),'source_table_candidates':len(manifest),'delivered_tables':len(kept),'excluded_table_files':len(reject),'text_pages':textpages,'pages':len(doc),'figures':9,'new_figures':4,'prohibited_term_count':fulltext.count('敏感性'),'new_statistical_tests':False,'visual_review':'PENDING'}
(q/'BUILD_VALIDATION.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(summary)
