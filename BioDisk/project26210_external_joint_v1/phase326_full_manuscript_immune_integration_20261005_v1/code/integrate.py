from pathlib import Path
import json,copy,re,hashlib,shutil
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
BASE=R/'phase273_clinical_suppfig_wording_20260928_v1'
IMM=R/'phase325_immune_results_mechanism_revision_20261005_v1'
P=R/'phase326_full_manuscript_immune_integration_20261005_v1'
pdf=BASE/'report/26210_manuscript_clinical_suppfig_v1.pdf'
assert hashlib.sha256(pdf.read_bytes()).hexdigest()=='25e580c0b0f2e56567d6f431012fc425c50188c832d06fefb5e7004048311112'
P.mkdir(exist_ok=False)
for x in ['code','inputs','report','figures','qc','logs','delivery','root_record_backups']:(P/x).mkdir()
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:shutil.copy2(R/name,P/'root_record_backups'/name)
base=json.loads((BASE/'report/report_materialized.json').read_text())
imm=json.loads((IMM/'inputs/DOCUMENT_CONTENT.json').read_text())
for i,x in enumerate(base):x['_base_index']=i
for i,x in enumerate(imm):x['_immune_index']=i
# Local immune references are mapped to the complete manuscript reference list.
refmap={1:23,2:24,3:25,4:19,5:26,6:27,7:28,8:29,9:20,10:30,11:31,12:32,13:33,14:34,15:35}
def refs(t):
 def f(m):
  vals=[]
  for z in m[1].replace('–','-').split(','):
   if '-' in z:
    a,b=map(int,z.split('-'));vals.extend(range(a,b+1))
   else:vals.append(int(z))
  vals=[refmap[v] for v in vals]
  return '['+','.join(map(str,vals))+']'
 return re.sub(r'\[(\d+(?:[-–]\d+)?(?:,\d+(?:[-–]\d+)?)*)\]',f,t)
def p(t,**kw):return dict(type='p',text=t,**kw)
def h(t):return dict(type='h2',text=t,_new_heading=True)
def cv(x):
 y=copy.deepcopy(x)
 if x['type']=='h':y['type']='h1' if x['level']==1 else 'h2';y.pop('level',None)
 if 'text' in y:y['text']=refs(y['text'])
 if x['type']=='table':
  n=int(re.match(r'表免疫(\d+)',x['title'])[1]);y['caption']=f'[[ITABLE_{n}]] '+re.sub(r'^表免疫\d+\s*','',x['title'])
  y['_immune_table_id']=n;y['_col_weights']=y.pop('widths',None);y.pop('title',None)
 return y
def group(title):
 i=next(i for i,x in enumerate(imm) if x.get('text')==title)
 j=next((j for j in range(i+1,len(imm)) if imm[j]['type']=='h'),len(imm))
 return [cv(x) for x in imm[i+1:j]]
def narrative(title):return [x for x in group(title) if x['type']=='p']
def tables(title):return [x for x in group(title) if x['type']=='table']
before={};after={}
def add_before(i,items):before.setdefault(i,[]).extend(items)
def add_after(i,items):after.setdefault(i,[]).extend(items)
# Background in the original scientific-question section; no second independent report.
add_before(11,[
p('TYMP相关代谢状态与巨噬细胞功能之间的联系，是理解肾癌代谢微环境的重要研究问题。TYMP参与核苷代谢；多胺合成伴随MTA生成，而SPP1富集巨噬细胞具有脂质处理和细胞内降解相关特征。本研究关注这些代谢与细胞状态之间的对应，并进一步评价其免疫功能和组织分布。'),
p('已有结直肠癌研究显示，SPP1阳性巨噬细胞与FAP阳性成纤维细胞及基质环境相联系，部分干预模型支持SPP1参与T细胞扩增调节；肾癌研究也报道了SPP1相关巨噬特征与免疫治疗结局的联系[24,25,19]。据此，本研究分别评价巨噬细胞自身的免疫程序、同患者T细胞读数、空间基质及免疫分布，并以公开干预研究提供功能参照。')])
# Full data description retained, consolidated into the original data-source section.
add_after(18,group('数据来源与独立生物学样本'))
# Detailed immune methods remain in Methods, positioned after original macrophage methods.
for title in ['巨噬细胞成员与免疫指标定义','患者内关联与患者层统计汇总','相互条件调整与程序基因敏感性分析','跨细胞与组织免疫关联']:
 add_before(49,[h(title),*group(title)])
space_methods=narrative('空间免疫分布与外部干预资料')
assert len(space_methods)==3
add_before(68,[h('空间免疫界面与基质程序分析'),*space_methods[:2]])
add_before(76,[h('公开干预研究与源数据分析'),space_methods[2]])
newstat=narrative('多重检验与结果报告')[0]['text']
newstat=newstat[newstat.index('巨噬免疫发现的'):]
add_after(81,[p(newstat)])
# Four immune result sections between macrophage state and clinical outcomes.
r_titles=['TYMP与SPP1在巨噬细胞内的正相关得到独立队列重复','TYMP与SPP1对应不同的免疫功能程序','两基因均与糖酵解相关而总体炎症关联较不稳定','组织SPP1与免疫程序升高相关但未显示一致的T细胞抑制关系']
for title in r_titles:
 g=group(title)
 if title==r_titles[0]:
  g.insert(0,p('在SPP1富集状态及其功能特征得到定位后，进一步比较整体巨噬细胞和C4内部的TYMP相关表达。该分析将状态间比较延伸为同患者内的连续基因关联。'))
 if title==r_titles[1]:
  g.append(p('上述结果与前述C4相对其他巨噬的抗原呈递及干扰素程序降低，分别回答不同问题：C4比较描述两个细胞群的平均状态差异，本节描述患者内TYMP或SPP1表达与功能程序的连续变化，且免疫程序拓展采用独立的Obradovic队列。因而，两类结果共同呈现群间状态与群内表达联系，不能互相替代。'))
 add_before(173,[h(title),*g])
# Clinical results are already in the original manuscript; add only complementary interpretation/data.
clinical=narrative('SPP1相关评分与部分免疫治疗结局相关')
add_after(178,[
p('与OS关联不同，文献五基因评分的PFS高低组比较及CM025的PFS治疗交互未获支持。JAVELIN联合治疗组中，SPP1单基因与较短PFS为名义关联，临床调整HR=1.218，P=0.027445、FDR=0.10978；五基因评分未获得同样结果。C4的独立OS及相对治疗获益交互也未获相应支持（附表S19、S20）。'),
*clinical[-2:],*tables('SPP1相关评分与部分免疫治疗结局相关')])
# Spatial expression/perturbation evidence beside spatial programs, and interface evidence beside T-cell distributions.
sp=narrative('SPP1与基质相关并有部分外部干预证据支持免疫作用')
assert len(sp)==6
add_before(193,[h('SPP1与基质程序相联系且外部干预支持免疫调节作用'),*sp[:3],sp[4],*tables('SPP1与基质相关并有部分外部干预证据支持免疫作用')])
add_after(195,[sp[3],sp[5]])
# Discussion integrates mechanism and important boundaries, rather than appending another discussion.
add_after(231,[h('TYMP与SPP1共同表达对应不同的免疫特征'),*group('TYMP和SPP1共同升高不代表免疫作用相同'),
p('C4的群间低抗原呈递特征与TYMP的细胞内正向免疫程序关联并不矛盾。前者比较一种SPP1富集状态与其他巨噬，后者追踪同患者巨噬内部的表达共变化。这一区分将C4的状态定位、整体巨噬的连续关联和外部队列验证连接起来，避免将不同分析对象压缩为同一免疫功能结论。'),
h('糖酵解提供共同代谢背景而免疫程序具有方向差异'),*group('糖酵解是共同特征但总体炎症不宜作为主要机制')])
disc_space=narrative('SPP1下游实验支持免疫作用但肾癌屏障尚未被证实')
add_after(235,[disc_space[0]])
add_after(253,[
p('免疫程序与组成。TYMP和SPP1对干扰素γ及抗原呈递程序的关联在条件模型及去共享基因后保留；IFNα去重叠后未达到名义显著，总体炎症依赖细胞和程序定义，糖酵解在联合删除少数高影响基因后明显减弱。组织SPP1与细胞毒及功能障碍信号同时升高，而细胞类型限定的巨噬—CD8分析未显示一致抑制关系。这些结果支持具体免疫转录特征，不支持普遍炎症增强或所有SPP1相关巨噬均抑制T细胞的概括。')])
add_after(265,[
p('Obradovic的相互控制、基因删除和免疫程序拓展使用同一批11名患者，属于同队列内的机制定位和稳健性检查，不构成多次独立验证。巨噬细胞内TYMP关联与肿瘤来源TYMP处理效应分别解释。不同癌种或宿主Spp1干预的结果提供功能背景，不能替代肾癌中的TYMP干预，也不等同于已证明TYMP经MTA、SPP1、基质与T细胞构成完整因果路径。')])
# Append reproducibility index while keeping all historical statistical tables.
add_before(367,[
h('免疫程序与空间功能证据的完整统计'),
p('新增免疫部分保留原始患者资格、程序成员、相互条件模型、去共享基因与逐基因敏感性、跨细胞端点、空间界面及公开干预源表。完整数值采用各自原始检验家族；巨噬免疫发现224项、基因对24项、条件20项及去重叠36项分别记录。所有患者方向和不可估计端点均保留于配套免疫统计与来源资料。'),
p('正文区分本项目肾癌数据、公开肿瘤模型干预及文献引用。公开CRC源表提供Spp1缺失后CAF和CD8变化；原始流程的完整重现、数据源层复算与肾癌外部对应分别登记。文献评分、SPP1单基因、单细胞C4及旧空间参考状态保持原定义。')])
# Explicit edits needed for coherent synthesis; numerical facts are preserved.
changes={}
def replace(i,t):changes[i]=t
replace(0,'肾癌TYMP相关代谢与SPP1巨噬细胞的免疫及空间特征')
replace(1,'肾癌单细胞、多组学、免疫、临床与空间转录组联合研究　2026-10-05')
replace(3,base[3]['text']+' 免疫拓展另使用11名未治疗肾癌患者的单细胞数据及44名供者的组织转录组，并结合公开Spp1干预资料。')
add_after(5,[p('巨噬细胞内TYMP与SPP1正相关，但两者对应的免疫程序不同。独立11名患者中，TYMP与IFNγ、MHC-I和吞噬程序正相关，SPP1与IFNγ及MHC-II负相关；五条关系均11/11同向，P=0.0009766、原发现家族FDR=0.003716，控制另一基因及去共享基因后保持。两者共同关联糖酵解，提示共同代谢背景与免疫功能差异并存。')])
replace(8,'这些结果表明，TYMP相关恶性细胞状态伴随核苷、蛋氨酸及多胺代谢改变，SPP1富集巨噬具有脂质处理和溶酶体特征，并与相关代谢区域共同分布。TYMP与SPP1在巨噬细胞内对应不同的抗原处理及免疫响应程序；SPP1与基质的联系和公开干预中的CD8变化提供组织及功能参照。MTA保留为候选代谢联系，完整的TYMP经SPP1形成免疫屏障的因果路径尚未建立。')
replace(11,base[11]['text'].replace('巨噬细胞状态、临床结局及其空间分布','巨噬细胞状态及免疫程序、临床结局及其空间分布'))
replace(77,base[77]['text'].replace('独立生物学样本量以患者为单位','肾癌观察性分析以独立患者为生物学样本单位；公开干预按独立动物或独立混样池计数').replace('患者内相关和配对差异汇总后，采用','原恶性细胞和髓系状态分析的患者内相关及配对差异汇总后，采用'))
replace(81,base[81]['text']+' Obradovic患者内免疫关联采用相应方法小节定义的双侧精确符号检验，分家族BH结果与原患者配对Wilcoxon分析分别保留。')
replace(193,'SPP1相关区域保留T细胞相对分布及CD8界面差异')
replace(241,'这一工作模型将代谢候选、细胞状态、免疫转录特征及组织位置联系起来。巨噬细胞内TYMP与SPP1的共同变化并不要求免疫作用同向；两者与糖酵解的共同联系可以与干扰素和抗原处理程序的差异并存。SPP1相关巨噬与基质程序相联系，外部Spp1干预又支持部分免疫作用，为组织层面的机制解释提供参照。模型中的各联系来自不同测量层级，分段证据不等同于已经证实的连续代谢介导路径。')
replace(248,base[248]['text']+' 本文新增的巨噬细胞内TYMP—SPP1正相关描述同一细胞谱系内的共变，并不改变肿瘤来源TYMP与巨噬SPP1之间跨细胞联系的上述结果。')
replace(251,base[251]['text'].replace('成纤维及ECM程序的邻近检验也未获得一致正向支持，当前数据不支持普遍的SPP1-CAF屏障模式。','既有映射丰度与成纤维及ECM程序的邻近检验未获得一致正向支持。新增的公开标记表达程序分析支持SPP1状态与基质相关，但它使用不同于旧映射丰度邻近检验的读出；两类结果分别保留，不据程序正相关宣称普遍的SPP1-CAF屏障。')+' 固定界面内侧CD8信号较低与局部SPP1的对应并不稳定，因而免疫空间分区不能直接归因于SPP1特异的物理阻挡。')
replace(267,'本研究确立了TYMP相关恶性上皮状态的核苷、蛋氨酸及部分多胺代谢特征，识别出SPP1富集且具有胆固醇处理和溶酶体程序的巨噬状态，并在多来源空间组织中观察到代谢背景与SPP1相关区域的共同分布。配对代谢组和实验方向将MTA保留为候选代谢联系。巨噬细胞内TYMP与SPP1正相关，但分别对应不同的干扰素、抗原呈递及吞噬相关程序，共同糖酵解背景与免疫差异并存。SPP1与基质程序、部分临床结局及外部干预中的CD8变化进一步提供组织和功能对应。上述结果支持肿瘤代谢与巨噬免疫功能适应之间的分层联系，尚未证明TYMP通过SPP1形成阻挡T细胞的连续因果机制。')
replace(269,base[269]['text'].replace('rho 为患者内 Spearman 相关的跨患者中位数','rho 的估计对象随端点确定：细胞内分析报告患者内相关的跨患者中位数，组织及临床相关采用相应样本层系数；具体以表题及方法为准'))
replace(334,base[334]['text'].replace('图注说明','补充图S2展示SPP1相关评分与C4状态的临床关联。新增免疫细胞内及空间功能结果以正文统计表呈现。图注说明'))
# Materialize the full source.
out=[]
for i,x in enumerate(base):
 out.extend(copy.deepcopy(before.get(i,[])))
 y=copy.deepcopy(x)
 if i in changes:y['text']=changes[i]
 out.append(y)
 out.extend(copy.deepcopy(after.get(i,[])))
# Append unique immune references before supplement figures; retain the original 22 reference identities.
newrefs=[]
for x in imm:
 t=x.get('text','')
 m=re.match(r'^\[(\d+)\]',t)
 if m and refmap[int(m[1])]>=23:newrefs.append(cv(x))
j=next(i for i,x in enumerate(out) if x.get('text')=='补充图')
out[j:j]=sorted(newrefs,key=lambda x:int(re.match(r'^\[(\d+)\]',x['text'])[1]))
# Renumber all subsection headings in their parent chapter.
chapter=None;cnt={};sectionmap={}
for x in out:
 if x['type']=='h1':
  m=re.match(r'^([1-4]) ',x['text']);chapter=m[1] if m else None
 if x['type']=='h2' and chapter:
  cnt[chapter]=cnt.get(chapter,0)+1;num=chapter+'.'+str(cnt[chapter])
  oldnum=re.match(r'^\d+\.\d+\s+',x['text'])
  if oldnum:sectionmap[oldnum[0].strip()]=num
  title=re.sub(r'^\d+\.\d+\s+','',x['text'])
  x['text']=num+' '+title
# Determine main-table numbering by actual reading order.
tmap={};k=0
for x in out:
 if x['type']=='table' and not x['caption'].startswith('附表'):
  k+=1
  if x.get('_immune_table_id'):key='I'+str(x['_immune_table_id'])
  else:key='B'+re.match(r'^表\s*(\d+)',x['caption'])[1]
  tmap[key]=k
def renumtext(t):
 t=re.sub(r'\[\[ITABLE_(\d+)\]\]',lambda m:'表'+str(tmap['I'+m[1]]),t)
 t=re.sub(r'(?<!附)表\s*([1-6])(?!\d)',lambda m:'表'+str(tmap['B'+m[1]]),t) if '[[ITABLE_' not in t else t
 return t
# Apply base references only to base text; immune captions use their own table identities.
for x in out:
 for field in ['text','caption']:
  if field not in x:continue
  t=x[field]
  if '[[ITABLE_' in t:
   t=re.sub(r'\[\[ITABLE_(\d+)\]\]',lambda m:'表'+str(tmap['I'+m[1]]),t)
  elif '_base_index' in x:
   t=re.sub(r'(?<!附)表\s*([1-6])(?!\d)',lambda m:'表'+str(tmap['B'+m[1]]),t)
  # Existing Figure4 caption is the only prose link to numbered subsections.
  t=t.replace('见2.11和3.8节','见2.11和'+sectionmap['3.8']+'节')
  x[field]=t
# A figure source never changes; copies support stable native-vector merge.
for f in (BASE/'figures').glob('*.pdf'):shutil.copy2(f,P/'figures'/f.name)
(P/'report/report_materialized.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
renderer=(BASE/'code/render_report.py').read_text().replace(BASE.name,P.name).replace('26210_manuscript_clinical_suppfig_v1','26210_manuscript_immune_integrated_v2')
renderer=renderer.replace("if x['caption'].startswith('附表 S1B'):","if x.get('_col_weights'):weights=x['_col_weights']\n  if x['caption'].startswith('附表 S1B'):")
renderer=renderer.replace("'SPP1与C4临床关联 含补充图S1-S2'","'TYMP相关代谢 SPP1巨噬免疫与空间联系 含补充图S1-S2'")
# Supplementary data-source notes now remain attached to their tables.
renderer=renderer.replace("story += ([] if len(rows)<=10 else [CondPageBreak(80)])+[cap,tt,Spacer(1,6)]","story += ([] if len(rows)<=10 else [CondPageBreak(80)])+[cap,tt,Spacer(1,6)]\n  if x.get('_immune_table_id') and x.get('source'):story.append(para(x['source'],'note'))")
(P/'code/render_report.py').write_text(renderer)
# Proof that every original table and every newly reviewed immune table is retained verbatim.
bt=[x for x in base if x['type']=='table'];it=[x for x in imm if x['type']=='table'];nt=[x for x in out if x['type']=='table']
assert len(nt)==41,(len(nt),len(bt),len(it))
for x in bt+it:
 assert any(x['headers']==y['headers'] and x['rows']==y['rows'] for y in nt),x
assert cnt=={'1':25,'2':26,'3':12,'4':16},cnt
assert len(newrefs)==13
assert len({x['_base_index'] for x in out if '_base_index' in x})==len(base)
sources=[pdf,BASE/'report/report_materialized.json',IMM/'inputs/DOCUMENT_CONTENT.json',IMM/'delivery/26210_免疫章节_结果机制解释_v5.docx',*list((BASE/'figures').glob('*.pdf'))]
(P/'inputs/SOURCE_MANIFEST.json').write_text(json.dumps([{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sources],ensure_ascii=False,indent=2))
(P/'inputs/USER_PDF_IDENTITY.json').write_text(json.dumps({'user_path':r'D:\电脑管家迁移文件\xwechat_files\wxid_913cjq0lvmfd32_3e57\msg\file\2026-09\26210_manuscript_clinical_suppfig_v1.pdf','user_sha256':'25e580c0b0f2e56567d6f431012fc425c50188c832d06fefb5e7004048311112','server_pdf':str(pdf),'exact_match':True},ensure_ascii=False,indent=2))
(P/'qc/CONTENT_MERGE_VALIDATION.json').write_text(json.dumps({'status':'PASS','all_34_original_tables_preserved':True,'all_7_immune_tables_preserved':True,'all_original_blocks_retained':True,'section_counts':cnt,'section_map':sectionmap,'main_table_map':tmap,'reference_map':refmap,'new_references':13,'new_statistics':False,'source_figures':7},ensure_ascii=False,indent=2))
(P/'qc/EDITORIAL_CHANGE_LOG.json').write_text(json.dumps({'base_text_changes':[{'index':i,'before':base[i]['text'],'after':t} for i,t in changes.items()],'integration_choices':['Full immune methods integrated by analysis type; clinical method duplicate absorbed by original more detailed method','Four immune result sections inserted after original2.12 and before clinical; clinical duplicate effects not repeated','New ECM/source-intervention result placed beside spatial results; CD8 interface evidence merged with T-cell spatial section','Discussion distinguishes C4 between-state effects from TYMP within-cell effects, old abundance adjacency from expression-program correlation, and tumour source from myeloid source','Abstract, scientific question, working model and conclusion updated without new statistical analyses']},ensure_ascii=False,indent=2))
print('MERGED',len(out),'blocks',cnt,'41 tables',len(newrefs)+22,'references',tmap)
