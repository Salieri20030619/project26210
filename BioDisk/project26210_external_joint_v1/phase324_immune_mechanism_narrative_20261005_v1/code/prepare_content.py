from pathlib import Path
import json,copy,re,shutil,hashlib
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
S=R/'phase323_immune_manuscript_adaptation_20261005_v1'
P=R/'phase324_immune_mechanism_narrative_20261005_v1'
P.mkdir(exist_ok=False)
for d in ['code','inputs','delivery','qc','logs','render','root_record_backups']:(P/d).mkdir()
for n in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:shutil.copy2(R/n,P/'root_record_backups'/n)
old=json.loads((S/'inputs/DOCUMENT_CONTENT.json').read_text())
main=[]
def p(x):main.append({'type':'p','text':x})
def h(x,l=2):main.append({'type':'h','level':l,'text':x})
main.append({'type':'title','text':'TYMP与SPP1相关巨噬细胞的代谢及免疫特征'})
h('研究背景',1)
p('前述研究从肾癌细胞的TYMP相关代谢变化出发，发现核苷、蛋氨酸和多胺代谢与SPP1富集巨噬细胞之间存在联系。TYMP编码胸苷磷酸化酶，参与胸苷向胸腺嘧啶和脱氧核糖磷酸的转化。[1] 本项目识别出的SPP1富集C4巨噬细胞具有胆固醇处理和溶酶体特征。研究由此延伸到一个具体问题：这些代谢与细胞状态的变化，对应巨噬细胞怎样的免疫活动，又与周围T细胞和基质有何联系。')
p('这一问题具有明确的组织背景。已有结直肠癌研究显示，SPP1阳性巨噬细胞与FAP阳性成纤维细胞相联系；在部分肿瘤模型中，SPP1还参与胶原沉积及T细胞扩增调节。肾癌研究也报道了SPP1相关巨噬细胞与CD8 T细胞功能和免疫治疗结局的关系。[2–4] 因此，除关注巨噬细胞数量，还需要分析其承担的抗原处理、免疫信号响应和基质相关活动。')
p('本部分围绕TYMP与SPP1的共同变化及功能差异，连接前文的代谢发现与后续免疫、空间和临床观察。C4用于定位SPP1富集状态，整体巨噬细胞用于判断这些联系是否具有更广泛的细胞背景；肿瘤来源TYMP与巨噬细胞内TYMP分别解释。')
h('研究方法',1)
p('利用前述肾癌单细胞数据、独立的11名未治疗患者数据，以及组织、空间和临床队列，比较TYMP、SPP1表达与免疫及代谢相关基因程序。先在每名患者内评价巨噬细胞关联，再汇总不同患者的方向，并通过更严格的细胞筛选、控制另一基因和去除共享基因检查结果。随后比较同患者T细胞状态、SPP1相关区域的基质与免疫分布，并结合公开干预实验及免疫治疗随访资料。统计以独立患者或动物为单位；完整方法、效应量、P值、FDR和敏感性结果见配套技术附录。')
h('主要发现',1)
h('TYMP与SPP1及巨噬细胞输出相关基因共同变化')
p('TYMP与SPP1在巨噬细胞内呈正相关，并在独立肾癌队列中得到重复。原集合中29名患者有24名同向；独立队列严格筛选后，11名患者中有10名同向。原集合的患者内相关中位数为0.0491，P=0.000546、FDR=0.002184；独立队列为0.0470，P=FDR=0.011719。方向在多数患者中一致，细胞内关联幅度较小。')
p('TYMP还与TGFB1及部分溶酶体相关输出基因相联系。原集合中TYMP与TGFB1在26/31名患者中同向，相关中位数为0.0526，P=0.000192、FDR=0.000769。与前文C4的溶酶体和脂质处理特征结合，这些结果把TYMP的髓系联系具体定位到SPP1相关细胞状态及其伴随的功能过程。')
h('TYMP与SPP1对应不同的免疫响应和抗原处理特征')
p('两基因虽然共同变化，但各自对应的免疫程序不同。TYMP较高时，干扰素γ响应、MHC-I抗原呈递和吞噬相关基因表达较高；SPP1较高时，干扰素γ响应及MHC-II抗原呈递相关基因表达较低。独立队列的11名患者在这五条关系上均呈一致方向，原发现分析均达到校正后显著。TYMP三条关系的相关中位数依次约为0.189、0.205和0.194，SPP1两条关系约为−0.088和−0.135；P均为0.0009766，FDR均约为0.003716。')
p('控制另一基因的表达，并排除程序间共享的基因后，上述五条关系仍然保留。由此可见，TYMP与SPP1各自反映了巨噬细胞内不同的免疫转录特征。两者的共同表达，不足以概括细胞的全部免疫状态。')
h('不同免疫特征伴随共同的糖酵解相关变化')
p('TYMP与SPP1均与糖酵解相关基因程序正相关，11名患者中各有10名同向；相关中位数分别为0.1065和0.1043，P均为0.011719、FDR均为0.025485。两基因的免疫程序不同，却可伴随相同方向的代谢变化。')
p('TYMP与总体炎症程序在宽松细胞集合中也达到校正后显著，但对细胞筛选和程序成员更敏感；TYMP与TNFα/NFκB程序则呈负向关联。这表明与TYMP相联系的免疫变化具有具体过程差异，不能统一概括为炎症增强。相关敏感性结果及原统计保留在技术附录。')
h('SPP1相关状态与基质环境相联系')
p('肾癌空间数据表明，SPP1相关巨噬状态与成纤维细胞及细胞外基质程序相联系。去除共享及状态定义基因后，SPP1状态与基质程序的相关中位数在TJ-RCC的12名患者中为0.195，在GSE175540的23名患者中为0.164，分别为P=0.006348、FDR=0.006925和P=FDR=0.000488。该结果将SPP1相关细胞状态与其所在的组织环境连接起来。')
p('外部干预研究提供了相应的功能依据。在公开的结直肠癌宿主Spp1敲除模型中，与对照相比，癌相关成纤维细胞数量减少44.6%，CD8 T细胞数量增加53.0%，CD8中增殖标志Ki67阳性比例由30.56%升至50.25%。CAF数量比较为6只对照与7只敲除动物，模型P=0.000485、两模型复核家族FDR=0.000970；Ki67每组5只动物，混合模型估计增加19.59个百分点，P=FDR=0.000986。[3] 这些变化支持SPP1在该模型中参与基质和免疫调节。')
h('组织免疫信号和临床关联呈现情境差异')
p('在44名供者的组织转录组中，SPP1同时与细胞毒和功能障碍相关程序正相关，相关系数分别约为0.671和0.492，P分别约为0.000010和0.000790，FDR分别约为0.000020和0.000790。但在细胞类型限定的分析中，巨噬SPP1与CD8功能程序没有显示一致的抑制关系。组织层面的免疫活跃程度与单一细胞类型的功能状态，呈现了不同的信息。')
p('临床队列中，SPP1相关文献评分与nivolumab治疗后的较短总生存期相关，临床调整后每增加1个标准差的风险比为1.270，P=0.016110、FDR=0.064440，属于名义显著。不同治疗方案、随访终点及C4评分的方向并不完全一致。临床结果为SPP1相关免疫状态提供了结局层面的对应，完整模型和交互结果保留于技术附录。')
h('机制解释与研究意义',1)
h('TYMP相关特征连接免疫信号响应与抗原处理')
p('干扰素γ是免疫细胞之间传递的信号，其下游基因参与细胞防御和抗原处理。MHC-I和MHC-II可理解为细胞向T细胞展示蛋白片段的两套系统，分别主要对应CD8和CD4 T细胞识别；吞噬则涉及细胞对外界颗粒或细胞成分的摄取。[5–7] 因此，TYMP同时关联干扰素响应、MHC-I和吞噬程序，指向了一组免疫信号响应与物质处理相伴随的转录特征。')
p('SPP1与干扰素γ及MHC-II的负向联系，则显示它对应另一种免疫状态。共同表达的TYMP与SPP1，可以分别关联细胞内不同的功能活动。这一差异也解释了为什么不能仅凭SPP1富集，或者TYMP与SPP1正相关，就把整群巨噬细胞视为具有单一功能。')
h('代谢和溶酶体特征为免疫状态提供细胞背景')
p('前文的核苷、多胺及胆固醇处理发现，与本部分的糖酵解、吞噬和抗原处理结果在研究问题上相衔接：它们共同描述了肿瘤代谢背景及巨噬细胞的物质处理和免疫转录状态。糖酵解提供了两基因共同关联的代谢特征，但它没有规定免疫变化的方向；TYMP和SPP1可以在这一共同背景下对应不同的免疫程序。')
p('当前可提出的工作模型是：TYMP相关代谢变化与SPP1相关巨噬状态相联系，巨噬细胞在这一环境中呈现不同的抗原处理和免疫响应特征，SPP1相关状态又与成纤维细胞及基质环境对应。这是由分层证据形成的机制假说。其中，巨噬细胞内的表达联系和空间基质联系来自本项目；SPP1的部分功能作用来自外部干预实验。TYMP经MTA直接调节SPP1，再影响T细胞的连续因果路径尚未建立。')
h('基质重塑为SPP1相关免疫作用提供组织层面的解释')
p('成纤维细胞和细胞外基质构成了免疫细胞活动的组织环境。SPP1相关巨噬与这些成分的联系，使其研究意义超出了巨噬细胞本身。已有实验中，Spp1缺失同时伴随CAF减少和CD8数量、增殖增加，说明基质状态与T细胞变化可以在同一干预背景下出现。[3] 这一证据为SPP1参与肿瘤免疫环境提供了功能依据，也为解释肾癌中的空间关联提供了参考。')
p('本项目已经观察到肿瘤表达区内侧CD8信号较低，但这种内外差没有稳定地由局部SPP1解释；巨噬SPP1与CD8功能程序也未显示一致的抑制关系。因此，当前数据支持SPP1相关基质联系和免疫空间分区，但不支持将两者合并表述为本项目已证实的SPP1特异性物理屏障。这一结果保留了组织关联本身，也界定了机制假说中尚未被直接证明的环节。')
h('结果的适用范围')
p('主要结果来自RNA表达，所描述的是功能相关的转录特征。抗原呈递、吞噬和糖酵解程序的相关变化，不能直接等同于实际抗原呈递能力、吞噬量或代谢通量的改变。去除少数影响较大的基因后，糖酵解效应减弱，尤其是TYMP；总体炎症和干扰素α结果更依赖细胞及基因集合。因而，机制解释应重点采用复核后仍保留的干扰素γ、抗原处理和吞噬相关证据。')
p('外部干预证明的是相应癌种和模型中的SPP1作用，不能代替肾癌中的TYMP干预。巨噬细胞内TYMP关联也不能代替肿瘤来源TYMP的处理效应。临床关联在队列与终点之间存在差异，尚不构成可用于选药的预测工具。总体而言，本研究将TYMP与SPP1的联系从表达共变推进到具体的免疫功能特征，并将这些特征与代谢背景和基质环境对应起来，为TYMP参与SPP1相关过程提供来源明确的生物学线索。')
h('参考文献',1)
refs=[
'[1] MedlinePlus Genetics. TYMP gene. https://medlineplus.gov/genetics/gene/tymp/',
'[2] Qi等. Single-cell and spatial analysis reveal interaction of FAP+ fibroblasts and SPP1+ macrophages in colorectal cancer. Nature Communications. 2022;13:1742. https://doi.org/10.1038/s41467-022-29366-6',
'[3] Henriques等. TGF-β builds a dual immune barrier in colorectal cancer by impairing T cell recruitment and instructing immunosuppressive SPP1+ macrophages. Nature Genetics. 2025;57:3050–3065. https://doi.org/10.1038/s41588-025-02380-2',
'[4] Jiang等. SPP1+ Tumor-Associated Macrophages Drive Immunotherapy Resistance via CD8+ T-cell Dysfunction in Clear-Cell Renal Cell Carcinoma. Cancer Immunology Research. 2025;13:1533–1546. https://doi.org/10.1158/2326-6066.CIR-24-1146',
'[5] Reactome. Expression of IFNG-stimulated genes. https://reactome.org/content/detail/R-HSA-1031716',
'[6] Reactome. Antigen Presentation: Folding, assembly and peptide loading of class I MHC. https://reactome.org/content/detail/R-HSA-983170',
'[7] Reactome. MHC class II antigen presentation. https://reactome.org/content/detail/R-HSA-2132295'
]
for t in refs:p(t)
# Technical appendix preserves every detailed method, original result paragraph and all seven tables.
start=next(i for i,b in enumerate(old) if b.get('text')=='研究方法')
end=next(i for i,b in enumerate(old) if b.get('text')=='讨论')
ref=next(i for i,b in enumerate(old) if b.get('text')=='参考文献')
appendix=[{'type':'title','text':'免疫章节详细方法与统计证据'},{'type':'p','text':'本附录保留全文适配稿v3的详细方法、全部结果数值和七张统计及来源表。与机制叙述版正文配套使用；正文采用独立的参考文献编号，本附录沿用原资料的局部编号。'}]+copy.deepcopy(old[start:end])+copy.deepcopy(old[ref:])
(P/'inputs/MAIN_CONTENT.json').write_text(json.dumps(main,ensure_ascii=False,indent=2))
(P/'inputs/APPENDIX_CONTENT.json').write_text(json.dumps(appendix,ensure_ascii=False,indent=2))
builder=(S/'code/build_docx.py').read_text()
builder=builder.replace("P=Path(__file__).resolve().parents[1];", "P=Path(__file__).resolve().parents[1];")
for key,filename in [('MAIN','26210_免疫章节_机制叙述_v4.docx'),('APPENDIX','26210_免疫章节_详细方法与统计附录_v4.docx')]:
 code=builder.replace("inputs/DOCUMENT_CONTENT.json",f"inputs/{key}_CONTENT.json")
 code=code.replace("26210_免疫章节_全文适配_v3.docx",filename)
 (P/f'code/build_{key.lower()}.py').write_text(code)
shutil.copy2(S/'code/render_docx.py',P/'code/render_docx.py')
sources=[S/'inputs/DOCUMENT_CONTENT.json',S/'delivery/26210_免疫章节_全文适配_v3.docx']
(P/'inputs/SOURCE_MANIFEST.json').write_text(json.dumps([{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sources],ensure_ascii=False,indent=2))
checks={'all_seven_tables_preserved': [b for b in old if b['type']=='table']==[b for b in appendix if b['type']=='table'],'all_detailed_methods_results_preserved':appendix[2:2+end-start]==old[start:end],'main_has_no_statistical_tables':not any(b['type']=='table' for b in main),'main_methods_one_paragraph':main[next(i for i,b in enumerate(main) if b.get('text')=='研究方法')+2]['text']=='主要发现'}
assert all(checks.values()),checks
(P/'qc/CONTENT_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':checks,'no_new_models':True,'main_mechanism_is_interpretation_not_new_experiment':True},indent=2))
(P/'inputs/MECHANISTIC_SOURCES.json').write_text(json.dumps({'references':refs,'scope':'Added general biological explanations and evidence-labelled working model; no newly inferred causal or clinical claim','retrieved':'2026-10-05'},ensure_ascii=False,indent=2))
print('PREPARED',len(main),len(appendix))

