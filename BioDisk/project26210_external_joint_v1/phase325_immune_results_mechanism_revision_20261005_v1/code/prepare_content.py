from pathlib import Path
import copy,json,shutil,hashlib,re
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
S=R/'phase323_immune_manuscript_adaptation_20261005_v1'
P=R/'phase325_immune_results_mechanism_revision_20261005_v1'
P.mkdir(exist_ok=False)
for d in ['code','inputs','delivery','qc','logs','render','root_record_backups']:(P/d).mkdir()
for n in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:shutil.copy2(R/n,P/'root_record_backups'/n)
old=json.loads((S/'inputs/DOCUMENT_CONTENT.json').read_text())
start=next(i for i,b in enumerate(old) if b.get('text')=='分析结果')
end=next(i for i,b in enumerate(old) if b.get('text')=='讨论')
bounds=[i for i in range(start+1,end) if old[i]['type']=='h']+[end]
add=[
[
'上述联系把TYMP的髓系作用线索从“与SPP1共同表达”推进到具体的细胞过程。CTSB编码的组织蛋白酶B参与溶酶体内蛋白降解，TGFB1编码的TGF-β1参与细胞间信号调节。[10,11] TYMP与这些基因的关联，分别对应细胞内物质处理和细胞间调节两个方面，与前文C4的溶酶体及脂质处理特征形成了功能上的衔接。',
'结合这些结果，TYMP相关巨噬细胞的生物学特征更适合从SPP1状态、溶酶体处理和调节信号共同理解。这里的机制线索来自基因表达的协同变化；TGFB1和CTSB是联系TYMP与巨噬功能过程的具体候选，并非本项目已经证实的分泌效应或TYMP下游介质。'
],
[
'干扰素γ反应反映细胞对免疫信号的转录响应，其中包含多种防御及抗原处理相关基因。MHC-I和MHC-II则是细胞向T细胞展示蛋白片段的两套系统，分别主要供CD8和CD4 T细胞识别。[12–14] 因而，TYMP与干扰素γ、MHC-I及吞噬程序同向，将其相关状态具体定位到免疫信号响应、抗原处理和物质摄取等过程。',
'SPP1与干扰素γ及MHC-II的负向联系，说明共同表达的两基因并不对应完全相同的功能活动。TYMP较高的细胞同时呈现较强的免疫响应及处理相关特征，SPP1较高的细胞则呈现不同的抗原呈递转录状态。这种差异存在于巨噬细胞内部，可帮助解释为何TYMP与SPP1共同升高，却不能被统一贴上“免疫激活”或“免疫抑制”的标签。',
'相互控制和去共享基因后仍保留这些关系，说明上述差异不完全来自另一基因的共同变化，也不完全来自重复计算相同程序成员。这里观察到的是功能相关的转录特征，而非干扰素分泌量、抗原呈递能力或吞噬量的直接测定。'
],
[
'糖酵解是葡萄糖分解并形成ATP及代谢中间产物的过程。[15] 两基因与糖酵解程序的共同联系，提供了理解其免疫差异的代谢背景：巨噬细胞可在相近的代谢转录背景下，呈现不同的干扰素响应和抗原处理特征。因此，共同的糖酵解联系并不要求TYMP和SPP1具有相同的免疫功能方向。',
'这也为前文TYMP相关核苷、蛋氨酸和多胺代谢发现提供了免疫层面的对应，但尚不能将不同数据层的关联直接串成代谢反应链。少数基因删除后效应衰减，提示糖酵解联系部分集中在特定程序成员；总体炎症与TNFα/NFκB的方向又不相同。机制解释应落实到具体代谢与免疫过程，而不将其概括为全局代谢增强或普遍炎症激活。'
],
[
'组织表达是多种细胞共同贡献的信号。SPP1较高的组织同时出现较高的细胞毒和功能障碍相关程序，说明两类转录信号可以在同一组织中共存；它们未必来自同一批T细胞，也不能仅凭组织均值判断某个细胞已经失去功能。',
'限定巨噬来源与CD8接收端后，未观察到一致的抑制关系。因此，组织SPP1关联更适合作为肿瘤免疫环境的整体特征，不能直接替代“巨噬SPP1抑制CD8”的细胞间作用证据。这个层级差异解释了组织免疫信号升高与细胞类型限定结果并不一致的原因范围，并保留了两类观察各自的信息。'
],
[
'这些结果将SPP1相关状态与细胞外环境连接起来。成纤维细胞参与基质形成，基质的组成和结构影响免疫细胞所处的组织环境。肾癌中的程序共分布，以及外部敲除模型中CAF减少、CD8数量和增殖增加，分别提供了组织关联和干预层面的支持。[2,3] 因此，SPP1相关免疫作用可从巨噬—成纤维细胞—基质相互联系的角度理解。',
'这一路径的功能意义不只在于T细胞数量，还涉及T细胞能否到达相应区域、接触肿瘤细胞并维持增殖。已有文献为这些问题提供了机制依据；本项目的空间数据目前直接支持SPP1与基质程序联系，以及CD8分布不均。局部SPP1与CD8内外差未形成稳定联系，意味着从基质关联到特异性屏障之间的环节仍应按其真实证据层级解释。'
],
[
'临床关联把细胞和组织层面的发现延伸到治疗结局。SPP1文献评分代表一组SPP1相关巨噬特征；C4评分描述本项目定义的细胞状态，两者不是同一个指标，也不是细胞数量的直接测量。它们对应不同结局或时间模式，说明不同髓系状态在临床上不能被简化为同一种“耐药巨噬细胞”。',
'CM025中的治疗交互为SPP1相关状态与相对治疗优势之间的联系提供了探索方向；但这种关联不能单独解释其通过抗原呈递、基质还是其他过程影响治疗。结合前述结果，临床证据更适合支持髓系免疫状态具有结局相关性，具体作用过程由细胞、空间及干预证据共同解释。'
]]
blocks=copy.deepcopy(old[:start+1])
for j in range(6):
 section=copy.deepcopy(old[bounds[j]:bounds[j+1]])
 # Keep scientific observations and interpretations together, then present their numeric evidence.
 blocks.extend([b for b in section if b['type']!='table'])
 blocks.extend({'type':'p','text':x} for x in add[j])
 blocks.extend([b for b in section if b['type']=='table'])
blocks.extend(copy.deepcopy(old[end:]))
refs=[
'[10] NCBI Gene. CTSB cathepsin B. https://www.ncbi.nlm.nih.gov/gene/1508',
'[11] NCBI Gene. TGFB1 transforming growth factor beta 1. https://www.ncbi.nlm.nih.gov/gene/7040',
'[12] Reactome. Expression of IFNG-stimulated genes. https://reactome.org/content/detail/R-HSA-1031716',
'[13] Reactome. Antigen Presentation: Folding, assembly and peptide loading of class I MHC. https://reactome.org/content/detail/R-HSA-983170',
'[14] Reactome. MHC class II antigen presentation. https://reactome.org/content/detail/R-HSA-2132295',
'[15] Reactome. Glycolysis. https://reactome.org/content/detail/R-HSA-70171']
blocks.extend({'type':'p','text':x} for x in refs)
(P/'inputs/DOCUMENT_CONTENT.json').write_text(json.dumps(blocks,ensure_ascii=False,indent=2))
builder=(S/'code/build_docx.py').read_text().replace('26210_免疫章节_全文适配_v3.docx','26210_免疫章节_结果机制解释_v5.docx')
(P/'code/build_docx.py').write_text(builder)
shutil.copy2(S/'code/render_docx.py',P/'code/render_docx.py')
checks={
'full_methods_exactly_unchanged':blocks[:start]==old[:start],
'all_original_result_paragraphs_preserved':all(b in blocks for b in old[start:end] if b['type']=='p'),
'all_seven_tables_exact':all(b in blocks for b in old if b['type']=='table') and sum(b['type']=='table' for b in blocks)==7,
'discussion_unchanged':blocks[blocks.index(old[end]):blocks.index(old[end])+len(old[end:])]==old[end:],
'six_results_receive_explanation':len(add)==6 and all(len(x)>=2 for x in add),
'new_interpretation_paragraphs':sum(map(len,add))==13
}
assert all(checks.values()),checks
(P/'qc/CONTENT_VALIDATION.json').write_text(json.dumps({'status':'PASS','checks':checks,'no_new_models':True,'rejected_draft':'Phase324 compressed-method draft not delivered or adopted','scope':'Only results expanded with mechanisms; detailed methods remain in same document'},ensure_ascii=False,indent=2))
(P/'inputs/SOURCE_MANIFEST.json').write_text(json.dumps([{'path':str(S/'inputs/DOCUMENT_CONTENT.json'),'sha256':hashlib.sha256((S/'inputs/DOCUMENT_CONTENT.json').read_bytes()).hexdigest()}],indent=2))
(P/'inputs/MECHANISM_SOURCE_NOTES.json').write_text(json.dumps({'new_references':refs,'new_text':add,'evidence_type':'functional explanation and explicitly labelled interpretation, no new causal finding'},ensure_ascii=False,indent=2))
print('PREPARED',len(blocks),'13 interpretation paragraphs; complete methods retained')

