from pathlib import Path
import json,hashlib,datetime,shutil,zipfile,html,os,csv,io
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent;D=P/'delivery';D.mkdir(exist_ok=True)
v=json.loads((P/'qc/FINAL_ASSOCIATION_VALIDATION.json').read_text());assert v['status']=='PASS';assert (P/'logs/validate_associations.exit').read_text().strip()=='0'
w=pd.read_csv(P/'tables/WITHIN_MACRO_SUMMARY.tsv',sep='\t');e=pd.read_csv(P/'tables/ALL_PATIENT_ASSOCIATIONS.tsv',sep='\t');counts=pd.read_csv(P/'tables/LINEAGE_COUNTS_BY_PATIENT.tsv',sep='\t');individual=pd.read_csv(P/'tables/WITHIN_MACRO_PATIENT_EFFECTS.tsv',sep='\t');oq=pd.read_csv(P/'tables/OBRADOVIC_QC_BY_LIBRARY.tsv',sep='\t')
core=w[w.dataset.eq('OBRADOVIC')&~w.model.str.endswith('hypoxia')].iloc[0];env=w[w.dataset.eq('OBRADOVIC')&w.model.str.endswith('hypoxia')].iloc[0];assert core.patients==11 and core.positive==10 and env.positive==10
main=e[e.dataset.eq('OBRADOVIC')&e.model.eq('prespecified_adjusted')];immune=main[main.stage.eq(3)&main.source.isin(['Macrophage:TYMP','Macrophage:SPP1'])];assert len(immune)==4 and immune.P.ge(.05).all()
def row(src,tgt):return main[main.source.eq(src)&main.target.eq(tgt)].iloc[0]
up=row('Tumour_like:TYMP','Macrophage:SPP1');mature=sum(counts[counts.dataset.eq('OBRADOVIC')].Macrophage);cd8=sum(counts[counts.dataset.eq('OBRADOVIC')].CD8)
fmt=lambda x:'—' if pd.isna(x) else f'{x:.6f}'
it=[]
for r in immune.itertuples():it.append(f'{r.source} → {r.target}：n={r.patients}，rho={r.rho:.6f}，P={r.P:.6f}，FDR={r.FDR:.6f}。')
paired=counts[counts.dataset.eq('OBRADOVIC')][['patient','Macrophage','CD8','Tumour_like','Fibro_strict_candidate','Stromal_comarker']]
st=[]
for src in ['Tumour_like:TYMP','Macrophage:TYMP','Macrophage:SPP1']:
 rr=[row(src,'Stromal_comarker:'+ep) for ep in ['collagen_maturation','focal_adhesion','actomyosin']];st.append(src+'：胶原成熟 '+fmt(rr[0].rho)+'；黏着斑 '+fmt(rr[1].rho)+'；收缩装置 '+fmt(rr[2].rho)+'（均n=8，仅描述）。')
txt=f'''26210 免疫单细胞队列验证结果
Phase305　2026年10月2日

结论
新增Obradovic未治疗ccRCC免疫富集队列显示：11位独立患者中10位的成熟巨噬细胞内TYMP与SPP1呈正关联。控制测序深度、复杂度、文库和髓系亚群后，中位rho={core.median_rho:.6f}，双侧方向检验P={core.P:.8f}、FDR={core.FDR:.8f}；进一步控制干扰素与缺氧程序后，中位rho={env.median_rho:.6f}，10/11同向，P={env.P:.8f}、FDR={env.FDR:.8f}。两模型在预设的2项检验家族内均达到FDR<0.05。
这是对TYMP–SPP1巨噬细胞状态联系的新增独立队列支持。此前Phase280的C4和整体巨噬细胞正面结果继续保留；本轮没有将新11例与旧29例直接合并重检验。
下游患者级分析中，巨噬TYMP/SPP1与CD8功能的四个主要关联均未达到原始P<0.05。肿瘤样TYMP到巨噬SPP1、纯成纤维接收端及基质到免疫功能未形成一致连接，因此本轮支持上游细胞内联系，尚未证实完整的TYMP驱动免疫屏障机制。

一　免疫富集队列使主验证扩大到11位合格患者
Obradovic队列：官方原始矩阵归档2,289,775,883字节，按用户要求直连后台下载，官方SHA256与ZIP CRC通过。使用11位未治疗ccRCC患者的13个肿瘤CD45+文库，另保留8位患者的8个肿瘤CD45−文库作肿瘤及基质来源配套。B/C各有两个免疫文库，均合并到各自患者，不增加n。
质控后共有164,610细胞，其中CD45+分选文库{int(oq[oq.selection.eq('CD45+')].singlets.sum()):,}细胞。经髓系细分与谱系标记复核，主分析保留{mature:,}个成熟巨噬细胞及{cd8:,}个CD8 T细胞；11位患者全部达到相应细胞数及TYMP/SPP1检出门槛。
HCRN GU16-260：用户提供的Data_S1.zip为463,255,827字节，本地和服务器SHA256均与官方一致；外层和19个内层ZIP均通过CRC，矩阵维度、唯一条码及基因测量核对通过。作者临床表逐一匹配后，17位单细胞患者中13位为ccRCC、3位乳头状、1位未分类；本轮主对象只纳入13位ccRCC患者的14个文库，质控后60,911细胞。
HCRN的ccRCC部分包含8份治疗前和6份进展时取材，260-1095前后配对属于同一位患者。文件名PR代表进展取材时点，临床疗效字段中的PR才代表部分缓解，两者分别记录。实际达到成熟巨噬>=100细胞的只有治疗前1例、进展时5例；因此不能把13位总患者直接写成免疫关联分析n=13，也不因CAF不足而排除符合免疫问题资格的样本。

Obradovic患者细胞覆盖
{paired.to_string(index=False)}

二　TYMP–SPP1正关联在独立未治疗队列得到验证
主要模型：n=11，10正1负，中位rho={core.median_rho:.6f}，范围{core.min_rho:.6f}至{core.max_rho:.6f}；P={core.P:.8f}，FDR={core.FDR:.8f}。
环境调整模型：n=11，10正1负，中位rho={env.median_rho:.6f}，范围{env.min_rho:.6f}至{env.max_rho:.6f}；P={env.P:.8f}，FDR={env.FDR:.8f}。同一数据的两种调整用于稳健性评价，不算两次独立复制。PatientB在两个模型中均为负，保留在完整患者表内。
HCRN治疗前仅1例可评估，两个模型均为正，不做推断检验。进展时5例中，基础模型3正2负，中位rho=0.045129；环境调整后5/5正，中位rho=0.046933，作为小样本方向补充，P/FDR留空。按供者最早可得时点保留一份样本的敏感性汇总也单列在表中，不重复计算患者。
对于肿瘤来源，Obradovic的8位CD45−配套患者中，肿瘤样细胞TYMP与巨噬SPP1的调整后rho={up.rho:.6f}，没有形成明确正向联系。细胞内TYMP–SPP1联系不能替代肿瘤细胞来源的跨细胞证据。

三　巨噬TYMP与部分基质过程同向 但SPP1未形成一致基质连接
Obradovic有8位患者具备分选后的基质配套数据。主要基质细胞具有成纤维与周细胞/平滑肌共标记，单列为Stromal_comarker；严格成纤维候选只有1位患者达到20细胞，不能用于跨患者纯CAF验证。
控制作者pT类别后，各来源与共标记基质三程序的效应为：
{chr(10).join(st)}
其中巨噬TYMP与黏着斑程序rho=0.736321，提供基质过程的探索性方向线索；巨噬SPP1与胶原成熟同向，但与黏着斑、收缩装置反向，不能将这些结果串成一致的SPP1基质执行链。共标记基质与壁细胞池相互重叠，完整壁细胞结果保留在统计表中，不视为独立重复。

四　11位患者的CD8功能分析未达到名义显著
主要模型控制分选文库的技术批次与作者pT1b/pT3a类别，得到：
{chr(10).join(it)}
因此，当前数据不支持巨噬TYMP/SPP1表达与上述CD8功能读出存在显著患者级关联。未调整及调整结果均完整保留；例如巨噬TYMP与CD8功能障碍的未调整rho=0.554545、P=0.076652，调整后降为0.132997、P=0.733023，未达到预设原始P<0.05筛选标准。
在8例基质配套患者中，共标记基质胶原成熟与CD8细胞毒程序rho=-0.228916、与功能障碍程序rho=-0.380235，仅作描述。HCRN的治疗前/进展时功能联系分别报告，因有效患者不足不作推断。这些结果尚未显示从基质变化到免疫功能受损的一致后果。

五　与免疫屏障主线的关系
此次扩大免疫队列后，明确补强的是“TYMP与SPP1相关巨噬细胞状态相联系”这一上游环节。免疫屏障主线还需要同患者的基质执行变化及抗肿瘤免疫受损相连接。本轮已经依次检验这些关系，下游没有得到同等强度的支持。
建议将本轮主结论写为：“独立免疫富集队列确认TYMP与SPP1在成熟巨噬细胞内的正向表达联系；其与基质重塑和CD8功能之间的患者级连接尚未建立。”保留已成立的上游结果，同时保留会改变整链判断的关键下游结果。

六　方法
所有样本保留原始计数及来源映射。QC使用UMI>=500、检出基因200–8000、线粒体比例<=20%；每文库采用Scanpy原生Scrublet、期望doublet率0.06、种子305。使用CP10K/log1p、2000个HVG、30主成分、15邻居及Leiden分群；TYMP/SPP1不用于身份HVG选择。髓系另行细分，将成熟巨噬、单核、树突、骨吸收相关及明显混合群分开。
成熟巨噬需属于经标记复核的髓系亚群，至少检出2个巨噬支持标记，并排除多重T/上皮/内皮标记；CD8需属于T/NK候选群，至少检出2个CD3D/CD3E/TRAC及1个CD8A/CD8B，排除多重髓系/上皮/内皮标记。Obradovic免疫主分析仅使用CD45+分选，肿瘤与基质配套使用CD45−分选。身份规则在目标关联计算前冻结，旧C4及既有对象不改。
患者内TYMP–SPP1计算秩残差相关，控制UMI、基因数、文库及髓系亚群；环境模型加入既有IFN及缺氧面板，移除TYMP/SPP1。每患者>=100巨噬、两目标各>=10细胞检出、剩余维度>=8。独立患者n>=10时用双侧精确二项方向检验；同队列/同层内两模型构成2项BH家族。只以患者为推断单位，不使用细胞级P。
跨细胞室读出按患者和治疗层合并原始计数，再计算log1p(CP10K)及固定基因面板均值。每来源/接收池保持原文库基因测量资格；单基因至少5细胞检出，程序至少80%基因测得且至少50%面板基因在该类细胞中检出率>=5%。阈值为巨噬100、CD8 20、肿瘤样50、严格成纤维候选/共标记基质20、壁细胞50。
跨患者报告未调整及预设调整后的Spearman/偏Spearman相关。Obradovic控制技术批次与作者pT类别，HCRN按治疗时点分层。n>=10且自由度>=5时使用双侧t近似P；n<10仅报效应。多重校正在同队列、同治疗层、同来源与接收池的同类终点及两模型内进行：基质3×2，免疫2×2，上游1×2。保留原P和FDR，以原P<0.05筛选名义信号，FDR<0.05才称校正后显著。

七　讨论
本轮患者内正关联的中位效应较小，但方向一致性通过患者级检验，不能写成强效应或因果调控。排序后的原始表达关系仍不同于蛋白功能、细胞间信号方向或扰动效应。肿瘤样细胞本轮未新增CNV确认；基质共标记不等于已证明纯CAF或谱系转分化；RNA程序也不直接测量胶原排列、T细胞迁移或物理排斥。
FACS富集队列用于细胞状态与患者配对表达，不能将分选后的细胞比例当作原组织丰度。HCRN的治疗与取材背景不同于未治疗原发肾肿瘤；治疗前1例、进展时5例的可评估规模只能提供方向补充。原有29例整体巨噬证据及15例C4证据保持独立呈现，本轮没有拼成40例同质队列，也没有重选阈值追求阳性。

八　验收及文件
下载和两队列预处理正式退出0。首次关联脚本在HCRN零巨噬样本遇到空数组问题，增加长度判断后重跑，阈值/模型/检验家族未变，原失败代码与日志完整存档。
从逐文库原始QC计数独立重算表达和目标细胞向量，并以独立正交投影复核偏相关、组合数复核方向检验、SciPy复核BH；{v['checks']}项检查全部通过。统计表包含全部240个跨患者模型组合及8项患者内汇总，涵盖可评估、描述性与不可评估结果。
完整患者表、逐基因覆盖、原P/FDR和检验家族位于tables；冻结规则位于resources；验收位于qc。原始矩阵与细胞对象继续保存在服务器Phase305，不随轻量结果包重复复制。无新增生物学图件，不覆盖既有客户Word与正式文章。

来源
Obradovic官方数据：https://data.mendeley.com/datasets/nc9bc8dn4m/1
作者代码及患者注释：https://github.com/Aleksobrad/single-cell-rcc-pipeline
HCRN官方数据：https://datadryad.org/dataset/doi:10.5061/dryad.dv41ns28r
HCRN原文：https://pubmed.ncbi.nlm.nih.gov/39992403/
HCRN临床表为出版社Figshare Table S1及Table S5，下载入口及MD5保存在metadata/figshare_28918534.json和figshare_28918492.json。
'''
stem='26210_Phase305_免疫队列验证结果';(D/(stem+'.txt')).write_text(txt,encoding='utf-8-sig');body=''.join('<p>'+html.escape(s).replace('\n','<br>')+'</p>' for s in txt.split('\n\n'));(D/(stem+'.html')).write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>26210 免疫队列验证结果</title><style>body{max-width:1000px;margin:40px auto;padding:0 24px;font:16px/1.8 "Microsoft YaHei",sans-serif;color:#20303c}p{margin:1.3em 0}table{border-collapse:collapse;font-size:13px}td,th{padding:6px;border:1px solid #bbc8d0}@media print{body{margin:0;font-size:11pt}}</style><body>'+body+'<hr><p>患者覆盖表</p>'+paired.to_html(index=False)+'</body></html>',encoding='utf-8')
for folder in ['tables','resources','qc','code']:(D/folder).mkdir(exist_ok=True)
for n in ['OBRADOVIC_SAMPLE_PLAN.tsv','HCRN_SAMPLE_PLAN.tsv','OBRADOVIC_QC_BY_LIBRARY.tsv','HCRN_QC_BY_LIBRARY.tsv','LINEAGE_COUNTS_BY_PATIENT.tsv','PATIENT_MEASURES.tsv','PROGRAM_GENE_COVERAGE.tsv','WITHIN_MACRO_PATIENT_EFFECTS.tsv','WITHIN_MACRO_SUMMARY.tsv','ALL_PATIENT_ASSOCIATIONS.tsv','ASSOCIATION_PATIENT_PAIRS.tsv','STEP1_PATIENT_ASSOCIATIONS.tsv','STEP2_PATIENT_ASSOCIATIONS.tsv','STEP3_PATIENT_ASSOCIATIONS.tsv']:shutil.copy2(P/'tables'/n,D/'tables'/n)
for n in ['LINEAGE_IDENTITY_FROZEN.json','PATIENT_ASSOCIATION_PROTOCOL.json','MYELOID_REFINEMENT_FROZEN.json','OBRADOVIC_QC_PROTOCOL.json','HCRN_QC_PROTOCOL.json']:shutil.copy2(P/'resources'/n,D/'resources'/n)
for n in ['FINAL_ASSOCIATION_VALIDATION.json','ASSOCIATION_VALIDATION_CHECKS.json','HCRN_ARCHIVE_INPUT_VALIDATION.json','HCRN_LOCAL_TRANSFER_VERIFIED.json','OBRADOVIC_DOWNLOAD_VERIFIED.json','EMPTY_INPUT_REPAIR.json']:shutil.copy2(P/'qc'/n,D/'qc'/n)
for n in ['refine_myeloid.py','qualify_patient_inputs.py','run_patient_associations.py','validate_associations.py','deliver_results.py']:shutil.copy2(P/'code'/n,D/'code'/n)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
manifest={str(f.relative_to(D)):sha(f) for f in D.rglob('*') if f.is_file() and f.name!='SHA256.json'};(D/'SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));archive=P/(stem+'.zip')
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in D.rglob('*'):
  if f.is_file():z.write(f,f.relative_to(D))
(P/'qc/DELIVERY_VALIDATED.json').write_text(json.dumps({'status':'PASS','zip_bytes':archive.stat().st_size,'zip_sha256':sha(archive),'files':len(manifest),'statistical_checks':v['checks']},indent=2));print('DELIVERY_READY',archive.stat().st_size,sha(archive))
