from pathlib import Path
import json,datetime,hashlib,shutil,os,zipfile,html
import pandas as pd
P=Path(__file__).resolve().parents[1];R=P.parent;D=P/'delivery';D.mkdir(exist_ok=True)
assert json.loads((P/'qc/FINAL_VALIDATION.json').read_text())['status']=='PASS'
assert (P/'logs/validation.exit').read_text().strip()=='0'
edges=pd.read_csv(P/'tables/ALL_SEQUENTIAL_EDGES.tsv',sep='\t');summ=pd.read_csv(P/'tables/STEP1_WITHIN_MACRO_SUMMARY.tsv',sep='\t');cnt=pd.read_csv(P/'tables/LINEAGE_COUNTS_BY_PATIENT.tsv',sep='\t');w=pd.read_csv(P/'tables/STEP1_WITHIN_MACRO_PATIENT.tsv',sep='\t')
assert len(edges)==120 and int(cnt.cells.sum())==81535
s0=summ[summ.model.eq('depth_complexity')].iloc[0];s1=summ[summ.model.eq('depth_complexity_environment')].iloc[0]
assert s0.n_patients==4 and s0.positive==4 and s1.positive==4
q=cnt[cnt.dataset.eq('GSE210038')&cnt.main]
cols=['sample','Macrophage','RCC_like','Fibro_strict','Fibro_transition','CD8']
selected=q[cols].rename(columns={'sample':'文库','Macrophage':'巨噬细胞','RCC_like':'肿瘤样细胞','Fibro_strict':'严格成纤维候选','Fibro_transition':'成纤维/壁细胞共标记','CD8':'CD8 T细胞'})
txt=f'''26210｜新单细胞选择与 TYMP–SPP1 免疫屏障主线首轮分析
Phase303 · 2026-10-02

结论
本轮新增下载并处理 GSE210038、GSE224630、GSE222703，共15个肿瘤文库、81,535个质控后细胞；其中 GSE210038 的5例未治疗患者最适合先做患者配对的方向复核。它覆盖肿瘤样细胞、巨噬细胞、带成纤维及壁细胞共标记的基质细胞、CD8 T细胞，但严格成纤维细胞配对仍不足。
新数据再次观察到巨噬细胞内部 TYMP–SPP1 正关联：4/4合格未治疗患者同向，测序因素调整后的中位偏相关为 {s0.median_rho:.6f}，加入干扰素和缺氧因素后为 {s1.median_rho:.6f}。这是新增的小样本方向支持。
同患者三段分析已依次完成。跨患者的肿瘤样细胞 TYMP–巨噬细胞 SPP1、SPP1–胶原成熟、胶原成熟–CD8功能没有形成一致的屏障增强链条。既往 TYMP–SPP1 正面结果继续保留；本轮没有证实肿瘤 TYMP 经 SPP1 驱动基质屏障。

一、选什么数据，为什么
GSE210038：7例肾肿瘤，T1既往Nivolumab+Ipilimumab、T3既往Sunitinib，主分析为T2/T4/T5/T6/T7共5例未治疗患者。质控后44,283细胞，主集合31,655细胞。5例均有至少100巨噬细胞、50肿瘤样细胞、20 CD8 T细胞及20共标记基质细胞，可做分段配对。严格成纤维候选达到20细胞的仅1例，因此不能用它完成纯CAF主线。
GSE224630：5例未治疗患者；冻存RCC4和RCC5血栓是重复/不同部位，不增加患者数。22,147细胞中仅35个巨噬细胞，各患者0–18个；CD8细胞也不足。适合补充肿瘤—基质描述，不能承担完整巨噬—基质—免疫链。
GSE222703：3例患者，15,105个肿瘤细胞。作者注释与本轮标记核对后，成熟巨噬细胞每患者61/116/44个，严格成纤维候选19/3/1个；整链配对不足。可以补充3例肿瘤样TYMP—CD8功能描述。治疗史在GEO样本中未报告；分选富集设计不用于推断原组织细胞比例。
另外，GSE202813为骨转移/骨髓背景，其整合资源还含既往已用来源，不作为本轮原发肾肿瘤主输入。GSE171306的左右肿瘤样本与供者独立性需进一步核清，本轮未纳入。GSE159115属于既有分析来源，不当作新增独立队列。
更匹配基质问题的 Ma 等2025年队列 HRA007082 包含22例未治疗RCC、其中19例ccRCC，发现队列同时含原发瘤和免疫/基质细胞；其原始数据为受控访问，本轮未取得，也未向作者或数据委员会发送申请。公开代码不等于已获得细胞计数矩阵。

GSE210038主集合实际细胞覆盖
{selected.to_string(index=False)}

二、巨噬细胞内 TYMP–SPP1 正关联获得新增方向支持
按患者先控制测序深度与基因检出复杂度，T2/T5/T6/T7共4例均为正关联；进一步控制IFN与缺氧程序后4例仍为正。T4的SPP1仅在9个巨噬细胞检出，未达到患者内关联所需的10个检出细胞，因此不进入此项；它符合患者汇总表达所需的5个检出细胞，可进入另行定义的跨细胞配对分析。两种分析的n不同有明确资格原因。
既往Phase280结果继续成立：C4中15例患者13例正向，中位rho=0.134872，P=0.007385、FDR=0.029541；环境调整15/15正向，中位rho=0.156478，P=0.00006104、FDR=0.0003662。整体成熟巨噬细胞29例中24例正向，中位rho=0.049133，P=0.0005461、FDR=0.002184。以上为既往结果，不与本轮4例直接合并或重新校正。
与此同时，本轮5例患者的“肿瘤样细胞TYMP—巨噬细胞SPP1”rho=-0.20，当前数据不支持预期的正向跨细胞联系。细胞内部关联与肿瘤来源关联是两个问题，不能相互替代。

三、SPP1与部分基质装置同向，但胶原成熟没有接上
在GSE210038的5例同患者共标记基质细胞中，巨噬SPP1与黏着斑程序rho=0.60、与收缩装置程序rho=0.50，给出部分基质过程的方向线索；与胶原成熟程序rho=-0.30。
肿瘤样TYMP与这三个程序分别为：胶原成熟-0.70、黏着斑-0.60、收缩装置-0.90；巨噬TYMP分别为-0.50、0.00、-0.60。不能只保留SPP1对应的两个正方向并宣布整条基质链成立。
GSE224630的5例肿瘤样TYMP对应基质程序分别为0.00、-0.80、-0.80，亦未提供一致正向补充。严格成纤维接收池在各队列均没有至少3例配对，状态为“不可评估”，不是生物学阴性。壁细胞池完整结果另列统计表；它与共标记池有细胞重叠，不能算独立重复。

四、CD8功能端未显示一致的屏障后果
GSE210038的5例中，巨噬SPP1与CD8功能障碍程序rho=0.60，但与CD8细胞毒程序rho=0.30，并未伴随预期的细胞毒程序下降。巨噬TYMP对应两项为0.60和-0.30，肿瘤样TYMP对应两项为0.00和-0.30。
共标记基质胶原成熟与CD8细胞毒程序rho=0.20、与功能障碍程序rho=0.10，未提供清晰的免疫功能受损方向。GSE222703的3例肿瘤样TYMP与CD8细胞毒/功能障碍分别为-0.50/0.50，作为小样本描述保留。
与“免疫屏障”的直接关系是：若TYMP相关巨噬状态参与屏障，应在同患者观察到基质执行程序变化，并伴随抗肿瘤免疫功能受损。本轮补强了最上游的巨噬细胞内联系，尚未把中间基质和末端免疫功能连接成一致证据。

五、方法
以患者为独立单位，保留原始下载、作者注释、计数对象及全部日志。单细胞QC使用UMI>=500、检出基因200–8000、线粒体比例<=20%，每文库Scanpy Scrublet自动阈值、期望doublet率0.06、随机种子303。归一化CP10K后log1p，以每研究2000个HVG、30主成分、15邻居、Leiden分群；TYMP/SPP1不参与HVG身份分群。
GSE222703使用作者原始计数与细胞注释，另外两研究使用标记核查后的簇身份。成纤维候选需至少2个core及1个support基因检出，排除多重免疫/内皮/上皮标记；有周细胞/平滑肌共标记者单列。CD8依据T细胞及CD8谱系标记或作者CD8T注释确定，均不以TYMP/SPP1结果筛细胞。
每患者每类细胞合并原始计数、转为log1p(CP10K)，程序为冻结基因面板等权均值。单基因要求原文库存在该基因且至少5细胞检出；程序要求原文库测量不少于80%面板基因，至少50%面板基因在该类细胞中检出率达到5%。身份门槛分别为巨噬100、肿瘤样50、严格成纤维/共标记基质20、壁细胞50、CD8 20。
患者内TYMP–SPP1为秩残差相关，依次控制深度/复杂度及干扰素/缺氧程序；各研究分别汇总方向及中位效应。跨患者各段采用同文库配对Spearman相关，至少3例且两端有变化才报告rho。本轮各研究有效n均少于10，按项目既定规则仅报效应与方向，不计算患者推断P/FDR，也不使用细胞级P代替患者证据。共预设120个研究×端点组合，所有有效和不可评估结果均交付。
胶原成熟9基因、黏着斑9基因、收缩装置8基因沿用Phase292冻结面板；CD8细胞毒5基因和功能障碍8基因沿用Phase213。详细基因、资格与逐基因读出见随附资源和统计表。

六、讨论与后续定位
新增数据适合“小样本方向复核”，尚不足以承担纯CAF屏障的完整验证。肿瘤样细胞身份来自转录标记，本轮没有新增CNV验证；成纤维/壁细胞共标记不等于已经证实谱系转分化；程序表达反映相关转录活动，不等于实测胶原排列、细胞迁移或物理阻挡。单细胞数据用于细胞来源、状态与功能联系，空间排斥仍需独立空间读出。
本轮4例正方向来自既往Phase280五研究以外的新来源，但GSE210038与已用GLMF空间数据属于同一原始研究，不能再称其为独立跨组学验证。5例跨患者相关易受单例和临床异质性影响，不称统计显著，不据此推断因果或临床价值。
本轮分析结论应写为：TYMP与SPP1巨噬状态的细胞内联系获得新增一致方向支持；肿瘤来源、纯成纤维接收端和免疫屏障功能后果仍需更多匹配患者数据。不能因后两段未接通而删去上游已成立的正面证据，也不能用上游阳性代替整链证明。
下一次扩样应优先取得同患者、未治疗原发ccRCC、成熟巨噬与明确成纤维/CD8均有足够测量的完整对象；HRA007082是适配但受控的候选。本轮不通过放宽细胞身份、混入骨转移或把技术样本当患者来补n；不重复对当前5例换阈值寻找正向结果。

七、验收与交付
分析按TYMP–SPP1→基质过程→免疫功能执行，三个完成标记及全部正式程序退出码已核验。独立验收从原文库重算表达、逐项配对和全部相关，1273项通过。下载归档、SHA256、原始输入、脚本、失败日志及修复记录全部保存在服务器Phase303。原C4、旧RDS/H5AD、既有统计和客户报告保持。
完整统计表在tables/，冻结规则在resources/，必要验收在qc/；空P/FDR代表本轮仅描述或不可评估，不能解释成P=1。

官方数据与原始研究
GSE210038：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210038
GSE224630：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224630；原文 https://pubmed.ncbi.nlm.nih.gov/36801089/
GSE222703：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222703
GSE202813：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202813
GSE171306：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171306
Ma等，FAP+ fibroblasts orchestrate tumor microenvironment remodeling in renal cell carcinoma with tumor thrombus：https://www.nature.com/articles/s41467-025-64447-2 （Data availability明确受控访问）
'''
(D/'26210_Phase303_单细胞筛选与主线首轮结果.txt').write_text(txt,encoding='utf-8-sig')
# Readable standalone report; no new biological figures.
body=''.join('<p>'+html.escape(p).replace('\n','<br>')+'</p>' for p in txt.split('\n\n'))
ht='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>26210 · 单细胞主线首轮结果</title><style>body{max-width:1050px;margin:40px auto;padding:0 24px;font:16px/1.85 "Microsoft YaHei",sans-serif;color:#20303c}p{margin:1.3em 0}table{border-collapse:collapse;font-size:13px}td,th{padding:6px;border:1px solid #bbc8d0}a{color:#165b91}@media print{body{margin:0;font-size:11pt}}</style><body>'+body+'<hr><p>附：GSE210038 主分析患者细胞覆盖</p>'+selected.to_html(index=False)+'<p>完整结果：<a href="tables/ALL_SEQUENTIAL_EDGES.tsv">所有120项端点</a>；<a href="tables/ALL_MATCHED_PATIENT_PAIRS.tsv">患者配对数值</a>。</p></body></html>'
(D/'26210_Phase303_单细胞筛选与主线首轮结果.html').write_text(ht,encoding='utf-8')
for folder in ['tables','resources','qc','code']:(D/folder).mkdir(exist_ok=True)
table_names=['SAMPLE_PLAN.tsv','QC_BY_LIBRARY.tsv','LINEAGE_COUNTS_BY_PATIENT.tsv','PATIENT_LINEAGE_MEASURES.tsv','FUNCTION_GENE_COVERAGE.tsv','LINEAGE_IDENTITY_GENE_COVERAGE.tsv','STEP1_WITHIN_MACRO_PATIENT.tsv','STEP1_WITHIN_MACRO_SUMMARY.tsv','STEP1_TUMOR_MACRO.tsv','STEP2_STROMA.tsv','STEP3_IMMUNE.tsv','ALL_SEQUENTIAL_EDGES.tsv','ALL_MATCHED_PATIENT_PAIRS.tsv']
for n in table_names:shutil.copy2(P/'tables'/n,D/'tables'/n)
for f in (P/'resources').glob('*.json'):shutil.copy2(f,D/'resources'/f.name)
for n in ['FINAL_VALIDATION.json','VALIDATION_CHECKS.json','STEP1_COMPLETE.json','STEP2_COMPLETE.json','STEP3_COMPLETE.json']:shutil.copy2(P/'qc'/n,D/'qc'/n)
for f in (P/'code').glob('*.py'):shutil.copy2(f,D/'code'/f.name)
ledger=[{'accession':'GSE210038','status':'Selected direction replication','main_patients':5,'cells_total':44283,'limitation':'Only one patient with >=20 strict fibro candidates; stromal co-marker pool distinct from pure CAF'},{'accession':'GSE224630','status':'Tumour-stroma secondary only','main_patients':5,'cells_total':22147,'limitation':'Only 35 macrophages and <20 CD8 per patient'},{'accession':'GSE222703','status':'Tumour-CD8 secondary only','main_patients':3,'cells_total':15105,'limitation':'Insufficient mature macrophage-fibro pairing; therapy unknown; sorted enrichment'},{'accession':'GSE202813','status':'Not selected primary','limitation':'Bone metastasis/marrow context; integrated resources overlap prior studies'},{'accession':'GSE159115','status':'Already used in Phase280','limitation':'Not new independent source'},{'accession':'GSE171306','status':'Metadata qualification unresolved; not acquired','limitation':'Left/right tumour samples; independent donor mapping unresolved'},{'accession':'HRA007082','status':'Promising, controlled-access; not acquired','limitation':'Original study includes 19 ccRCC donors among 22 RCC; no raw count object obtained; no access request sent'}]
(P/'resources/SOURCE_SELECTION_FINAL.json').write_text(json.dumps(ledger,ensure_ascii=False,indent=2));shutil.copy2(P/'resources/SOURCE_SELECTION_FINAL.json',D/'resources/SOURCE_SELECTION_FINAL.json')
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
manifest={str(f.relative_to(D)):sha(f) for f in D.rglob('*') if f.is_file() and f.name!='SHA256.json'}
(D/'SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
zip_path=P/'26210_Phase303_单细胞筛选与主线首轮结果.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
 for f in D.rglob('*'):
  if f.is_file():z.write(f,f.relative_to(D))
(P/'qc/DELIVERY_VALIDATION.json').write_text(json.dumps({'files':len(manifest),'zip_sha256':sha(zip_path),'zip_bytes':zip_path.stat().st_size,'report_edges':len(edges),'checks':1273},indent=2))
print('DELIVERY_READY',zip_path,zip_path.stat().st_size,sha(zip_path))
