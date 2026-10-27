from pathlib import Path
import json,datetime,hashlib,zipfile,shutil,html,os,csv,io,importlib.metadata
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent;D=P/'delivery';now=datetime.datetime.now().astimezone();stamp=now.isoformat();day=now.strftime('%Y-%m-%d')
v=json.loads((P/'qc/PHASE307_FINAL_VALIDATION.json').read_text());assert v['status']=='PASS';assert (P/'logs/final_verify.exit').read_text().strip()=='0'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
dl=json.loads((P/'qc/OT1_DOWNLOAD_COMPLETE.json').read_text());e=pd.read_csv(P/'tables/CD44_PSEUDOBULK_ASSOCIATIONS_V2.tsv',sep='\t');loc=pd.read_csv(P/'tables/GUIDE_LOCALIZATION_EFFECTS.tsv',sep='\t');lq=pd.read_csv(P/'tables/GUIDE_LOCALIZATION_QUALIFICATION.tsv',sep='\t');spq=pd.read_csv(P/'tables/SECTION_QUALIFICATION.tsv',sep='\t');panels=json.loads((P/'resources/CD44_LOSS_PANELS_FROZEN.json').read_text());pure=spq[spq.assignment.eq('single_target')];lr=loc[loc.assignment.eq('single_target')&loc.endpoint.eq('distance_Spp1_detected_M')];main=e[e.dataset.eq('OBRADOVIC')]
assert len(main)==8 and main.P.ge(.05).all();assert pure.eligible.sum()==0;assert lr.eligible.sum()==4
def get(src,pr,mod):return main[main.source.eq(src)&main.program.eq(pr)&main.model.eq(mod)].iloc[0]
a=get('SPP1','Coculture','batch_stage_adjusted');b=get('SPP1','Coculture','batch_stage_CD8depth_adjusted')
source_files=[P/'downloads/RCC_SourceData.xlsx',P/'downloads/SPAC_mmc6.xlsx',P/'downloads/SPAC_mmc7.xlsx',P/'metadata/OT1_md5.xlsx',P/'resources/HOM_MouseHumanSequence.rpt']
sources={'date':stamp,'RCC_mIF_paper':'https://www.nature.com/articles/s41467-026-76858-w','RCC_source_url':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-76858-w/MediaObjects/41467_2026_76858_MOESM10_ESM.xlsx','RCC_blockade_NCBI_check':'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=SRR32360655&rettype=full&retmode=xml','SPAC_paper':'https://doi.org/10.1016/j.cell.2026.04.049','SPAC_processed_data':'https://spac.pku-genomics.org/','SPAC_code':'https://github.com/zenglab-pku/SPAC-seq','SPAC_bulk_DE':'https://ars.els-cdn.com/content/image/1-s2.0-S0092867426005167-mmc6.xlsx','SPAC_neighborhood_DE':'https://ars.els-cdn.com/content/image/1-s2.0-S0092867426005167-mmc7.xlsx','orthology':'https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt','file_hashes':{str(f.relative_to(P)):sha(f) for f in source_files},'single_file_download_rule':'all24SPACfiles<2GB,proxy;largest1239751801bytes','access_limits':'No private access bypass or messages to authors. NCBI public metadata corroborates unresolved RCC pairing; mIF Fig5C lacksIDs; cell segmented OT1 mapping not in fetched processed bin objects.'}
(P/'resources/SOURCES_AND_HASHES.json').write_text(json.dumps(sources,indent=2));(P/'resources/SOFTWARE_VERSIONS.json').write_text(json.dumps({n:importlib.metadata.version(n) for n in ['numpy','pandas','scipy','anndata','scikit-learn','h5py','openpyxl','requests']},indent=2))
audit=[{'candidate':'CIR-24-1146_RCC_SPP1_blockade','available':'NCBI/SRA publicmetadata;OMIXcountmatrixlisted','decision':'NOT_PAIRED_ANALYSIS_READY','reason':'Same-numberB/SPP1age conflicts;NCBIB1normal_tissue;noexplicit treatment/donor map;yearconflict;39GBFASTQnotdownloaded','next_required':'14-column countmatrix and explicit sample-treatment-donor table, ideally individualflow values'}, {'candidate':'NatCommun76858_ccRCC_mIF','available':'Officialsourcexlsx,11ccRCC+10tRCC','decision':'NO_CROSS_PANEL_PATIENT_JOIN','reason':'Sheet1G hasID/densCD8/densCD8PD1/histology;5C onlyperCD163p/perCD163pSPP1p/perCD163pMRC1p/histology, noIDorhiddencolumn;rowordercannotestablishpairing','next_required':'5C-to1GpatientIDmapping;denominator definitions;ideally percellintensities/coordinates'}, {'candidate':'SPAC_OT1_spatial_Cd44','available':'All6sections RNA/guide/annotation/mapping;24officialMD5matchedfiles','decision':'IMPLEMENTED_NO_CONSISTENT_SUPPORT','reason':'PureTrecipient0/6qualified;guide-location4/6qualified,1positive2negative1zero;strict2UMIsensitivity0/6qualified','next_required':'Segmented-cell sgRNA mapping and independentmouse/sectionIDs for literalT-cell spatialcomparison'}, {'candidate':'SPAC_CD44_KO_purified_CD8_bulk','available':'Officialmmc6,Coculture/Staticgeneeffects/rawP/FDR;MGI1to1orthology','decision':'IMPLEMENTED_RCC_PATIENT_VALIDATION','reason':'Coculture402/Static694orthologgeneprograms;Obradovicn11/n10;all8source-modelrawP>=0.05;co-culturerhodirectionflipswithdepth','next_required':'ExternalhumanSPP1-specificperturbation with donor-labelledreceiverRNA/functionalreadouts;notadditionalfilters on existingRCCdonors'}]
pd.DataFrame(audit).to_csv(P/'tables/CANDIDATE_ACTION_AUDIT.tsv',sep='\t',index=False)
for name in ['1G','5C']:pd.read_excel(P/'downloads/RCC_SourceData.xlsx',sheet_name=name).to_csv(P/f'tables/RCC_mIF_SOURCE_{name}.tsv',sep='\t',index=False)
sptext='\n'.join(f"{r['sample']}：CD44信号位置{int(r['Cd44_sites'])}，非靶向对照{int(r['NTC_sites'])}，资格{'通过' if r['eligible'] else '未通过'}。" for r in lq[lq.assignment.eq('single_target')].to_dict('records'))
lrtext='\n'.join(f'{r["sample"]}：距离中位数比={r["raw_distance_ratio"]:.6f}，标准化距离差={r["effect"]:.6f}。' for r in lr[lr.eligible].to_dict('records'))
etext='\n'.join(f'{r.source}—{r.program}，{r.model}：n={r.patients}，rho={r.rho:.6f}，P={r.P:.6f}，FDR={r.FDR:.6f}。' for r in main.itertuples())
hc=e[e.dataset.eq('HCRN')].groupby(['stratum','program']).patients.first().to_dict()
txt=f'''26210 空间扰动与CD44实验响应验证
Phase307｜{day}

本轮结论
已将新找到的空间CRISPR数据和纯化CD8的真实CD44敲除响应落实为分析。本轮未得到稳健的新增下游阳性：空间扰动信号的距离变化不一致，外部CD44敲除程序在肾癌CD8中的患者级关联均未达到原始P<0.05。可继续保留既有TYMP—SPP1上游阳性及Phase306外部功能证据；目前不能升级为本项目已验证SPP1驱动CD8受抑或物理屏障。

一、已取得完整空间扰动输入并完成两层复核
Cell2026 SPAC-seq研究提供第4、7、10天各2张切片。本轮下载了全部表达矩阵、sgRNA矩阵、细胞类型注释和坐标映射，共24文件、{dl['bytes']:,}字节（{dl['bytes']/1e9:.2f}GB）。最大单文件1.24GB，全部按用户规则走代理；官方MD5逐文件一致。首轮TLS断连后复用完整文件补齐，失败记录保留。
表达与guide对象的条码后缀和行数不同，所有RNA格点均按同切片唯一精确坐标匹配到guide；未按行号、去后缀或近邻猜配。下载对象是空间bin而非核分割单细胞。官方校验表注明10μm/50μm，映射CSV字段却写20μm/100μm；本轮使用原存储坐标及无量纲距离比，不虚构微米比例。

二、纯T接收端测量不足，区域定位没有一致方向
最初固定的纯T接收端要求：作者T注释、Cd3d/Cd3e/Trac至少一项检出、RNA计数>=20，并且仅检出一个目标gene的sgRNA。每组至少20个位置。6张切片中的合格Cd44位置只有0–12个，因此原96项空间/表达端点均记录为不可评估，没有用格点数补出生物学P值。
在确认上述资格失败后，另行冻结“sgRNA信号区域定位”问题：实验中的sgRNA来源为导入的OT1，但一个格点的RNA类型不等于该guide来自其中哪一细胞。区域分析只比较guide信号位置，不再解释为单个T细胞或CD8内在功能；原纯T门槛及失败结果保持。
{sptext}
可评估的4张切片中，距离SPP1检出巨噬区域1张增加、2张略减、1张持平：
{lrtext}
调整总体巨噬接近程度后仍没有一致方向。每目标至少2个guideUMI且占全部guide>=80%的固定敏感性均未达到每组20个位置。第10天一张切片无Cd44检出，另一张无单一Cd44目标位置，未补零或计入有效比较。当前数据不支持将该区域复核写成稳定的“Cd44扰动后远离SPP1巨噬”。

三、纯化CD8的CD44敲除响应已在肾癌患者级验证
从出版社原始补充表6取得两种条件的CD44_KO_vs_Ctrl差异表达结果：与MC38-OVA肿瘤细胞共培养，以及静态培养。原始基因P、FDR、效应量均保存。使用官方MGI同源类中唯一小鼠/人一对一关系，未靠大小写猜同源。
计算前固定作者原始P<0.05、|log2FC|>=0.5及baseMean>=50；排除CD44、SPP1、TYMP和CD3/CD8身份基因。共培养程序402基因（KO上调210、下调192），静态程序694基因（上调454、下调240）。患者CD8原始计数汇总后，计算上调组平均log1pCP10K减去下调组均值，得到CD44缺失方向程序。
Obradovic共培养程序11位、静态程序10位未治疗患者合格。巨噬SPP1—共培养CD44缺失程序在批次/分期调整后rho={a.rho:.6f}、P={a.P:.6f}、FDR={a.FDR:.6f}；加CD8深度后rho={b.rho:.6f}、P={b.P:.6f}、FDR={b.FDR:.6f}。方向随深度控制改变。巨噬TYMP相关模型和静态培养程序也没有原P<0.05的结果。全部8项主要模型见附录，未按显著性删选。
HCRN各时点按患者分别测量，资格数为{hc}；有效样本不足10例的模型只保留方向和效应，无P/FDR，不混合基线与进展或重复供者。

四、患者汇总的测量资格修订完整留档
最初沿用小型RNA标记程序的“每方向至少50%基因在>=5%单细胞检出”门槛。对402/694个bulk实验响应基因，这使所有患者不合格，未产生可评估关联。覆盖度审计显示大部分基因已被测量，但在稀疏单细胞中不到5%的细胞检出。
随后在任何有效关联计算前，单列患者伪bulk v2：沿用既有单基因测量支持标准（至少5个CD8细胞检出），要求每方向至少50%基因满足；仍要求CD8>=100且每方向>=80%基因在全部贡献文库中被测量。没有更换基因集、效应权重、方向、来源门槛或模型。静态程序的一位患者仍未通过门槛。原v1、覆盖审计、修订方案和v2完整保留，不能把v2描述为事先从未修订的方案。

五、两条肾癌直接数据入口仍缺患者映射
既有CIR-24-1146的SPP1阻断队列：本轮向NCBI原站交叉核查，B1被标为normal tissue；此前同编号B/SPP1年龄及采集年份冲突没有解决，缺少明确处理组和供者对应。公开原始测序并不自动使它具备配对分析资格，未下载39GB后按编号强配。
新NatCommun76858的蛋白验证包含11例ccRCC与10例tRCC，原xlsx已取得。Sheet1G提供患者ID及CD8/PD1-CD8密度，Sheet5C提供CD163/SPP1/MRC1比例但无患者ID，也无隐藏ID列。两表行顺序不同且没有可核实映射，因此未拼接相关性；两种组织学也未混成21例ccRCC。
当前真正能够重新开放这两条分析的资料是：干预计数矩阵及样本—处理—供者配对表；多重染色5C与1G的患者映射及比例分母定义，最好同时有逐细胞坐标/标记强度。SPAC空间再进一步则需要带独立小鼠ID的切片对应和分割细胞—sgRNA对应表。未向作者发信或提交访问申请。

六、统计方法与解释边界
空间所有效应按切片描述；独立小鼠对应未核实，不能写n=6只独立小鼠，也不计算格点P。主要区域端点为最近SPP1检出巨噬距离的log1p(距离/坐标网格步长)中位数差（Cd44减非靶向），并保留原距离中位数比、总体巨噬、SPP1未检出巨噬、成纤维距离和同格点比例。两种分配规则及全部72项区域组合均保存。
患者程序关联使用秩偏相关；Obradovic控制作者批次和pT，敏感性加CD8深度。n>=10、剩余自由度>=5才给双侧t近似P；按每队列/时点/来源的两程序×两模型四项家族BH，保留原P及FDR。HCRN分别报告各时点，细胞数不是样本量。
CD44存在多个配体，本轮程序是CD44缺失实验响应，不是SPP1特异扰动签名、耗竭定论或治疗预测。小鼠到人一对一同源映射也不能消除培养、物种和细胞状态差异。空间sgRNA位置的多数信号仅1个UMI，混合格点、环境RNA和吞噬/残留转录本仍会影响其细胞归属。最近距离端点也不同于原论文其他共定位算法和湿实验，当前复核结果不替代或推翻作者所有实验结论。
本轮在既有Obradovic/HCRN患者上进行外部实验引导的探索，不是新的独立人类验证队列。Phase305/306主要阴性和反向结果继续保留；上游患者内TYMP—SPP1正关联保持，不反向否定该事实。

七、验收与交付
正式下载、空间分析、独立空间验收、区域定位、两版本响应分析及最终验收均成功退出。24个原始输入MD5一致，最终共{v['checks']}项检查通过，包括独立坐标字典、BallTree与cKDTree距离一致性、合并对象与逐文库基因计数/检出一致性、患者资格、评分及P/FDR。验收脚本首轮把允许的缺失基因误判为必须全部存在，已修正为核验缺失状态和评分排除；没有改变分析数据、门槛或统计结果。
本轮原始矩阵及稽核缓存保存在sgs10，轻量包包含报告、全部统计/资格表、冻结规则、源表、代码与校验。没有修改旧C4、既有对象、Phase302客户Word或正式文章。无新增生物学图件。

附录：Obradovic完整关联结果
{etext}

来源
SPAC-seq原文：https://doi.org/10.1016/j.cell.2026.04.049
作者处理后数据：https://spac.pku-genomics.org/
作者分析代码：https://github.com/zenglab-pku/SPAC-seq
肾癌多重染色来源：https://www.nature.com/articles/s41467-026-76858-w
MGI同源基因：https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
具体源表URL与完整SHA256见resources/SOURCES_AND_HASHES.json；公开分享的临时签名链接不随交付包散发。
'''
stem='26210_Phase307_空间扰动与CD44响应验证';(D/(stem+'.txt')).write_text(txt,encoding='utf-8-sig');body=''.join('<p>'+html.escape(s).replace('\n','<br>')+'</p>' for s in txt.split('\n\n'));(D/(stem+'.html')).write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>26210 下游实施结果</title><style>body{max-width:1050px;margin:40px auto;padding:0 24px;font:16px/1.9 "Microsoft YaHei",sans-serif;color:#263442}p{margin:1.4em 0}@media print{body{font-size:11pt;margin:0}}</style><body>'+body+'</body></html>',encoding='utf-8')
for folder in ['tables','resources','qc','code']:(D/folder).mkdir(exist_ok=True)
for f in (P/'tables').glob('*.tsv'):shutil.copy2(f,D/'tables'/f.name)
for f in (P/'resources').iterdir():
 if f.is_file():shutil.copy2(f,D/'resources'/f.name)
qc_names=['OFFICIAL_MD5_VALIDATION.json','INPUT_IDENTITY.json','ANALYSIS_COMPLETE.json','FINAL_VALIDATION.json','LOCALIZATION_COMPLETE.json','CD44_RESPONSE_COMPLETE.json','CD44_PSEUDOBULK_VALIDATION_V2.json','PHASE307_FINAL_VALIDATION.json','VALIDATOR_MISSING_GENE_REPAIR.json']
for name in qc_names:shutil.copy2(P/'qc'/name,D/'qc'/name)
for f in (P/'code').glob('*.py'):shutil.copy2(f,D/'code'/f.name)
(D/'source_data').mkdir(exist_ok=True)
for f in source_files:
 if f.suffix=='.xlsx':shutil.copy2(f,D/'source_data'/f.name)
manifest={str(f.relative_to(D)):sha(f) for f in D.rglob('*') if f.is_file() and f.name!='SHA256.json'};(D/'SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));archive=P/(stem+'.zip')
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in D.rglob('*'):
  if f.is_file():z.write(f,f.relative_to(D))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
(P/'qc/DELIVERY_VALIDATED.json').write_text(json.dumps({'status':'PASS','bytes':archive.stat().st_size,'sha256':sha(archive),'files':len(manifest),'checks':v['checks']},indent=2))
message=f'''## Phase307终态（{day}）：空间CRISPR与CD44真实实验响应已实施
用户“继续查找下一步实施”。新Cell2026 SPAC-seq六张Day4/7/10切片24文件共{dl['bytes']}字节，单文件均<2GB代理，官方MD5全部一致。RNA/guide条码不同，按精确唯一坐标对齐，全部RNA匹配。纯T+单guide分析Cd44仅0–12位置/片，6片不合>=20资格，96端点不可评估。另冻结sgRNA区域定位（不是纯T细胞）：4片合格，SPP1巨噬距离比Day4rep1/rep2/Day7rep1/rep2为1.000/1.10086/0.98005/0.98985；1正2负1平。Day10两片不合格，>=2UMI/80%敏感性全不合；72完整组合，无生物/格点P，无一致新支持。
取得出版社mmc6纯化CD8 CD44KO共培养/静态实验DE表，按原P<.05、|log2FC|>=.5、baseMean>=50及MGI唯一1对1映射，排除目标和身份基因。固定402(210up/192down)/694(454up/240down)基因程序。患者伪bulk旧5%单细胞检出v1全不合，在任何可评估关联前按患者测量单位单列v2（>=50%方向基因至少5细胞检出，保留>=80%测量/CD8>=100）；基因集/评分/模型不变，原v1与审计保留。Obradovic共培养n11/静态n10；SPP1共培养rho-0.488355/P.182251/FDR.243002，加深度rho+0.532710/P.174047/FDR.243002；全部8主要模型原P>=.05。HCRN各时点描述，不混供者。CD44缺失程序不是SPP1特异响应，旧Phase305/306阴性及上游阳性保持。
肾癌直接入口：NCBI B1标normal tissue加深旧干预标签冲突；NatCommun76858源表11ccRCC/10tRCC已取，1G有ID而5C无ID/隐藏列，不能按行拼；未跑伪配对。作者处理空间是bin，独立小鼠映射未定，不称6鼠。新的空间端点不一致不等于推翻作者其他实验。下一有效输入是带处理/供者映射的RCCSPP1计数/流式、mIF跨面板ID/逐细胞坐标、带鼠ID的分割sgRNA表；未联系作者/发访问申请。
权威目录phase307_spatial_perturbation_and_RCC_protein_20261002_v1；独立原始对象与逐文库复算/几何/资格/统计总{v['checks']}检查PASS。报告delivery/{stem}.txt和html；包{stem}.zip。TLS重试及验收missinggene断言修复均存档，不改生物统计结果。全部授权首轮实施完成，无分析后台。旧C4/对象/正式文章/Phase302Word未覆盖，无新图。当前未新增稳健下游阳性，不为阳性继续改阈值。

'''
backup=P/'resources'/('root_before_completion_'+now.strftime('%Y%m%d_%H%M%S'));backup.mkdir()
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=R/name;old=f.read_bytes();(backup/name).write_bytes(old);assert f.read_bytes()==old;t=f.with_name(f.name+'.phase307.tmp');t.write_bytes(message.encode()+old);os.replace(t,f)
f=R/'RESULTS_MANIFEST.tsv';old=f.read_bytes();(backup/f.name).write_bytes(old);fields=next(csv.reader(io.StringIO(old.decode('utf-8-sig')),delimiter='\t'));buf=io.StringIO();wr=csv.DictWriter(buf,fieldnames=fields,delimiter='\t',lineterminator='\n')
for ident,path,typ in [('phase307_report',D/(stem+'.txt'),'report'),('phase307_spatial',P/'tables/GUIDE_LOCALIZATION_EFFECTS.tsv','spatial_perturbation'),('phase307_cd44_response',P/'tables/CD44_PSEUDOBULK_ASSOCIATIONS_V2.tsv','patient_statistics'),('phase307_validation',P/'qc/PHASE307_FINAL_VALIDATION.json','validation'),('phase307_delivery',archive,'archive')]:wr.writerow({'result_id':ident,'stage':'Phase307','artifact_type':typ,'absolute_path':str(path),'generated_at':stamp,'validation_status':'PASS','validation_evidence':str(P/'qc/PHASE307_FINAL_VALIDATION.json'),'notes':'New spatial perturbation and externalCD44 response executed; no consistent supportive signal; full negative and qualification results preserved'})
assert f.read_bytes()==old;t=f.with_name(f.name+'.phase307.tmp');t.write_bytes(old+(b'' if old.endswith(b'\n') else b'\n')+buf.getvalue().encode());os.replace(t,f)
(P/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_VALIDATED','time':stamp,'active_analysis_jobs':[],'checks':v['checks'],'new_robust_downstream_support':False,'archive_sha256':sha(archive)},indent=2));print('DELIVERY_COMPLETE',archive.name,archive.stat().st_size,sha(archive),len(manifest),flush=True)
