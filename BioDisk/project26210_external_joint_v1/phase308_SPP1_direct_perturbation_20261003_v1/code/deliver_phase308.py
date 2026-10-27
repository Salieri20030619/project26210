from pathlib import Path
import datetime, hashlib, html, json, os, shutil, zipfile, importlib.metadata
import pandas as pd

P=Path(__file__).resolve().parents[1]; R=P.parent
now=datetime.datetime.now().astimezone().isoformat()
validation=json.loads((P/'qc/FINAL_VALIDATION.json').read_text())
assert validation['status']=='PASS' and validation['checks']==2183
for name in ['prepare_inputs','qc_cluster','refine_identity','analyze_pools','validate_pools']:
    assert (P/'logs'/f'{name}.exit').read_text().strip()=='0',name
effects=pd.read_csv(P/'tables/PERTURBATION_EFFECTS_ALL.tsv',sep='\t')
pools=pd.read_csv(P/'tables/POOL_ENDPOINTS_ALL.tsv',sep='\t')
counts=pd.read_csv(P/'tables/LINEAGE_COUNTS_BY_POOL.tsv',sep='\t')
assert len(effects)==35 and len(pools)==210 and counts.QC_cells.sum()==77087
D=P/'delivery'; D.mkdir(exist_ok=True)
title='26210_Phase308_SPP1直接扰动与基质执行证据'
package=D/title; package.mkdir(exist_ok=True)

sources=[
 {'source':'GSE274609 / PDAC 宿主 Spp1 缺失单细胞','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274609','article':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12479095/','availability':'已下载并完成原始计数重分析；696801280 bytes；3混样池/组','role':'本阶段实施的外部真实遗传扰动'},
 {'source':'JITC 2026 巨噬/髓系 Spp1 条件缺失','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12829364/','article':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12829364/','availability':'正文数据声明为合理请求提供，未找到可直接下载的原始计数 accession；未联系作者','role':'文献机制支持，不能记为本项目已验证数据；其 CXCL9/10 上升与本阶段结果不同'},
 {'source':'LUAD POSTN CAF–SPP1 巨噬细胞互作','url':'https://link.springer.com/article/10.1186/s12967-026-08746-2','article':'https://link.springer.com/article/10.1186/s12967-026-08746-2','availability':'公开 GEO 为观察性队列；功能实验原始值按请求提供，未取得干预矩阵','role':'SPP1/CD44、基质与 CD8 迁移实验的文献背景'},
 {'source':'GSE316918 / 滑膜 OPN 刺激成纤维细胞','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316918','article':'https://pubmed.ncbi.nlm.nih.gov/41832953/','availability':'已取元数据，未下载计数矩阵','role':'非肿瘤炎症消退情境，保留来源审计，不作为肿瘤屏障的同向验证'},
 {'source':'小鼠 Tymp 官方基因标识','url':'https://www.ncbi.nlm.nih.gov/gene/72962','article':'https://www.informatics.jax.org/marker/MGI%3A1920212','availability':'官方 symbol、别名和 ENSMUSG00000022615 均核对；六库 feature 表均不含对应项','role':'未测量审计；不将缺失特征记为零表达'}
]
pd.DataFrame(sources).to_csv(P/'tables/EXTERNAL_SOURCE_AUDIT.tsv',sep='\t',index=False)
def effect(lin,ep):
    z=effects[(effects.lineage==lin)&(effects.endpoint==ep)]
    assert len(z)==1
    return z.iloc[0]
names={'collagen_maturation':'胶原成熟','focal_adhesion':'黏着斑','actomyosin':'收缩装置','myCAF_markers':'肌成纤维标志','steady_state_markers':'稳态成纤维标志'}
matrix=[]
for ep,name in names.items():
    z=effect('Fibroblast',ep)
    matrix.append(f"{name}\t{z.KO_median:.6f}\t{z.WT_median:.6f}\t{z.median_difference:+.6f}\t{int(z.negative_comparisons)}/9 降低；{int(z.positive_comparisons)}/9 升高")
matrix_text='端点\tKO 中位数\tWT 中位数\tKO−WT\t全部跨组比较方向\n'+'\n'.join(matrix)
report=f'''26210 项目 Phase308：宿主 SPP1 缺失支持成纤维细胞基质执行状态改变
生成时间：{now}
状态：分析及独立原始计数验收完成；2183 项检查通过。

一、结论：补到了基质端的外部真实扰动证据
在 GSE274609 小鼠胰腺癌单细胞数据中，宿主 Spp1 缺失后，严格鉴定的成纤维细胞黏着斑、收缩装置和肌成纤维标志程序均降低，三份 KO 混样池均低于三份 WT 混样池。胶原成熟程序的组中位数也降低，9 次跨组比较中 7 次同向。这是新的基质端探索性扰动证据。
与此同时，成纤维细胞捕获比例增加。因此支持的是 SPP1 缺失伴随成纤维细胞状态改变，不能表述为“成纤维细胞减少导致屏障消失”。免疫端的 T/CD8 捕获比例存在池间不一致，巨噬细胞 CXCL9/10 程序降低；本阶段不支持稳定的 CD8 招募恢复机制。

二、基质执行程序降低，细胞捕获比例增加
{matrix_text}
分数单位均为固定基因集的平均 log1p(CP10K)，以每份混样池、每种谱系的原始计数进行伪 bulk 汇总；差值不是百分比，也不是蛋白或组织收缩的倍数。每组有效独立测量单位为 3 个混样池，不计算 P 值或 FDR。9 次比较只是 3×3 个跨组方向比较，不能当作 n=9。
成纤维细胞捕获比例中位数 KO 3.3075%、WT 0.9838%，为 3.362 倍；三份 KO 均高于三份 WT。成纤维细胞 Spp1 分数降低（中位差 −0.394031，9/9 次比较降低）。肌成纤维和稳态标志为预先记录的次要读数，未用来重新定义细胞群。

三、免疫组成改善不一致，CXCL9/10 方向相反
T 细胞占全部 QC 细胞的中位比例为 KO 3.6358%、WT 1.8769%（1.937 倍）；KO1、KO2 较高，KO3 较低。CD8 占全部 QC 细胞的中位比例为 0.9174% 对 0.7209%（1.272 倍），6/9 次跨组比较升高。改用免疫细胞作分母，CD8 中位比例为 2.1047% 对 1.9174%（1.098 倍），只有 5/9 次跨组比较升高。
CD8 功能程序要求每池至少 100 个合格 CD8 细胞；KO3 为 51、WT2 为 65，故每组只有 2 池合格，未达到预设完整 3 对 3 比较条件，功能端点记为不可完整评估。所有逐池数值与资格原因均保留，没有删除不合格池后重新声明完整比较。
巨噬细胞 CXCL9/10 分数 KO 0.537052、WT 0.723616，中位差 −0.186564，三份 KO 均低于三份 WT。这与另一项髓系 Spp1 条件缺失文献中的 CXCL9/10 上升不同，必须保留。巨噬细胞 Spp1 分数也降低（中位差 −1.271763，9/9 次比较降低）。

四、方法：按混样池进行推断单位控制与原始计数汇总
实验为 WT 6620c1 胰腺癌细胞原位植入 WT 或全身宿主 Spp1 KO 小鼠。每组 6 个肿瘤，每 2 个肿瘤混合建库，因此实际为 3 个独立混样池/组；KO1 与 WT1 等同号文库不配对。六个文库为 GSM8453641–GSM8453646。
原始归档 696801280 字节，通过代理下载，SHA256 为 169f9c3a6d5dc22a02438fce95f9a7570e4c0f7c36375fab5d33390797fde380。六个内部归档和压缩计数文件均完成 gzip CRC 检验；上述 SHA256 是本地完整性记录，不宣称为官方发布的校验值。
每库 QC：至少 200 个基因、500 UMI，至多 6000 个基因，线粒体比例不超过 20%，去除该库 UMI 高于 99.5 分位的细胞；Scrublet 预期双细胞率 0.06，按实际预测去除 3 个细胞。最终保留 77087 个细胞、32285 个特征。QC 细胞总数并非统计样本量。
初始 log1p(CP10K)、按文库选 3000 个高变基因，排除 Spp1/Tymp，PCA40、邻居20、Leiden0.7。对混有腺泡背景或低复杂度信号的 1/4/6/7 簇进行事先记录的身份细化：2000 高变基因、PCA35、邻居20、Leiden0.8，额外排除线粒体、核糖体和指定预测基因。最终 32 个群，保留全部 QC 细胞。
在计算 KO/WT 效应前固定群身份。仅 0 簇且满足至少两个 Col1a1/Col1a2/Col3a1/Dcn/Lum/Pdgfra、同时排除跨谱系信号的细胞纳入成纤维细胞。R4/R8 的低复杂度或混合基质群未归为纯 CAF。未建立可靠的周细胞或间皮群，对应端点记为不可评估，不解释为生物学缺失。巨噬、T、CD8 及免疫分母均按固定簇与多标志规则确定，完整规则和代码随包提供。
基质三程序沿用 Phase292 冻结基因集，通过 MGI 一对一同源映射；胶原成熟 9 基因、黏着斑 9 基因、收缩装置 8 基因。次要 myCAF 为 Postn/Acta2/Lrrc15，稳态为 Pi16/Dpt。CD8 细胞毒性为 Gzmb/Prf1/Nkg7/Ctsw；功能障碍为 Pdcd1/Ctla4/Lag3/Havcr2/Tigit/Tox/Cxcl13/Entpd1。小鼠无对应人 GNLY 的直接一对一同源，已记录。
伪 bulk 汇总每池每谱系原始计数，按该谱系全部计数归一到 CP10K，log1p 后对预设基因等权平均。巨噬、CD8、上皮样细胞至少 100 个/池，成纤维、周细胞、间皮样细胞至少 20 个/池；要求至少 80% 面板基因在 feature 表中测量。
在任何组间端点效应计算前，修订了从观察性流程继承的“至少 5 个检出细胞”资格门槛：本次干预分析保留实际测量到的零值，资格由测量覆盖和谱系细胞数决定，检出细胞数仅保留为表达支持描述。原始协议、修订说明及时间均留档，未根据效应修改基因集、谱系、预期方向或样本数门槛。
完整比较要求各组全部 3 池均合格，报告中位数、范围、中位差、组成比例的中位数比及全部跨组方向；不作细胞级显著性检验。全部 210 条逐池端点、35 条组间端点以及原始 P/FDR 空值列均保留。

五、讨论：它补强了免疫屏障主线的哪一段
免疫屏障主线包含基质形成/维持与免疫细胞进入/功能两个不同环节。这里新增加的结果支持前一环节：宿主 SPP1 缺失伴随成纤维细胞黏着斑、收缩和肌成纤维程序一致下降，为既有 TYMP–SPP1 正向关系及基质成熟关联提供外部功能背景。此前人肾癌上游 TYMP–SPP1 阳性结果继续保留，不被本次免疫端的不一致覆盖。
本实验是胰腺癌的全身宿主 Spp1 缺失，作者强调抗原呈递型 CAF 是其主要基质 SPP1 来源；不能把全部效应归因于巨噬细胞特异性 SPP1，也不能直接外推为肾癌或 TYMP 的处理效应。各组只有 3 个混样池，无法恢复池内单个肿瘤变异，结果定位为探索性。
所有六库的特征表均不含 Tymp 官方标识、别名或 ENSMUSG00000022615，本阶段不能评价 TYMP 表达变化；缺少该特征不等于表达为零。Spp1 RNA 残余也不等于基因型无效：实验使用 WT 肿瘤细胞，且转录读数不能直接替代功能蛋白鉴定。
转录程序下降没有直接测量胶原几何、组织刚度或 T 细胞跨越屏障的行为；单细胞捕获比例不等于绝对浸润量。CD8 组成异质、功能端点未满足完整资格、巨噬 CXCL9/10 反向，故当前数据不支持“基质改变已经稳定恢复免疫进入/功能”，也未闭合 TYMP→SPP1→肾癌免疫屏障的完整因果链。Phase305–307 的关键阴性、反向结果继续保留。

六、可用于客户报告的结果段落
“为检验 SPP1 与基质屏障维持的关系，我们进一步重分析了宿主 Spp1 缺失的小鼠肿瘤单细胞数据。Spp1 缺失后，成纤维细胞的黏着斑、收缩装置及肌成纤维标志程序在三份独立混样池中均低于对照，胶原成熟程序总体也降低。这一结果为 SPP1 参与成纤维细胞基质执行状态提供了外部探索性扰动证据，并与本项目已发现的 TYMP–SPP1 正向关系及基质成熟关联相衔接。免疫细胞组成的变化仍存在样本异质性，尚不能将该基质状态变化等同于免疫排斥的解除。”
该段为本阶段结果稿；未覆盖已有客户 Word、正式文章或旧分析对象。

七、数据入口审计与后续优先级
本次实际实施的是 GSE274609，其余入口的可用性见 EXTERNAL_SOURCE_AUDIT.tsv。JITC 髓系 Spp1 条件缺失和 LUAD SPP1/CD44 功能实验提供文献背景，但功能原始数据按请求提供，未当作已取得的新队列。滑膜 OPN 刺激属于不同炎症情境，未为追求同向结果继续展开。
后续若继续实施，优先取得具备个体/处理映射的肾癌直接 SPP1 干预计数或流式，其次是可跨面板配对的 mIF 个体 ID/逐细胞坐标。单纯增加同批患者相关分析不会替代这些证据。本阶段不以等待作者答复作为交付条件，未对外联系或提交申请。

八、复核与文件
独立校验直接读取六个原始文库计数，重新汇总每池每谱系基因计数和分母，重新计算归一化、端点中位数、全部跨组方向、资格判定及样本映射，2183 项检查通过；各分析进程退出码为 0。计数/H5AD 和下载归档保留在服务器权威目录，轻量交付包含报告、全量端点表、身份/QC表、冻结协议、来源审计、脚本、校验记录和逐文件 SHA256。
服务器权威目录：{P}

参考来源
'''+ '\n'.join(f"- {x['source']}\n  {x['url']}\n  {x['article']}\n  {x['availability']}" for x in sources)+'\n'
(package/(title+'.txt')).write_text(report,encoding='utf-8')
table=effects[(effects.lineage=='Fibroblast')&effects.endpoint.isin(names)][['endpoint','KO_median','WT_median','median_difference','negative_comparisons','positive_comparisons']].copy()
table['endpoint']=table.endpoint.map(names)
sections=report.split('\n\n')
body=''.join('<section><pre>'+html.escape(s)+'</pre></section>' for s in sections)
document='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>'+html.escape(title)+'</title><style>body{font:17px/1.8 system-ui,"Microsoft YaHei",sans-serif;color:#173042;background:#f4f7fa;margin:0}main{max-width:1100px;margin:28px auto;padding:32px;background:white}h1{font-size:28px}section{border-bottom:1px solid #e2e8ef;padding:12px 0}pre{font:inherit;white-space:pre-wrap;overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:8px;border:1px solid #ccd7e0;text-align:right}th:first-child,td:first-child{text-align:left}@media print{body{background:white}main{margin:0;padding:0}}</style><main><h1>SPP1 缺失后的基质执行状态</h1><p>Phase308 · 外部真实扰动 · 3 个混样池/组 · 2183 项检查通过</p>'+table.to_html(index=False,float_format=lambda x:f'{x:.6f}',border=0)+body+'</main></html>'
(package/(title+'.html')).write_text(document,encoding='utf-8')
for directory in ['tables','resources','qc','code']:
    target=package/directory;target.mkdir(exist_ok=True)
    for f in (P/directory).iterdir():
        if f.is_file() and not f.name.endswith('.tsv.gz') and f.name not in ['DELIVERY_VALIDATED.json','LOCAL_DELIVERY_VERIFIED.json']:
            shutil.copy2(f,target/f.name)
(package/'logs').mkdir(exist_ok=True)
for f in (P/'logs').glob('*.exit'):
    if f.stem!='deliver_phase308':shutil.copy2(f,package/'logs'/f.name)
(package/'metadata').mkdir(exist_ok=True)
shutil.copy2(P/'metadata/FETCH_MANIFEST.json',package/'metadata/FETCH_MANIFEST.json')
versions={name:importlib.metadata.version(name) for name in ['scanpy','anndata','numpy','pandas','scipy','scikit-learn','h5py']}
versions['scrublet_implementation']='scanpy.pp.scrublet; bundled in scanpy, no separate distribution metadata'
(package/'qc/SOFTWARE_VERSIONS.json').write_text(json.dumps(versions,indent=2))
records=[{'path':str(f.relative_to(package)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted(package.rglob('*')) if f.is_file() and f.name!='SHA256_MANIFEST.json']
(package/'SHA256_MANIFEST.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
archive=D/(title+'.zip')
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for f in sorted(package.rglob('*')):
        if f.is_file():z.write(f,f.relative_to(package))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for row in records:assert hashlib.sha256(z.read(row['path'])).hexdigest()==row['sha256']
delivery={'status':'PASS','created':now,'archive':str(archive),'archive_name':archive.name,'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'manifest_files':len(records),'archive_files':len(records)+1,'archive_CRC':'PASS','all_packaged_sha256':'PASS','analysis_validation':validation,'reports':[str(package/(title+ext)) for ext in ['.txt','.html']]}
(P/'qc/DELIVERY_VALIDATED.json').write_text(json.dumps(delivery,ensure_ascii=False,indent=2))
status={'status':'COMPLETE_VALIDATED','time':now,'new_external_matrix_directional_support':True,'full_TYMP_RCC_immune_barrier_causal_chain':False,'sample_unit':'3 pooled libraries/arm,2tumors/pool','analysis_jobs':'none','delivery':delivery,'local_delivery':'pending verification'}
(P/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
message=f'''## Phase308终态（2026-10-03）：直接宿主SPP1扰动补强基质执行状态
GSE274609已实施：WT6620c1胰腺癌植入WT或全身宿主Spp1KO，3混样池/组、2肿瘤/池。非巨噬特异、非RCC、非TYMP扰动。77,087 QC细胞；固定身份后按池伪bulk，全部210池端点/35组端点保留，无P/FDR。
新探索性阳性：严格成纤维细胞黏着斑Δ−0.079298、收缩装置−0.182098、myCAF−0.236029，全部3KO低于全部3WT（9/9跨组比较）。胶原成熟Δ−0.151457，7/9降低；稳态标志+0.136479，7/9升高。成纤维捕获比例3.362倍，不能写为CAF减少。
免疫主要矛盾保留：T/CD8全QC比例中位数1.937/1.272倍，但KO3反向；CD8/免疫仅1.098倍、5/9正。CD8功能因KO3/WT2不足100细胞而每组仅2池，不满足完整3对3。巨噬CXCL9/10Δ−0.186564且9/9降低，不支持此模型中趋化恢复。六库无Tymp特征（symbol/别名/Ensembl核对），不可把缺测当零或称TYMP处理效应。既有TYMP–SPP1上游阳性及Phase305–307阴性全部保持。
原观察性检出门槛在计算任何效应前修订为保留已测量零值、按测量覆盖/细胞数评估资格，初始协议/修订/时间保留；不依效应调面板/身份/样本门槛。JITC条件KO/LUAD功能数据仅按请求提供，GSE316918滑膜情境仅元数据审计；未联系作者。
权威目录：{P.name}。原始六库独立复算2183检查PASS，各分析退出0；交付ZIP逐文件SHA256与CRC通过。包：delivery/{archive.name}。本地副本由qc/LOCAL_DELIVERY_VERIFIED.json记录；此处仅确认服务器交付已完成。无分析后台，无新图，旧Word/正式文章/对象不覆盖。
下一优先输入：具备个体/处理映射的RCC直接SPP1干预计数或流式，其次跨面板可配对mIF个体ID/逐细胞坐标。当前新增的是外部基质状态支持，完整免疫屏障因果链未闭合。

'''
backup=P/'resources'/('root_before_final_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'));backup.mkdir()
updates={}
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
    old=(R/name).read_bytes();updates[name]=(old,message.encode()+old)
mf=R/'RESULTS_MANIFEST.tsv';old=mf.read_bytes()
entries=[('phase308_report','report',package/(title+'.html'),'external SPP1 perturbation; exploratory fibroblast execution-state support; immune heterogeneity retained'),('phase308_effects','table',P/'tables/PERTURBATION_EFFECTS_ALL.tsv','35 endpoints;3pools/arm;no P/FDR'),('phase308_validation','validation',P/'qc/FINAL_VALIDATION.json','2183 original-count independent checks'),('phase308_delivery','archive',archive,'ZIP CRC and all file SHA256 validated')]
rows=''.join('\t'.join([rid,'Phase308',kind,str(path),now,'PASS',str(P/'qc/FINAL_VALIDATION.json'),notes])+'\n' for rid,kind,path,notes in entries)
assert b'phase308_delivery\t' not in old
updates[mf.name]=(old,old+(b'' if old.endswith(b'\n') else b'\n')+rows.encode())
for name,(old,new) in updates.items():
    f=R/name;(backup/name).write_bytes(old);assert f.read_bytes()==old
    temporary=f.with_name(f.name+'.phase308.final.tmp');temporary.write_bytes(new)
    assert f.read_bytes()==old;os.replace(temporary,f)
print(json.dumps(delivery,ensure_ascii=False,indent=2),flush=True)
