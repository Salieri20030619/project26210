from pathlib import Path
import os,json,datetime,hashlib,shutil,zipfile,html,importlib.metadata
import pandas as pd
P=Path(__file__).resolve().parents[1];R=P.parent
v=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());assert v['status']=='PASS'
for s in ['acquire','identity_audit','analyze_spatial','validate_spatial']:assert (P/'logs'/f'{s}.exit').read_text().strip()=='0'
now=datetime.datetime.now().astimezone().isoformat();title='26210_Phase309_肾癌逐细胞空间排斥检验';D=P/'delivery';pkg=D/title;pkg.mkdir(exist_ok=True)
e=pd.read_csv(P/'tables/SPATIAL_ASSOCIATIONS_ALL.tsv',sep='\t');b=pd.read_csv(P/'tables/SPATIAL_BLOCK_SUMMARY.tsv',sep='\t');c=pd.read_csv(P/'tables/FIBRO_NEIGHBORHOOD_CONTRASTS.tsv',sep='\t');n=json.loads((P/'qc/IDENTITY_AUDIT_COMPLETE.json').read_text())['counts']
display=e.pivot(index=['target','fibro_radius_um','CD8_cells'],columns='model',values='rho').reset_index();display=display[['target','fibro_radius_um','CD8_cells','unadjusted','density_vascular_depth_xy']]
tabletext=display.to_csv(sep='\t',index=False,float_format='%.6f')
sources=[
 {'source':'10x RCC新版 RNA/蛋白同细胞空间数据','url':'https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-ccrcc','status':'已实施','sample_unit':'1供者、1组织切片；690322检测细胞','details':'2026-08-19页面；XOA4.0；实际矩阵477RNA/27蛋白；TYMP/SPP1未测量。坐标、作者细胞分组和原始矩阵共100995534字节，全部代理下载并核SHA；与旧465534细胞版本分开，不当独立患者复制'},
 {'source':'Proseg 四患者RCC Xenium','url':'https://figshare.com/articles/dataset/Xenium_Renal_Cell_Carcinoma/25685961','status':'元数据已取；未启动矩阵下载','sample_unit':'作者声明4患者；公开10个同名transcripts.csv.gz文件','details':'API未提供本次能够核实的文件到患者映射，未取得作者已注释矩阵；不能把10文件当10患者或猜配。7个<2GB？以逐文件表为准，不在此手写计数'},
 {'source':'Nicheverse 93患者/173样本RCC图谱','url':'https://nicheverse.org/','status':'公开样本目录已取，完整逐细胞入口未取得','sample_unit':'论文93患者、173样本，不能把173当患者数','details':'gallery_data.js有173样本元数据和预览链接；viewer配置指向https://data.nicheverse.org/v1，索引请求两次SSL失败。公开可读页面不等于全量计数/坐标已取。旧core21只有7824x366计数与坐标、sample_id，无作者细胞注释且缺TYMP/SPP1；本轮未对其重新猜标签'},
 {'source':'Proseg已发表空间邻近方法','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12285883/','status':'方法背景','sample_unit':'原文4患者','details':'原文以局部细胞组成归一的随机游走评估邻近；本轮使用直接距离加预设密度调整，未声称复现其随机游走算法'}
]
sources[1]['details']='API未提供本次能够核实的文件到患者映射，未取得作者已注释矩阵；不能把10文件当10患者或猜配。实际大小/官方MD5逐文件保存，未下载或重建未配对矩阵。'
pd.DataFrame(sources).to_csv(P/'tables/SOURCE_FEASIBILITY_AUDIT.tsv',sep='\t',index=False)
j=json.loads((P/'metadata/RCC_figshare.json').read_text());pd.DataFrame(j['files']).to_csv(P/'tables/PROSEG_PUBLIC_FILE_MANIFEST.tsv',sep='\t',index=False)
report=f'''26210项目 Phase309：肾癌逐细胞空间排斥检验
完成时间：{now}

一、结论：观察到基质候选丰富邻域的CD8空间分离，密度调整后的关联较弱
本阶段已从解离单细胞的总量/程序比较推进到真实组织中的逐细胞坐标分析。在一个肾透明细胞癌切片中，成纤维细胞候选较多的邻域，CD8到作者标注癌细胞的距离更远。主尺度下，低成纤维邻域的距离中位数为37.33μm，高成纤维邻域为63.23μm；20μm内存在癌细胞的CD8比例分别为29.06%和9.79%。这是实际坐标上的空间分离现象。
加入局部癌细胞密度、总体细胞密度、血管相关细胞密度、RNA计数和位置趋势后，主关联rho由0.205639降至0.022175。分区后236个合格空间块中126个正向、110个负向，块内rho中位数0.012441。故本案例提供描述性的基质候选—CD8癌细胞接近关系，当前数据不支持强而一致、独立于这些局部组织因素的关联，未确认基质的物理阻挡作用。
有效生物学样本量为1个供者。细胞数和空间块数均不作独立生物学重复，全阶段不计算P值、FDR或细胞/区域推断性置信区间。

二、直接空间端点与全部预设模型
主问题：在癌细胞周围200μm内的CD8中，50μm邻域成纤维候选数量是否与到最近癌细胞的距离正相关？固定半径计数与对应密度成正比，表内保存原计数以便复算。
主参考为作者明确标注的三类癌细胞；敏感性另要求这些癌细胞的CA9 RNA至少1计数，并把成纤维邻域半径改为100μm。CA9参考对应的合格CD8集合不同，作为预设测量定义敏感性保留，不称相同细胞集合的独立复制。
{tabletext}
Tumor：全部合格作者癌细胞；Tumor_CA9：其CA9检出子集。unadjusted为未调整Spearman；density_vascular_depth_xy为预设协变量调整后的秩残差相关。正rho表示邻域成纤维候选越多、CD8到癌细胞越远。
癌细胞参考有32949个，CA9限定参考18314个；29118个合格CD8中，28003个进入主200μm分析区，12195个进入CA9限定分析区。全部CD8到主癌细胞参考的距离中位数为57.63μm；15.41%的CD8在20μm内存在该参考癌细胞。
主50μm成纤维邻域的下四分位界为1个、上四分位界为5个。低组7028个CD8、高组9446个CD8；局部100μm癌细胞数中位数分别为4与2，总QC细胞数中位数为276与266。两组的癌细胞密度本来就不同，不能将未经调整的37.33对63.23μm直接解释为处理效应或阻挡强度。

三、区域方向并不一致
固定500μm网格、每块至少30个合格CD8，并使用与全切片相同的预设调整。全部区域含不可评估原因均保留。
{b.to_csv(sep='\t',index=False,float_format='%.6f')}
236或121是同一切片的技术空间分区数，不是患者数，也不是新增验证队列。CA9限定50μm模型中58块正、63块负，中位rho−0.010522；100μm模型60正、61负，中位−0.002975。不同区域的异质性必须与全切片正向rho一起呈现。

四、方法：先核细胞身份，再固定几何与调整变量
来源为10x Genomics 2026-08-19发布的单供者FFPE ccRCC RNA/蛋白空间数据，Xenium Onboard Analysis4.0。下载cells.csv.gz、cell_feature_matrix.h5和cell_groups.csv，文件均小于2GB，经代理获取并记录字节数、SHA256及HTTP ETag。坐标和矩阵按唯一cell_id精确对齐；477个Gene Expression、27个Protein Expression分别识别，控制探针不混作基因。
基础QC：RNA计数至少10、坐标有限、细胞面积为正、nucleus_count等于1。690322个输入细胞中609525个通过。该规则用于控制低信息和多核分割混合，在任何空间效应之前固定。
作者分组文件包含359085个细胞、18类注释，主要为免疫和肿瘤亚群；未注释细胞保留为未注释，不直接改名CAF。
CD8要求位于预设作者T细胞组、CD8A/B至少一个RNA标记检出，同时CD3D/CD3E/TRAC至少一个检出，排除至少两个髓系/上皮/内皮标记或任一B标记，得到29118个。
成纤维候选要求原作者未注释、固定DPT/COL5A2/PDGFRA/PCOLCE/MFAP5/OGN/SFRP2/THBS2/TNC/FBLN1中至少两个RNA标记检出，同时排除T、髓系、上皮、内皮、壁细胞的多标记及B标记，得到48348个。候选群PDGFRA检出48.72%、COL5A2检出34.82%、FBLN1检出41.18%、TNC检出67.43%。保留“候选”命名，不将所有ACTA2阳性或壁细胞称为纯CAF。
癌细胞限定作者Hypoxic Tumor、Differentiated Proximal-Like Tumor、Growth Signaling-Enriched Tumor三类，并排除跨谱系多标记，得到32949个。三个原作者癌细胞群在基础QC后CA9检出率为34.42%–65.88%；CA9限定敏感性使用其检出子集。未使用CD8位置或空间效应选择癌细胞标签。
身份审计检查RNA和蛋白分开标记；正式CD8、基质和癌细胞纳入以冻结的RNA/作者标签规则为准，不把蛋白原始计数大于0当作可靠阳性阈值。
所有空间测量使用细胞质心坐标，单位μm。KD-tree计算每个CD8最近癌细胞距离、50/100μm成纤维候选数，以及100μm范围癌细胞、其他QC细胞、内皮和壁细胞数。总细胞邻域分母扣除CD8自身。辅助内皮、壁细胞按各自至少两个RNA标记及跨谱系排除计算，与成纤维候选不互换。
主模型对暴露、距离和连续协变量取秩，通过最小二乘去除局部癌细胞数、总体QC细胞数、内皮数、壁细胞数、CD8 RNA计数、x、y趋势后计算残差相关；未调整模型另存。常数协变量省略并记录，设计秩不满或残差自由度小于10时不估计。网格原点为全部QC坐标最小值，500μm宽度与30细胞资格在效应计算前固定。
冻结记录为resources/IDENTITY_FROZEN.json与SPATIAL_PROTOCOL_FROZEN.json。后者于2026-10-03T16:01:55+08:00冻结；首次空间分析完成于16:02:38。没有根据方向改半径、阈值或候选定义。

五、讨论：本阶段补到了哪一段，尚不能解释哪一段
本阶段直接测到了“CD8与癌细胞在组织中相隔多远”，比解离单细胞的总CD8比例更贴近空间排斥问题。基质候选丰富邻域与较远距离的原始方向是真实可描述结果，应保留；调整后较弱并不抹去这一空间分布事实。
密度调整回答的是在局部癌细胞数量、组织细胞密度及其他预设变量相近时是否仍有剩余关联，与原始空间分布是不同问题。调整后减弱提示两者与局部组织构成存在明显重叠，不能仅凭减弱就断言密度是唯一原因，亦不能排除密度本身处于某种生物过程路径中。这里没有处理实验或时间序列，无法把这些解释分开。
距离测量对象是被作者明确标注且通过资格的癌细胞质心；部分真实癌细胞可能未被该注释/门槛捕获。它不是细胞膜接触、病理确认肿瘤巢边界、胶原纤维几何、组织刚度或真实穿越速度。未绘制新的病理边界，也未把低癌细胞邻域自行标成肿瘤外间质。因此本轮尚不能给出严格病理“巢内/间质”面积归一CD8密度比。
RNA面板缺少TYMP、SPP1、COL1A1、COL3A1、POSTN和FAP；这些缺测没有填零、插补或改用别名冒充。此数据可以检验基质细胞—CD8空间接近这一末端，但不能直接验证TYMP或SPP1的局部来源，也不能重现Phase308全基质执行面板。477基因与旧465534细胞版本分别记录，不视为两个独立患者。
Phase308宿主Spp1缺失后成纤维基质程序下降，以及既有人肾癌TYMP–SPP1上游阳性继续保留。本阶段与它们是不同数据的分段证据，不能拼成同患者、同处理的完整因果链。当前新增的是直接空间测量及描述性方向，未新增稳健的功能屏障确认。

六、更大或多患者数据入口的实际核查
Proseg发表的四患者RCC数据：已取得官方Figshare清单，10个同名转录本文件，尚缺本次可核实的文件到患者映射和已注释细胞对象，未将文件数当患者数。暂未下载大型转录本或臆配样本。
Nicheverse预印本93患者/173样本：已取公开图谱目录，确有173条RCC样本预览元数据；逐细胞浏览器配置指向data.nicheverse.org/v1，两次读取索引均出现SSL连接失败。因此没有取得全队列坐标/注释/患者映射，不能称大队列已分析。项目原有core21教学对象7824细胞只有sample_id、计数和坐标，缺作者注释，未凭预览图猜标签。
这些入口的状态完整记录于SOURCE_FEASIBILITY_AUDIT.tsv；公开网页、预览图与可复算原始数据明确区分。没有联系作者、提交受控访问申请或要求客户新增湿实验。
后续真正的增量应优先为有患者映射、癌细胞/基质/CD8注释及坐标的多患者肾癌对象，或同切片病理区室/基质结构数据。现有单切片结果不再通过换半径或择区域寻找阳性。

七、独立验收与交付
{v['checks']}项独立检查通过：重新从H5按另一方向存储重建全部标记、QC与身份；全部CD8用BallTree复算邻域计数和最近距离，另取固定64个细胞用暴力欧氏距离复核；用QR投影重算8个全切片模型和全部空间块；复算高低邻域对比、资格、样本量和P/FDR空值。各正式阶段退出码均为0。
权威目录：{P}
轻量包包含本报告、所有汇总/区域结果、身份标记审计、完整协议、代码、来源清单与校验记录。69万细胞原始矩阵、全量身份与逐CD8坐标明细保留在服务器，交付中提供其SHA256位置清单；旧对象、正式文章、Phase302客户Word及Phase308交付未覆盖。本轮无新图。

来源链接
'''+ '\n'.join(f"{s['source']}：{s['url']}" for s in sources)+'\n'
(pkg/(title+'.txt')).write_text(report,encoding='utf-8')
body=''.join('<section><pre>'+html.escape(part)+'</pre></section>' for part in report.split('\n\n'))
(pkg/(title+'.html')).write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>'+title+'</title><style>body{font:17px/1.8 system-ui,"Microsoft YaHei",sans-serif;color:#203545;background:#f3f6f8}main{max-width:1100px;margin:24px auto;padding:28px;background:white}pre{font:inherit;white-space:pre-wrap;overflow-wrap:anywhere}section{padding:12px 0;border-bottom:1px solid #ddd}table{border-collapse:collapse;font-size:14px;width:100%}td,th{padding:8px;border:1px solid #ccd7df}h1{font-size:27px}</style><main><h1>肾癌逐细胞空间检验</h1><p>1个供者 · 实际细胞坐标 · 不作细胞级推断</p>'+display.to_html(index=False,float_format=lambda x:f'{x:.6f}',border=0)+body+'</main></html>',encoding='utf-8')
server_files=[]
for f in [*list((P/'downloads').iterdir()),*list((P/'tables').glob('*.gz'))]:
 if f.is_file():server_files.append({'absolute_path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(P/'qc/SERVER_LARGE_ASSET_MANIFEST.json').write_text(json.dumps(server_files,indent=2))
(P/'qc/SOFTWARE_VERSIONS.json').write_text(json.dumps({name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','scikit-learn','h5py']},indent=2))
for directory in ['tables','resources','qc','code']:
 (pkg/directory).mkdir(exist_ok=True)
 for f in (P/directory).iterdir():
  if f.is_file() and not f.name.endswith('.tsv.gz') and f.name not in ['DELIVERY_VALIDATED.json','LOCAL_DELIVERY_VERIFIED.json']:shutil.copy2(f,pkg/directory/f.name)
(pkg/'logs').mkdir(exist_ok=True)
for f in (P/'logs').glob('*.exit'):
 if f.stem!='deliver':shutil.copy2(f,pkg/'logs'/f.name)
records=[{'path':str(f.relative_to(pkg)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted(pkg.rglob('*')) if f.is_file() and f.name!='SHA256_MANIFEST.json']
(pkg/'SHA256_MANIFEST.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
archive=D/(title+'.zip')
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(pkg.rglob('*')):
  if f.is_file():z.write(f,f.relative_to(pkg))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for row in records:assert hashlib.sha256(z.read(row['path'])).hexdigest()==row['sha256']
delivery={'status':'PASS','created':now,'archive':str(archive),'archive_name':archive.name,'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'manifest_files':len(records),'archive_files':len(records)+1,'CRC':'PASS','all_file_SHA256':'PASS','analysis_validation':v}
(P/'qc/DELIVERY_VALIDATED.json').write_text(json.dumps(delivery,ensure_ascii=False,indent=2))
(P/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_VALIDATED','time':now,'analysis_jobs':'none','new_cellular_spatial_case':True,'new_robust_barrier_support':False,'effective_donors':1,'local_delivery':'pending','delivery':delivery},ensure_ascii=False,indent=2))
msg=f'''## Phase309终态（2026-10-03）：肾癌逐细胞癌细胞接近检验完成
用户“继续尝试下一步”。新10x RCC XOA4.0同切片RNA/蛋白输入100995534bytes代理下载，1供者690322细胞；原始矩阵477RNA/27蛋白，TYMP/SPP1未测量，不能称其直接验证。QC609525，作者T组+CD8/T RNA及跨谱系排除29118 CD8，未注释多标记成纤维候选48348，作者明确癌细胞+排除32949（CA9限定18314）。身份/几何在效应前冻结，原对象只读。
空间末端：主200μm癌周区28003CD8，50μm成纤维候选高/低邻域癌细胞距离中位63.23/37.33μm，20μm癌接近比例9.79%/29.06%。原始rho0.205639；局部癌/总细胞/内皮/壁细胞密度+RNA+x/y调整后rho0.022175。100μm调整0.030228；CA9限定50/100调整0.091320/0.083405。236块126正110负，中位0.012441；CA9块58正63负，中位−0.010522。分离现象保留，当前无强一致的条件关联；密度调整与总空间分布为不同估计，不据此断言唯一混杂或因果中介。
1供者，无细胞/区域P或FDR。质心距离不等于膜接触/病理肿瘤巢边界/胶原结构/穿越实验。全部8模型及区域/不可评估完整保留，不继续换尺度挑区域。旧TYMP–SPP1上游阳性、Phase308 Spp1KO基质程序下降和此前空间/免疫阴性均保持。
来源增量：Proseg4患者Figshare25685961 10同名transcripts文件，未核清文件-患者映射，不猜配/未下大文件。Nicheverse目录173样本（论文93患者）已取，data.nicheverse.org/v1/index.json两次SSL失败，未取得全队列坐标。旧core21无作者标签，不据预览猜注释。下一有效输入是带患者映射的肾癌坐标/注释或同切片病理区室与基质结构，不再同一切片换阈值。
权威目录{P.name}；原H5身份与全CD8 BallTree/64点暴力几何、QR模型和全部区域等{v['checks']}项检查PASS，正式阶段exit0。报告/完整轻量包delivery/{archive.name}；ZIP CRC和逐文件SHA通过，本地交付见qc/LOCAL_DELIVERY_VERIFIED.json。无计算后台，无新图，未覆盖Word/正式文章/旧对象，未联系作者。

'''
backup=P/'resources'/('root_before_final_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'));backup.mkdir()
updates={}
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 old=(R/name).read_bytes();updates[name]=(old,msg.encode()+old)
name='RESULTS_MANIFEST.tsv';old=(R/name).read_bytes();assert b'phase309_delivery\t' not in old
entries=[('phase309_report','report',pkg/(title+'.html'),'Single-donor spatial case; raw separation, adjusted weak and regional heterogeneity retained'),('phase309_spatial','table',P/'tables/SPATIAL_ASSOCIATIONS_ALL.tsv','8 predeclared models;no biological P/FDR'),('phase309_validation','validation',P/'qc/FINAL_VALIDATION.json',f"{v['checks']} original-count and independent geometry checks"),('phase309_delivery','archive',archive,'CRC/all-file SHA256 validated')]
rows=''.join('\t'.join([rid,'Phase309',kind,str(path),now,'PASS',str(P/'qc/FINAL_VALIDATION.json'),notes])+'\n' for rid,kind,path,notes in entries);updates[name]=(old,old+(b'' if old.endswith(b'\n') else b'\n')+rows.encode())
for name,(old,new) in updates.items():
 f=R/name;(backup/name).write_bytes(old);assert f.read_bytes()==old;t=f.with_name(f.name+'.p309.final.tmp');t.write_bytes(new);assert f.read_bytes()==old;os.replace(t,f)
print(json.dumps(delivery,ensure_ascii=False,indent=2),flush=True)
