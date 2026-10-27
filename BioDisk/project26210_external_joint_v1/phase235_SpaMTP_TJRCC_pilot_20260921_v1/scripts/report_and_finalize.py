from pathlib import Path
import json,hashlib,datetime,shutil,zipfile,csv
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent
def fmt(x):
 if pd.isna(x):return '—'
 if isinstance(x,(float,np.floating)):return f'{x:.5g}'
 return str(x).replace('|','/')
def table(d):
 return '\n'.join(['| '+' | '.join(d.columns)+' |','| '+' | '.join(['---']*len(d.columns))+' |']+['| '+' | '.join(fmt(v) for v in row)+' |' for row in d.itertuples(index=False,name=None)])
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 now=datetime.datetime.now().astimezone().isoformat()
 validations={}
 for name in ['RAW_SPECTRA_VALIDATION','SPAMTP_MODULE_VALIDATION','SERIALIZED_OBJECT_VALIDATION','PILOT_DATA_VALIDATION','CONDITIONAL_ASSOCIATION_VALIDATION']:
  o=json.loads((P/(name+'.json')).read_text());assert o['status']=='PASS';validations[name]=o
 qc=pd.read_csv(P/'tables/SPECTRAL_TECHNICAL_QC.tsv',sep='\t')
 hits=pd.read_csv(P/'tables/TARGETED_MASS_CANDIDATES_PASS5PPM.tsv',sep='\t')
 hits=hits.merge(qc[['sample','polarity','nonzero_pixels']],on=['sample','polarity'],validate='many_to_one')
 hits['positive_percent_within_nonzero_pixels']=hits.detected_pixels/hits.nonzero_pixels*100
 hits.to_csv(P/'tables/CANDIDATE_COVERAGE_NONZERO_PIXELS.tsv',sep='\t',index=False)
 pair=pd.read_csv(P/'tables/PAIRING_FEASIBILITY.tsv',sep='\t')
 cond=pd.read_csv(P/'tables/CONDITIONAL_DIRECTION_STABILITY.tsv',sep='\t')
 modules=pd.read_csv(P/'tables/SPAMTP_MODULE_QC.tsv',sep='\t')
 chunks=[
 '# SpaMTP与TJ-RCC空间代谢组首轮试点\n\nPhase235｜2026-09-21｜数据与官方模块可运行；正式跨模态配准尚未冻结\n',
 '''## 结论

已实际开展SpaMTP试点，而非仅作方案评估。通过本机7897代理定向获取R114_T、Y12_T两例完整正/负离子模式，完成原始导出谱解析、候选质量核查、SpaMTP官方导入/TIC归一化/空间GraphPCA与4区表示。4份谱共64222扫描像素，其中41022个非零像素；生物学样本量始终为2例，未计算新的P/FDR。

两例均有与亚精胺、精胺精确质量相符的正离子候选峰；当前作者注释表仅对应532个负离子特征，因此正离子候选仍是质量相容线索。MTA、SAM、SAH等在本轮公开导出峰表、固定加合物及5ppm规则下未出现匹配峰，不等于原始组织中没有这些分子。

SM与ST组织轮廓最佳重合Dice为0.891–0.953，但近似同样合理的朝向不唯一。进一步保留全部11个候选几何映射后，多数分子—RNA关系改变方向；只有个别病例/端点保持同向。当前结果可支持继续核查多胺空间线索，不能直接作为已确证的TYMP—MTA—SPP1空间验证。

## Methods

### 数据取得与试点选择
原研究Figshare24599295缓存清单与Zenodo8063124公开镜像共同提供数据身份：spatial_metabolomics.zip总347579559字节、原MD5 b06916ac6ce28deca12ad1fb45ffc637；空间代谢注释表78880字节、MD5 edadc8a02aebc8024cdb96baa29977ac。注释表整文件MD5通过。数据包目录含75条文件记录，覆盖12个肿瘤样本和2个正常组织样本，各具正、负模式的imzML/ibd。

按完整双极性文件配套及压缩大小选择最小的两个肿瘤样本R114_T、Y12_T，选择先于查看原始候选峰。整包直连约20KB/s，随后使用本机127.0.0.1:7897代理下载两个字节范围，分别10922725和14052170字节，再上传已授权sgs10。提取9文件均通过ZIP条目CRC32、长度校验；4个imzML/ibd配对进一步通过文件内SHA-1、UUID校验，保存全部SHA256。未下载完整ZIP，因此不声称整包MD5已通过；原直连片段保留。

### 谱数据与质量候选
公开文件为Cardinal2.8.0导出的continuous imzML/ibd；负离子共同532峰，正离子共同767峰。按XML声明的64-bit m/z与32-bit intensity及外部偏移读取，核验谱数、统一轴、像素坐标唯一、非负有限强度。所有TIC>0像素进入技术对象，零强度像素仅保留在原始矩阵/QC中。负离子轴与作者532行注释逐行一致（仅小数位舍入差）。

计算前锁定MTA、SAM、SAH、甲硫氨酸、亚精胺、精胺、腐胺、鸟氨酸、胸苷、胸腺嘧啶、尿苷、尿嘧啶、腺苷及公式。正离子检验M+H、M+Na、M+K、M+NH4、M+2H（双电荷）；负离子检验M-H、M+Cl、M+FA-H、M-H2O-H；绝对质量误差≤5ppm。完整未匹配最近峰同样保留，不增加误差窗追候选。名称来自预列分子，仅表示质量相容，未将质量匹配提升为结构鉴定。

### SpaMTP实际执行
固定官方代码提交323511b500e00bfe054365c735a974e4c59a6b81，版本1.1.0.9000；6份原始R源文件均经Git blob SHA-1及本地SHA256核对。执行方式是加载未修改的官方源模块，使用既有Seurat5.4.0/SeuratObject5.3.0环境；没有安装完整SpaMTP包，也没有改共享R环境。

实际运行ReadSM_mtx、NormalizeSMData、RunSpatialGraphPCA、GetKmeanClusters。先按原像素建立SPM assay/FOV并作TIC归一化。技术空间表示固定4×4原始SM格点分箱，每箱取非零像素强度均值；再按TIC归一化、非恒定峰标准化、GraphPCA lambda0.5/20成分/platform ST默认4邻居，k-means固定4区/nstart20/seed262235。产生4个像素对象与4个分箱空间对象，共8个新版本RDS。区域编号仅为各切片内部技术区域，跨样本不直接对应。无图。

官方导入矩阵与输入逐元素一致，TIC输出与独立公式一致，8个RDS重读有效；GraphPCA输出有限、图对称且无自环。默认ST邻接包含少量跨1格以上边（最大约2.24个分箱格点），本轮图用于技术表示，没有将其改称经过组织空白约束的精确空间传播模型。

### 配准可行性与条件性对应
核对同名ST输入、原in_tissue标记、H&E及坐标；R114_T为1869个QC spot，Y12_T为1661个QC spot。SM文件像素尺度均写1且无可靠物理单位，先用格点处理，不直接当μm。组织轮廓以SM非零TIC和独立ST in_tissue构建，不使用TYMP/SPP1或候选峰强度优化。检查8种旋转/翻转及固定平移网格，记录最佳Dice和朝向分差。

进一步在查看新跨模态效应前固定：保留各朝向最优平移且Dice距全局最佳≤0.03的全部几何假设；归一化轮廓空间最近SM像素距离≤1.5个SM间距才纳入。以既有ST4×4分块的算术均值作粗区域描述，每块≥10配对spot、至少8块、至少50配对spot且覆盖≥50%QC spot；各候选峰对TYMP RNA、MTA转录背景、SPP1绝对q05、巨噬内SPP1比例及去SPP1程序计算Spearman效应。11个候选映射、145条效应记录中140可估计，全部保留。仅报告每病例的映射效应范围与方向变化，不把映射数量/像素数当生物重复，范围不是置信区间。
''',
 '## Results 1：两例正负模式及原始文件完整可用\n\n'+table(qc[['sample','polarity','pixels','features','nonzero_pixels','zero_pixels','uuid_verified','ibd_sha1_verified']]),
 '## Results 2：多胺与核苷候选峰提供可继续核查的实测强度入口\n\n亚精胺相容峰m/z146.165169，误差−0.0377ppm；精胺相容峰203.223203，误差0.8845ppm。两例均有检出，Y12_T的非零扫描像素中覆盖较广，R114_T更局限。鸟氨酸相容峰133.097268较稀疏。负离子279.038980与尿苷/伪尿苷等异构候选相关，不能唯一指定尿苷；111.020015对应尿嘧啶等作者多候选注释。\n\n'+table(hits[['sample','polarity','molecule','observed_mz','ppm','detected_pixels','nonzero_pixels','positive_percent_within_nonzero_pixels']])+'\n\nMTA、SAM、SAH、甲硫氨酸、腐胺、胸苷、胸腺嘧啶及腺苷在本轮固定质量检查中未匹配；完整目标/加合物及最近非匹配峰见TARGETED_MASS_CANDIDATES_ALL.tsv。该结论限定于公开导出的有限峰表。',
 '## Results 3：SpaMTP技术模块与对象均通过运行检查\n\n'+table(modules[['label','binned','pixels_or_bins','features','normalization_max_error','graph_edges','max_neighbor_distance_grid','graph_components']])+'\n\n各4区多胺候选的原始强度/TIC值和检出比例已保存于CANDIDATE_REGION_DESCRIPTIVE.tsv。区域由同一SM峰矩阵获得，因此此表为描述，不计算“独立区域富集”P值，不把技术区域命名为已确定肿瘤/巨噬区域。',
 '## Results 4：轮廓匹配支持继续配准，但朝向仍有歧义\n\n'+table(pair[['sample','polarity','ST_QC_spots','SM_nonzero_pixels','best_mask_dice','orientation_margin','best_flip','best_rotation90','mapping_accepted']])+'\n\n四份谱最佳形态候选均为不翻转、90°旋转，但次优朝向与最佳仅差0.010–0.023。现有SM包没有独立形态图像或现成SM→ST变换，不能由高Dice独自确立对应坐标。',
 '## Results 5：条件性区域对应中出现个别方向线索，多数关系对配准敏感\n\nR114_T亚精胺候选峰与巨噬内SPP1比例在5个保留朝向均为正，rho范围0.1192–0.5235。Y12_T精胺候选峰与TYMP RNA、SPP1绝对q05及去SPP1程序在2个保留朝向均为正，其rho范围分别0.0294–0.2059、0.0176–0.4088和0.0941–0.3206。\n\n这些分别涉及不同分子、不同端点的一名患者，不能拼成两例复现同一联系。其他主要关系，包括亚精胺—TYMP/MTA转录背景及多项SPP1关联，随候选配准可变号。所以下一步优先解决配准，而不是从中选择给出正向结果的朝向。\n\n'+table(cond[cond.molecule.isin(['Spermidine','Spermine'])][['sample','molecule','RNA_target','retained_mappings','estimable_mappings','rho_min','rho_max','positive_mappings','negative_mappings','direction_stable']]),
 '''## Discussion与下一步

本轮确定了真实空间代谢数据入口、可运行的SpaMTP核心模块，以及与前面多胺/核苷主线相关的质量候选峰。第四部分可以继续推进区域层对应，不需要设置完整因果中介或临床显著性门槛。当前尚未获得可靠的MTA/SAM空间分子身份和唯一SM–ST坐标；已有Phase234的RNA区域联系仍作为原有证据保留。

下一优先任务为：寻找/核对SM切片形态图、原配准坐标或可独立识别的组织轮廓标志，冻结两例坐标后再进行相同候选的区域对应。可保留亚精胺、精胺候选作为多胺分支，同时核对异构物/加合物；MTA未进入公开峰表不应被解释为组织不存在MTA，也不能用精胺替代MTA实测结论。

若无法确定唯一朝向，条件性方向范围仍可作为内部可行性资料；不能在正文只挑最佳正相关。2患者仅作病例描述，像素/分箱/候选配准不构成生物重复。SPP1空间变量仍为旧映射SPP1-state，不能改名新C4。未进行新P/FDR、因果/治疗模型或绘图。

## 文件与验收

- results/SpaMTP_objects：8个新RDS；results/matrices：4份原强度与坐标NPZ及对应输入CSV。
- metadata：原作者清单、注释表身份、固定官方提交、原始源模块、预定候选与技术参数。
- tables：谱QC、全部m/z、全部候选/未匹配、技术区域值、配准候选、全部条件性效应与方向稳定性。
- audit：9文件CRC/SHA256、4个ibd内部SHA-1/UUID、8对象读回、原注释轴、ST条码坐标和333旧RDS保留检查。
- FINAL_VALIDATION.json：本轮最终验收；完整源代码在scripts。
- 本轮9个已提取原文件完整，347.6MB整包未完整下载；旧直连片段仍在downloads。其DOWNLOAD_EXIT_STATUS=1对应主动中止的整包路线，由PILOT_DOWNLOAD_VERIFIED.json与本轮完成记录说明后续成功的定向下载。

## 资料来源

[原研究与空间代谢组方法](https://www.nature.com/articles/s41588-024-01662-5)；[原作者数据仓库](https://zenodo.org/records/8063124)；[SpaMTP原论文](https://www.nature.com/articles/s41592-026-03140-8)；[固定代码提交](https://github.com/GenomicsMachineLearning/SpaMTP/tree/323511b500e00bfe054365c735a974e4c59a6b81)。
质量公式核对：[MTA](https://pubchem.ncbi.nlm.nih.gov/compound/methylthioadenosine)、[SAM](https://pubchem.ncbi.nlm.nih.gov/compound/S-adenosylmethionine)、[SAH](https://pubchem.ncbi.nlm.nih.gov/compound/Adenosylhomocysteine)。质量相容并非独立结构鉴定。
''']
 (P/'PILOT_RESULTS_CN.md').write_text('\n\n'.join(chunks),encoding='utf-8')
 summary='''# SpaMTP首轮试点摘要
Phase235｜2026-09-21

已通过本机7897代理取得R114_T、Y12_T两例正负模式完整数据，9个文件完成校验并上传sgs10；原慢速整包路线已停止，片段保留。

实际跑通SpaMTP官方源模块：导入、TIC归一化、空间GraphPCA及4区技术表示，得到8个独立新对象，重读检查通过。使用锁定官方源文件，没有声称完整SpaMTP包已安装。2患者、4份谱、64222扫描像素，其中41022非零；不作患者推断P值，不生成图。

可继续的线索：两例均有亚精胺和精胺精确质量相容峰；Y12_T信号覆盖较广。MTA/SAM在本次有限公开峰表和固定5ppm规则下未匹配，不能解释为组织中不存在。

配准轮廓最佳Dice约0.89–0.95，但朝向存在歧义。保留11个形态相容映射后，R114_T的亚精胺候选—巨噬内SPP1比例、Y12_T的精胺候选—TYMP/SPP1绝对量有病例级同向线索；多数其他关系会随配准变号。不同病例/分子/端点不能合并成一致复现。

下一步优先核对SM–ST共同组织标志与配准，再把多胺候选和SPP1区域对应接入主线。完整正负范围、未匹配项、数据资格与源脚本见PILOT_RESULTS_CN.md及附件。
'''
 (P/'SUMMARY_CN.md').write_text(summary,encoding='utf-8')
 final=dict(status='PASS',time=now,scope='technical SpaMTP source-module pilot and conditional alignment exploration complete;validated cross-modal mapping not complete',samples=2,raw_files=9,acquisitions=4,scan_pixels=64222,nonzero_pixels=41022,Seurat_objects=8,mapping_hypotheses=11,conditional_effects=145,estimable_effects=140,new_P_values=0,new_figures=0,old_RDS_unchanged=333,full_archive_verified=False,selected_entries_verified=True,full_SpaMTP_package_installed=False,official_source_modules_executed=True,validation_sources=list(validations))
 (P/'FINAL_VALIDATION.json').write_text(json.dumps(final,ensure_ascii=False,indent=2))
 status=dict(status='TECHNICAL_PILOT_COMPLETE',time=now,next_step='resolve cross-modal registration using independent tissue landmarks;do not select orientation by molecular correlation',formal_crossmodal_alignment_frozen=False,new_figures=0)
 (P/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2));(P/'EXIT_STATUS').write_text('0\n');(P/'COMPLETE_VERIFIED').write_text(now+'\n')
 # Exact hashes for every generated object; no need to transfer raw matrices to local review package.
 model_files=[dict(path=str(f.relative_to(P)),bytes=f.stat().st_size,sha256=sha(f)) for f in (P/'results/SpaMTP_objects').glob('*.rds')]
 (P/'audit/MODEL_SHA256.json').write_text(json.dumps(model_files,indent=2))
 files=[P/'PILOT_RESULTS_CN.md',P/'SUMMARY_CN.md',P/'FINAL_VALIDATION.json',P/'STATUS.json',P/'EXIT_STATUS',P/'COMPLETE_VERIFIED']
 files += [P/(k+'.json') for k in validations]
 for sub in ['tables','scripts','audit','environment']:
  files.extend(f for f in (P/sub).rglob('*') if f.is_file() and '__pycache__' not in str(f) and 'python_lib' not in str(f))
 files.extend(P/'metadata'/n for n in ['TECHNICAL_PILOT_FREEZE.json','SPAMTP_EXECUTION_FREEZE.json','PAIRING_DIAGNOSTIC_PLAN.json','CONDITIONAL_ASSOCIATION_FREEZE.json','SAMPLE_ARCHIVE_INVENTORY.json','ZIP_DIRECTORY.json'])
 manifest=[dict(path=str(f.relative_to(P)),bytes=f.stat().st_size,sha256=sha(f)) for f in files]
 (P/'OUTPUT_SHA256.json').write_text(json.dumps(manifest,indent=2));files.append(P/'OUTPUT_SHA256.json')
 delivery=P/'delivery'
 for n in ['SUMMARY_CN.md','PILOT_RESULTS_CN.md']:shutil.copy2(P/n,delivery/n)
 archive=delivery/'Phase235_SpaMTP_pilot_statistics.zip'
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for f in files:z.write(f,f.relative_to(P).as_posix())
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None
  for f in manifest:assert hashlib.sha256(z.read(f['path'])).hexdigest()==f['sha256']
 dh=[dict(name=f.name,bytes=f.stat().st_size,sha256=sha(f)) for f in delivery.iterdir() if f.name!='DELIVERY_SHA256.json']
 (delivery/'DELIVERY_SHA256.json').write_text(json.dumps(dh,indent=2))
 record=f"""## {now} Phase235：SpaMTP–TJ-RCC首轮技术与条件配准试点完成
用户明确开始尝试，随后要求慢速下载使用代理。已通过本机7897代理定向取得R114_T/Y12_T两例正负模式（按配套完整及压缩量最小选，不按效果），9文件CRC32/SHA256、4对ibd内部SHA-1/UUID通过。14样本75文件总目录已盘点；整包347.6MB未完整下载，原7.47MB片段保留，慢速整包任务已停止；已有其他下载任务未操作。
固定SpaMTP官方提交323511b500e00bfe054365c735a974e4c59a6b81，实际执行未修改官方ReadSM_mtx/NormalizeSMData/RunSpatialGraphPCA/GetKmeanClusters模块，非完整包安装。2患者4谱64222像素/41022非零，8个新Seurat对象通过重读。官方源哈希、原始数据至对象/归一化一致，旧333RDS及ST保持，无图/新P。
负532/正767公共导出峰；两例亚精胺m/z146.165169、精胺203.223203有5ppm质量候选，非结构确证；MTA/SAM等固定候选未匹配，不等于组织中不存在。正离子未有作者532行负离子注释表的对应注释，不伪造分子身份。
组织轮廓最佳Dice0.891–0.953但朝向分差仅0.010–0.023，尚无独立变换/组织标志确认。按形态最佳差<=0.03保留11候选映射，4x4粗区域145效应140可估计，仅病例方向。R114亚精胺候选—SPP1巨噬内比例5映射均正；Y12精胺候选—TYMP/绝对SPP1/程序2映射均正，但其他主要关系可变号，不能拼成统一两例验证或挑正向朝向。
权威{P.name}/PILOT_RESULTS_CN.md；完整源脚本、表、8对象及未匹配/反向留档；FINAL_VALIDATION PASS，EXIT_STATUS0。当前断点：首轮技术及条件性探索已完成；下一优先为独立组织形态/配准依据核对后冻结共同坐标，才形成正式跨模态区域对应。不要把技术4区直接命名巨噬/肿瘤，不把旧SPP1-state叫C4；MTA和原Phase234结果保留。
"""
 for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
  f=R/name;b=P/'root_record_backups'/name
  if not b.exists():shutil.copy2(f,b)
  old=f.read_text()
  if 'Phase235：SpaMTP–TJ-RCC首轮' not in old:f.write_text(record+'\n\n'+old)
 with (R/'RESULTS_MANIFEST.tsv').open('a',newline='') as f:
  w=csv.writer(f,delimiter='\t')
  for name in ['PILOT_RESULTS_CN.md','SUMMARY_CN.md','tables/CONDITIONAL_DIRECTION_STABILITY.tsv','FINAL_VALIDATION.json','delivery/Phase235_SpaMTP_pilot_statistics.zip']:
   w.writerow(['phase235_'+Path(name).stem,'phase235_SpaMTP_pilot','report' if name.endswith('.md') else 'artifact',str(P/name),now,'PASS',str(P/'FINAL_VALIDATION.json'),'technical pilot;conditional geometry not validated mapping;no new P/figures'])
 print(json.dumps(dict(final=final,delivery=dh),ensure_ascii=False),flush=True)
if __name__=='__main__':main()
