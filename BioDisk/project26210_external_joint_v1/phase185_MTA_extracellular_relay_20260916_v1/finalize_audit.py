from pathlib import Path
import json,hashlib,shutil,datetime,xml.etree.ElementTree as E,csv
root=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
r=root/'phase185_MTA_extracellular_relay_20260916_v1'
report="# Phase185：MTA胞外中转机制与SPP1端点核查\n\n## 一、结论\n\n本轮找到有直接培养实验支持的候选解释：**MTA经成纤维细胞处理后改变胞外环境，进而改变髓系细胞状态。** 这为“增加中间环节”提供了依据，但尚未证明该过程提高巨噬细胞SPP1，不能将其画成已打通的链条。\n\n继续保留两条候选路线：\n\n1. MTA—A2B/STAT3功能依赖—巨噬细胞响应：继承前期分段证据，STAT3转录代理的公共阴性结果不变。\n2. MTA—成纤维细胞代谢/分泌变化—巨噬细胞响应：本轮补充原始文献依据，SPP1接收端待检验。\n\nIL10不再作为通用必经中介；PRMT5本轮不升格为正向中介。免疫功能和空间组织仍为两个独立、暂缓的下游终点。\n\n## 二、直接MTA与条件培养液产生不同响应\n\n[Valera等，PNAS 2023，原始研究](https://doi.org/10.1073/pnas.2311674120) 的Fig.2、4、5分别涵盖代谢检测、RAW264.7和THP-1条件培养/共培养实验。\n\n- 成纤维细胞消耗补充MTA，并改变胞外腺嘌呤、次黄嘌呤水平；有LC–MS支持。\n- RAW264.7中，直接MTA条件主要表现为NOS2上升、ARG1变化小；加入MTA的成纤维—MTAP缺失肿瘤共培养来源培养液则表现为ARG1上升、NOS2下降。\n- THP-1中，纯MTA与包含成纤维细胞的培养条件也呈现不同的TNF/CD163方向。\n- 本文列出的髓系qPCR端点是ARG1、NOS2、TNF、CD163；正文XML未出现SPP1/osteopontin读出。不能用这些标志物代替SPP1，也不能直接认定免疫抑制功能已经被这些qPCR证明。\n- 培养体系同时包含代谢变化和其他分泌变化，尚不能从这组实验指定某个单一代谢物或细胞因子承担全部作用。不能由腺嘌呤/次黄嘌呤结果推成“腺苷—A2B”已被测得。\n\n**解释：** 胞外MTA浓度及其作用对象受到局部摄取和处理影响。该文支持继续检验环境依赖机制；它并不证明ccRCC中已经存在同样中转，更未验证TYMP上游。\n\n## 三、公开数据资格\n\n[GSE249082官方记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249082) 是乳腺成纤维细胞HBF，4个对照文库、4个MTA处理文库。处理为10 µM、24小时；不是巨噬细胞转录组，也不是8名患者。\n\n逐样本元数据已归档。注意其treatment_protocol字段在对照中也复制了统一处理说明；分组须结合title及characteristics中的Untreated/MTA标签，不得仅用该字段重分类。重复标为rep1–4，供体和配对关系未确认，不将编号相同视为配对，不将文库数直接当独立供体数。\n\n作者另外公开了[Zenodo原始分析数据包](https://zenodo.org/records/8130216)，网页列出PNAS_Data.zip约30.5 MB及MD5。此次服务器整包下载长时间未产生可用内容，已停止；网页预览也未成功。因此：\n\n- 已检查论文正文及GEO样本元数据；\n- **尚未完成补充PDF和ZIP内部所有端点的检查**；\n- “正文未报告SPP1”不能扩写成“所有原始文件确定没有SPP1”；\n- 不从不完整下载中生成效应量，不将该包登记为已通过完整性校验。\n\n## 四、另外两类命中不补上该箭头\n\n### PRMT5：保留方向约束\n\n[JITC 2025原始研究](https://doi.org/10.1136/jitc-2024-011299) 中，髓系Prmt5缺失/抑制降低M2相关标志及脂质代谢程序，涉及STAT6–PPARγ；这不是直接MTA处理实验，也未在正文报告SPP1端点。若借“MTA抑制PRMT5”补成正向促SPP1链，至少需要额外的方向及功能证据，不能因C4具有脂质程序就直接接上。该文数据声明为依合理请求获得；列出的GSE5099/GSE154763是复用数据，不登记为新PRMT5扰动数据。\n\n### MTAP缺失与SPP1：先识别细胞来源\n\n[Cells 2020原始研究](https://doi.org/10.3390/cells9020492) 在U251肿瘤细胞NanoString结果中列出SPP1；这既不是髓系细胞，也不是直接MTA暴露。SPP1未进入正文所列的校正P与绝对log2FC≥2同时满足的九基因名单；本轮未取得其补表精确效应，不推断方向、显著性或MTA依赖。该命中不作为巨噬SPP1支持。\n\n## 五、下一步的有限验证方案\n\n先完成中转数据资格及预设小面板核查，不扩大到全基因搜中介。\n\n1. 恢复Zenodo时使用明确总时限/可续传方式，先看目录和qPCR端点；若确有SPP1测量，再核查独立实验重复、处理时序、归一化及原始数值。\n2. GSE249082仅用于识别“MTA改变成纤维输出”的候选环节。预先固定HGF、TGFB1/2/3、IL6、CSF1、WNT5A、PTGS2/PTGES、SPP1；后两类分别注意合成酶表达不等于PGE2分泌、成纤维SPP1不等于巨噬SPP1。同时记录MTAP及生存/整体表达质量。\n3. 重复独立性未确认前只报告组间方向、幅度、分布，不作推断性P/FDR，不训练用于筛选本项目患者的响应签名。\n4. 仅对确有MTA响应、且存在直接髓系受体/阻断文献的候选，检查本项目来源与受体可测性。HGF等预设名称只是待检验面板，不提前称中介。\n5. 若没有候选同时满足“外部扰动响应—髓系作用依据—项目可测性”，停止该分支。保留MTA候选及已成立的STAT3—SPP1关联，不再换评分追求阳性。\n\n## 六、对正文的当前影响\n\n可加入机制讨论：MTA相关髓系响应具有环境依赖性，基质细胞处理可作为待检验的中转机制。暂不新增“成纤维细胞介导MTA诱导SPP1”的Results结论。\n\n前期TYMP相关代谢背景、整体巨噬STAT3—SPP1关联、C4状态定位均保留；公共MTA—STAT3代理不支持稳定正向联系、成熟人巨噬IL10/STAT3扰动的反向结果也完整保留。本轮无新患者统计、无新图、未重跑旧模型。\n"
(r/'MTA_EXTRACELLULAR_RELAY_AUDIT_CN.md').write_text(report,encoding='utf-8')
now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
checks={}
for pmc in ['PMC10756296','PMC11973769','PMC7072758']:
 x=E.parse(r/'sources'/f'{pmc}.xml')
 checks[pmc+'_parse']=True
 checks[pmc+'_identity']=any(e.text==pmc for e in x.iter('article-id'))
s=(r/'sources/GSE249082_samples.soft').read_text()
checks['GEO_eight_samples']=s.count('^SAMPLE = ')==8
checks['GEO_four_untreated']=s.count('!Sample_characteristics_ch1 = treatment: Untreated')==4
checks['GEO_four_MTA']=s.count('!Sample_characteristics_ch1 = treatment: 10 µm MTA for 24 hours')==4
assert all(checks.values()),checks
status={'time':now,'scope':'literature and sample metadata audit; no new inferential analysis','checks':checks,'zip_status':'incomplete retrieval terminated; integrity not validated','zenodo_api_status':'HTTP 504','supplement_pdf':'not inspected','active_jobs':[]}
(r/'VALIDATION.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
with open(r/'SOURCE_INVENTORY.tsv','w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['file','bytes','sha256','status'])
 for p in sorted((r/'sources').iterdir()):
  w.writerow([p.name,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest(),'UNVALIDATED_PARTIAL' if p.suffix=='.zip' else 'ARCHIVED'])
backup=r/'project_record_backups';backup.mkdir(exist_ok=True)
note=f"""## Phase185 MTA胞外中转原始文献与数据资格核查（{now}）
PNAS2023原文显示纯MTA与成纤维/肿瘤条件培养背景产生不同髓系标志响应，支持中转候选；ARG1/NOS2和TNF/CD163不是SPP1。GSE249082为HBF成纤维4+4库、10µM24h，不是髓系/8供体；样本已归档，配对与供体独立性未确认。
三篇XML及8样本资格检查通过；Zenodo整包下载未完成已停止，API504、补充PDF未查，不能声称全部原始端点均无SPP1。PRMT5干预非MTA且方向不支持直接填正向中介；U251的SPP1属于肿瘤细胞。
报告：{r}/MTA_EXTRACELLULAR_RELAY_AUDIT_CN.md
下一断点：有限恢复原始端点清单；预设成纤维输出小面板及重复资格，再决定是否有资格检查项目接收端。无新统计/图/签名，旧正负结果保留；免疫空间继续分开暂缓。无活跃任务。
"""
for name in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']:
 p=root/name
 shutil.copy2(p,backup/name)
 if name!='RESULTS_MANIFEST.tsv':p.write_text(note+'\n'+p.read_text(),encoding='utf-8')
 else:
  with open(p,'a') as f:f.write('\n'+'\t'.join(['phase185_MTA_extracellular_relay','phase185','report',str(r/'MTA_EXTRACELLULAR_RELAY_AUDIT_CN.md'),now,'PARTIAL_SOURCE_ACCESS','9 metadata/XML checks passed; ZIP/SI not inspected','Literature audit complete within accessed sources; no new causal/SPP1 validation'])+'\n')
d=r/'delivery';d.mkdir(exist_ok=True)
for name in ['MTA_EXTRACELLULAR_RELAY_AUDIT_CN.md','VALIDATION.json','SOURCE_INVENTORY.tsv']:shutil.copy2(r/name,d/name)
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in d.iterdir()}
(d/'SHA256.json').write_text(json.dumps(hashes,indent=2))
print(json.dumps({'checks':checks,'delivery_hashes':hashes},ensure_ascii=False,indent=2))
