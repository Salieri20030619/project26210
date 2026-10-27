from pathlib import Path
import json,hashlib,datetime,subprocess,os,csv,io,shutil
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent
def write(p,s):
 with p.open('x',encoding='utf-8') as f:f.write(s)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
with open(P/'logs/memory_gate.log','a') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,check=True)
s1=json.loads((P/'step1_receiver_network/COMPLETION.json').read_text())
s2=json.loads((P/'step2_immune_phenotype/COMPLETION.json').read_text())
s3=json.loads((P/'step3_reference_alignment/CONCORDANCE_COMPLETE.json').read_text())
ref=pd.read_csv(P/'step3_reference_alignment/REFERENCE_GENE_CONCORDANCE.tsv',sep='\t')
sel=ref.source_P.lt(.05);ss=sel&~ref.definition_gene
s3.update(source_P_lt05_genes=int(sel.sum()),source_P_lt05_concordant=int(ref.loc[sel,'median_concordant'].sum()),source_P_lt05_nondefinition_genes=int(ss.sum()),source_P_lt05_nondefinition_concordant=int(ref.loc[ss,'median_concordant'].sum()),tabula_status='TARGET_STATES_ABSENT_FROM_DOWNLOADED_MODEL',model_labels=45,model_features=6436,website_access='curl connection timeout after10001ms',status='COMPUTED_CONCORDANCE_AND_CLOSED_STATE_TRANSFER_INELIGIBLE')
write(P/'step3_reference_alignment/COMPLETION.json',json.dumps(s3,indent=2))
ref.assign(source_nominal_filter=sel,source_nominal_nondefinition_filter=ss).to_csv(P/'step3_reference_alignment/SOURCE_P_FILTER_SENSITIVITY.tsv',sep='\t',index=False)
partial=P/'resources/cd8ippred/R/sysdata.rda'
if partial.exists():partial.rename(partial.with_suffix('.rda.partial'))
source_audit=[
('NicheNet','https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md','official gene-set/background/ligand workflow; exact3human prior MD5 match'),
('NicheNet prior','https://zenodo.org/records/7074291','official MD5 listing; API timeout resolved using verified web checksums'),
('RCC histology','https://www.nature.com/articles/s41379-021-00864-0','manual IT/PT; actual lymphocyte image quantities; supplementary PDF incomplete not used'),
('CD8 classifier','https://github.com/bedapub/cd8ippred','commit9d226f57e710dde49c70188de6a0413ed81abf83;12filesverified;raw countVST required;model sysdata incomplete and reference distribution not downloaded'),
('TabulaTIME','https://github.com/wanglabtongji/TabulaTiME','commita07d177a7c3be0fec2d96821e156ff43855a733f;metadata45broadclasses not targetstates'),
('TabulaTIME fulltext','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12643937/fullTextXML','complete parsed original article XML saved'),
('Tomko2018','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6113240/fullTextXML','complete parsed original article;Table1 genuine27genes25up2down;distinct24alignmentcorrelated proteins'),
('DDR1','https://pmc.ncbi.nlm.nih.gov/articles/PMC8839149/','27gene reference traced toTomko2018;exactDDR1direction/ssGSEAimplementation not reconstructed')]
pd.DataFrame(source_audit,columns=['resource','source_url','qualification']).to_csv(P/'resources/SOURCE_AUDIT_FINAL.tsv',sep='\t',index=False)
write(P/'logs/NETWORK_ATTEMPTS_RECORDED.txt','Zenodo API expression-preparation process stalled after all expression inputs written; confirmed own process command then SIGTERM; exit143 retained. Official web MD5 all3 matched, input/score checks pass.\nRCC supplemental curl45s:exit28,2376767bytes. Resume60s:exit28,5577441of49352980remaining bytes. Total expected51729747bytes;pdfinfo missingxref;retained.partial.\nTabulaTIME website curl: connection timeout10001ms.\nCD8 sysdata: original fetch audit preserves timeout;partial renamed. No prediction generated from incomplete input.\n')
ranks=pd.read_csv(P/'step1_receiver_network/LIGAND_ACTIVITIES_ALL.tsv',sep='\t')
dim=pd.read_csv(P/'step1_receiver_network/NETWORK_DIMENSIONS.tsv',sep='\t').iloc[0]
table='| 配体 | 全部配体排名 | 巨噬表达配体排名 | Pearson模型分值 |\n|---|---:|---:|---:|\n'
for name in ['TGFB1','IL1B','IL1A','SPP1']:
 z=ranks[ranks.test_ligand.eq(name)].iloc[0]
 table+=f'| {name} | {int(z.rank_all)}/797 | {int(z.rank_macrophage)}/190 | {z.pearson:.6f} |\n'
report=f"""# Phase295：文献方法顺序尝试完成记录

当前权威：{R}。计划在计算前记录于CURRENT_STATE_AND_SEQUENTIAL_PLAN.md。三项按序尝试已产生终态；“完成”指本轮试验与资格核查完成，不代表TYMP到免疫屏障因果链已证明。

## 接收靶基因网络优先指向TGFB1候选输出
7位独立患者、4个研究，同患者同文库的原始Pericyte/SMC类别内，冻结活化状态与非活化状态各至少20细胞。一患者仅留一文库，按较小状态细胞数最大、总细胞数及字典序选择，不按表达结果选择。完整资格、来源表达资格、所有基因效应均保存。
接收背景{int(dim.background)}基因，描述性目标集1,479基因；预设中位log2(CPM+1)差≥0.25且≥70%患者正向，排除Phase286身份及Phase289全部定义/程序基因。靶集来自全基因变化，没有只挑胶原或SPP1下游。
全部受体合格配体797个，巨噬表达合格配体190个，完整排序如下文件所载：step1_receiver_network/LIGAND_ACTIVITIES_ALL.tsv。关键候选：

{table}
TGFB1在当前接收状态模型中优先于SPP1。探索性解释应区分“SPP1相关巨噬状态的候选输出”与“SPP1分子直接驱动”。C4表达标记仅说明配体表达资格，不证明C4专属性、SPP1状态特异性或分泌蛋白量。全排名及候选靶连接同时保留，不删除IL1/SPP1低排序结果。
官方NicheNet人类v2三个先验文件全部MD5匹配；R函数Pearson用独立二元目标向量逐配体重算误差<1e-10。32项输入/结果检查通过，网络/验证脚本退出0。预处理中的Zenodo API请求超时后定点终止，退出143保留，官方网页校验与已写输入核查完成恢复；不是未处理的统计失败。

## 免疫表型分类器的原始方法已定位，当前输入不合格
已取得cd8ippred作者代码和三个类别系数文件等12个校验通过的资源。原始代码要求raw counts经DESeq2大小因子和训练集固定离散函数转换，并比较训练参考分布。既有TCGA是batch-normalized RSEM，CPTAC/RC18/RC20是TPM；单细胞或Visium计数的测量分布也不能直接替代bulk输入。
代码实际n_out_of_range≥3将预测置NA，注释描述>6，差异如实记录。训练模型元数据下载不全、参考分布未下载；未声称模型完整可运行，未用平均基因分数代替模型。
RCC论文的病理法要求人工IT/PT区域和淋巴细胞图像定量。其补充PDF约51.7MB，两次限时下载均不完整，PDF结构检查失败，已标为partial，未提取患者标签。本轮分类预测0，终态INPUT_NOT_QUALIFIED，不能将“未分类”写成无排斥。
重新开启需取得相容bulk原始计数并完成参考校准，或完整公开病理标签/有效图像分区。这是本轮输入缺口，不要求客户先做新实验。

## 原始ECM参考提供方向对应，但不能证明纤维几何
已从Tomko2018原文表1完整提取27个基因、原始方向和P值。25个上调、LAMB2/OGN下调。原文这27个是normal与aligned-IDC的差异/趋势蛋白集合；与文中19正相关+5负相关的24个排列相关蛋白不同。DDR1引用27基因的确切计分方向未完全重建，因此本轮只做逐基因对应，不宣称复现DDR1评分。
全部27基因在现有矩阵中有基因记录；21/27中位方向对应，21个在至少5/7患者中对应。去除7个已有定义/程序基因后，15/20对应。OGN、VTN、FBN2未通过接收表达资格，记录为低表达/不可稳定评价；TGFBI、ACAN、AGRN中位反向，其中AGRN为7/7反向，均完整保留。
原文FMOD、BGN、TGFBI的P≥0.05，不能称为原文名义显著。按原P<0.05另保留筛选敏感性：{s3['source_P_lt05_concordant']}/{s3['source_P_lt05_genes']}方向对应；再排除定义基因为{s3['source_P_lt05_nondefinition_concordant']}/{s3['source_P_lt05_nondefinition_genes']}。完整27集合仅用于忠实对照发表参考，不作为自造显著基因集。
TabulaTIME官方预训练元数据有6,436特征、45个宽泛细胞类别，包含Fibroblasts、Myofibroblasts、Pericytes、SMC、Mono/Macro，但没有Fibro_CTHRC1或Macro_SLPI状态标签。网页连接本轮超时；现有下载模型不具备目标状态转移资格。没有将壁细胞改为CAF，也没有改C4成员。

## 证据链和使用边界
保留此前已确立结果：巨噬胞内TYMP—SPP1患者级相关，巨噬TYMP—间质胶原成熟的13患者关联，独立4患者方向对应。肿瘤TYMP到巨噬SPP1、稳定T/CD8进入减少和完整屏障链仍无新增支持。本轮没有重新检验这些旧假设。
本轮新证据来自相同7患者接收状态数据；外部文献特征独立于项目提出，但不构成独立生物队列验证。NicheNet是先验匹配的模型排名，无患者处理效应P值；n7不作推断检验、不补P/FDR。状态及聚类来源仍有选择依赖，即使排除了直接定义基因也不能消除全部循环性。RNA对应不等于蛋白变化、胶原排列、屏障功能或因果中介。
下一轮值得优先提出的独立问题是：TYMP相关巨噬状态是否伴随TGFB1输出变化，以及该输出是否与独立接收端靶反应一致。不能将本轮TGFB1排名直接写成TYMP→TGFB1→屏障已成立。

## 文件与验收
计划：CURRENT_STATE_AND_SEQUENTIAL_PLAN.md。
第一步：step1_receiver_network/，含冻结协议、患者资格、全基因对照、全部配体排名、受体/靶连接、独立检查和COMPLETION.json。
第二步：step2_immune_phenotype/，含方法输入资格和关闭原因。
第三步：step3_reference_alignment/，含原文27基因、全部患者逐基因结果、原文P筛选敏感性、模型标签和COMPLETION.json。
资源出处、Git版本/哈希、异常与退出状态：resources/、logs/；最终验收：qc/FINAL_VALIDATION.json。
所有新记录仅服务器保存，旧结果、对象、边界和文章未修改；未创建新图。
"""
write(P/'delivery/Phase295_文献方法顺序尝试终态报告.md',report)
now=datetime.datetime.now().astimezone().isoformat()
entry=f"""## Phase295终态（2026-09-30）：三项文献方法已按序尝试
详见{P.name}/delivery/Phase295_文献方法顺序尝试终态报告.md。
1. NicheNet真实网络完成：7患者/1,479目标/797全部合格配体/190巨噬表达配体；TGFB1全排名1，SPP1巨噬排名162。模型探索性，不等于TYMP或SPP1因果。32项核查与独立Pearson通过；API终止143经官方MD5校验与输入验证恢复。
2. 92基因CD8作者代码12文件验证；现有标准化RSEM/TPM不满足raw-count/训练分布要求；RCC病理补充PDF未完整；关闭于输入资格，无分类结果、不自造评分。
3. Tomko原文27基因完整对应21/27，排除定义基因为15/20；原文P<0.05再排除定义为{s3['source_P_lt05_nondefinition_concordant']}/{s3['source_P_lt05_nondefinition_genes']}。27差异蛋白≠24排列相关蛋白。RNA方向对应不等于实测排列。TabulaTIME已下模型45宽类不含目标CTHRC1/SLPI状态，未做错误转移。
当前状态：本轮试验完成，功能屏障证据缺口保留；无后台。所有来源/C4/对象/边界/旧结果保持，服务器为唯一进度依据。以下为历史启动记录。

"""
names=['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']
snap=P/'root_record_backups/before_final';snap.mkdir()
old={n:(R/n).read_bytes() for n in names}
for n,b in old.items():(snap/n).write_bytes(b)
artifacts=[
('phase295_trial_report','report',P/'delivery/Phase295_文献方法顺序尝试终态报告.md','COMPLETE'),
('phase295_nichenet_all_ligands','table',P/'step1_receiver_network/LIGAND_ACTIVITIES_ALL.tsv','VALIDATED_DESCRIPTIVE'),
('phase295_nichenet_terminal','validation',P/'step1_receiver_network/COMPLETION.json','PASS'),
('phase295_cd8_eligibility','qualification',P/'step2_immune_phenotype/COMPLETION.json','INPUT_NOT_QUALIFIED'),
('phase295_reference_concordance','table',P/'step3_reference_alignment/REFERENCE_GENE_CONCORDANCE.tsv','VALIDATED_DESCRIPTIVE'),
('phase295_reference_terminal','validation',P/'step3_reference_alignment/COMPLETION.json','COMPLETE'),
('phase295_source_audit','provenance',P/'resources/SOURCE_AUDIT_FINAL.tsv','RECORDED')]
changes={}
for n in names:
 if n.endswith('.md'):changes[n]=entry.encode()+old[n]
 else:
  fields=old[n].decode('utf-8-sig').splitlines()[0].split('\t');out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
  for rid,typ,path,status in artifacts:
   assert path.is_file()
   w.writerow(dict(result_id=rid,stage='Phase295',artifact_type=typ,absolute_path=str(path),generated_at=now,validation_status=status,validation_evidence=str(P/'qc/FINAL_VALIDATION.json'),notes='SHA256='+sha(path)))
  changes[n]=old[n]+(b'' if old[n].endswith(b'\n') else b'\n')+out.getvalue().encode()
for n in names:assert (R/n).read_bytes()==old[n],'Concurrent record change '+n
for n in names:
 tmp=R/(n+'.phase295final.tmp')
 with tmp.open('xb') as f:f.write(changes[n])
 shutil.copymode(R/n,tmp);os.replace(tmp,R/n)
status=json.loads((P/'qc/STATUS.json').read_text());status.update(status='COMPLETE_WITH_QUALIFICATION_CLOSURES',completed_at=now,step1=s1['status'],step2=s2['status'],step3=s3['status'],background_jobs=0)
(P/'qc/STATUS.json').write_text(json.dumps(status,indent=2))
# Final integrity checks against immutable startup protocol and terminal exit markers.
checks={}
checks['plan_unchanged']=sha(P/'CURRENT_STATE_AND_SEQUENTIAL_PLAN.md')==status['plan_sha256']
checks['step1_R_exit0']=(P/'logs/step1_network.exit').read_text().strip()=='0'
checks['step1_validation_exit0']=(P/'logs/step1_validate.exit').read_text().strip()=='0'
checks['step3_exit0']=(P/'logs/step3_alignment.exit').read_text().strip()=='0'
checks['step2_no_prediction']=s2['predictions']==0
checks['27_genes_and7patients']=len(ref)==27 and ref.patient_n.eq(7).all()
checks['all_four_records_written']=all((R/n).read_bytes()==changes[n] for n in names)
checks['histories_preserved']=all(changes[n].endswith(old[n]) for n in names if n.endswith('.md'))
checks['manifest7_unique']=all((R/'RESULTS_MANIFEST.tsv').read_text().count(rid+'\t')==1 for rid,_,_,_ in artifacts)
checks['partial_explicit']=(P/'resources/RCC2021_supplement.pdf.partial').exists() and not (P/'resources/RCC2021_supplement.pdf').exists()
assert all(checks.values()),checks
write(P/'qc/FINAL_VALIDATION.json',json.dumps(dict(status='PASS',at=now,checks=checks,checks_passed=len(checks),artifacts=[dict(path=str(p),sha256=sha(p)) for _,_,p,_ in artifacts],root_hashes={n:sha(R/n) for n in names}),indent=2))
print(json.dumps(dict(status=status['status'],final_checks=len(checks),s3=s3,report=str(P/'delivery/Phase295_文献方法顺序尝试终态报告.md')),ensure_ascii=False,indent=2))
