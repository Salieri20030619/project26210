from pathlib import Path
import json,csv,io,os,hashlib,shutil,datetime,subprocess
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def dump(n,o):(P/n).write_text(json.dumps(o,ensure_ascii=False,indent=2))
q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(P/'logs/finalize_memory.log').write_text(q.stdout+q.stderr);assert q.returncode==0
for name in ['extract','qualify','sourcepair','external_stats','raw_verify','public_inspect','remote_zip','10x_protein','perturbation','terminal_evidence_v3']:
 assert (P/'logs'/f'{name}.exit').read_text().strip()=='0',name
components=['RAW_EXTRACTION_VALIDATION.json','EXTERNAL_STATS_VALIDATED.json','PERTURBATION_VALIDATION.json','TERMINAL_EVIDENCE_VALIDATION.json']
validation=[json.loads((P/'qc'/f).read_text()) for f in components]
assert all(v['status']=='PASS' for v in validation)
nchecks=sum(v['checks'] for v in validation);assert nchecks==716
protein=json.loads((P/'inputs/10x_protein_panel.json').read_text())['payload']['panel']
pt=pd.DataFrame(protein['targets']);pt.to_csv(P/'tables/TENX_OLD_PROTEIN_TARGETS.tsv',sep='\t',index=False)
text=' '.join(pt.fillna('').astype(str).to_numpy().ravel()).lower()
assert all(g not in text for g in ['tymp','spp1','osteopontin','thymidine phosphorylase'])
dump('qc/PROTEIN_TARGET_QUALIFICATION.json',dict(status='PASS_TECHNICAL_NOT_ELIGIBLE_TARGET',declared_targets=protein['num_targets'],records=len(pt),TYMP_measured=False,SPP1_measured=False))
dump('qc/NEW_10X_ACCESS_LIMIT.json',dict(url='https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-ccrcc',server_HTTP=429,windows_HTTP=429,local_file_created=False,conclusion='newer_690322_cell_dataset_panel_not_obtained; do_not_infer_its_panel_from_old465534_dataset'))
pa=pd.read_csv(P/'tables/EXTERNAL_SOURCEPAIR_FREEZE.tsv',sep='\t')
base=pa[(pa.dataset=='GSE264586')&(pa.stratum=='Baseline')&(pa.receiver=='strict_fibro')&(pa.source=='macro_TYMP')&(pa.endpoint=='collagen_maturation')].copy()
tissue=pd.read_csv(P/'tables/GSE264586_SAMPLE_TISSUE_ALIGNMENT.tsv',sep='\t')
base=base.merge(tissue,left_on='libraries',right_on='library',validate='one_to_one')
base.to_csv(P/'tables/INDEPENDENT_COLLAGEN_DIRECTION_PATIENTS.tsv',sep='\t',index=False)
assert len(base)==4
a=pd.read_csv(P/'tables/EXTERNAL_ASSOCIATIONS.tsv',sep='\t')
row=a[(a.dataset=='GSE264586')&(a.stratum=='Baseline')&(a.receiver=='strict_fibro')&(a.source=='macro_TYMP')&(a.endpoint=='collagen_maturation')&(a.model=='unadjusted')].iloc[0]
assert row.n==4 and abs(row.rho-1)<1e-12 and pd.isna(row.P)
# Unchanged input scope and panel validation; source objects were independently hashed by R.
for rec in json.loads((P/'qc/EXTERNAL_INPUT_MANIFEST.json').read_text()):assert sha(rec['path'])==rec['sha256']
for name,h in json.loads((P/'qc/SOURCEPAIR_FREEZE.json').read_text())['files'].items():assert sha(P/name)==h
stages=[
dict(route='independent_stromal',status='COMPLETE_DESCRIPTIVE',evidence='2 cohorts/202924 cells audited;strict candidates1062+1348;max source-paired stratum n4;no inferential P',adds='independent small-sample macroTYMP-collagen directional correspondence',does_not_add='tumorTYMP causality;pure CAF identity;T-cell exclusion'),
dict(route='nicheverse_RCC_spatial',status='CLOSED_TARGET_UNMEASURED',evidence='7824x366 real public RCC core;official Git blob verified;TYMP/SPP1 absent',adds='definitive panel qualification',does_not_add='direct spatial target validation;do not impute targets'),
dict(route='10x_old_RNA_protein',status='CLOSED_TARGET_UNMEASURED',evidence='official ZIP directory and panels range-extracted;465534 cells/405 unique RNA genes;RNA/protein both lack TYMP/SPP1',adds='version-resolved target qualification without36GB download',does_not_add='TYMP/SPP1 spatial validation'),
dict(route='10x_new_690322',status='CLOSED_THIS_RUN_ACCESS_LIMIT',evidence='public page known;metadata download HTTP429 from server and Windows;panel unresolved',adds='explicit acquisition boundary',does_not_add='absence of genes cannot be inferred from old panel'),
dict(route='GSE204809',status='COMPLETE_DESCRIPTIVE_DESIGN_LIMIT',evidence='all3 XLSX/all14samples verified;date and TPI aliased;designrank4/5;source donor independence unconfirmed',adds='endothelial context and complete processed-count evidence',does_not_add='unconfounded directTYMP drug effect;SPP1 induction;macrophage-CAF barrier'),
dict(route='Meyer_IMC',status='CLOSED_ACCESS_RESTRICTED',evidence='Phase292 Zenodo official access restricted;no new access grant',adds='none',does_not_add='no executable independent cohort'),
dict(route='full93patient_Xenium',status='CLOSED_UNACQUIRED_TARGET_SCOPE',evidence='fullcohort not acquired;available real366panel core lacksboth targets;fullstudy access unresolved',adds='public exemplar assessed',does_not_add='93 new validated patients'),
dict(route='GSE202374_existing',status='CAPACITY_AUDIT_ONLY',evidence='GEO6samples radiotherapy vs control;2/4groups insufficient for n>=10;outside two-cohort frozen computational scope',adds='documented alternative low-n source',does_not_add='not analyzed or independent replication claimed')]
pd.DataFrame(stages).to_csv(P/'tables/ROUTE_TERMINAL_STATUS.tsv',sep='\t',index=False)
report="""# Phase293：独立小样本支持巨噬TYMP—胶原成熟方向，直接免疫屏障证据未增加
2026-09-30；Phase292后续执行终态。用户要求继续到完成，本轮完成可执行的独立输入分析、公共扰动描述、空间面板实证资格及统一验收；受限或不具测量资格的路线明确关闭，不以缺失数据写成生物学阴性。

## 1. 新增独立患者的方向支持
GSE264586的基线4名同文库合格患者中，巨噬TYMP与严格成纤维候选胶原成熟程序Spearman rho=1.0。这是小样本描述性排序一致，不计算P/FDR，不称校正显著。肿瘤样上皮—成纤维同文库仅2名，不足以评价预设关联。
4例分别为RCC18（肺）、RCC19（腹膜）、RCC20（淋巴结）、RCC24（肺），均有GEO ccRCC疾病字段，不能称原发肾肿瘤独立队列。19份外部文库与GEO标题一对一匹配，C2/EOT只作明确时间点别名。全部原始患者分数、细胞数和组织来源见INDEPENDENT_COLLAGEN_DIRECTION_PATIENTS.tsv及EXTERNAL_LIBRARY_INPUTS.tsv。
三个端点全部保留：基线严格候选巨噬TYMP—胶原成熟/黏着斑/收缩装置rho分别1.0/0.4/0；巨噬SPP1对应0.8/0.8/0.6，同样不作推断。不能将这些分段排序拼成TYMP→SPP1→CAF中介或T细胞迁移效应。
GSE264586后续时间点和GSE294109各合格治疗层仅1名，保留患者值，相关不可估计；不合并时间点/治疗层来扩充n。

## 2. 身份资格先于关联
完整读取GSE264586 72422细胞、GSE294109 130502细胞，共202924细胞。旧Fibroblast标签分别3135/5736；按预设核心/支持标记、竞争谱系及强壁细胞排除得到1062/1348个严格成纤维候选。严格候选是标记支持的本轮资格层，不是新确证CAF标签；原RDS与冻结标签没有改写。
三室同时达标资格在原样记录之外，关联前补登记“来源对”资格：巨噬—间质无需额外足量肿瘤室，肿瘤样—间质无需额外足量巨噬室。仍要求同文库来源>=100巨噬或>=50肿瘤样细胞，严格接收>=20/原池>=50。不能把不同来源对的患者拼成同病例链。
两个研究、时间点和治疗层分开，患者是统计单位。三个面板完全继承Phase292，未更改基因或优化分组。阴性、不可估计和广义原池敏感性全部保留。
Phase292的13患者巨噬TYMP—间质胶原成熟结果（rho0.8096/P0.00895/FDR0.02685；深度敏感FDR0.03159）原样保留。本轮4患者独立方向一致增加了有限外部支持；既有肿瘤TYMP三个功能端点不支持的结论未改变。

## 3. 高分辨率空间路线完成实证资格核查
### 作者公开RCC示例
从作者Nicheverse仓库取得真实rcc_tma_core21 H5AD，7824细胞×366基因，原整数计数和微米坐标通过检查，文件与官方Git blob SHA一致。面板没有TYMP或SPP1；有POSTN/DCN/CD68/CD163/CD8A等不能替代缺失靶点，因此不进行TYMP/SPP1补值、代理分群或直接屏障检验。
[作者数据说明](https://nicheverse.org/tutorials/platform_xenium)，[作者仓库](https://github.com/digvijayky/nicheverse)。

### 10x旧版RNA+蛋白
通过[BioImage Archive S-BIAD2470](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD2470)官方记录定位10x原始公开ZIP，文件总长36149509228字节；仅HTTP范围读取目录和experiment/gene_panel/protein_panel，未下载36GB影像。
实际版本为465534细胞、405个唯一RNA靶标（430条基因面板记录存在重复），不是较新690322细胞数据。RNA及蛋白面板均没有TYMP/SPP1/osteopontin/thymidine phosphorylase，因此该路线不具直接靶点验证资格。完整面板和范围请求证据已保存。
较新690322细胞页面在服务器和Windows均返回HTTP429。本轮未取得其面板，不依据旧版推断新版也缺靶标；作为本轮访问受限终态保留，未称已完成该数据的生物学验证。没有创建本地文件。
93患者预印本全队列没有取得；公开366面板示例完成核查，不能据一个芯片称93例已验证。Meyer IMC仍按官方Restricted记录排除可执行路线，未申请受控访问或联系作者。

## 4. 公共扰动完成计数表分析，但不支持独立TYMP抑制验证
GSE204809三份官方XLSX均37247个基因，合计14文库；ZIP检查、逐列GSM映射、独立XML重读全部计数通过。作者处理计数含小数，原值保留，不四舍五入伪装整数原始计数。
两套HUVEC无/有TPI文库分别带RNA201106SA/RNA210902SA日期标签，TPI与日期标签完全重合；设计矩阵[截距,RIG,TPI,日期,RIG×TPI]秩4/5，药物主效应不能与日期项分离。这证明设计别名，不证明实际批次效应的大小。供者/培养独立性仍未确认，无处理推断P/FDR。
同表RIG相对control的TYMP log2(CPM+1)均值差：HMVEC无TPI +2.6523、HUVEC无TPI +3.7799、HUVEC含TPI +4.8716。这是RIG响应，不等于TPI增加TYMP；酶抑制也不必降低其RNA。
SPP1在HMVEC无TPI及HUVEC无TPI表内全零，标记不可评估；含TPI六样本总处理计数仅12，RIG对比差-0.042884，不能据此宣称SPP1稳定上调/抑制或屏障参与。ICAM1/VCAM1及CXCL9/10/11等固定端点与三个转录装置程序全部保留；这些是内皮表达背景，不能改称成纤维功能。
全基因CPM/上四分位方向、三个程序、样本映射和跨日期RIG响应差均归档；响应差的差只作描述，不当独立TYMP处理效应。
[官方GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204809)，[原始研究](https://www.nature.com/articles/s41419-023-05821-0)。

## 5. 验收、资源及故障
R_4.5.2/Seurat5.4.0顺序读取两对象；Python batch_bench进行资格、描述和校验。每步骤50GiB门禁通过，运行时可用内存约371GiB，swap几乎用满已记录；无训练、无依赖升级。
原Seurat全部16639768个面板计数、总UMI、基因数、严格候选ID及输入哈希独立核验12项通过；患者资格/聚合/相关复算595项；XLSX输入与设计35项；全基因效应及面板/来源74项，共716项通过。全基因效应最大误差3.6e-15。另完成蛋白面板靶标资格检查。
两次末轮验证器失败分别是内部NPZ字符串读取开关和GEO时间点缩写，已按原数据修正；失败脚本/日志/退出码完整保留，阈值和效应未改变。网络部分初次70秒超时后小文件续传成功，官方Git大小/hash校验通过；新版10x429不是生物学阴性。
无新图、未改旧C4/H5AD/RDS/空间边界/Phase273正文；进度和结果仅服务器。

## 6. 最终证据定位
本轮增加：独立ccRCC小样本中巨噬TYMP与成纤维候选胶原成熟的方向对应。
本轮没有增加：肿瘤TYMP直接驱动SPP1—CAF屏障、独立TYMP处理效应、靶点定位到同片屏障、T细胞跨界进入减少。
因此正文可保留“TYMP相关髓系—间质基质组织联系”的探索性补充，不能升级为已确认免疫屏障机制。没有通过新阈值寻找阳性，也没有用无法测量当阴性。
本轮三条路线均已得到终态：能完成的描述计算/核验完成，不具目标测量或访问资格者实证关闭。未来若取得新的目标覆盖空间/扰动数据，是新的数据增量任务；不把它挂为本轮持续后台任务。
"""
(P/'delivery/Phase293_独立证据补强终态报告.md').write_text(report)
(P/'delivery/README_REPRODUCE.md').write_text("""# Phase293 reproduction
Authority: /mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase293_independent_barrier_followthrough_20260930_v1
Review PLAN.md and both frozen protocols/addendum before any rerun. Existing outputs read-only; reruns require a new version.
Gate: /mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh && <analysis command>.
R: /mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla code/extract_external.R ; code/verify_external_raw.R
Python: /mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
Sequence: inventory.py; acquire_public.py/resume_small_inputs.py; external_qualify.py; freeze_source_pairs.py; external_analyze_validate.py; inspect_public_inputs.py; inspect_10x_remote_zip.py/inspect_10x_protein.py; perturbation.py; terminal_evidence_checks_v3.py; finalize.py.
Version1/2 terminal check failures preserved. No source RDS/H5AD or inherited labels modified.
Final statistics: EXTERNAL_ASSOCIATIONS.tsv (all n<10, P/FDR NA); PERTURBATION_TARGETS_QUALIFIED.tsv (descriptive only).
Final route status: ROUTE_TERMINAL_STATUS.tsv. Restricted or unmeasured sources are not biological nulls.
No local progress records, new figures, or manuscript edits.
""")
selected=['PLAN.md','EXTERNAL_IDENTITY_AND_STATISTICS_PROTOCOL.md','SOURCE_PAIR_ELIGIBILITY_ADDENDUM.md','PERTURBATION_PROTOCOL.md']
selected+=['tables/'+f for f in ['EXTERNAL_ASSOCIATIONS.tsv','EXTERNAL_LIBRARY_INPUTS.tsv','EXTERNAL_SOURCEPAIR_FREEZE.tsv','IDENTITY_CAPACITY.tsv','INDEPENDENT_COLLAGEN_DIRECTION_PATIENTS.tsv','PERTURBATION_TARGETS_QUALIFIED.tsv','PERTURBATION_SAMPLE_MAP.tsv','PERTURBATION_PROGRAM_SCORES.tsv','NICHEVERSE_PANEL_QUALIFICATION.tsv','TENX_OLD_RNA_PANEL_QUALIFICATION.tsv','TENX_OLD_PROTEIN_TARGETS.tsv','ROUTE_TERMINAL_STATUS.tsv','GSE264586_SAMPLE_TISSUE_ALIGNMENT.tsv']]
selected+=['qc/'+f for f in components+['PROTEIN_TARGET_QUALIFICATION.json','NEW_10X_ACCESS_LIMIT.json']]
for name in selected:shutil.copy2(P/name,P/'delivery'/Path(name).name)
files=sorted(f for f in (P/'delivery').iterdir() if f.is_file())
with (P/'delivery/DELIVERY_MANIFEST.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['file','bytes','sha256'],delimiter='\t');w.writeheader()
 for file in files:w.writerow(dict(file=file.name,bytes=file.stat().st_size,sha256=sha(file)))
now=datetime.datetime.now().astimezone().isoformat()
dump('qc/COMPLETION.json',dict(status='COMPLETE_AVAILABLE_ANALYSES_AND_ROUTE_DISPOSITION',at=now,technical_checks=716,additional_protein_panel_check='PASS',independent_spatial_biological_validation='NOT_PERFORMED_TARGET_UNMEASURED_OR_ACCESS_LIMIT',independent_tumorTYMP_perturbation='NOT_ESTABLISHED',new_support='4_patient_descriptive_macroTYMP_collagen_direction',all_executable_analysis_exits=0,no_background_tasks=True,no_local_progress_records=True))
entry=f"""## Phase293（2026-09-30）：Phase292后续独立证据补强执行终态
用户要求继续直到完成。权威目录{P.name}。报告delivery/Phase293_独立证据补强终态报告.md；三路线已计算/验收或给出资格关闭状态，不挂后台。
现存GSE264586/GSE294109共202924细胞原计数独立核查；3135/5736旧Fibroblast中固定标记通过1062/1348严格候选，原注释/RDS未改。源对资格在关联前补登记，三室资格原样保留。同文库并分时间/治疗后最大n4；所有外部P/FDR空，不拼接样本量。
GSE264586基线严格候选4患者巨噬TYMP→胶原成熟rho1.0，黏着斑.4、收缩0；SPP1对应.8/.8/.6，均描述。4例RCC18肺/RCC19腹膜/RCC20淋巴结/RCC24肺，GEO19库逐一匹配；不能称原发肾肿瘤或人群显著复制。肿瘤样TYMP配对仅2例，不作相关推断。GSE294109各合格治疗层仅1例。Phase292次级阳性及肿瘤TYMP阴性完整保留。
空间：作者真实RCC core21 7824x366原计数/坐标/Git blob通过，但TYMP/SPP1均不在面板。通过S-BIAD2470官方镜像范围读取10x36GB ZIP目录与RNA/蛋白面板，确证旧465534细胞/405唯一RNA版本RNA和蛋白都缺TYMP/SPP1；未下载大影像、不插补。新版690322细胞页面服务器及Windows均429，面板未取得，不将旧缺失迁移新版。Meyer IMC Restricted、93患者全队列未取得单列限制。
GSE204809三表37247基因x14库实际分析/XML重读；TPI与文库日期标签重合，设计rank4/5，供者独立性不清，仅描述，无处理P/FDR。RIG对比TYMP+2.6523/+3.7799/+4.8716（HMVEC无药/HUVEC无药/HUVEC含药）；非药物效应。SPP1无药两表全零不可评估，含药总counts12、差-.042884，不支持强SPP1输出。原小数计数保留。
技术验收12+595+35+74=716项通过，另蛋白资格通过；全基因效应误差<3.6e-15。内部NPZ读取/GEO缩写两个验证器失败修复，旧日志/exit保留。所有正式分析exit0。记录仅服务器，无新图、不改Phase273稿/C4/旧对象。
最终：新增独立小样本髓系—间质胶原成熟方向支持；未建立肿瘤TYMP直接SPP1屏障、无混杂TYMP处理效应或T进入减少。后续若新数据到位属于新增任务，不继续同批阈值搜索。

"""
updates={}
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
 target=R/name;old=target.read_bytes()
 with (P/'root_record_backups'/name).open('xb') as f:f.write(old)
 if name!='RESULTS_MANIFEST.tsv':updates[name]=entry.encode()+old
 else:
  fields=old.decode('utf-8-sig').splitlines()[0].split('\t');out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
  for rid,kind,rel,status,note in [
   ('phase293_protocol','protocol','PLAN.md','FROZEN','three routes;new inputs and old source identity preserved'),
   ('phase293_external_stromal','table','tables/EXTERNAL_ASSOCIATIONS.tsv','PASS_DESCRIPTIVE','max n4;no P/FDR;4patient macroTYMP-collagen rho1;not causal'),
   ('phase293_spatial_qualification','audit','tables/ROUTE_TERMINAL_STATUS.tsv','COMPLETE_WITH_LIMITS','two measured panels lack TYMP/SPP1;new10x429;no independent biological spatial validation'),
   ('phase293_perturbation','table','tables/PERTURBATION_TARGETS_QUALIFIED.tsv','PASS_DESCRIPTIVE','TPI/date alias;no inferential P;two SPP1 contrasts allzero'),
   ('phase293_validation','qc','qc/COMPLETION.json','PASS','716 checks plus protein panel eligibility'),
   ('phase293_report','report','delivery/Phase293_独立证据补强终态报告.md','PASS','all executable branches finished;explicit unmeasured/access-limited routes')]:
   w.writerow(dict(result_id=rid,stage='Phase293',artifact_type=kind,absolute_path=str(P/rel),generated_at=now,validation_status=status,validation_evidence=str(P/'qc/COMPLETION.json'),notes=note+';SHA256='+sha(P/rel)))
  updates[name]=old+(b'' if old.endswith(b'\n') else b'\n')+out.getvalue().encode()
for name,new in updates.items():
 target=R/name;assert target.read_bytes()==(P/'root_record_backups'/name).read_bytes()
 tmp=R/(name+'.phase293.tmp')
 with tmp.open('xb') as f:f.write(new)
 shutil.copymode(target,tmp);os.replace(tmp,target)
dump('qc/ROOT_RECORD_UPDATES.json',[dict(path=str(R/n),sha256=sha(R/n),backup=str(P/'root_record_backups'/n)) for n in updates])
# Only claim completed after all artifacts and record writes are checked.
for row in csv.DictReader((P/'delivery/DELIVERY_MANIFEST.tsv').open(),delimiter='\t'):
 file=P/'delivery'/row['file'];assert file.stat().st_size==int(row['bytes']) and sha(file)==row['sha256']
(P/'RUN_COMPLETE').write_text(now+'\n')
print(json.dumps(dict(status='COMPLETE',checks=716,delivery_files=len(files)+1,root_records=len(updates),directory=str(P)),ensure_ascii=False))
