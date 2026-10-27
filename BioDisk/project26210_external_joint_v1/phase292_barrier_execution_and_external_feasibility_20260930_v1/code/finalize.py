from pathlib import Path
import json,csv,io,os,shutil,hashlib,datetime,re,subprocess
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
for name in ['analyze','stats_v2','validate_v4','geo_sample_metadata']:
 assert (P/'logs'/f'{name}.exit').read_text().strip()=='0'
v=json.loads((P/'qc/INDEPENDENT_VALIDATION_v2.json').read_text());assert v['status']=='PASS' and v['checks']==364
d=pd.read_csv(P/'tables/PATIENT_INPUTS.tsv',sep='\t');a=pd.read_csv(P/'tables/ASSOCIATIONS_v2.tsv',sep='\t')
leave=pd.read_csv(P/'tables/LEAVE_OUT_DESCRIPTIVE_v2.tsv',sep='\t')
lg=pd.read_csv(P/'tables/LEAVE_GENE_DESCRIPTIVE_v2.tsv',sep='\t')
w=leave[(leave.source=='macro_TYMP')&(leave.endpoint=='collagen_maturation')&(leave.model=='study_plus_perivascular_fraction')&(leave.omission_type=='patient')]
g=lg[(lg.source=='macro_TYMP')&(lg.endpoint=='collagen_maturation')]
samples=json.loads((P/'resources/GSE204809_SAMPLE_DESIGN.json').read_text());assert len(samples)==14
design=[]
for x in samples:
 rec=dict(sample=x['sample'],title=';'.join(x['title']))
 for s in x['characteristics_ch1']:
  if ': ' in s:
   k,val=s.split(': ',1);rec[k]=val
 mm=re.search(r'RNA(\d+)SA',rec['title'])
 rec['library_batch_label']=mm.group(1) if mm else 'UNSPECIFIED'
 rec['donor_pairing']='not_declared_in_sample_record'
 design.append(rec)
des=pd.DataFrame(design);des.to_csv(P/'tables/GSE204809_DESIGN_AUDIT.tsv',sep='\t',index=False)
hu=des[des['cell line']=='HUVEC'];assert len(hu)==10
cross=pd.crosstab(hu.library_batch_label,hu.treatment)
cross.to_csv(P/'tables/GSE204809_BATCH_TREATMENT.tsv',sep='\t')
zen=json.loads((P/'resources/PUBLIC_zenodo_imc.txt').read_text())
zenstatus=zen.get('access',{});dump('resources/ZENODO_ACCESS_EXTRACT.json',zenstatus)
external=[
dict(candidate='GSE171542',kind='TYMP protein stimulation',status='ALREADY_ANALYZED_PHASE280',evidence='3 groups x3 samples; donor pairing unresolved; prior SPP1/C4 decreases retained',next_action='No duplicate analysis or claim of new validation',url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171542'),
dict(candidate='GSE204809',kind='Tipiracil endothelial perturbation',status='PUBLIC_METADATA_CONFIRMED_DESIGN_LIMITED',evidence='14 samples; HUVEC untreated-TPI arms carry RNA201106SA vs RNA210902SA batch labels; donor IDs unresolved; three XLSX count links',next_action='If used, inspect processed matrix and restrict to context/response contrast; do not label clean TYMP inhibition or tumor-to-CAF experiment',url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204809'),
dict(candidate='HCC TP overexpression 2019',kind='direct tumor perturbation literature',status='NO_REPLICATED_PUBLIC_MATRIX_CONFIRMED',evidence='Methods calls samples technical duplicates; no usable accession confirmed in inspected article; no inference from technical replication',next_action='Mechanistic context only pending independent biological replicates and actual matrix',url='https://www.nature.com/articles/s41419-018-1282-6'),
dict(candidate='10x Xenium RCC RNA+protein updated 2026-08-19',kind='single-cell spatial multiomics',status='CASE_ONLY_PANEL_PENDING',evidence='Official current page: 1 donor, 690322 detected cells; not older 465534-cell record; server metadata request HTTP429; TYMP/SPP1/panel coverage not yet confirmed',next_action='Prioritize small panel/file-list audit at next access opportunity; no inference at patient level from one donor; dedup older versions',url='https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-ccrcc'),
dict(candidate='Meyer RCC IMC myeloid Zenodo18701672',kind='multipatient spatial protein',status='FILES_RESTRICTED',evidence='Zenodo page and API confirm restricted file access; public GitHub says data will open upon publication; three panels on consecutive sections',next_action='Exclude from immediately executable data plan; do not treat cross-panel cells as same-section contacts',url='https://zenodo.org/records/18701672'),
dict(candidate='Developmental reversion ccRCC Xenium 2026 preprint',kind='multipatient single-cell spatial',status='ACCESSION_AND_PANEL_UNRESOLVED',evidence='Primary indexed preprint reports 173 specimens/93 patients; full-text access failed; deposited matrices, panel and independence not confirmed',next_action='Metadata lead only; do not count as acquired independent cohort',url='https://www.biorxiv.org/content/10.64898/2026.08.05.743137v1.full'),
dict(candidate='APOE-POSTN ccRCC paper spatial cohorts',kind='existing Visium cohorts',status='DUPLICATE_EXISTING_DATA',evidence='GSE175540 and GSE210041 already used; GLMF equals GSE210041',next_action='Reference only, no new independent replication',url='https://pmc.ncbi.nlm.nih.gov/articles/PMC13450211/')]
pd.DataFrame(external).to_csv(P/'tables/EXTERNAL_FEASIBILITY.tsv',sep='\t',index=False)
queries=['site.ncbi.nlm.nih.gov/geo TYMP knockdown','TYMP GSE silencing','thymidine phosphorylase RNA-seq knockdown','Thymidine phosphorylase facilitates retinoic GSE','renal clear cell carcinoma multiplex imaging Xenium CODEX public dataset','A spatial multi-omic portrait of survival outcome renal','Developmental reversion underlies resistance data availability']
dump('resources/SEARCH_LOG.json',dict(date='2026-09-30',queries=queries,scope='bounded primary-source feasibility review; not exhaustive absence proof',no_author_contact=True,no_large_dataset_download=True))
nom=a[(a.model=='study_only')&(a.source=='macro_TYMP')&(a.endpoint=='collagen_maturation')].iloc[0]
dep=a[(a.model.str.endswith('receiver_depth'))&(a.source=='macro_TYMP')&(a.endpoint=='collagen_maturation')].iloc[0]
report=f"""# Phase292：巨噬TYMP与间质胶原成熟的关联获得补充支持，肿瘤TYMP直接屏障联系仍未建立

## 已完成的工作
2026-09-30。先冻结方案、三面板及输入哈希，再执行同文库患者分析；实际计数来自Phase287，身份沿Phase289，来源伪bulk沿Phase157。核查公共扰动与高分辨率空间来源。原H5AD、C4、边界、旧结果及Phase273正文未改，无新图。
细胞总输入19073，净化间质池15438；16名同文库候选按原恶性>=50/巨噬>=100/间质>=50门槛纳入13名。患者是统计单位，不是细胞数。研究组成：GSE156632 2、GSE159115 1、GSE178481 1、GSE207493 9。
## 三个执行装置端点
冻结胶原成熟/交联9基因、黏着斑9基因、收缩装置8基因。移除TYMP/SPP1及定义性身份标记；MYL9因身份重叠排除。面板均为研究者预设转录装置代理，不是功能实测或已验证屏障signature。
胶原面板与Phase289 ECM组织共享LOX/LOXL2/PLOD2，因此本轮是同数据功能细化，非独立验证。完整基因、移除记录及覆盖见PANELS_FROZEN.json、PANEL_EXCLUSIONS.json和PANEL_COVERAGE.tsv。
## 巨噬TYMP—胶原成熟的联系通过家族内FDR
研究调整偏Spearman rho={nom.rho:.6f}，P={nom.P:.5f}，BH FDR={nom.FDR:.5f}（3个预设端点家族）。加入接收池每细胞UMI深度后rho={dep.rho:.6f}，P={dep.P:.5f}，FDR={dep.FDR:.5f}。
逐一移除患者后rho范围{w.rho.min():.6f}至{w.rho.max():.6f}，13/13正向；去除胶原面板任一基因后rho范围{g.rho.min():.6f}至{g.rho.max():.6f}。这些是同数据描述性稳定性检查，没有新增P或独立复制。
移除GSE207493后只剩4名患者且模型自由度不足，不宣称跨研究稳健复现。巨噬TYMP是本轮次级来源，不能自动替代客户关注的肿瘤TYMP。
## 肿瘤TYMP三个端点均未获得名义显著支持
研究调整：胶原成熟rho=-0.166250、P=0.76382、FDR=0.76382；黏着斑rho=0.233862、P=0.49165、FDR=0.737475；收缩装置rho=0.459789、P=0.14174、FDR=0.42522。加入深度后亦均P>=0.05。
巨噬SPP1—胶原成熟rho=0.666619、P=0.13567、FDR=0.361695，当前数据不支持该预设关联；不可把巨噬TYMP阳性拼成TYMP—SPP1—CAF串联或免疫进入因果。
## 接收细胞身份与常数协变量
13名合格患者的配对间质池原广义注释均为Pericyte/SMC，比例恒1。这是原注释口径，不等于否定Phase289的活化/过渡候选状态；也不能改称纯CAF。
组成常数项与截距冗余。v1机械满秩门槛把含此项模型列为不可估计，原表保留。v2按补充方案仅移除无变异常数项：组成调整模型等价研究调整；组成+深度等价研究+深度，不是独立验证，也不能声称已排除可变壁细胞/CAF混杂。非恒定共线性仍阻断。无阈值、人群、面板、效应方向变化。
## 外部证据可行性
1. GSE171542已由Phase280分析，保留外源TYMP刺激SPP1/C4反向结果，不重复宣称新验证。
2. [GSE204809](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204809)：官方14份样本和3个处理计数表入口已核对。HUVEC无TPI组使用RNA201106SA批次标签，含TPI组使用RNA210902SA；同一批内分别有对照/RIG-I对比。未声明供者配对，不能直接把跨批差异解释为独立TYMP抑制效应。只保留内皮情境候选，不证明肿瘤—巨噬—间质作用。[原始研究](https://www.nature.com/articles/s41419-023-05821-0)
3. [2019肝癌TP研究](https://www.nature.com/articles/s41419-018-1282-6)：方法为技术双份，尚未确认可用公开矩阵及独立重复；不能作为推断性公共扰动验证。
4. [10x更新版RCC Xenium](https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-ccrcc)：官方页面为1供者，690322检测细胞，2026-08-19版本；不能套用文献旧465534细胞版本。仅可考虑高分辨率案例，基因面板和版本去重待完成。本轮服务器请求HTTP429，未取得可核验的TYMP/SPP1面板。
5. [Meyer IMC](https://zenodo.org/records/18701672)：记录公开但文件Restricted，不能列为当前可执行外部验证。代码库公开不等于数据开放，且三面板为连续切片，不能直接拼成同片接触。
6. 93患者Xenium预印本仍为候选线索，入口/面板/独立性未核清，不登记为已取得队列。现有GSE175540/GSE210041重复来源继续排除。
## 验收及故障
独立计数重聚合、来源表达、患者资格、所有27行模型、P/FDR、等价参数化与输入SHA共364项通过；rho最大误差5.56e-16，P误差0。27行中组成模型与基础模型重复，只有18个不同模型组合，不按重复行计支持次数。
原独立验证器float32重算与分析float64精度不一致；初次诊断未完整定位，后确认并改用float64。另补齐原冻结的检测资格规则。失败脚本/日志及修正说明保留，分析数值、验收容差未改变。
分析、v2统计、v4终验、GEO样本核查均exit0。内存门禁PASS；SwapFree仅76kB已记录。环境batch_bench，版本见resources/SOFTWARE.json。没有训练或依赖修改。
## 后续计划
优先1：完善可直接使用的独立ccRCC成纤维/间质输入资格，区分纯CAF、壁细胞及过渡状态，先锁定身份与患者匹配再检验；不得为延续阳性重标现有池。
优先2：继续高分辨率数据的小型面板/文件目录核查，满足TYMP/SPP1/巨噬/间质/CD8标记与细胞定位后才冻结案例分析；多患者独立数据优先。
优先3：公共TYMP遗传/药理扰动继续要求处理与批次可区分、重复可核查。GSE204809仅可作为内皮情境补充候选；不为构造阳性直接合并两个批次。
当前完成Phase292首轮功能分析及可行性筛查；后续外部功能/空间验证尚未完成，无后台任务。本轮不启动大文件下载，不等待客户新增实验，不改既有文章。
"""
(P/'delivery/Phase292_计划执行与首轮结果.md').write_text(report)
(P/'delivery/README_REPRODUCE.md').write_text("""# Reproduction
Authority: /mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase292_barrier_execution_and_external_feasibility_20260930_v1
Existing outputs must be preserved. Review PLAN_AND_PROTOCOL.md and MODEL_PARAMETERIZATION_ADDENDUM_v2.md.
Environment: /mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
Threads: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
All analysis stages require /mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh && python <script>.
Executed scripts in order: code/preflight.py, code/analyze.py, code/validate.py (failed), code/validate_v2.py (failed), code/validate_v3.py (passed original models), code/stats_v2.py, code/validate_v4.py (passed final models), code/public_metadata.py, code/geo_sample_metadata.py, code/finalize.py.
Final numerical table: tables/ASSOCIATIONS_v2.tsv; final validation: qc/INDEPENDENT_VALIDATION_v2.json.
Historical ASSOCIATIONS.tsv and validator failures are audit records, not deleted.
No new figures or manuscript edits. No local progress documents.
""")
for rel in ['PLAN_AND_PROTOCOL.md','MODEL_PARAMETERIZATION_ADDENDUM_v2.md','PANELS_FROZEN.json','tables/ASSOCIATIONS_v2.tsv','tables/PATIENT_INPUTS.tsv','tables/EXTERNAL_FEASIBILITY.tsv','tables/GSE204809_DESIGN_AUDIT.tsv','tables/GSE204809_BATCH_TREATMENT.tsv','tables/LEAVE_OUT_DESCRIPTIVE_v2.tsv','tables/LEAVE_GENE_DESCRIPTIVE_v2.tsv','qc/INDEPENDENT_VALIDATION_v2.json']:
 shutil.copy2(P/rel,P/'delivery'/Path(rel).name)
now=datetime.datetime.now().astimezone().isoformat()
dump('qc/COMPLETION.json',dict(status='COMPLETE_FIRST_EXECUTION_AND_FEASIBILITY',at=now,analysis_exit=0,validation_exit=0,checks=364,unique_patient_models=18,reported_model_rows=27,external_validation='not_yet_performed',no_running_analysis=True,no_local_progress_records=True))
files=[x for x in (P/'delivery').iterdir() if x.is_file()]
pd.DataFrame([dict(file=x.name,bytes=x.stat().st_size,sha256=sha(x)) for x in sorted(files)]).to_csv(P/'delivery/DELIVERY_MANIFEST.tsv',sep='\t',index=False)
entry=f"""## Phase292（2026-09-30）：屏障执行功能首轮与外部数据资格核查完成
用户授权制定计划然后推进。权威目录{P.name}，冻结方案PLAN_AND_PROTOCOL.md及三面板；当前正式表ASSOCIATIONS_v2.tsv、报告delivery/Phase292_计划执行与首轮结果.md。
16候选同文库患者按旧门槛纳入13（GSE207493占9）。主来源肿瘤TYMP→胶原成熟/黏着斑/收缩装置均P>=.05。次级巨噬TYMP→胶原成熟rho.809604/P.00895/FDR.02685；深度敏感rho.817427/P.01053/FDR.03159。留一患者13/13正（.736122–.884011），删除GSE207493剩4且df不足；不称跨研究复制。
接收池原广义注释全部Pericyte/SMC，不能称纯CAF，也不因此否定既有活化/过渡状态。组成比例常数1；v1设计冗余拒绝模型，v2仅删常数协变量以等价参数化，不新增假设。组成模型等价study-only，不作额外独立证据；原27行和所有失败记录保留，18个不同模型组合。原验证器float32精度问题已改float64，364/364通过，rho误差<5.6e-16/P误差0。
外部：GSE171542已做，不重复；新核对GSE204809官方14样本、3计数表入口，但TPI和非TPI与RNA201106SA/RNA210902SA批次混杂、供者配对未声明，只保留内皮情境候选。2019肝癌TP技术双份不作推断。10x更新版RCC Xenium为1供者/690322细胞，面板待核、服务器HTTP429；不能套旧465534版本。Meyer IMC Zenodo18701672文件Restricted；93患者Xenium预印本入口/面板未核清。无新独立外部验证，不自动启动大文件下载。
下一优先：独立合格成纤维/间质输入、独立高分辨率面板与病例资格；先来源和身份再锁端点，保持肿瘤TYMP主线，巨噬阳性仅补充。原C4、H5AD、空间边界、Phase273正文不变，无新图。分析/统计/终验退出0；无后台。记录仅服务器。

"""
# Preserve existing root records and atomically prepend current state.
updates={}
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
 path=R/name;old=path.read_bytes();backup=P/'root_record_backups'/name
 with backup.open('xb') as f:f.write(old)
 if name!='RESULTS_MANIFEST.tsv':
  updates[name]=entry.encode()+old
 else:
  header=old.decode('utf-8-sig').splitlines()[0].split('\t');assert header==['result_id','stage','artifact_type','absolute_path','generated_at','validation_status','validation_evidence','notes']
  buf=io.StringIO();writer=csv.DictWriter(buf,fieldnames=header,delimiter='\t',lineterminator='\n')
  for key,kind,rel,status,notes in [
   ('phase292_plan','protocol','PLAN_AND_PROTOCOL.md','FROZEN','functional refinement; three endpoints; no new spatial thresholds'),
   ('phase292_models_v2','table','tables/ASSOCIATIONS_v2.tsv','PASS','13 patients;27 rows/18 distinct models;constant composition redundant;macro TYMP-collagen FDR support;malignant TYMP unsupported'),
   ('phase292_independent_validation','qc','qc/INDEPENDENT_VALIDATION_v2.json','PASS','364/364;counts and P/FDR replay;inputs unchanged'),
   ('phase292_external_feasibility','audit','tables/EXTERNAL_FEASIBILITY.tsv','REVIEWED_WITH_LIMITS','7 candidates;GSE204809 batch caveat;Xenium panel unresolved;IMC restricted;no new external biological validation'),
   ('phase292_report','report','delivery/Phase292_计划执行与首轮结果.md','PASS','first execution complete;future external work pending')]:
   file=P/rel
   writer.writerow(dict(result_id=key,stage='Phase292',artifact_type=kind,absolute_path=str(file),generated_at=now,validation_status=status,validation_evidence=str(P/'qc/INDEPENDENT_VALIDATION_v2.json'),notes=notes+';SHA256='+sha(file)))
  updates[name]=old+(b'' if old.endswith(b'\n') else b'\n')+buf.getvalue().encode()
for name,new in updates.items():
 target=R/name;backup=P/'root_record_backups'/name
 assert target.read_bytes()==backup.read_bytes()
 tmp=R/(name+'.phase292.tmp')
 with tmp.open('xb') as f:f.write(new)
 shutil.copymode(target,tmp);os.replace(tmp,target)
dump('qc/ROOT_RECORD_UPDATES.json',[dict(path=str(R/n),sha256=sha(R/n),backup_sha256=sha(P/'root_record_backups'/n)) for n in updates])
print(json.dumps(dict(status='FINALIZED',delivery_files=len(files)+1,checks=364,leave_gene_rho=[float(g.rho.min()),float(g.rho.max())],GSE204809_groups=des.groupby(['cell line','treatment','library_batch_label']).size().to_dict()),default=str))
