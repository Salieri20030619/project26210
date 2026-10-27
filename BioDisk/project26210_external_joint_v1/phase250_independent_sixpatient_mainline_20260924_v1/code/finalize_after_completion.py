from pathlib import Path
import os,sys,json,datetime,hashlib,subprocess,shutil,fcntl,csv,ctypes,select,signal
import numpy as np,pandas as pd,h5py
from scipy import sparse,stats
ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
A=ROOT/'phase248_independent_cases_mainline_20260924_v1';S=ROOT/'phase249_independent_spatial_sources_20260924_v1';P=ROOT/'phase250_independent_sixpatient_mainline_20260924_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2))
def gate(stage):
 r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(P/'logs/validation_memory.log','a') as h:h.write(stage+'\n'+r.stdout+r.stderr)
 assert r.returncode==0
# Native completion event, no polling or repeated status checks.
pid=int((S/'RECOVERY_PID').read_text());cmd=Path('/proc')/str(pid)/'cmdline'
if cmd.exists() and b'recover_and_analyze.py' in cmd.read_bytes():
 fd=ctypes.CDLL(None,use_errno=True).syscall(434,pid,0)
 if fd>=0:select.select([fd],[],[]);os.close(fd)
assert (S/'DOWNLOAD_EXIT_STATUS').read_text().strip()=='0','unresolved download failure'
assert (P/'EXIT_STATUS').read_text().strip()=='0','analysis failed'
assert json.loads((P/'FINAL_VALIDATION.json').read_text())['status'].startswith('PASS')
gate('independent_validation')
registry=pd.read_csv(P/'tables/CASE_IDENTITY_FROZEN.tsv',sep='\t')
plan=json.loads((P/'PLAN_BEFORE_ANALYSIS.json').read_text())
effects=pd.read_csv(P/'tables/ALL_SECTION_EFFECTS.tsv',sep='\t')
qc=pd.read_csv(P/'tables/INPUT_QC.tsv',sep='\t').set_index('sample')
score_checks=[];rho_checks=[]
for r in registry.itertuples():
 gate('case_'+r.sample)
 d=pd.read_csv(P/'inputs'/f'{r.sample}_MAINLINE_INPUTS.tsv.gz',sep='\t')
 assert len(d)==qc.loc[r.sample,'QC_spots'] and d.barcode.is_unique
 with h5py.File(S/'downloads'/f'{r.sample}_filtered_feature_bc_matrix.h5','r') as h:
  a=h['matrix'];m=sparse.csc_matrix((a['data'][:],a['indices'][:],a['indptr'][:]),shape=tuple(a['shape'][:])).tocsr().astype(float)
  genes=np.array([v.decode() for v in a['features/name'][:]]);bc=np.array([v.decode() for v in a['barcodes'][:]])
 ids=pd.Index(bc).get_indexer(d.barcode);assert (ids>=0).all();m=m[:,ids].T.tocsr();tot=np.asarray(m.sum(1)).ravel()
 for name,gs in plan['programs'].items():
  if not d[name].notna().all():continue
  used=[g for g in gs if (genes==g).any()]
  raw=np.column_stack([np.asarray(m[:,np.flatnonzero(genes==g)].sum(1)).ravel() for g in used])
  z=np.nan_to_num(stats.zscore(np.log1p(raw/tot[:,None]*10000),axis=0,ddof=0),nan=0).mean(1)
  error=float(np.max(np.abs(z-d[name])));assert error<1e-10
  score_checks.append(dict(sample=r.sample,program=name,max_difference=error))
 for row in effects[(effects['sample']==r.sample)&(effects.model=='A0')&effects.scale.isin(['whole_tissue','block4x4'])].itertuples():
  if pd.isna(row.effect):continue
  if row.scale=='whole_tissue':g=d[[row.exposure,row.target]].dropna()
  else:
   g=d.groupby('block4')[[row.exposure,row.target]].median();g=g[d.groupby('block4').size()>=10].dropna()
  rr=stats.spearmanr(g.iloc[:,0],g.iloc[:,1]).statistic;assert abs(rr-row.effect)<1e-12
  rho_checks.append(dict(sample=r.sample,exposure=row.exposure,target=row.target,scale=row.scale,rho=rr,status='PASS'))
pd.DataFrame(score_checks).to_csv(P/'tables/INDEPENDENT_PROGRAM_CHECKS.tsv',sep='\t',index=False)
pd.DataFrame(rho_checks).to_csv(P/'tables/INDEPENDENT_RHO_CHECKS.tsv',sep='\t',index=False)
# Patient aggregation independently rebuilt, keeping adjacent controls separate.
pat=pd.read_csv(P/'tables/PATIENT_EFFECTS.tsv',sep='\t');sm=pd.read_csv(P/'tables/PATIENT_DIRECTION_SUMMARY.tsv',sep='\t')
keys=['tissue','role','exposure','target','scale','model']
for row in pat.itertuples():
 mask=effects.patient.eq(row.patient)
 for key in keys:mask &= effects[key].eq(getattr(row,key))
 v=effects.loc[mask,'effect'].dropna().to_numpy()
 if len(v):
  z=float(np.tanh(np.mean(np.arctanh(np.clip(v,-.999999999,.999999999)))));assert abs(z-row.effect)<1e-12
 assert len(v)==row.n_sections
assert registry.patient.nunique()==6 and registry[registry.tissue=='tumor'].patient.nunique()==6
paired=pd.read_csv(P/'tables/PAIRED_TUMOR_ADJACENT_DESCRIPTIVE.tsv',sep='\t')
assert paired.patient.nunique()==2 and np.allclose(paired.tumor-paired.adjacent,paired.tumor_minus_adjacent)
assert not any(c in effects.columns for c in ['P','p','FDR'])
dump(P/'INDEPENDENT_VALIDATION.json',dict(status='PASS',at=now(),program_score_checks=len(score_checks),independent_rho_checks=len(rho_checks),patient_rows_checked=len(pat),paired_donors=2,n_inference_tests=0))
# Primary results across two independent external sources; no pooled patient test.
first=pd.read_csv(A/'tables/COHORT_DIRECTION_COMPARISON.tsv',sep='\t')
first=first[(first.dataset=='GSE294006')&(first.role=='mainline')&(first.scale=='block4x4')].copy()
second=sm[(sm.tissue=='tumor')&(sm.role=='mainline')&(sm.scale=='block4x4')].copy();second['dataset']='E-MTAB-12767'
cols=['dataset','exposure','target','n_patients','positive','negative','median_rho','min_rho','max_rho']
pd.concat([first[cols],second[cols]],ignore_index=True).to_csv(P/'tables/TWO_INDEPENDENT_COHORTS_MAINLINE.tsv',sep='\t',index=False)
# Quantify scale/adjustment direction differences without choosing a favored model.
external_whole=sm[(sm.tissue=='tumor')&(sm.role=='mainline')&(sm.scale=='whole_tissue')]
external_whole.to_csv(P/'tables/WHOLE_TISSUE_AND_ADJUSTMENT_ALL.tsv',sep='\t',index=False)
cd8=sm[(sm.tissue=='tumor')&(sm.role=='CD8_small_supplement')]
cd8.to_csv(P/'tables/CD8_SMALL_SUPPLEMENT.tsv',sep='\t',index=False)
lines=['# 独立病例对照优先：本轮交付','',f'完成时间：{now()}','',
'本轮主线为TYMP相关代谢转录背景与SPP1信号的独立空间病例对照。GSE294006三登记病例与E-MTAB-12767六登记患者分别报告；不合并成统一九患者推断，不把切片数当患者数。CD8仅作现成或固定小型补充。','',
'|独立来源|关系|有效患者|正/负|rho中位数|范围|','|---|---|---:|---:|---:|---|']
for table in [first,second]:
 for r in table.itertuples():lines.append(f'|{r.dataset}|{r.exposure}—{r.target}|{r.n_patients}|{r.positive}/{r.negative}|{r.median_rho:.4f}|{r.min_rho:.4f}至{r.max_rho:.4f}|')
lines+=['',
'以上为固定4x4主尺度；全组织、400um邻域与深度/谱系调整均完整保留。每个独立来源少于10患者，本轮只给方向、效应与病例一致性，不新增P/FDR。GSE294006两例TYMP块中位数变化不足，保留不可估计；E-MTAB六患者的多张肿瘤切片按患者合并，两张旁组织只作两患者配对对照。','',
'## 配对旁组织对照','配对信号采用mean log1p CP10k，未比较片内Z分数均值。完整两患者每个变量的肿瘤减旁组织效应见tables/PAIRED_TUMOR_ADJACENT_DESCRIPTIVE.tsv，不能推广为人群差异。','',
'## CD8少量补充']
for r in cd8.itertuples():lines.append(f'- 六患者来源，{r.model}：有效{r.n_patients}，{r.positive}正/{r.negative}负，rho中位数{r.median_rho:.4f}。')
lines+=['三病例CD8结果继续引用Phase247，不新增阈值或扩大CD8分析；不能凭程序或混合spot声称物理排斥。','',
'## 来源资格与边界','12患者Visium/CosMx研究所列Zenodo16833780返回404，文献参考短链返回403；本轮未取得该数据，未计入独立验证。六患者E-MTAB官方38文件全部通过大小与SHA256登记、H5结构和坐标核验。独立性依据作者个体字段和研究来源，未做跨队列基因型鉴定。','代谢程序是转录背景，不是实测MTA；混合spot的TYMP/SPP1或去SPP1程序不是C4/细胞特异因果链。全尺度反向、调整减弱及未达资格端点均保留。','',
'详细报告：Phase248/REPORT_CN.md；Phase250/REPORT_CN.md。脚本、病例登记、输入完整性、程序和逐病例/患者统计均在服务器版本目录保存。旧数据、RDS及已交付文稿未修改；无新图、无cell2location重训。']
(P/'SESSION_SUMMARY_CN.md').write_text('\n'.join(lines))
# Author-object backup is not necessary once every official source matrix is validated.
author_status='NOT_USED'
if (S/'AUTHOR_OBJECT_VALIDATION.json').exists():author_status='DOWNLOADED_OFFICIAL_MD5_VERIFIED_NOT_USED_IN_PRIMARY_ANALYSIS'
elif (S/'AUTHOR_OBJECT_PID').exists():
 apid=int((S/'AUTHOR_OBJECT_PID').read_text());base=Path('/proc')/str(apid)
 if base.exists() and b'download_author_object.py' in (base/'cmdline').read_bytes():
  child=[int(v) for v in (base/'task'/str(apid)/'children').read_text().split()]
  ids=[apid]+[v for v in child if b'ccRCC_spatial_EU_ImmuneUpload.rds' in (Path('/proc')/str(v)/'cmdline').read_bytes()]
  fds=[]
  for v in ids:
   fd=ctypes.CDLL(None,use_errno=True).syscall(434,v,0)
   if fd>=0:fds.append(fd)
  for v in ids:
   try:os.kill(v,signal.SIGTERM)
   except ProcessLookupError:pass
  for fd in fds:assert select.select([fd],[],[],10)[0];os.close(fd)
  author_status='OPTIONAL_BACKUP_STOPPED_SOURCE_MATRICES_COMPLETE_PARTIAL_PRESERVED'
dump(S/'AUTHOR_BACKUP_DISPOSITION.json',dict(at=now(),status=author_status,reason='complete official source matrices and coordinates used for all mainline analyses; no need to substitute immune-subset author object',partial_bytes=(S/'downloads/ccRCC_spatial_EU_ImmuneUpload.rds.part').stat().st_size if (S/'downloads/ccRCC_spatial_EU_ImmuneUpload.rds.part').exists() else 0))
(S/'REPORT_CN.md').write_text('# 独立空间来源资格与下载\n\nE-MTAB-12767官方SDRF确认六个个体、九张捕获阵列，七张肿瘤和两张旁组织。38个官方文件114321307字节全部完成并逐项大小/SHA256登记；读取与坐标验收见Phase250。\n\n论文12患者数据入口Zenodo16833780本次404、作者短链403，未取得数据、不计入验证。另一Zenodo7619249作者处理对象仅作备用，实际状态见AUTHOR_BACKUP_DISPOSITION.json，未用于本轮主分析。\n\n本轮下载SSL超时首轮未通过；单连接持续会话恢复后全部通过，原失败日志和断点记录保留。\n')
dump(S/'FINAL_VALIDATION.json',dict(status='PASS_SOURCE_DOWNLOAD_AND_ELIGIBILITY',at=now(),official_files=38,total_bytes=114321307,individuals=6,author_backup=author_status,inaccessible_mRCC_dataset_not_used=True))
dump(S/'STATUS.json',dict(status='COMPLETE_VERIFIED',at=now()));(S/'EXIT_STATUS').write_text('0\n');(S/'RUN_FINISHED').write_text(now())
# Output hash records and root progress updates.
for p in [A,S,P]:
 files=[]
 for f in sorted(p.rglob('*')):
  if f.is_file() and ('code'in f.parts or 'tables'in f.parts or f.parent==p) and f.name!='OUTPUT_MANIFEST.json':files.append(dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)))
 dump(p/'OUTPUT_MANIFEST.json',files)
stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
entry=f"""## {now()} Phase248–250 独立病例主线对照完成，CD8降为少量补充
用户最新优先级：先补独立病例对照，CD8有数据就少量补充。本轮不以CD8为扩展主导。
Phase248：GSE294006三病例固定TYMP/代谢背景/SPP1五端点、两功能补充，主4x4及全组织/400um/调整完整保存。代谢背景—SPP1 RNA3/3正，中位rho0.558824；去SPP1状态2正1负。两例TYMP区域变化不足，主尺度仅1例可评；不放宽资格。旧38片按同RNA定义做方向比较，不覆盖q05/C4结果。独立验证12程序/35相关通过。
Phase249：发现12患者研究列出的Zenodo16833780为404、作者短链403，未取得或计入数据。实际新增六患者E-MTAB-12767官方38文件114321307bytes全部验收；九捕获阵列，七肿瘤+两旁组织，患者按官方individual字段冻结。首轮SSL失败留档，持续会话恢复完成。
Phase250：六患者主线分析完成；HG1两肿瘤阵列按患者合并，旁组织两配对患者单列；LG2同阵列肿瘤/瘤栓不拆新患者。主尺度及所有尺度/调整、两功能补充、CD8固定两口径均保留。各来源<10患者，只给方向/效应，不作P/FDR、不合并成九患者推断。完整数值见Phase250/tables/TWO_INDEPENDENT_COHORTS_MAINLINE.tsv及SESSION_SUMMARY_CN.md。
三个阶段均技术验收通过，外部六患者独立复核{len(score_checks)}程序/{len(rho_checks)}相关及患者合并。作者RDS备用状态{author_status}；正式分析使用完整官方矩阵，未用免疫子集替代整片。原数据/RDS/已交付文稿不改，无新图/模型重训。独立病例是RNA空间对应，不称实测MTA、C4定位或因果链；CD8仅辅助。
"""
(P/'root_record_backups').mkdir(exist_ok=True)
with open(ROOT/'phase248250_record_update.lock','a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']:shutil.copy2(ROOT/name,P/'root_record_backups'/f'{stamp}_{name}')
 for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
  f=ROOT/name;tmp=f.with_name(f.name+'.phase248250.tmp');tmp.write_text(entry+'\n\n'+f.read_text());os.replace(tmp,f)
 f=ROOT/'AGENTS.md';f.write_text(f.read_text()+'\n\n## 2026-09-24 用户优先独立病例主线对照\n用户明确“优先补充独立病例对照，CD8有的话就补一点”。以TYMP—代谢转录背景—SPP1相关状态的病例对照为先，CD8仅作已有或顺带可得的小型补充，不作为扩大分析主导。保持独立患者口径、各来源分开、反向与资格不足结果完整；不为名义显著挑选样本或端点。\n')
 f=ROOT/'RESULTS_MANIFEST.tsv'
 with open(f) as h:headers=h.readline().rstrip('\n').split('\t')
 rows=[]
 for stage,p,note in [('Phase248',A,'3 external mainline RNA cases; sparse TYMP eligibility retained'),('Phase249',S,'6 individuals verified,38 official processed files; unavailable12patient source excluded'),('Phase250',P,'6 independent patient mainline cases and2 paired adjacent controls;CD8 small supplement; no inferential P')]:
  r=dict(result_id=stage.lower()+'_independent_cases_v1',stage=stage,artifact_type='report',absolute_path=str(p/'REPORT_CN.md'),generated_at=now(),validation_status='PASS_DESCRIPTIVE' if stage!='Phase249' else 'PASS',validation_evidence=str(p/'FINAL_VALIDATION.json'),notes=note);rows.append({h:r.get(h,'') for h in headers})
 with open(f,'a') as h:csv.DictWriter(h,fieldnames=headers,delimiter='\t').writerows(rows)
 fcntl.flock(lock,fcntl.LOCK_UN)
dump(P/'RECORD_UPDATE.json',dict(status='PASS',at=now(),root_records_updated=5,mainline_priority_recorded=True,local_proxy_tunnel_pid_pending_close=3196))
print('ALL_FINALIZED',len(score_checks),len(rho_checks));print(second[cols].to_string(index=False));print('PAIRED',paired.to_string(index=False))

