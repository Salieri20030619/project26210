from pathlib import Path
import json,datetime,hashlib,subprocess,shutil,os,fcntl,csv
import numpy as np,pandas as pd
from scipy import sparse,stats
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra
ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
A=ROOT/'phase245_spatial_validation_sources_20260924_v1';B=ROOT/'phase246_spatial_marker_structure_20260924_v1';C=ROOT/'phase247_external_spatial_cases_20260924_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2))
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def gate(stage):
 r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with open(C/'logs/final_memory.log','a') as h:h.write(stage+'\n'+r.stdout+r.stderr)
 assert r.returncode==0
gate('external_independent_verification')
assert (C/'EXIT_STATUS').read_text().strip()=='0'
qc=pd.read_csv(C/'tables/EXTERNAL_INPUT_QC.tsv',sep='\t');effects=pd.read_csv(C/'tables/EXTERNAL_CASE_EFFECTS.tsv',sep='\t')
checks=[]
for r in qc.itertuples():
 gate('validate_'+r.sample)
 s=pd.read_csv(C/'tables'/f'{r.sample}_SPOT_MARKERS.tsv.gz',sep='\t');assert len(s)==r.QC_spots and s.barcode.is_unique
 uv=s[['x_um','y_um']].to_numpy();pairs=cKDTree(uv).query_pairs(100.01,output_type='ndarray')
 G=sparse.coo_matrix((np.full(2*len(pairs),100.),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(s),len(s))).tocsr()
 rr=[];cc=[]
 for start in range(0,len(s),256):
  dst=dijkstra(G,directed=False,indices=np.arange(start,min(len(s),start+256)),limit=400)
  a,b=np.where((dst>0)&np.isfinite(dst));rr.extend(a+start);cc.extend(b)
 H=sparse.coo_matrix((np.ones(len(rr)),(rr,cc)),shape=G.shape).tocsr();assert H.diagonal().sum()==0
 assert all(0<=i<len(s) for i in rr)
 # Full source eligibility is independently reconstructed from raw counts only for macrophage markers.
 import gzip
 from scipy.io import mmread
 prefix=r.accession+'_'+r.sample+'_ST.'
 features=pd.read_csv(C/'extracted'/(prefix+'features.tsv.gz'),sep='\t',header=None)
 barcodes=pd.read_csv(C/'extracted'/(prefix+'barcodes.tsv.gz'),sep='\t',header=None)[0]
 m=mmread(C/'extracted'/(prefix+'matrix.mtx.gz')).tocsr()
 ids=pd.Index(barcodes).get_indexer(s.barcode);assert (ids>=0).all()
 gi=np.flatnonzero(features[1].isin(['CD68','CSF1R','FCER1G','AIF1']))
 source=np.asarray((m[gi,:][:,ids]>0).sum(0)).ravel()>=2
 n=np.asarray(H.sum(1)).ravel();W=sparse.diags(1/np.maximum(n,1))@H;source &= n>=3
 assert source.sum()==r.macro_marker_source_spots
 for target in ['CD8lineage','Endothelial','Fibroblast']:
  row=effects[(effects['sample']==r.sample)&(effects.endpoint=='SPP1_RNA__'+target+'_400um')&(effects.model=='A0')].iloc[0]
  if pd.notna(row.effect):
   rrho=stats.spearmanr(s.SPP1_RNA.to_numpy()[source],(W@s[target].to_numpy())[source]).statistic
   assert abs(rrho-row.effect)<1e-12
   checks.append(dict(sample=r.sample,target=target,rho_saved=row.effect,rho_independent=rrho,status='PASS'))
pd.DataFrame(checks).to_csv(C/'tables/INDEPENDENT_EXTERNAL_RHO_CHECKS.tsv',sep='\t',index=False)
assert len(effects)==27 and not any(k in effects.columns for k in ['P','p','FDR'])
dump(C/'INDEPENDENT_VALIDATION.json',dict(status='PASS',at=now(),source_markers_recomputed_cases=3,neighborhood_rho_checks=len(checks),spots_checked=int(qc.QC_spots.sum()),no_inferential_P=True))
note='新增三病例未重复出稳定CD8分离。sc06T未通过CD8谱系标记资格；sc24T/sc27T未调整400um相关rho分别0.02445/0.00294，调整后−0.2252/0.04757；到CD8表达富集区的距离差分别−100/0um，没有一致远离方向。成纤维表达邻域未调整在三例均负（−0.1404、−0.6285、−0.3314），调整后−0.0675、−0.2982、0.0013，不支持将其写为普遍SPP1—CAF邻近增强。三例仅描述，不计算P/FDR。'
f=C/'REPORT_CN.md';old=f.read_text();shutil.copy2(f,C/'REPORT_CN.before_final_interpretation.md');f.write_text(old.replace('\n## Methods','\n\n'+note+'\n\n## Methods'))
# Correct paper-to-dataset attribution using the actually retrieved XML; preserves correction history.
f=A/'tables/EXTERNAL_VALIDATION_QUALIFICATION.tsv';d=pd.read_csv(f,sep='\t')
d.loc[d.source=='PMC11687132','reason']='retrieved paper accession GSE175540 is already in the current spatial cohort; selected published mIF figures do not provide our patient-level raw segmented-cell validation'
d.to_csv(f,sep='\t',index=False)
f=A/'REPORT_CN.md';old=f.read_text();shutil.copy2(f,A/'REPORT_CN.before_accession_correction.md')
f.write_text(old.replace('旧SPP1-CD44论文GSE210041与GLMF存在重复','旧SPP1-CD44论文所列GSE175540已在当前空间队列中'))
f=A/'code/review_sources.py';code=f.read_text();shutil.copy2(f,A/'code/review_sources.before_accession_correction.py.txt')
code=code.replace("'; GSE210041 already contains GLMF cases;", "'; this accession is already in the current spatial cohort;")
code=code.replace("GSE210041 already contains GLMF cases","GSE175540 already exists in the current spatial cohort")
code=code.replace('旧SPP1-CD44论文GSE210041与GLMF存在重复','旧SPP1-CD44论文所列GSE175540已在当前空间队列中')
f.write_text(code)
dump(A/'ACCESSION_CORRECTION.json',{'at':now(),'paper':'PMC11687132','actual_retrieved_accession':'GSE175540','prior_report_prose':'erroneous GSE210041 attribution corrected; old prose preserved','qualification_unchanged':'already used spatial dataset, not independent new raw validation'})
# Final per-item completion map prevents technical completion being confused with biological proof.
items=[
 ['SPP1/macrophage and receiver identity','RNA_CONCORDANCE_COMPLETED','38 sections; independent marker sets concordant with mapping','independent protein staining/segmented cell identity not obtained'],
 ['CD8 spatial structure','EXPLORATORY_SIGNAL_WITH_FAILED_EXTERNAL_REPLICATION','GSE175540 distance16 patients +300um P0.03516 FDR0.05273; unadjusted neighborhood20 patients rho-0.2814 FDR0.03545','adjusted association not significant; TJ and new3 cases do not reproduce consistent separation; no cell-exclusion claim'],
 ['vascular/CAF/pathology boundary','COMPUTATIONAL_CHECK_COMPLETED_UNSUPPORTED_OR_UNVALIDATED','fixed400um lineage marker and graph distance tests completed, negatives preserved','no consistent positive CAF/vascular neighborhood; pathology-reviewed boundaries absent'],
 ['shared-component independence','NOMINAL_SUPPORT_AFTER_REMOVING_SPP1','TJ12:10negative/2positive at both grids; rho-0.4897/-0.3873;P0.03857 FDR0.05786','direct T abundance9/12 positive,P0.146; relative balance not absolute cell loss'],
 ['MSI anatomy registration','SOURCE_AUDIT_COMPLETED_MATERIAL_MISSING','all75 archive members and author source inspected','independent optical/author transform absent in accessible materials'],
 ['metabolite structural identity','SOURCE_AUDIT_COMPLETED_MATERIAL_MISSING','28 imzML CV inspected; mslevel0; precursor/fragment metadata0;19ibh caches','no MS/MS/standard match; spermine/spermidine remain mass candidates'],
 ['independent cohort','THREE_NEW_CASES_ANALYZED_NO_STABLE_CD8_REPLICATION','GSE294006:27files,3cases,11684QC spots; all directions exported','small descriptive cohort; source RNA proxy differs from old V3; no treatment effect'],
 ['C4 spatial mapping','PREVIOUS_QUALIFICATION_FAILURE_RETAINED','existing leave-study mixture recovery false positive16.67%,worst-study error25pp','no new qualified reference; no C4 identity upgrade or cell2location retraining']]
pd.DataFrame(items,columns=['validation_item','status','evidence','remaining_gap']).to_csv(C/'tables/ALL_VALIDATION_ITEMS_STATUS.tsv',sep='\t',index=False)
lines=['# 26210空间补充验证交付','',f'日期：{now()}','',
'已完成现有材料可执行的空间补充及独立三病例对照；技术完成不等于所有生物学验证通过。','',
'1. 原38片的谱系标记一致性、固定400um邻域与表达区域距离分析完成；以患者为单位分来源统计。','2. 去除T-cell balance参照中的SPP1成分后，TJ-RCC两尺度均10/12负向，rho−0.4897/−0.3873，P0.03857、FDR0.05786，为名义支持；直接T细胞丰度未下降。','3. GSE175540的CD8表达区域距离差中位+300um（16患者，12正3负1零），P0.03516/FDR0.05273；未调整400um邻域rho−0.2814/FDR0.03545，调整后P0.2632。', '4. '+note,
'5. 75个MSI归档成员、28份imzML和19个解析缓存完成源核查，未找到可确证峰结构或独立解剖配准的材料。','6. C4原混合恢复门槛仍未通过；不放宽门槛、不重训cell2location。独立染色/病理审核材料仍缺。','',
'完整资格与剩余缺口：tables/ALL_VALIDATION_ITEMS_STATUS.tsv。全部负向、零效应、不可估计和校正未通过结果留在报告及统计表。','',
f'来源审查：{A}/REPORT_CN.md',f'计算补充：{B}/REPORT_CN.md',f'外部病例：{C}/REPORT_CN.md','',
'旧RDS、原始数据、Phase243/244结果和已交付文稿未改；没有新图。公开三病例下载已完成，本轮专用SSH代理PID16072已关闭。','',
'环境：part4_advanced_py_v1（Python3.11.15/numpy1.26.4/pandas2.3.3/scipy1.17.1）。正式代码分别保存于各阶段code/，逐阶段/样本调用50GiB门禁。']
(C/'SESSION_SUMMARY_CN.md').write_text('\n'.join(lines))
# Final status and manifests.
for p in [A,B,C]:
 assert (p/'EXIT_STATUS').read_text().strip()=='0'
 assert json.loads((p/'FINAL_VALIDATION.json').read_text())['status'].startswith('PASS')
 artifacts=[]
 for f in sorted(p.rglob('*')):
  if f.is_file() and ('code' in f.parts or 'tables' in f.parts or f.parent==p) and f.name not in ['OUTPUT_MANIFEST.json']:
   artifacts.append(dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)))
 dump(p/'OUTPUT_MANIFEST.json',artifacts)
stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
entry=f"""## {now()} Phase245–247 空间验证补充完成，未解决项明确保留
用户授权补充空间身份、结构、去共享成分、MSI配准/峰身份及独立来源验证。完成38片谱系RNA一致性、固定400um及固定20%表达区域最短图距离；不把RNA程序当独立细胞身份。
主要新增：TJ12去掉T-balance参照SPP1后两尺度10/12负，rho-0.4897/-0.3873,P0.038574,FDR0.057861；直接T绝对量9/12正,P0.146，不能说T细胞绝对减少。GSE175540 CD8区域距离16患者12正3负1零、中位+300um,P0.035156,FDR0.052734；未调整邻域20患者rho-0.28137,FDR0.035454，完整调整P0.263176。CAF与普遍正向血管结构未获补强。
新增外部GSE294006下载111564800bytes/27文件验收通过，sc06/sc24/sc27三登记治疗后病例共11684QC spots，仅描述。CD8一例资格不足、另两例A0 rho0.02445/0.00294、距离-100/0um，没有稳定分离重复；所有27条效应保留，无P/FDR。
MSI源审查75成员：28imzML/28ibd/19ibh解析缓存；mslevel字段0且无前体/碎片CV；未取得独立MSI光学/作者变换或结构确证。SPP1独立染色/病理边界、C4合格映射仍缺；原Phase220资格失败保留，不重训cell2location。外部IMC原数据未公开，S-BIAD1661面板无所需SPP1/CD8；旧论文所列GSE175540已在队列，不能再算独立。
目录：{A}；{B}；{C}。计算和来源核查验收通过，原38片独立复核55项患者汇总/112项相关/28条旧相关；外部独立复核{len(checks)}项相关及原始源资格。全缺口见Phase247/tables/ALL_VALIDATION_ITEMS_STATUS.tsv。旧原始数据/RDS/文稿不改，无新图。本轮分析、下载均结束，专用代理PID16072已关闭。后续可基于明确方向组织文字；独立染色、病理、MS/MS/作者配准材料未取得前不标为验证完成。
"""
(C/'root_record_backups').mkdir(exist_ok=True)
with open(ROOT/'phase245247_record_update.lock','a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:shutil.copy2(ROOT/name,C/'root_record_backups'/f'{stamp}_{name}')
 for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
  f=ROOT/name;tmp=f.with_name(f.name+'.phase245247.tmp');tmp.write_text(entry+'\n\n'+f.read_text());os.replace(tmp,f)
 f=ROOT/'RESULTS_MANIFEST.tsv'
 with open(f) as h:header=h.readline().rstrip('\n').split('\t')
 rows=[]
 for stage,p,status,note in [('Phase245',A,'PASS_WITH_EVIDENCE_GAPS','MSI source and independent validation material review; no anatomy/structural confirmation'),('Phase246',B,'PASS_EXPLORATORY','38sections marker concordance, shared-component removal, CD8 expression-region distance; adjustment failure retained'),('Phase247',C,'PASS_DESCRIPTIVE_EXTERNAL_CASES','3new treated ccRCC cases, no stable CD8 separation; full negative and eligibility outcomes retained')]:
  v=dict(result_id=stage.lower()+'_spatial_validation_v1',stage=stage,artifact_type='report',absolute_path=str(p/'REPORT_CN.md'),generated_at=now(),validation_status=status,validation_evidence=str(p/'FINAL_VALIDATION.json'),notes=note)
  rows.append({h:v.get(h,'') for h in header})
 with open(f,'a') as h:csv.DictWriter(h,fieldnames=header,delimiter='\t').writerows(rows)
 fcntl.flock(lock,fcntl.LOCK_UN)
dump(C/'RECORD_UPDATE.json',dict(status='PASS',at=now(),files=['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv'],old_inputs_preserved=True,new_figures=0,remaining_material_gaps_explicit=True))
print('FINALIZED',len(checks),'external correlations')

