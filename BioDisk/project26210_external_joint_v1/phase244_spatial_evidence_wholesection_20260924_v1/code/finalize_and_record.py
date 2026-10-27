from pathlib import Path
import sys,json,datetime,hashlib,subprocess,shutil,csv,fcntl,os
import numpy as np,pandas as pd
from scipy.stats import spearmanr,false_discovery_control
ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
A=ROOT/'phase243_donor_verified_immune_20260924_v1';B=ROOT/'phase244_spatial_evidence_wholesection_20260924_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,allow_nan=False))
def write_atomic(f,text):
    tmp=f.with_name(f.name+'.phase243244.tmp');tmp.write_text(text);os.replace(tmp,f)
gate=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(B/'logs/final_validation_memory.log').write_text(gate.stdout+gate.stderr)
assert gate.returncode==0
for p in [A,B]:
    assert (p/'EXIT_STATUS').read_text().strip()=='0'
    assert (p/'RUN_FINISHED').is_file()
    assert json.loads((p/'STATUS.json').read_text())['status']=='COMPLETE_VERIFIED'
    assert json.loads((p/'FINAL_VALIDATION.json').read_text())['status']=='PASS'
# Independent validation uses scipy rank correlations and its BH implementation.
ia=pd.read_csv(A/'tables/IMMUNE_ASSOCIATIONS_ALL.tsv',sep='\t')
for spec,g in ia.groupby('specification'):
    d=pd.read_csv(A/'tables'/f'{spec}_INPUT.tsv',sep='\t')
    assert len(d)==44 and d.run.is_unique
    for r in g.itertuples():
        assert abs(spearmanr(d[r.exposure],d[r.endpoint]).statistic-r.rho)<1e-12
        assert r.ci_low<=r.ci_high and r.loo_min<=r.loo_max and 0<r.p<=1
    for fam,s in g.groupby('family'):
        assert len(s)==2
        assert np.max(np.abs(false_discovery_control(s.p.to_numpy())-s.FDR.to_numpy()))<1e-12
freeze=json.loads((A/'IDENTITY_FREEZE.json').read_text())
assert len(freeze['selected'])==44 and 'SRR32389767' not in freeze['selected']
pairs=pd.read_csv(A/'tables/IDENTITY_PAIRS_REPRODUCED.tsv',sep='\t')
selected=pairs[pairs.run_a.isin(freeze['selected'])&pairs.run_b.isin(freeze['selected'])]
assert len(selected)==946 and selected.classification.eq('DISTINCT_SUPPORTED').all()
d=pd.read_csv(B/'tables/WHOLE_SECTION_PATIENT_INPUT.tsv',sep='\t')
wb=pd.read_csv(B/'tables/WHOLE_SECTION_ASSOCIATIONS_ALL.tsv',sep='\t')
for spec,g in wb.groupby('specification'):
    assert len(g)==6
    suffix='mean_TIC' if spec=='main_pixelwise_mean' else 'pooled_TIC'
    for r in g.itertuples():
        xx=d[r.molecule+'_'+suffix];yy=d[r.endpoint]
        assert len(d)==r.n_patients==12 and abs(spearmanr(xx,yy).statistic-r.rho)<1e-12
    assert np.max(np.abs(false_discovery_control(g.p.to_numpy())-g.FDR.to_numpy()))<1e-12
mq=pd.read_csv(B/'tables/WHOLE_SECTION_MASS_QC.tsv',sep='\t')
assert len(mq)==24 and mq.features.nunique()==1 and mq.features.iloc[0]==767 and mq.candidate_present.all()
for rec in json.loads((B/'INPUT_MANIFEST_BEFORE.json').read_text()):
    f=Path(rec['path']);assert sha(f)==rec['sha256'] and f.stat().st_mtime_ns==rec['mtime_ns']
for rec in json.loads((A/'INPUT_MANIFEST.json').read_text()):
    assert sha(rec['path'])==rec['sha256']
# Make consequential alternative-normalization result explicit in report prose, alongside main results.
note='整片主分析中，精胺质量候选—SPP1 RNA为正向名义显著：12患者rho=0.6224，P=0.03468，六项家族FDR=0.2081；留一患者rho范围0.5364–0.7273均为正。总强度之比的归一化对照方向仍为正，但P=0.05605/FDR=0.3363。因此属于依赖分析定义的探索性候选线索，不能称校正后显著、跨归一化均显著或局部共定位。TYMP及巨噬内SPP1比例未同步支持；这条线索不组成TYMP—精胺—SPP1因果链。'
f=B/'REPORT_CN.md';old=f.read_text();assert note not in old
shutil.copy2(f,B/'REPORT_CN.before_final_interpretation.md')
f.write_text(old.replace('\n## Methods','\n\n'+note+'\n\n## Methods'))
immune_note='SPP1与细胞毒、功能障碍程序均为正向且通过各自家族FDR；五基因定义方向相同。替换重复供者代表后SPP1—细胞毒rho=0.6712、SPP1—功能障碍rho=0.5081，均保持校正后支持。结果显示组织免疫程序共变化，不能将细胞毒程序正向结果改写为细胞毒受抑制。'
f=A/'REPORT_CN.md';old=f.read_text();shutil.copy2(f,A/'REPORT_CN.before_final_interpretation.md')
f.write_text(old.replace('\n## Methods','\n\n'+immune_note+'\n\n## Methods'))
summary=['# 26210样本核实、免疫及空间接续结果','',f'完成时间：{now()}','',
'样本：60份文库1770对关系复算一致。1份低信息观察文库暂缓，1对同供者去重后44名合格供者；7对候选干预供者成立，但处理标签未确认。被暂缓文库未假称解决。','',
'|分析|独立n|效应rho|双侧P|FDR|定位|','|---|---:|---:|---:|---:|---|',
'|SPP1—细胞毒程序|44|0.670754|0.0000099999|0.0000199998|家族校正后支持，正向|',
'|SPP1—功能障碍程序|44|0.491896|0.000789992|0.000789992|家族校正后支持，正向|',
'|整片精胺候选—SPP1 RNA|12|0.622378|0.0346797|0.208078|事后探索，名义显著|','',
immune_note,'',note,'',
'空间原局部6主终点FDR=1的结果保留。Phase234转录背景—SPP1丰度校正后结果及TJ-RCC固定400um SPP1—ITGAV/ITGB1名义关联也保留；完整正负统计未删。12正谱8组有多个朝向、6组对掩模不一致，独立解剖确认0组；正式局部对应优先补MSI自身光学图/作者变换/独立形态标志。','',
'本轮没有新图或修改旧RDS、原始数据和已交付文稿。所有新任务已结束，EXIT_STATUS=0且FINAL_VALIDATION=PASS。新增6个整片端点和6个归一化对照均完整报告，不因达到或未达到名义显著继续扩展筛选。','',
f'免疫报告：{A}/REPORT_CN.md',f'空间报告：{B}/REPORT_CN.md','',
'实际环境：/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python（Python3.11.15、numpy1.26.4、pandas2.3.3、scipy1.17.1）。',
f'复现命令（新版本路径，原版结果只读）：python {A}/code/run_analysis.py；python {B}/code/run_analysis.py。脚本内部逐阶段/样本调用50GiB内存门禁；应复制到新版本目录再运行，保留现有结果。']
(B/'SESSION_SUMMARY_CN.md').write_text('\n'.join(summary))
v={'status':'PASS','at':now(),'immune_independent_rho_checks':8,'spatial_independent_rho_checks':12,'scipy_BH_families_checked':8,'selected_donor_pairs_distinct':946,'all_12_positive_spectra_feature_count':767,'original_input_hashes_unchanged':True,'reports_include_normalization_counterevidence':True,'new_figures':0}
dump(B/'INDEPENDENT_FINAL_VALIDATION.json',v)
for p in [A,B]:
    artifacts=[]
    for f in sorted(p.rglob('*')):
        if f.is_file() and ('tables' in f.parts or 'code' in f.parts or f.parent==p) and f.name not in ['OUTPUT_MANIFEST.json','RUNNER_PID','RECORD_UPDATE.json']:
            artifacts.append({'path':str(f),'bytes':f.stat().st_size,'sha256':sha(f)})
    dump(p/'OUTPUT_MANIFEST.json',artifacts)
stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
entry=f"""## {now()} Phase243/244 样本关系复核、免疫执行及空间整片探索完成
用户授权先核实样本关系、启动免疫、再改进空间；名义显著如实分级，不以P阈值挑样本、配准或模型。
Phase243：60文库1770对SNP复算一致。观察46文库中SRR32389767覆盖534位点低于原1000门槛暂缓；SRR32389752/55同供者去重，冻结44供者，946对均DISTINCT_SUPPORTED。原45对未决保留，不改原Phase240门槛。14候选支持7对供者，但处理标签/来源字段冲突未解，不做干预比较。
免疫：SPP1—细胞毒rho0.670754/P9.9999e-6/FDR1.99998e-5；SPP1—功能障碍rho0.491896/P0.000789992/FDR0.000789992。五基因定义与重复代表替换全部保留。双正向组织关联，不称细胞毒受抑制；原Phase213阴性保留。
Phase244：免配准整片12患者的精胺质量候选—SPP1 RNA rho0.622378/P0.0346797/FDR0.208078，为新事后探索名义信号；总强度归一化对照rho0.573427/P0.0560494/FDR0.336297。TYMP/巨噬内比例未同步支持，不能当局部共定位或完整因果链。完整六主端点和六归一化对照均保存；原Phase242局部主要阴性不变。旧Phase234空间1001项符号P独立复算一致，已有转录背景阳性与配受体名义线索保留。
新版本：{A}；{B}。两任务EXIT_STATUS0/FINAL_VALIDATION PASS；独立复核20项rho、8校正家族、原输入哈希通过。无新图，原始数据/旧RDS/既有文稿不改。所有本轮任务已结束；后续空间局部验证优先取得独立形态/作者配准依据，不按分子方向选朝向。
"""
with open(ROOT/'phase243244_record_update.lock','a') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX)
    for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:
        shutil.copy2(ROOT/name,B/'root_record_backups'/f'{stamp}_{name}')
    for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
        f=ROOT/name;write_atomic(f,entry+'\n\n'+f.read_text())
    f=ROOT/'RESULTS_MANIFEST.tsv'
    with open(f) as h:header=h.readline().rstrip('\n').split('\t')
    newrows=[]
    for rid,stage,path,note in [
        ('phase243_donor_verified_immune_v1','Phase243',A/'REPORT_CN.md','44 donor-qualified observations; four fixed immune tests plus duplicate sensitivity; unresolved low-information and treatment identities retained'),
        ('phase244_wholesection_exploratory_v1','Phase244',B/'REPORT_CN.md','12 patients; exploratory spermine-candidate versus SPP1 RNA nominal P0.03468/FDR0.2081; alternative normalization P0.05605 retained; not local colocalization'),
        ('phase243244_independent_validation','Phase243/244',B/'INDEPENDENT_FINAL_VALIDATION.json','20 rho checks, 8 BH families, 946 distinct retained donor pairs; inputs preserved')]:
        row={h:'' for h in header}
        vals={'result_id':rid,'stage':stage,'artifact_type':'validation' if path.suffix=='.json' else 'report','absolute_path':str(path),'generated_at':now(),'validation_status':'PASS','validation_evidence':str(B/'INDEPENDENT_FINAL_VALIDATION.json'),'notes':note}
        row.update({k:v for k,v in vals.items() if k in row});newrows.append(row)
    with open(f,'a') as h:csv.DictWriter(h,fieldnames=header,delimiter='\t').writerows(newrows)
    fcntl.flock(lock,fcntl.LOCK_UN)
dump(B/'RECORD_UPDATE.json',{'status':'PASS','at':now(),'updated':['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv'],'backup_directory':str(B/'root_record_backups'),'no_local_progress_records':True})
print(json.dumps(v));print('RECORDS_UPDATED')

