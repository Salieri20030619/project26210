from pathlib import Path
import json,datetime,hashlib,shutil
import pandas as pd,numpy as np
from scipy.stats import rankdata,t,spearmanr,false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase154_same_library_TYMP_C4_test_20260915_v1'
d=pd.read_csv(O/'SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t')
res=pd.read_csv(O/'TYMP_C4_AND_SPP1_RESULTS.tsv',sep='\t')
cells=pd.read_csv(O/'SELECTED_CELL_TWO_GENE_COUNTS.tsv.gz',sep='\t',dtype={'candidate':str})
frozen=set(pd.read_csv(R/'09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz',sep='\t').cell_id)
checks=[]
def ck(name,v):checks.append(dict(check=name,passed=bool(v)))
ck('16_unique_patients_libraries',len(d)==16 and d.patient.is_unique and d.sample_id.is_unique)
ck('raw_counts_nonnegative_integer',(cells[['TYMP_counts','SPP1_counts']].to_numpy()>=0).all() and np.equal(cells[['TYMP_counts','SPP1_counts']],np.floor(cells[['TYMP_counts','SPP1_counts']])).all().all())
ck('library_size_identity',np.allclose(cells.total_counts_raw,cells.recomputed_library_counts,rtol=0,atol=1e-5))
for _,row in d.iterrows():
 z=cells[cells.sample_id==row.sample_id];a=z[z.lineage=='Malignant'];b=z[z.lineage=='Macrophage']
 ck(row.sample_id+'_identity',set(z.curated_patient_id)=={row.patient} and len(a)==row.malignant_cells and len(b)==row.macrophage_cells)
 count=sum(b.cell_id.isin(frozen));ck(row.sample_id+'_frozen_C4',count==row.C4_cells and np.isclose(count/len(b),row.C4_fraction))
for i,row in res[res.adjustment=='raw'].iterrows():
 x=d[row.exposure].rank();y=d[row.target].rank();r=x.corr(y);pv=2*t.sf(abs(r)*np.sqrt((len(d)-2)/(1-r*r)),len(d)-2)
 ck('independent_spearman_'+str(i),np.isclose(r,row.rho) and np.isclose(pv,row.P))
for key,g in res.groupby(['exposure','adjustment']):
 ck('BH_'+str(key),np.allclose(false_discovery_control(g.P),g.FDR_three_endpoints))
# Reproduce study-preserving permutation results using explicit rank residual projection.
rng=np.random.default_rng(262154);targets=['C4_fraction','SPP1_pseudobulk','SPP1_detection_fraction'];C=pd.get_dummies(d.study,dtype=float).to_numpy();projection=C@np.linalg.pinv(C)
groups=[np.asarray(v) for v in d.groupby('study').groups.values()]
for exposure in ['TYMP_pseudobulk','TYMP_mean_log1pCP10K']:
 x=rankdata(d[exposure]);Y=np.column_stack([rankdata(d[k]) for k in targets]);x=x-projection@x;Y=Y-projection@Y;den=np.sqrt(np.sum(x*x)*np.sum(Y*Y,axis=0));rho=x@Y/den;hits=np.zeros(3,int)
 for _ in range(9999):
  xp=x.copy()
  for ix in groups:xp[ix]=rng.permutation(xp[ix])
  rp=xp@Y/den;hits+=(np.abs(rp)>=np.abs(rho)-1e-14)
 for j,target in enumerate(targets):
  row=res[(res.exposure==exposure)&(res.target==target)&(res.adjustment=='within_study_rank_permutation')].iloc[0]
  ck(exposure+target+'_study_permutation',np.isclose(rho[j],row.rho) and np.isclose((hits[j]+1)/10000,row.P))
# Planned technical read-depth diagnostic, no gate on significance.
C=np.column_stack([np.ones(len(d)),rankdata(np.log1p(d.malignant_total_counts)),rankdata(np.log1p(d.macrophage_total_counts))]);project=C@np.linalg.pinv(C);x=rankdata(d.TYMP_pseudobulk);xr=x-project@x;rows=[]
for target in targets:
 y=rankdata(d[target]);yr=y-project@y;r=np.corrcoef(xr,yr)[0,1];df=len(d)-np.linalg.matrix_rank(C)-1;p=2*t.sf(abs(r)*np.sqrt(df/(1-r*r)),df)
 rows.append(dict(target=target,n_patients=len(d),rho=r,P=p,adjustment='rank_control_malignant_and_macrophage_total_counts'))
tech=pd.DataFrame(rows);tech['FDR_three_endpoints']=false_discovery_control(tech.P);tech.to_csv(O/'TECHNICAL_COUNT_ADJUSTED_DIAGNOSTIC.tsv',sep='\t',index=False)
loo=pd.read_csv(O/'LEAVE_ONE_PATIENT_STUDY_DIRECTIONS.tsv',sep='\t');c4=loo[loo.target=='C4_fraction']
for i,row in c4.iterrows():
 col=row.omission_type;z=d[d[col]!=row.omitted];rr=z.TYMP_pseudobulk.rank().corr(z.C4_fraction.rank());ck('leave_one_'+str(i),len(z)==row.n_patients and np.isclose(rr,row.rho))
ck('no_figures',not any(p.suffix.lower() in ['.png','.svg','.pdf','.jpg'] for p in O.rglob('*')))
assert all(x['passed'] for x in checks),[x for x in checks if not x['passed']]
(O/'INDEPENDENT_VALIDATION.json').write_text(json.dumps(dict(status='PASS',checks=len(checks),passed=len(checks),details=checks),indent=2))
def table(df):
 def fmt(v):
  if isinstance(v,float):return 'NA' if np.isnan(v) else f'{v:.6g}'
  return str(v)
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
primary=res[(res.exposure=='TYMP_pseudobulk')&(res.target=='C4_fraction')&(res.adjustment=='raw')].iloc[0]
study=res[(res.exposure=='TYMP_pseudobulk')&(res.target=='C4_fraction')&(res.adjustment!='raw')].iloc[0]
lopo=c4[c4.omission_type=='patient'];loso=c4[c4.omission_type=='study']
by=pd.read_csv(O/'STUDY_STRATIFIED_RESULTS.tsv',sep='\t');by=by[by.target=='C4_fraction'][['study','n_patients','rho','P','FDR_three_endpoints']]
report="""# 同肿瘤文库TYMP—C4直接检验

## TYMP与C4比例呈正向趋势，留一患者及留一研究方向均保留

在16名同肿瘤文库配对患者中，恶性细胞TYMP伪汇总表达与冻结C4占成熟巨噬比例的Spearman相关为rho=0.467647，P=0.067761。可将其作为本组数据中的探索性正向关联趋势，不能称名义显著。FDR=0.203283保留在统计表，不作为结果可用或继续工作的硬门槛。

去掉任一患者的16次分析全部保持正向，rho范围0.371429–0.585714。去掉任一研究的4次分析也全部保持正向，rho范围0.328671–0.714286；其中去掉GSE207493后只剩6名患者，相关系数仅作小样本方向描述。该方向不由单个患者或单个研究的存在所决定，不能据此称每个独立研究都已显著验证。

这项结果直接回答“肿瘤TYMP是否与新C4对应”，不借用旧SPP1-state或总体SPP1表达的显著性。

## 研究内调整仍为正向，证据力度较有限

研究内秩调整后的rho=0.302925，研究内9,999次置换P=0.2851、FDR=0.8553。原始与调整后方向一致，研究构成影响估计幅度；研究内统计证据较有限。两个仅1名患者的研究在去除研究均值后的系数中不提供研究内变异信息，实际有信息的研究为GSE156632（4人）和GSE207493（10人）。

{by}

使用恶性细胞逐细胞log1pCP10K均值作为替代表达指标，方向仍正，rho=0.364706、P=0.164868；研究内调整后rho=0.049454、P=0.8671。说明总体正方向保留，但同研究内强度对TYMP汇总方式敏感。不能把本结果与此前17患者z_TYMP_RNA结果当作同一端点的简单数值更新。

技术读数调整仅作补充诊断：
{tech}

## SPP1表达与C4比例分开报告

同一16患者、同一恶性TYMP伪汇总指标下：
{allprimary}

总体巨噬SPP1伪汇总表达的rho=0.105882、P=0.696334；SPP1检出比例rho=0.329412、P=0.212801。三个端点方向各自记录，不能相互转移显著性。当前C4端点的正向趋势与总体SPP1表达关联较弱可以同时成立。

## 当前可用于正文的表述

“在同肿瘤样本的恶性细胞与巨噬细胞配对分析中，TYMP表达与SPP1富集C4巨噬状态的相对丰度呈正向趋势（n=16，Spearman rho=0.468，P=0.0678）。该方向在留一患者及留一研究分析中均保留，提示TYMP较高的肿瘤伴随较高C4比例的候选联系。”

统计附件同时列出FDR、研究内调整和替代TYMP指标，保持探索性定位。本结果不要求达到很强的效应阈值，也不因FDR未通过而删除明确的方向；现有观察性证据不确定TYMP是否直接促进C4形成。

## 方法、输入与复核

16个文库对应16名独立患者：GSE156632 4名、GSE159115 1名、GSE178481 1名、GSE207493 10名。资格在表达检验前冻结：同Tumor文库至少20个正式恶性细胞、50个成熟巨噬。C4=0的样本原规则允许保留；实际这16个合格样本均有C4，不是先筛选C4阳性再分析。

恶性细胞使用原eligible_for_subclone_selection_v1；成熟巨噬使用新分群的12个成熟候选，排除Monocyte和C1边界群。C4成员固定，分母为同文库全部成熟巨噬。从原H5AD layers/counts分块读取所需23,018个细胞的TYMP/SPP1原始counts，并重算每细胞文库总counts；与total_counts_raw逐一一致。

每文库分别对恶性细胞、巨噬细胞求原始counts和总counts，转log2(CPM+1)。以16名患者计算Spearman；三桥接端点为C4比例、SPP1伪汇总表达、SPP1检出比例，每暴露/调整方式分别BH。没有细胞级推断或强制显著性门槛。研究内相关以rank去除研究均值后计算，研究内置换9,999次，seed262154。技术诊断控制恶性和巨噬伪汇总总counts的秩，单独保存，不替代主结果。

没有重训聚类、Monocle3、SCENIC或空间模型，没有生成图件。独立复核{ncheck}项通过，涵盖患者/文库身份、冻结C4成员、原始counts、伪汇总、独立Spearman计算、研究内置换、BH和留一结果。完整患者指标及统计表均保存于本目录。

原17患者探索结果、旧SPP1-state分析和既有阴性保持原口径。下一步可围绕当前正向趋势开展有细胞来源约束的空间/独立队列证据核查；本报告不把未完成的外部复现写为已完成。
"""
report=report.replace('{by}',table(by)).replace('{tech}',table(tech)).replace('{allprimary}',table(res[(res.exposure=='TYMP_pseudobulk')&(res.adjustment=='raw')][['target','n_patients','rho','P','FDR_three_endpoints']])).replace('{ncheck}',str(len(checks)))
(O/'TYMP_C4_DIRECT_TEST_RESULTS_CN.md').write_text(report,encoding='utf-8')
shutil.copy2(R/'scripts/40_same_library_TYMP_C4_test_20260915_v1.py',O/'01_compute.py')
B=O/'state_before_update';B.mkdir(exist_ok=False)
for p in [P/'CODEX_HANDOFF.md',P/'PROJECT_STATUS.md',P/'DECISIONS.md',P/'RESULTS_MANIFEST.tsv',R/'STATUS.json']:shutil.copy2(p,B/p.name)
note="""## 2026-09-15 Phase154 同文库TYMP—C4直接检验完成

用户明确不以FDR通过或很强效应作为硬门槛，先检验TYMP—C4。已完成16患者/16肿瘤文库原始counts伪汇总。主结果rho0.467647、P0.067761、FDR0.203283；探索性正向趋势，不能称名义显著。16次留一患者均正，rho0.371429–0.585714；4次留一研究均正，rho0.328671–0.714286（去GSE207493只剩n6，描述性）。
研究内调整rho0.302925、置换P0.2851/FDR0.8553；方向正但研究内证据有限。替代TYMP逐细胞均值rho0.364706、P0.164868；同研究调整0.049454/P0.8671。SPP1伪汇总与检出比例补充端点独立保留。
权威结果phase154_same_library_TYMP_C4_test_20260915_v1/TYMP_C4_DIRECT_TEST_RESULTS_CN.md；原始选中23018细胞的2基因counts及文库总量已逐一核对。C4=0样本规则保留，本次16合格文库均有C4，未按C4出现筛样本。不重训/不绘图，FDR完整保留但不作统计通过门槛。
下一步围绕TYMP—C4探索性正向趋势核查空间和独立资料；原17患者及旧SPP1阴性不覆盖。本文不把相关性升级为促进C4或直接诱导SPP1。
"""
for fn in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 p=P/fn;p.write_text(note+'\n'+p.read_text(),encoding='utf-8')
p=R/'STATUS.json';d0=json.loads(p.read_text());d0['remaining_analysis'].update(current_step='direct_TYMP_C4_test_complete',direct_C4_result=str(O/'TYMP_C4_DIRECT_TEST_RESULTS_CN.md'),FDR_is_not_action_gate=True);p.write_text(json.dumps(d0,ensure_ascii=False,indent=2))
with (P/'RESULTS_MANIFEST.tsv').open('a') as f:
 for name in ['TYMP_C4_DIRECT_TEST_RESULTS_CN.md','TYMP_C4_AND_SPP1_RESULTS.tsv','INDEPENDENT_VALIDATION.json']:
  f.write('\t'.join(['phase154_'+Path(name).stem,'phase154','same_library_bridge',str(O/name),datetime.datetime.now(datetime.timezone.utc).isoformat(),'COMPUTATION_REVIEW_COMPLETE','Exploratory positive direction; no significance gate','No figures'])+'\n')
status=dict(status='COMPLETE_DIRECT_TYMP_C4_TEST_AND_REVIEW',patients=16,validation_checks=len(checks),all_checks_pass=True,primary_rho=float(primary.rho),primary_P=float(primary.P),primary_FDR=float(primary.FDR_three_endpoints),FDR_required_to_proceed=False,figures_created=False,model_retrained=False)
(O/'STATUS.json').write_text(json.dumps(status,indent=2));print(json.dumps(status));print('TECHNICAL',tech.to_string(index=False));print('SHA256',hashlib.sha256((O/'TYMP_C4_DIRECT_TEST_RESULTS_CN.md').read_bytes()).hexdigest())

