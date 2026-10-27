from pathlib import Path
import json,sys,datetime,hashlib,subprocess,shutil
import numpy as np,pandas as pd,scipy
from scipy.stats import binomtest,false_discovery_control,spearmanr
P=Path(__file__).resolve().parents[1];ROOT=P.parent;OLD=ROOT/'phase234_part4_replan_execution_20260921_v1'
def now():return datetime.datetime.now().astimezone().isoformat()
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2))
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
r=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(P/'logs/validation_memory.log').write_text(r.stdout+r.stderr);assert r.returncode==0
assert (P/'EXIT_STATUS').read_text().strip()=='0' and json.loads((P/'FINAL_VALIDATION.json').read_text())['status']=='PASS'
checks=[];keys=['dataset','family','endpoint','model','effect_type']
for prefix in ['MARKER_STRUCTURE','COMPOSITE_RECHECK']:
 su=pd.read_csv(P/'tables'/f'{prefix}_COHORT_SUMMARY.tsv',sep='\t')
 pat=pd.read_csv(P/'tables'/f'{prefix}_PATIENT_EFFECTS.tsv',sep='\t')
 su['n_nonzero']=su.positive+su.negative;su['FDR_wilcoxon']=np.nan
 for row in su.itertuples():
  m=np.ones(len(pat),bool)
  for k in keys:m &= pat[k].eq(getattr(row,k)).to_numpy()
  g=pat[m];v=g.effect.dropna().to_numpy()
  assert g.patient.is_unique and len(v)==row.n
  if len(v):assert abs(float(np.median(v))-row.median)<1e-12
  if pd.notna(row.P):
   pp=float(binomtest(int((v>0).sum()),int((v!=0).sum()),.5).pvalue) if (v!=0).any() else 1.
   assert abs(pp-row.P)<1e-12
  else:assert len(v)<10
  checks.append(dict(module=prefix,dataset=row.dataset,endpoint=row.endpoint,model=row.model,n=row.n,nonzero=row.n_nonzero,status='PASS'))
 for key,g in su.groupby(['dataset','family']):
  q=false_discovery_control(g.P.fillna(1).to_numpy())
  assert np.allclose(q[g.P.notna()],g.loc[g.P.notna(),'FDR'])
  qw=false_discovery_control(g.P_wilcoxon.fillna(1).to_numpy())
  su.loc[g.index,'FDR_wilcoxon']=np.where(g.P_wilcoxon.notna(),qw,np.nan)
 su.to_csv(P/'tables'/f'{prefix}_COHORT_SUMMARY_WITH_EFFECTIVE_N.tsv',sep='\t',index=False)
# Marker score/mapping relations independently recomputed from exported spot tables.
secs=pd.read_csv(P/'tables/MARKER_STRUCTURE_SECTION_EFFECTS.tsv',sep='\t')
regs=pd.read_csv(OLD/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t');corrchecks=0
for r in regs.itertuples():
 d=pd.read_csv(OLD/'inputs'/f'{r.sample}_INPUTS.tsv.gz',sep='\t')
 s=pd.read_csv(P/'inputs'/f'{r.sample}_MARKERS.tsv.gz',sep='\t');assert np.array_equal(d.spot_id,s.spot_id)
 for target,col in [('Macrophage','total_macro_fraction'),('Tcell','Tcell_fraction'),('Endothelial','endothelial_fraction')]:
  row=secs[(secs['sample']==r.sample)&(secs.endpoint==col+'__'+target)].iloc[0]
  if pd.notna(row.effect):
   rr=spearmanr(d[col],s[target]).statistic;assert abs(rr-row.effect)<1e-12;corrchecks+=1
# Reproduction compares exact original composite correlations rather than claiming a new smaller legacy FDR.
orig=pd.read_csv('/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/tymp_vascular_niche_v1/tables/sample_spatial_block_tests_v1.tsv',sep='\t')
current=pd.read_csv(P/'tables/COMPOSITE_RECHECK_SECTION_EFFECTS.tsv',sep='\t')
count=0
for r in current[current.endpoint=='niche__T_cell_balance'].itertuples():
 g=orig[(orig['sample']==r.sample)&(orig.grid==r.model)&(orig.source=='spp1_vascular_niche')&(orig.target=='T_cell_balance')]
 assert len(g)==1 and abs(g.rho.iloc[0]-r.effect)<1e-12;count+=1
pd.DataFrame(checks).to_csv(P/'tables/INDEPENDENT_PATIENT_STATISTIC_CHECKS.tsv',sep='\t',index=False)
lines=['# 空间补充验证：去共享成分后保留T细胞相对平衡线索，CD8区域分离仍需独立确认','',
'## Results：去SPP1分母后仍有负向关系，直接T细胞丰度未降低','TJ-RCC原复合生态位—T-cell balance两尺度均12/12负向。移除参照中的SPP1丰度后，两尺度为10/12负向，中位rho分别−0.4897与−0.3873，主要符号P均0.03857，本轮六项重检家族FDR均0.05786，属于名义显著。原始结果的原家族FDR0.002197仍保留，不以本轮缩小家族重算值替换。','直接T细胞绝对丰度两尺度均9/12正向，中位rho0.2662/0.2673，P=0.1460。负向相对平衡不能解释为T细胞绝对量下降；去共享后仍有关系，也不能把原关联全部归因于数学构造。','',
'## Results：GSE175540出现CD8表达区域分离线索，调整后邻域关联减弱','固定400um、SPP1巨噬内比例与CD8谱系程序：GSE175540有效20患者中16负，中位rho−0.2814，P0.01182/FDR0.03545；加入来源和邻域深度、组成及来源位置目标程序后，rho−0.09162，P0.2632/FDR0.3948。','按固定高/低20%分层，GSE175540有效16患者中12正、3负、1零；SPP1高比例区域到CD8表达富集区的中位距离差+300um，P0.03516/FDR0.05273。距离是组织图上的表达区域距离。TJ-RCC同一距离端点仅5例合格，2正2负1零，中位差0，未形成重复方向；CD8邻域仅9例合格，按描述性输出不检验。','',
'## Results：既有映射与谱系RNA一致，CAF和普遍正向血管结构未获补强','巨噬、T细胞、内皮映射与各自谱系RNA在多数来源具有一致正向关系。SPP1巨噬内比例—SPP1 RNA在GSE175540为20/23正，rho0.2485/P0.0004883/FDR0.0004883；TJ-RCC10/12正，rho0.1784/P0.03857/FDR0.03857。','不共享SPP1的内皮标记400um邻域在GSE175540呈负相关：19/23负，rho−0.2442/P0.002599/FDR0.01560；完整调整后P0.09314。成纤维标记邻域和距离没有一致支持SPP1邻近增强。内皮距离中大量病例为0差值，完整有效n、正负零数保留，不能仅凭非零病例的符号P宣称普遍分离。','',
'## Methods','输入为原38张切片114640个QC spot，按原患者表汇总：GLMF两例、TJ-RCC12例、GSE17554023例（c_2双制备合并）。读取原始counts，按每spot总量归一到10000后log1p，片内基因Z分数等权聚合。完整标记名单、共同检出和基因完整性资格在PLAN_BEFORE_ANALYSIS.json冻结。','CD8目标要求CD8B与至少一个T谱系标记共同检出；所有标记集至少20个spot满足检出条件。复用100um第一圈空间图与400um可达邻域，排除自身，不跨组织分量。来源为成熟巨噬构成≥0.05。距离比较以来源SPP1状态高/低20%及目标谱系表达高20%冻结，目标至少20spot、两来源组各至少20可达spot。','每片一个效应；同患者相关系数Fisher-z等权合并、距离差算术均值。患者n≥10采用双侧精确符号检验，Wilcoxon伴随、10000次患者bootstrap；不足10只描述。所有正负零和非零有效n分别输出。按来源×模块家族BH，邻域全部3端点×2模型共同6项，区域距离3项，映射一致性4项；复合重检3定义×2尺度共同6项。','',
'## Discussion','这些新分析是已知历史结果基础上的锁定探索。谱系RNA与原映射来自同一表达数据，标记并非原模型训练留出，因此属于RNA一致性，不是独立蛋白或单细胞身份验证。旧SPP1-state不改名为C4；CD8/内皮/成纤维表达区域不等同已分割细胞、接触事件或病理界面。','距离线索与调整减弱、TJ-RCC资格不足同时保留。它支持进一步核对CD8区域分布，尚不能确认SPP1导致CD8排斥。原血管绝对丰度、转录背景及复合生态位结果保持各自定义，不用本轮阴性否定不同端点的既有阳性。','',
'## 复现与验收',f'独立核对{len(checks)}项患者汇总、{corrchecks}项映射标记相关和{count}条原复合相关；原始H5AD只读、无新图、无cell2location重训。代码code/run_analysis.py，环境part4_advanced_py_v1，完整结果见tables/。']
(P/'REPORT_CN.md').write_text('\n'.join(lines))
dump(P/'INDEPENDENT_VALIDATION.json',dict(status='PASS',at=now(),patient_summary_rows=len(checks),independent_spearman_checks=corrchecks,legacy_rho_reproduced=count,all_BH_checked=True,effective_nonzero_n_exported=True,new_figures=0))
dump(P/'ENVIRONMENT.json',dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
print('VALIDATED',len(checks),corrchecks,count)

