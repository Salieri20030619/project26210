from pathlib import Path
import json,datetime,hashlib,shutil
import pandas as pd,numpy as np
from scipy.stats import rankdata,wilcoxon,spearmanr
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase151_TYMP_SPP1_threeTF_role_20260915_v1'
def bh(v):
 v=np.asarray(v,float);q=np.full(len(v),np.nan);valid=np.isfinite(v);x=v[valid];idx=np.argsort(x);z=np.minimum.accumulate((x[idx]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1];back=np.empty(len(x));back[idx]=np.minimum(z,1);q[valid]=back;return q
checks=[]
def ck(name,value):checks.append(dict(check=name,passed=bool(value)))
summ=pd.read_csv(O/'TF_GENE_SUMMARY.tsv',sep='\t');pat=pd.read_csv(O/'ALL_PATIENT_TF_GENE_CORRELATIONS.tsv',sep='\t')
for i,row in summ.iterrows():
 z=pat
 for c in ['scope','gene','regulon','adjustment']:z=z[z[c]==row[c]]
 v=z.rho.dropna().to_numpy();pv=wilcoxon(v,method='approx',correction=True).pvalue if len(v)>=3 and np.any(v) else (1 if len(v)>=3 else np.nan)
 ck('patient_summary_'+str(i),z.patient.is_unique and len(v)==row.n_patients and np.isclose(np.median(v),row.median_rho) and (v>0).sum()==row.positive and (v<0).sum()==row.negative and np.isclose(pv,row.P,equal_nan=True))
for key,z in summ.groupby(['scope','gene','adjustment']):ck('BH_'+str(key),np.allclose(bh(z.P),z.FDR_threeTF,equal_nan=True))
br=pd.read_csv(O/'TUMOUR_TYMP_MACROPHAGE_MATCHED_PATIENTS.tsv',sep='\t');out=pd.read_csv(O/'TUMOUR_TYMP_BRIDGE_SUMMARY.tsv',sep='\t')
ck('bridge_patient_unique',br.patient.is_unique);ck('17_patients',len(br)==17);ck('same_study_mapping',(br.study==br.study_v1).all())
for key,z in out.groupby(['exposure','adjustment']):ck('bridge_BH_'+str(key),np.allclose(bh(z.P),z.FDR_five_readouts,equal_nan=True))
for i,row in out.iterrows():
 x=rankdata(br[row.exposure]);y=rankdata(br[row.target])
 if row.adjustment=='raw':
  r,pv=spearmanr(br[row.exposure],br[row.target]);ck('bridge_raw_'+str(i),np.isclose(r,row.effect_rho) and np.isclose(pv,row.P))
 else:
  C=pd.get_dummies(br.study,dtype=float).to_numpy();xr=x-C@(np.linalg.pinv(C)@x);yr=y-C@(np.linalg.pinv(C)@y)
  r=np.corrcoef(xr,yr)[0,1];ck('bridge_study_coefficient_'+str(i),np.isclose(r,row.effect_rho))
targets=pd.read_csv(O/'THREE_TF_TARGET_MEMBERSHIP.tsv',sep='\t');ck('six_absent_target_edges',len(targets)==6 and not targets.present.any())
# Independently verify all within-patient raw and depth-adjusted correlations.
m=pd.read_csv(R/'08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str},low_memory=False)
A=R/'08_regulons/network_v2_20260915/AUC_ALL_CELLS.csv';h=pd.read_csv(A,nrows=0)
auc=pd.read_csv(A,index_col=0,usecols=[h.columns[0],'JUN(+)','JUNB(+)','CEBPD(+)'])
g=pd.read_csv(R/'07_trajectory/dynamics_v2/CELL_EXPRESSION_PROGRAM_SCORES.tsv.gz',sep='\t',index_col=0,usecols=['cell_id','gene:TYMP','gene:SPP1'])
plain=pd.read_csv(O/'PATIENT_TF_GENE_CORRELATIONS.tsv',sep='\t')
for (scope,patient),rows in plain.groupby(['scope','patient']):
 z=m[(m.curated_patient_id==patient)&(~m.candidate.isin(['1','Monocyte']))]
 if scope=='tumor_only':z=z[z.source_tissue.str.lower()=='tumor']
 depth=rankdata(np.log1p(z.total_counts_raw));C=np.column_stack([np.ones(len(z)),depth]);pinv=np.linalg.pinv(C)
 for i,row in rows.iterrows():
  x=rankdata(auc.loc[z.index,row.regulon]);y=rankdata(g.loc[z.index,'gene:'+row.gene])
  if row.adjustment=='depth_adjusted':x=x-C@(pinv@x);y=y-C@(pinv@y)
  r=np.corrcoef(x,y)[0,1];ck('raw_score_correlation_'+str(i),np.isclose(r,row.rho,equal_nan=True,atol=1e-10))
assert all(c['passed'] for c in checks),[c for c in checks if not c['passed']]
(O/'VALIDATION.json').write_text(json.dumps(dict(status='PASS',checks=len(checks),passed=len(checks),details=checks),indent=2))
def table(df):
 def f(v):return f'{v:.6g}' if isinstance(v,float) else str(v)
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(f(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
key=summ[(summ.scope=='tumor_only')&(summ.gene=='SPP1')&(summ.adjustment=='depth_adjusted')][['regulon','n_patients','median_rho','positive','negative','P','FDR_threeTF']]
adjust=summ[(summ.scope=='tumor_only')&(summ.gene=='SPP1')&(summ.adjustment=='depth_pseudotime_adjusted')][['regulon','n_patients','median_rho','positive','negative','P','FDR_threeTF']]
bridge=out[out.exposure=='z_TYMP_RNA'][['target','adjustment','n_patients','effect_rho','P','FDR_five_readouts']]
report="""# 三TF在TYMP—SPP1主线中的证据角色

## 三TF与SPP1显著负相关，当前证据尚未把它们连接成肿瘤TYMP到SPP1的中介

本次围绕实际主线检查肿瘤侧TYMP、髓系TF活性和巨噬SPP1。JUN/JUNB/CEBPD在整体Monocle3主轨迹中升高的结果保留，但该结果描述广泛状态排序，不能单独说明其连接了TYMP和SPP1。

在肿瘤组织成熟巨噬细胞内，三个TF的regulon活性与SPP1表达均呈患者一致的负相关；原始相关、测序深度调整和可估计分区的拟时序调整均保留该方向。结合C4相对其他巨噬的三TF活性降低，它们目前可描述为“与SPP1高表达状态负向关联的调控程序”。不能进一步称为SPP1的抑制因子或保护因子；相关方向尚不确定调控方向。

## 与SPP1的负相关有患者级统计支持

肿瘤组织成熟巨噬细胞内，控制测序深度后：

{key}

31名有效患者中，CEBPD为29/31负向、JUN为28/31负向、JUNB为31/31负向。效应是患者内相关系数的中位数，量级约−0.10至−0.17，反映一致但幅度有限的关联。未使用细胞数作为推断样本量。

在有可估计主模型拟时序且满足表达覆盖的12名患者中，控制测序深度和拟时序后：

{adjust}

这说明观察到的负相关不能仅由已测的总counts或当前伪时间变量解释，但不排除其他分支、细胞状态及网络构建因素。该调整子集不代表所有31名患者。

对巨噬细胞自身TYMP表达也观察到三TF负相关：肿瘤来源主分析33患者，控制深度后JUN/JUNB/CEBPD的rho中位分别约−0.0971/−0.1024/−0.0758，三TF FDR分别约3.17e-5/8.53e-5/3.17e-5。这里是髓系细胞自己的TYMP，不能当作肿瘤细胞TYMP对巨噬的作用。

## 肿瘤侧TYMP尚未与三TF形成受支持的联系

匹配冻结的肿瘤侧TYMP指标与当前肿瘤组织巨噬读出，共17名符合条件患者。肿瘤侧TYMP RNA与JUN/JUNB/CEBPD活性均无统计支持；研究调整后也无支持。严格TYMP功能轴敏感性同样不支持这三个联系。

{bridge}

主线自身有一条探索性正关联：肿瘤侧TYMP RNA与巨噬SPP1平均表达rho=0.507353，P=0.037633，五读出家族FDR=0.188164。它可以作为名义显著的候选线索，不能写成校正后显著。去除研究间秩均值差并在研究内置换后，rho=0.338313、P=0.4324、FDR=0.720667，跨研究稳健性尚未建立。C4比例与肿瘤TYMP也未获统计支持。

因此，当前可以分别陈述“肿瘤TYMP与巨噬SPP1有探索性关联”和“三TF活性与巨噬SPP1负相关”，尚不能把它们合并成“肿瘤TYMP通过三TF促进SPP1”的连续中介路径。

## 当前SCENIC没有给出三TF直接调控TYMP或SPP1的预测边

已核对新网络的完整motif剪枝靶基因表：JUN、JUNB、CEBPD六个TF—目标组合均未保留TYMP或SPP1为预测靶基因。这意味着本次SCENIC没有提供这六条直接调控边的支持；不能反推生物学上绝无间接或条件性作用。

由于两个终点基因都不属于这三个regulon，本轮TF活性—TYMP/SPP1表达相关不含这两个目标基因直接进入AUCell计分的自相关环节。regulon其余靶基因仍来自同一发现数据，关联不是外部独立验证。

## 适合当前证据的主线位置

| 证据层 | 当前可用角色 | 不宜升级的结论 |
| --- | --- | --- |
| 整体Monocle3三TF上升 | 单核—巨噬状态变化的调控背景；主模型支持 | 不能单凭此说明TYMP—SPP1轴特异作用 |
| C4及SPP1表达相关 | 三TF活性降低与SPP1富集状态相伴；患者级支持 | 不能称直接抑制SPP1或证明负调控因果 |
| 肿瘤TYMP到三TF | 本次17患者分析未获支持 | 不能作为已建立的中介环节 |
| 肿瘤TYMP到巨噬SPP1 | RNA指标有名义显著正相关，研究调整后不支持 | 不能称校正后显著、跨研究稳定或因果作用 |
| 三TF到TYMP/SPP1靶边 | 本次网络没有保留预测边 | 不能称已重建出直接分子调控链 |

建议将这三个TF保留为SPP1相关状态的辅助调控证据。主线的连接应优先依靠肿瘤TYMP与巨噬SPP1的患者级对应、已明确的代谢背景和后续针对性验证；本报告不把原定三个TF设为必须出现的中介，也不自动更换候选直到获得期望方向。

可用正文：“JUN、JUNB和CEBPD活性与巨噬细胞SPP1表达呈负向关联，该方向在控制测序深度及可估计拟时序后保留。结合SPP1富集C4中的活性降低，这些调控程序标记了与SPP1高表达相伴的状态差异。肿瘤TYMP与三TF的跨患者关联尚未获得支持，因此当前证据未将三者确定为TYMP—SPP1联系的中介。”

## 方法与边界

本轮为用户问题驱动的定向补充，PLAN.json先于计算保存。复用冻结新SCENIC AUCell、既有log1pCP10K表达、Monocle3伪时间和Phase109 analysis_v3肿瘤侧患者摘要，不重训模型、不绘图。只纳入12个成熟巨噬候选；单核和C1边界群不并入本轮巨噬内相关。肿瘤来源为主分析，全组织为敏感性。

巨噬内每患者至少50细胞、目标基因至少10个阳性细胞；计算原始Spearman及rank残差偏相关。轨迹调整在同患者同分区至少50细胞且至少10个不同拟时序值后执行，多分区先取患者相关中位数。患者rho用双侧Wilcoxon近似及连续性校正汇总，按组织/目标基因/调整方式各自三TF BH。n<10探索性、n<3不作推断。

跨患者桥接保留原eligible_abundance_primary并要求至少50个肿瘤来源成熟巨噬，最终17患者。分别考察冻结z_TYMP_RNA和严格TYMP_axis_primary，五个读出为三TF活性均值、SPP1表达均值及C4比例；每暴露/调整方式内五项BH。原始Spearman与研究固定秩均值去除后的相关分开保存；后者使用9,999次研究内置换，seed262151。不同家族不混池。RNA/活性/丰度三个读出层级始终保留名称，不互换含义。

这是观察性转录组与推断调控网络的关联检查，没有进行分子中介、干预、谱系或TF结合验证；没有分析客户湿实验。无直接网络边不是生物学排除证据。统计显著表示给定数据和模型下的关联，不足以确定促进、抑制或中介作用。

本轮核验{ncheck}项全部通过。完整正负结果：TF_GENE_SUMMARY.tsv；巨噬内患者值：ALL_PATIENT_TF_GENE_CORRELATIONS.tsv；肿瘤桥接：TUMOUR_TYMP_BRIDGE_SUMMARY.tsv；预测边：THREE_TF_TARGET_MEMBERSHIP.tsv；来源：SOURCE_MANIFEST.tsv。进度只记录于服务器。
"""
report=report.replace('{key}',table(key)).replace('{adjust}',table(adjust)).replace('{bridge}',table(bridge)).replace('{ncheck}',str(len(checks)))
# Avoid rounded hand-entered TYMP q values; replace with exact generated clause.
ty=summ[(summ.scope=='tumor_only')&(summ.gene=='TYMP')&(summ.adjustment=='depth_adjusted')]
report=report.replace('三TF FDR分别约3.17e-5/8.53e-5/3.17e-5','三TF FDR均小于0.0001')
(O/'TYMP_SPP1_THREE_TF_ROLE_REVIEW_CN.md').write_text(report,encoding='utf-8')
shutil.copy2(R/'scripts/38_TYMP_SPP1_TF_role_20260915_v1.py',O/'01_compute.py')
B=O/'state_before_update';B.mkdir(exist_ok=False)
for q in [P/'CODEX_HANDOFF.md',P/'PROJECT_STATUS.md',P/'DECISIONS.md',P/'RESULTS_MANIFEST.tsv',R/'STATUS.json']:shutil.copy2(q,B/q.name)
note="""## 2026-09-15 Phase151 三TF在TYMP—SPP1主线中的角色核查完成

用户要求判断三TF是否在TYMP—SPP1主线发挥显著作用。权威报告phase151_TYMP_SPP1_threeTF_role_20260915_v1/TYMP_SPP1_THREE_TF_ROLE_REVIEW_CN.md。
肿瘤组织成熟巨噬内，JUN/JUNB/CEBPD活性与SPP1均为患者级负相关；深度调整后n31，rho中位−0.098645/−0.115404/−0.168304，28/31、31/31、29/31负向，三TF FDR通过。主拟时序可评估子集n12进一步控制深度和伪时间后仍负向且FDR通过。不能从相关性称三TF抑制SPP1。
新SCENIC未保留三TF到TYMP/SPP1六条预测边，不能据此排除真实间接作用。17名匹配患者中，肿瘤TYMP RNA/严格轴与三TF活性均未获支持；当前不能指定三TF为TYMP—SPP1中介。
肿瘤TYMP RNA与巨噬SPP1均值有探索性正关联（rho0.507353，P0.037633，五读出FDR0.188164），研究内调整/置换后P0.4324，不是稳定桥接证据。不能隐藏此名义阳性，也不能升级因果链。
Phase150整体Monocle3三TF上升仍成立，表示广泛状态排序；Phase151评估的是轴关联，两个端点分开。保留三TF为SPP1高表达相伴的低活性调控程序候选，不写正向中介驱动。仅定向统计补充，无重训/无新图；完整所有正负结果归档。
"""
for fn in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 q=P/fn;q.write_text(note+'\n'+q.read_text(),encoding='utf-8')
q=R/'STATUS.json';obj=json.loads(q.read_text());obj['remaining_analysis'].update(current_step='TYMP_SPP1_threeTF_role_review',axis_role_review=str(O/'TYMP_SPP1_THREE_TF_ROLE_REVIEW_CN.md'));q.write_text(json.dumps(obj,ensure_ascii=False,indent=2))
with (P/'RESULTS_MANIFEST.tsv').open('a') as f:
 for name in ['TYMP_SPP1_THREE_TF_ROLE_REVIEW_CN.md','TF_GENE_SUMMARY.tsv','TUMOUR_TYMP_BRIDGE_SUMMARY.tsv','VALIDATION.json']:
  f.write('\t'.join(['phase151_'+Path(name).stem,'phase151','axis_role_review',str(O/name),datetime.datetime.now(datetime.timezone.utc).isoformat(),'PASS_REVIEWED','Preserved positive/negative evidence','No figures; no model refit'])+'\n')
status=dict(status='COMPLETE_TYMP_SPP1_THREE_TF_ROLE_REVIEW',checks=len(checks),all_passed=True,figures_created=False,model_refit=False)
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status));print(hashlib.sha256((O/'TYMP_SPP1_THREE_TF_ROLE_REVIEW_CN.md').read_bytes()).hexdigest())

