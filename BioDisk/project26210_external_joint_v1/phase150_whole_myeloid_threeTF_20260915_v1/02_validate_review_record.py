from pathlib import Path
import json,datetime,hashlib,shutil
import pandas as pd,numpy as np
from scipy.stats import wilcoxon
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1';O=P/'phase150_whole_myeloid_threeTF_20260915_v1'
checks=[]
def ck(name,v):checks.append(dict(check=name,passed=bool(v)))
def bh(v):
 v=np.asarray(v,float);out=np.full(len(v),np.nan);valid=np.isfinite(v);a=v[valid];idx=np.argsort(a);q=np.minimum.accumulate((a[idx]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];back=np.empty(len(a));back[idx]=np.minimum(q,1);out[valid]=back;return out
m=pd.read_csv(R/'08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str},low_memory=False)
A=R/'08_regulons/network_v2_20260915/AUC_ALL_CELLS.csv';header=pd.read_csv(A,nrows=0);tf=['JUN(+)','JUNB(+)','CEBPD(+)'];a=pd.read_csv(A,index_col=0,usecols=[header.columns[0]]+tf)
pairs=pd.read_csv(O/'PATIENT_MACRO_MINUS_MONO.tsv',sep='\t')
for scope,md in [('all_tissues',m),('tumor_only',m[m.source_tissue.str.lower()=='tumor'])]:
 mm=md[md.candidate!='1'].copy();mm['lineage']=np.where(mm.candidate=='Monocyte','Monocyte','Macrophage')
 x=a.loc[mm.index].join(mm[['curated_patient_id','lineage']])
 means=x.groupby(['curated_patient_id','lineage'])[tf].mean()
 for _,row in pairs[pairs.scope==scope].iterrows():
  diff=means.loc[(row.patient,'Macrophage'),row.regulon]-means.loc[(row.patient,'Monocyte'),row.regulon]
  ck(scope+row.patient+row.regulon+'_independent_means',np.isclose(diff,row.effect,atol=1e-12))
d=pd.read_csv(O/'WHOLE_TRAJECTORY_TF_SUMMARY.tsv',sep='\t');s=pd.read_csv(O/'MACRO_MINUS_MONO_SUMMARY.tsv',sep='\t');pt=pd.read_csv(O/'PATIENT_TF_TRENDS.tsv',sep='\t')
for _,row in d.iterrows():
 g=pt
 for key in ['mode','scope','pt_version','population','regulon']:g=g[g[key]==row[key]]
 v=g[row.adjustment].dropna().to_numpy();pv=wilcoxon(v,method='approx',correction=True).pvalue if len(v)>=3 and np.any(v) else (1. if len(v)>=3 else np.nan)
 ck('trajectory_'+str(row.name),g.patient.is_unique and len(v)==row.n_patients and np.isclose(np.median(v),row.median_effect) and np.isclose(pv,row.P,equal_nan=True) and (v>0).sum()==row.positive)
for df,keys in [(s,['scope']),(d,['mode','scope','pt_version','population','adjustment'])]:
 for key,g in df.groupby(keys):ck(str(key)+'_BH',np.allclose(bh(g.P),g.FDR_threeTF,atol=1e-12,equal_nan=True))
t=pd.read_csv(R/'07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str,'partition':str},low_memory=False)
part=pd.read_csv(O/'PARTITION_PATIENT_TF_TRENDS.tsv.gz',sep='\t',dtype={'partition':str});part=part[(part['mode']=='primary')&(part.scope=='all_tissues')&(part.pt_version=='pt_primary')&(part.population=='all_partitions')]
nonc4=0
for (study,patient,partition),g in part.groupby(['study','patient','partition']):
 z=t[(t['mode']=='primary')&(t.study==study)&(t.curated_patient_id==patient)&(t.partition==partition)&(t.candidate!='1')];z=z[np.isfinite(z.pt_primary)]
 nonc4+=int(not (z.candidate=='4').any())
 for _,row in g.iterrows():
  r=pd.Series(z.pt_primary.to_numpy()).rank().corr(pd.Series(a.loc[z.index,row.regulon].to_numpy()).rank())
  ck('independent_rho_'+str(row.name),np.isclose(r,row.rho,equal_nan=True,atol=1e-12))
ck('non_C4_partitions_included',nonc4>0)
assert all(x['passed'] for x in checks),[x for x in checks if not x['passed']]
(O/'VALIDATION.json').write_text(json.dumps(dict(status='PASS',checks=len(checks),passed=len(checks),non_C4_patient_partitions_included=nonc4,details=checks),indent=2))
def table(df):
 def fmt(v):return f'{v:.6g}' if isinstance(v,float) else str(v)
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
primary=d[(d['mode']=='primary')&(d.scope=='all_tissues')&(d.pt_version=='pt_primary')&(d.population=='all_partitions')&(d.adjustment=='rho')]
cols=['regulon','n_patients','median_effect','positive','negative','P','FDR_threeTF']
text="""# 整体单核—巨噬体系的JUN/JUNB/CEBPD活性复核

## 整体Monocle3主轨迹支持三个TF活性上升

本次取消“分区必须含C4”的条件，使用全部可估计的单核—成熟巨噬分区。JUN、JUNB、CEBPD在31名有效患者中的活性与主模型拟时序均呈正向关联，三TF家族FDR均通过。效应为患者内Spearman相关的中位数，不是表达倍数；同患者多个分区先取相关中位数，患者才是统计单位。

{primary}

因此，对“整个单核—巨噬体系在Monocle3主模型中是否支持这三个TF随拟时序提升”的回答是支持。之前围绕C4得到的结论不能替代这项整体结果；Phase149中对“三TF整体不支持”的笼统概括须由本次明确范围的结果修正。C4内群间活性降低的原数值仍然成立，两项端点不同。

## 巨噬细胞整体相对单核，仅JUN明确升高

将全部12个冻结成熟巨噬候选合并，以同患者单核为对照，双方各至少20细胞，32名患者可评估：

{state}

JUN为25/32正向，AUCell差中位约0.005180，P=0.000395、FDR=0.001185。JUNB为20/32、CEBPD为19/32正向，但当前数据不支持其群间统一升高。限定肿瘤组织后，JUN仍有支持（21/28正，FDR=0.002550），另外两个不获支持。不能将“轨迹内相关为正”解释成“巨噬细胞整体均高于单核”。

## 共同上升依赖轨迹设置，尚不能写成稳健的三级联驱动

主模型控制测序深度后，三个TF仍为正向且通过三TF FDR。25%根位置敏感性也保留三者上升；50%根位置下只有JUN保留支持，JUNB的相关中位数转为轻微负值。去除应激/增殖群后，JUN保留名义正向关联（P=0.038806、FDR=0.116419），JUNB和CEBPD不获支持。DPT整体分析三者均无统计支持。

为更贴近单核—巨噬连续变化，本次另限制同患者同分区同时至少20个单核、20个成熟巨噬；13名患者可评估。JUN为11/13正向，rho中位0.380563、P=0.004166、FDR=0.012498；JUNB为10/13正向，rho中位0.292534、P=0.021098、FDR=0.031647；CEBPD为11/13正向、rho中位0.074009、P=FDR=0.069213。CEBPD有正向趋势，但当前统计不支持把三者都写成明确上升。该限定下JUNB深度调整后P=0.036032、FDR=0.050371，为名义显著。

完整结果：
{sensitivity}

## 报告建议措辞

“在整体单核—巨噬细胞的Monocle3主轨迹中，JUN、JUNB和CEBPD调控活性均随拟时序升高，且患者级三TF家族校正后获得支持。其中JUN同时表现为成熟巨噬细胞相对单核的活性升高；JUNB和CEBPD的动态关联对根位置和轨迹方法更敏感。”

可以将该主模型结果作为原方案的方向性支持。现有结果还不支持“这三个TF在所有分支稳定级联上扬并共同驱动TYMP+/SPP1+TAM终末分化”。提升指推断regulon活性，不等于三个TF的RNA表达均升高。

## 数据与方法

复用冻结新SCENIC的235-regulon AUCell及原10个Monocle3模型，本轮只分析原方案三个TF，不重新建模、不画图。使用7,685个原Monocyte和22,444个成熟巨噬细胞（39名出现患者）；2,969个C1边界候选不并入成熟巨噬，也不参与本轮主统计。图本身仍为既有完整输入构建。主模型中7,664个单核和8,307个成熟巨噬有有效拟时序，即成熟巨噬覆盖37.01%；无根分区保持缺失，不强行补时间。

每患者分区至少50细胞、10个不同拟时序值；分区相关先汇总为患者中位数。控制深度使用总counts的秩残差偏相关。状态与轨迹分别检验；每个组织范围/模型/根版本/分区规则/调整方式内三个TF分别BH，属于本轮整体范围的三TF家族，不能冒充此前全235个regulon的FDR。n<10为探索性，n<3不作检验。

方法完成与科学解释分开：整体主轨迹的正向关联成立，但横断面排序不等于真实谱系演变；全体分区包括单一谱系占主导的分区，混合分区敏感性单独呈现。不能将C4、全部巨噬与同分区动态混为同一端点。

本轮数值核验{checks}项通过，确认纳入{nonc4}个不含C4的患者分区，并独立复算群间均值差、主轨迹相关、患者汇总及三TF BH。全部表、来源哈希和脚本在本目录。
"""
state=s[s.scope=='all_tissues'][cols]
sens=d[(d.scope=='all_tissues')&(d.population=='all_partitions')&(d.adjustment=='rho')][['mode','pt_version']+cols]
text=text.replace('{primary}',table(primary[cols])).replace('{state}',table(state)).replace('{sensitivity}',table(sens)).replace('{checks}',str(len(checks))).replace('{nonc4}',str(nonc4))
assert '???' not in text
(O/'WHOLE_MYELOID_THREE_TF_REVIEW_CN.md').write_text(text,encoding='utf-8')
shutil.copy2(R/'scripts/37_whole_myeloid_TF_review_20260915_v1.py',O/'01_compute.py')
B=O/'state_before_update';B.mkdir(exist_ok=False)
for q in [P/'CODEX_HANDOFF.md',P/'PROJECT_STATUS.md',P/'DECISIONS.md',P/'RESULTS_MANIFEST.tsv',R/'STATUS.json']:shutil.copy2(q,B/q.name)
note="""## 2026-09-15 Phase150 整体单核—巨噬三TF范围纠正

用户指出原方案范围是Monocyte与Macrophage整体，而非仅C4。新增整体分析取消“必须含C4”分区筛选，复用新SCENIC与现有轨迹，无重训/无图。权威报告phase150_whole_myeloid_threeTF_20260915_v1/WHOLE_MYELOID_THREE_TF_REVIEW_CN.md。
整体Monocle3主模型31患者：JUN/JUNB/CEBPD活性与拟时序均正相关且三TF FDR通过，rho中位分别0.148627/0.292534/0.143427；明确修正此前将C4结果概括为整体不支持的表述。群间全部成熟巨噬比单核仅JUN支持升高（32患者25正，FDR0.001185）；另两TF无群间支持。主模型深度调整和25%根敏感性保留三TF；50%根、DPT、去应激增殖不支持三者一致。混合谱系分区n13只支持JUN/JUNB，CEBPD P0.069213。
因此可以写“整体Monocle3主轨迹支持三TF活性提升”，不能写跨方法稳健、C4特异或因果级联。主模型成熟巨噬伪时间覆盖37.01%。完整统计家族为本次整体三TF，不混用此前C4全regulon FDR。原Phase149的C4数值及标签对照保留，整体范围结论以Phase150为准。当前待用户审阅，继续禁止新图。
"""
for fn in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 q=P/fn;q.write_text(note+'\n'+q.read_text(),encoding='utf-8')
q=R/'STATUS.json';obj=json.loads(q.read_text());obj['remaining_analysis'].update(current_step='whole_myeloid_TF_user_review',whole_myeloid_review=str(O/'WHOLE_MYELOID_THREE_TF_REVIEW_CN.md'));q.write_text(json.dumps(obj,ensure_ascii=False,indent=2))
with (P/'RESULTS_MANIFEST.tsv').open('a') as f:
 for name in ['WHOLE_MYELOID_THREE_TF_REVIEW_CN.md','MACRO_MINUS_MONO_SUMMARY.tsv','WHOLE_TRAJECTORY_TF_SUMMARY.tsv','VALIDATION.json']:
  f.write('\t'.join(['phase150_'+Path(name).stem,'phase150','whole_scope_analysis',str(O/name),datetime.datetime.now(datetime.timezone.utc).isoformat(),'PASS_REVIEWED','No C4 selection; 3TF family','No figures; no model refit'])+'\n')
status=dict(status='COMPLETE_WHOLE_MYELOID_THREE_TF_REVIEW',checks=len(checks),all_passed=True,nonC4_partitions=nonc4,figures_created=False,model_refit=False)
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status));print('SHA256',hashlib.sha256((O/'WHOLE_MYELOID_THREE_TF_REVIEW_CN.md').read_bytes()).hexdigest())

