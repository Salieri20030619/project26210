from pathlib import Path
import json,datetime,hashlib,shutil,csv
import pandas as pd,numpy as np
R=Path.home()/'project26210_external_joint_v1';O=R/'phase214_literature_SPP1_signature_comparison_20260919_v1';P=R/'phase213_immune_restart_execution_20260919_v1';F=R/'phase211_immune_receiver_refreeze_20260919_v1/freeze_v2'
assert json.loads((O/'STATUS.json').read_text())['status']=='COMPUTE_COMPLETE_VALIDATED_REPORT_PENDING';assert all(json.loads((O/'VALIDATION.json').read_text()).values())
a=pd.read_csv(O/'LITERATURE_SIGNATURE_ASSOCIATIONS.tsv',sep='\t');old=pd.read_csv(P/'I2_main_statistics/PRIMARY_12_HYPOTHESES.tsv',sep='\t');old=old[old.predictor=='mature_SPP1_pseudobulk'];src=pd.read_csv(O/'FIVE_GENE_SOURCE_CONTRIBUTIONS.tsv',sep='\t')
s=src[src.source.isin(['C4','other_mature_macrophage'])].groupby(['study','patient','library','gene']).agg(mature_RNA_share=('share_of_gene_RNA','sum')).reset_index();s=s.groupby(['study','patient','gene']).mature_RNA_share.mean().reset_index();ss=s.groupby('gene').agg(patients=('patient','nunique'),median_mature_RNA_share=('mature_RNA_share','median')).reset_index();ss.to_csv(O/'SIGNATURE_GENE_SOURCE_SUMMARY.tsv',sep='\t',index=False)
comp=[]
for ep,g in a.groupby('endpoint'):
 base=old[old.endpoint==ep].iloc[0]
 for row in g.itertuples():comp.append(dict(endpoint=ep,n_patients=row.n_patients,new_predictor=row.predictor,new_rho=row.rho,new_P=row.p_value,new_FDR=row.FDR,original_single_SPP1_rho=base.rho,original_P=base.p_value,original_FDR=base.FDR,patient_and_receiver_sets_unchanged=True))
pd.DataFrame(comp).to_csv(O/'ORIGINAL_VS_LITERATURE_DEFINITION.tsv',sep='\t',index=False)
cn={'CD8_cytotoxicity':'CD8细胞毒','CD8_dysfunction':'CD8功能障碍相关程序','Treg_CD4':'Treg/CD4','NK_cytotoxicity':'NK细胞毒'};pc={'mature_macrophage_5gene_signature':'成熟巨噬内五基因','all_captured_cells_5gene_signature':'全细胞伪bulk五基因'}
lines=['| 定义 | 端点 | 患者n | 研究调整ρ | P | FDR |','|---|---|---:|---:|---:|---:|']
for r in a.itertuples():lines.append(f'| {pc[r.predictor]} | {cn[r.endpoint]} | {r.n_patients} | {r.rho:.3f} | {r.p_value:.5f} | {r.FDR:.5f} |')
report='''# 文献SPP1标准对照：固定五基因签名的探索性尝试

Phase214，2026-09-19。用户提出按文献中的SPP1标准尝试。本次为独立、版本化的文献引导探索性对照，原C4、接收端、患者名单、Phase213统计与报告均保留。

## 文献有明确可复用的五基因标准

[Jiang等，Cancer Immunology Research，2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485372/)，DOI 10.1158/2326-6066.CIR-24-1146，在Construction and validation of the SPP1+ TAM signature段规定：SPP1、FN1、CLEC5A、OLR1、IL1RN，取五基因log2(FPKM/TPM+1)表达的算术均值；已标准化的CheckMate数据直接平均，队列中位数用于高低分组。

原文另有蛋白标准：IHC/IF以CD68/SPP1双阳性定义SPP1+TAM，组织密度>12个/mm²为高组；CyTOF/流式使用CD45/CD68背景及样本SPP1中位水平。蛋白密度、蛋白门控和转录签名是不同定义。当前单细胞RNA数据不具备前两种测量，本次没有把RNA>0称为原文蛋白阳性，也未将C4重命名为文献SPP1+TAM。

## 方法：只改变SPP1表征，保持患者和接收端固定

在Phase211与Phase213的同一19个合格肿瘤文库中，抽取113469个捕获细胞的五基因原始计数及全基因总量，分别在完整成熟巨噬集合、全捕获细胞集合生成伪bulk。每基因计算log2(CPM+1)，再对固定五基因取算术均值；患者内对原端点的同一合格文库等权平均。接收端指标逐项复用Phase213，不重新定义终末功能障碍群、不添加新端点。

成熟巨噬内签名是来源限定对照；全捕获细胞签名是更接近bulk测量范围的伴随对照。两者均在计算前登记。保留连续变量的研究调整偏Spearman、99999次研究内双侧置换及5000次分层bootstrap。8项比较预先分为CD8功能4项、Treg组成2项、NK功能2项BH家族。中位数高低组仅提供描述表，无额外P值或最优截点搜索。

UMI伪bulk CPM不同于原文TPM/FPKM；全捕获细胞组成也不同于整块组织bulk。故本次是保留原五基因与聚合方式的RNA适配，不声称完整复现原队列实验。

## 结果：文献五基因标准在当前固定端点中仍未获得统计支持

'''+ '\n'.join(lines)+'''

8项预设比较P均>0.05，均未通过FDR。成熟巨噬内五基因与CD8细胞毒呈预期负方向（ρ=−0.191），但与CD8功能障碍和Treg/CD4也呈负方向，与NK细胞毒呈正方向；后面三项与预期免疫抑制方向相反。全细胞伪bulk结果同样没有形成统一的预期方向。

相对较强的两个效应也不能解释为支持：巨噬内签名—Treg/CD4为负向ρ=−0.526、P=0.08775、FDR=0.17550；全细胞签名—NK细胞毒为正向ρ=0.484、P=0.07595、FDR=0.15190。完整反向结果已保存，没有更换基因或门槛继续搜索。

巨噬内签名—Treg/CD4的分层百分位bootstrap区间不跨0，而预设置换P仍>0.05。该区间并非置换检验反演所得，小样本下两种重采样口径可不一致；两套结果均保留，以预先指定的P/FDR作为主要推断依据，不能挑选区间升级结论。

## 讨论：这次尝试回答了什么

本次说明：在现有冻结队列和平均RNA/组成端点下，将单基因SPP1改为该文献的五基因表征，尚未重现免疫抑制关联。它不否定文献的蛋白、组织及干预证据，也不能推断“SPP1没有免疫作用”。

文献的SPP1+巨噬密度、CD8终末功能障碍亚群及阻断实验，与本次患者平均RNA程序不完全等价。本轮有意只改变来源侧表征，以隔离定义差异；没有同时更换双方定义后挑选阳性。若后续核对文献的接收端，应先独立登记终末/早期状态、分母与有效患者，再开展新的探索性分支，原结果继续保留。

五基因均按原成员保留。SIGNATURE_GENE_SOURCE_SUMMARY.tsv提供其在同一捕获文库内成熟巨噬的RNA来源份额，不能把全细胞签名自动解释成SPP1+巨噬绝对密度。当前结果来自查看Phase213结果之后的新尝试，属于探索性对照，不是独立验证。

## 交付与验收

- LITERATURE_SIGNATURE_ASSOCIATIONS.tsv：全部8项效应、置换P、家族FDR、完整bootstrap区间和未经调整的对照。
- ORIGINAL_VS_LITERATURE_DEFINITION.tsv：同患者、同接收端的原单SPP1与五基因结果并排。
- FIVE_GENE_LIBRARY_SCORES.tsv、PATIENT_COMPARISON_INPUTS.tsv：完整五基因及文库/患者输入。
- MEDIAN_SPLIT_DESCRIPTIVE_ONLY.tsv：固定中位数分组描述，未生成额外显著性检验。
- PLAN_BEFORE_ANALYSIS.json、HYPOTHESIS_REGISTER.tsv、VALIDATION.json：预登记、完整检验家族与核验。

5/5基因存在；原始总计数匹配；所有患者/文库与接收端完全一致；QR独立效应、置换P、BH及五基因平均公式核验通过；原冻结SHA256保持不变。未改变C4或任何历史结果，无新图及临床/空间分析。
'''
(O/'LITERATURE_STANDARD_COMPARISON_CN.md').write_text(report,encoding='utf-8')
# Source identities and terminal state.
sourcepaths=[P/'I1_patient_inputs/PATIENT_HYPOTHESIS_INPUTS.tsv',P/'I2_main_statistics/PRIMARY_12_HYPOTHESES.tsv',P/'scripts/stats_shared.py',F/'FREEZE_SHA256.json'];(O/'INPUT_SHA256.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sourcepaths},indent=2))
seal=json.loads((F/'FREEZE_SHA256.json').read_text());assert all(hashlib.sha256((F/k).read_bytes()).hexdigest()==v for k,v in seal.items());status=json.loads((O/'STATUS.json').read_text());status['status']='COMPLETE_VALIDATED';status['completed_at']=datetime.datetime.now().astimezone().isoformat();status['report']=str(O/'LITERATURE_STANDARD_COMPARISON_CN.md');(O/'STATUS.json').write_text(json.dumps(status,indent=2))
backup=O/'root_record_backups';backup.mkdir(exist_ok=False)
for name in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md','RESULTS_MANIFEST.tsv']:shutil.copy2(R/name,backup/name)
note='''\n## 2026-09-19 Phase214：文献五基因SPP1标准对照完成\n用户提出使用文献SPP1标准尝试，已核实Jiang等CIR2025 DOI10.1158/2326-6066.CIR-24-1146的固定五基因SPP1/FN1/CLEC5A/OLR1/IL1RN及算术均值。原文蛋白双染/12个每mm²标准不能直接用于RNA。新建phase214_literature_SPP1_signature_comparison_20260919_v1，原C4/接收端/患者与Phase213不变。对同19文库分别计算成熟巨噬及全捕获细胞五基因伪bulk CPM适配，共8项固定比较，家族4/2/2。\n结果8项P均>0.05且FDR未支持。巨噬内四端点rho−0.191/−0.407/−0.526/+0.382，全细胞rho−0.240/−0.102/−0.216/+0.484。较强效应为Treg负向P0.08775/q0.1755及NK正向P0.07595/q0.1519，均与抑制预期相反，完整保留。Treg百分位bootstrap区间不跨零但置换P不支持，记录不同推断口径，不挑选升级结论。\n该任务是既有结果之后的文献引导探索性适配，不是独立复现。只有来源侧定义变化；未重定义终末CD8或修改受体端点。全部验证通过，报告LITERATURE_STANDARD_COMPARISON_CN.md；无运行任务，不继续调基因/阈值。\n'''
for name in ['PROJECT_STATUS.md','DECISIONS.md']:
 with (R/name).open('a') as f:f.write(note)
p=R/'CODEX_HANDOFF.md';p.write_text(note+'\n'+p.read_text())
with (R/'RESULTS_MANIFEST.tsv').open('a') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n')
 for name in ['LITERATURE_STANDARD_COMPARISON_CN.md','LITERATURE_SIGNATURE_ASSOCIATIONS.tsv','ORIGINAL_VS_LITERATURE_DEFINITION.tsv','PLAN_BEFORE_ANALYSIS.json','HYPOTHESIS_REGISTER.tsv','VALIDATION.json','STATUS.json']:
  w.writerow(['phase214_'+Path(name).stem,'phase214','literature_standard_exploratory',str(O/name),status['completed_at'],'COMPLETE_VALIDATED',str(O/'VALIDATION.json'),'Original results unchanged; all8 comparisons retained'])
artifacts={str(p.relative_to(O)):hashlib.sha256(p.read_bytes()).hexdigest() for p in O.rglob('*') if p.is_file() and 'root_record_backups' not in p.parts and p.name not in ['OUTPUT_SHA256.json','finalize.log']};(O/'OUTPUT_SHA256.json').write_text(json.dumps(artifacts,indent=2));print('COMPLETE_PHASE214',len(artifacts),'artifacts');print(ss.to_string(index=False))
