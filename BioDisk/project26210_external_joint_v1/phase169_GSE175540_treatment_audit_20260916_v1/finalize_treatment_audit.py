from pathlib import Path
import pandas as pd,json,hashlib,datetime,csv,shutil,re
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase169_GSE175540_treatment_audit_20260916_v1');R=O.parent
m=pd.read_csv(R/'phase168_reference_source_audit_20260916_v1/SPATIAL_26_SECTION_PATIENT_MAP.tsv',sep='\t',keep_default_na=False)
m=m[m.dataset=='GSE175540'].copy()
m['confirmed_drug_regimen']='';m['confirmed_treatment_before_sampling']='';m['confirmed_sampling_time_relative_to_therapy']=''
m['secondary_RI_NRI_label']=m.source_id.map(lambda x:'RI' if x.startswith('c_') else 'NRI')
m['secondary_label_source']='IJMS2024 25:3887 supplementary Table S1; not adopted as verified exposure'
m['treatment_status']='UNRESOLVED_NO_PRIMARY_PATIENT_CROSSWALK'
m.to_csv(O/'GSE175540_TREATMENT_AUDIT_24_SECTIONS.tsv',sep='\t',index=False)
assert len(m)==24 and m.patient_id.nunique()==23 and (m.secondary_RI_NRI_label=='RI').sum()==16
soft=(R/'phase165_A_input_confirmation_20260915_v1/metadata/GSE175540_family.soft').read_text()
selected=[]
for b in soft.split('^SAMPLE = ')[1:]:
 gsm=b.splitlines()[0].strip()
 if gsm in set(m.section):
  selected.append(dict(section=gsm,characteristics=[l for l in b.splitlines() if l.startswith('!Sample_characteristics')],treatment_characteristic_present=any(re.search(r'treat|therap|drug|nivolumab|ipilimumab',l,re.I) for l in b.splitlines() if l.startswith('!Sample_characteristics'))))
assert len(selected)==24 and not any(x['treatment_characteristic_present'] for x in selected)
(O/'GEO_CHARACTERISTICS_REVIEW.json').write_text(json.dumps(selected,indent=2))
report="""# GSE175540逐例治疗核查：未找到可核验的具体方案及取样时点

## 结论与处理
本轮重新检查GEO原始24条样本信息、原始Immunity研究已取得的三份补充表与补充图、Science Advances已有临床表、BIONIKK官方试验登记，并追溯一篇再分析的全文和补充表。没有找到可将23位GSE175540患者逐一对应到具体药物、治疗起始时间及组织取样时间的原始证据。正式字段继续缺失；不能把缺失填成未治疗，不能由BIONIKK身份推定取样前已用ICI。不继续无限扩展检索，本轮关闭该项资料核查，保留以后新证据补录入口。

## 各来源核查
1. GEO GSE175540：24条样本characteristics均无治疗字段，来源及疾病信息有记录。
2. 原始官方Table S1（mmc2.xlsx）：130例对象、pT和各实验的实施情况，已解决患者映射；无药物及取样相对用药时间列。
3. 原始官方Table S2（mmc3.xlsx）：BCR克隆及TLS状态；无用药对应。
4. 本轮新增官方Table S3（mmc4.xlsx）：TLS、浆细胞、Ig染色、凋亡、NK及巨噬细胞密度；无治疗列。完整表只读留档，没有重算作者实验。
5. Science Advances合并稿补充表2：GSE175540的Treatment与临床Stage为NA；另两例GSE210041治疗记录不能移用至本队列。
6. BIONIKK NCT02960906登记提供Nivolumab、Nivolumab+Ipilimumab、VEGFR-TKI的试验臂及入排标准，但没有c_2等研究ID至治疗臂的对照。登记要求转移性病程未接受既往系统治疗、分子分型在随机前完成；还规划治疗前及进展时采样，故不能据试验设计把每张GEO片一概标为治疗前或治疗后。也不能将“无既往mRCC系统治疗”扩写成终生从未接受任何治疗。
7. IJMS2024再分析（10.3390/ijms25073887）补充Table S1确有24切片的RI/NRI标签：全部16个c组切片为RI，8个a/b组切片为NRI。正文同时将c组说明为可能接受N、NI或TKI，未给逐例药物和时间。VEGFR-TKI不等同ICI，患者之后接受治疗也不证明取样时已暴露。因此该标签只存为二次文献注释，不写入confirmed treatment字段，不用于主分析的治疗分层。

## 对当前分析的影响
23位患者映射及c_2双制备结论保持。Phase167仍为上游24片/23患者、接收端20片/19患者。可继续患者层空间方向和效应分析，来源队列及制备方式保留为明确已知信息；治疗信息缺失统一放入研究边界，不据二次RI/NRI标签声称治疗效应、ICI耐药或治疗前后改变。无新关联检验、无分析图、未联系作者。

## 可追溯来源
- GEO：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175540
- 官方补充Table S1：https://ars.els-cdn.com/content/image/1-s2.0-S1074761322000814-mmc2.xlsx
- 官方补充Table S2：https://ars.els-cdn.com/content/image/1-s2.0-S1074761322000814-mmc3.xlsx
- 官方补充Table S3：https://ars.els-cdn.com/content/image/1-s2.0-S1074761322000814-mmc4.xlsx
- 官方试验：https://clinicaltrials.gov/study/NCT02960906
- 再分析全文：https://doi.org/10.3390/ijms25073887
- 再分析补充表：https://mdpi-res.com/d_attachment/ijms/ijms-25-03887/article_deploy/ijms-25-03887-s001.zip

## 可访问性与核验
官方原文网页/PDF本轮仍403，不能声称已经检查其完整主文图中所有患者标签；额外mmc5地址404。上述失败响应留档并明确不作为PDF/XLSX有效输入。原始三份补充表、试验登记、再分析全文XML及补充表已取得。再分析补充表第1页已逐行视觉核对。输出24条唯一GSM、23独立ID，与Phase168映射一致；全部confirmed字段留空，16/8二次文献标签另列。用户已明确找不到可接受，本轮未找到充分原始证据，不将线索强行补成事实。
"""
(O/'TREATMENT_AUDIT_CN.md').write_text(report,encoding='utf-8')
now=datetime.datetime.now().astimezone().isoformat()
status=dict(stage='Phase169',status='AUDIT_COMPLETE_TREATMENT_UNRESOLVED',sections=24,patients=23,confirmed_individual_drug_regimens=0,confirmed_sampling_relative_to_therapy=0,secondary_RI_sections=16,secondary_NRI_sections=8,secondary_labels_adopted_for_analysis=False,association_tests=0,new_analysis_figures=0,completed_at=now)
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
pd.DataFrame([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(O.rglob('*')) if p.is_file() and 'authority_backups' not in p.parts and p.name not in ['ARTIFACT_HASHES.tsv','FINAL_STATE_VALIDATION.json']]).to_csv(O/'ARTIFACT_HASHES.tsv',sep='\t',index=False)
B=O/'authority_backups';B.mkdir(exist_ok=False)
names=['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']
for n in names:shutil.copy2(R/n,B/n)
notice=f"""
## Phase169 GSE175540逐例治疗核查结束（{now}）
原始GEO及官方补充S1/S2/S3未找到23患者逐例药物及取样相对治疗时间；正式治疗字段保留缺失。BIONIKK试验提供臂/入排标准，未给c_ID映射。IJMS2024补充S1有c组16片RI、a/b8片NRI，正文混含N/NI/TKI且无时点，仅二次文献注释，不作验证过的ICI暴露或耐药分组。
报告：{O/'TREATMENT_AUDIT_CN.md'}
按用户“实在找不到也没办法”，本轮停止扩展该项检索，可继续患者层空间分析。Phase168患者映射有效：上游24片/23患者，接收20片/19患者；需先固定c_2双制备处理。无新关联/图。
"""
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 p=R/n;t=p.with_name(n+'.phase169.tmp');t.write_text(notice+p.read_text(),encoding='utf-8');t.replace(p)
p=R/'RESULTS_MANIFEST.tsv'
with p.open(newline='') as f:rows=list(csv.reader(f,delimiter='\t'))
assert len(rows[0])==8 and not any(x and x[0].startswith('phase169_') for x in rows[1:])
reg=['TREATMENT_AUDIT_CN.md','GSE175540_TREATMENT_AUDIT_24_SECTIONS.tsv','GEO_CHARACTERISTICS_REVIEW.json','original_mmc4.xlsx','MMC4_EXTRACT.json','bionikk_trial.json','reanalysis_fulltext.xml','reanalysis_supp.zip','FETCH_LOG.json','SECONDARY_FETCH_LOG.json','STATUS.json','ARTIFACT_HASHES.tsv','finalize_treatment_audit.py']
for n in reg:rows.append(['phase169_'+n,'Phase169','treatment_source_audit',str(O/n),now,'PASS','Audit complete; confirmed treatment unavailable','sha256='+sha(O/n)])
t=p.with_name(p.name+'.phase169.tmp')
with t.open('w',newline='') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerows(rows)
t.replace(p)
gate=dict(stage='Phase169',registered_paths_exist=all((O/n).exists() for n in reg),state_files_updated=all((R/n).read_text().startswith(notice) for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']),patient_mapping_retained=True,rows=24,confirmed_fields_blank=bool((m.confirmed_drug_regimen=='').all() and (m.confirmed_treatment_before_sampling=='').all()))
assert all(v for k,v in gate.items() if isinstance(v,bool))
(O/'FINAL_STATE_VALIDATION.json').write_text(json.dumps(gate,indent=2));print(json.dumps(status));print(json.dumps(gate))
