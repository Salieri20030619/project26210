from pathlib import Path
import json,hashlib,datetime,shutil,zipfile,html,os,csv,io
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent;D=P/'delivery';stamp=datetime.datetime.now().astimezone().isoformat()
v=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());assert v['status']=='PASS'
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
e=pd.read_csv(P/'tables/RECIPIENT_ASSOCIATIONS_ALL.tsv',sep='\t');sp=pd.read_csv(P/'tables/EXTERNAL_SPATIAL_EFFECTS.tsv',sep='\t');pe=pd.read_csv(P/'tables/EXTERNAL_STIMULATION_EFFECTS.tsv',sep='\t');ko=pd.read_csv(P/'tables/EXTERNAL_KO_EFFECTS.tsv',sep='\t');a=pd.read_csv(P/'tables/RCC_ALL_RUN_METADATA.tsv',sep='\t')
audit=[]
for n in range(1,8):
 b=a[a.alias.eq(f'B{n}')].iloc[0];s=a[a.alias.eq(f'SPP1_{n}')].iloc[0]
 audit.append({'numeric_suffix':n,'B_run':b.run_accession,'B_alias':b.alias,'B_age':b.age,'SPP1_run':s.run_accession,'SPP1_alias':s.alias,'SPP1_age':s.age,'age_matches':b.age==s.age,'pairing_status':'NOT_CONFIRMED_DO_NOT_PAIR_BY_SUFFIX','treatment_status':'alias_not_explicit_treatment_metadata'})
pd.DataFrame(audit).to_csv(P/'tables/RCC_CANDIDATE_PAIR_AUDIT.tsv',sep='\t',index=False)
assert not any(x['age_matches'] for x in audit)
source_audit={
 'RCC_paper':'https://pmc.ncbi.nlm.nih.gov/articles/PMC12485372/',
 'RCC_count_matrix':'OMIX010399-02:4.9MB,14samples/7patients per paper; controlled-access and unavailable as checked 2026-10-02',
 'RCC_clinical':'OMIX010400 controlled access; no clinical download',
 'RCC_SRA':'PRJNA1223824 PUBLIC,60runs,total177410289249compressed_fastq_bytes; candidate14 B1-B7/SPP1_1-SPP1_7 total39023903366bytes. All same-number pairs disagree in age; group treatment and donor pairing not explicitly annotated. Collection dates2024 differ from manuscript blockade collection2023. No FASTQ downloaded; no invented paired model.',
 'RCC_supplement':'Author TableS1-S9 obtained via publisher Figshare and verifiedMD5; cohort characteristics/antibodies, not individual intervention counts or flow values.',
 'Trehan_source':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-59529-0/MediaObjects/41467_2025_59529_MOESM4_ESM.xlsx',
 'Trehan_source_sha256':sha(P/'downloads/Trehan_SourceData.xlsx'),
 'Trehan_Nanostring':'Only normalized pathway scores(two values/arm) in source workbook, not785gene matrix or unambiguous donor mapping. GitHub rajuvee/TASrev2 tree contains analysisRscript/README, no expression matrix. No derived perturbation gene signature manufactured.',
 'Trehan_discrepancies':['Fig6L/M source14paired rows vs caption15: retain descriptive source values, no newP/FDR.','Fig5E source9weights/group vs caption5: do not pool rows or use tumor-weight endpoint.'],
 'access_actions':'Public metadata and small source files via authorized proxy. No access request sent; no attempt to circumvent controlled access; no new raw sequencing job.',
 'checked':stamp}
(P/'resources/SOURCE_ACCESS_AND_UNIT_AUDIT.json').write_text(json.dumps(source_audit,ensure_ascii=False,indent=2))
# Correct a wording-only counting typo without changing the frozen statistical family.
proto=P/'resources/DOWNSTREAM_PROTOCOL.json';orig=proto.read_bytes();j=json.loads(orig);j['external_spatial']=j['external_spatial'].replace('Eight distance hypotheses G-N/O excluding density F form9 endpoint family G-O','Nine distance hypotheses G-O excluding density F form a9 endpoint family')
if json.loads(orig)!=j:
 (P/'resources/DOWNSTREAM_PROTOCOL_before_wording_correction.json').write_bytes(orig);proto.write_text(json.dumps(j,ensure_ascii=False,indent=2));(P/'qc/PROTOCOL_WORDING_CORRECTION.json').write_text(json.dumps({'only_change':'Eight -> Nine distance hypotheses wording; planned family9, code, thresholds and results unchanged','original_sha256':hashlib.sha256(orig).hexdigest(),'corrected_sha256':sha(proto)},indent=2))
def one(src,ep,model='batch_stage_adjusted'):return e[e.dataset.eq('OBRADOVIC')&e.source.eq(src)&e.endpoint.eq(ep)&e.model.eq(model)].iloc[0]
g=sp[sp.panel.eq('6G')].iloc[0];k=sp[sp.panel.eq('6K')].iloc[0];ex=pe[pe.panel.eq('4I')].iloc[0];rev=one('SPP1','PDCD1_ENTPD1_double_detected_fraction');depth=one('SPP1','PDCD1_ENTPD1_double_detected_fraction','batch_stage_CD8depth_adjusted');tr=one('TYMP','PDCD1_ENTPD1_double_detected_fraction','batch_stage_CD8depth_adjusted')
main=e[e.dataset.eq('OBRADOVIC')].copy()
summary='\n'.join(f'{r.source} → {r.endpoint}；{r.model}；n={r.patients}，rho={r.rho:.6f}，P={r.P:.6f}，FDR={r.FDR:.6f}。' for r in main.itertuples())
text=f'''26210 下游功能与空间证据补强结果
Phase306｜2026年10月2日

结论
本轮补入了可复核的外部SPP1扰动与基质区域空间证据，使“SPP1与免疫屏障有什么关系”有了具体解释：SPP1可以影响CD8的效应状态，基质富集区域也可伴随巨噬与CD8更近的空间关系。现有肾癌单细胞的定向接收端复核没有新增稳健的预期阳性。因此，外部功能支持已加强，本项目“TYMP—SPP1—基质—CD8”完整链仍未建立。

一、肾癌SPP1阻断提供了直接的功能性文献支持
2025年Cancer Immunology Research原始研究在ccRCC新鲜肿瘤离体培养中比较SPP1阻断与同型对照。作者报告：CD8中的IFNγ和GZMB升高，CTLA4、CD39、TOX及PD-1/CD39双阳性比例降低；联合PD-1阻断还提高了CD8增殖和肿瘤细胞凋亡。主要流式面板原文每组n=16；转录组为7位患者、14份处理/对照样本。
这项工作可直接解释SPP1与抗肿瘤免疫受抑的关系。上述为作者实验及作者结论，本轮没有取得逐患者流式源值，未宣称独立复算其处理效应或临床疗效。

二、公开SPP1/CD44源数据支持功能抑制的部分可逆性
Trehan等Nature Communications 2025研究对纯化小鼠CD8设置IgG、重组Spp1、Spp1加抗CD44三组。本轮从官方源xlsx重新提取5个NanoString程序的全部数值。
耗竭评分组均值为IgG {ex.IgG:.3f}、Spp1 {ex.Spp1:.3f}、Spp1+抗CD44 {ex.Spp1_aCD44:.3f}。SPP1刺激增加{ex.stimulation_difference:.3f}评分单位；加入CD44阻断后降低{-ex.CD44_blockade_difference:.3f}，相当于抵消刺激所增加分值的{ex.attenuation_fraction:.1%}。阻断后仍高于对照，因此应写“部分回落”，不能写“完全恢复”。其他JAK-STAT、MAPK、NFκB及TLR程序的均值也呈刺激升高、阻断回落。
官方表每组只有2个评分值，缺少明确供者映射；原文图注同时写n=6及两个重复。这里仅报告源值、均值差和回落幅度，不制造独立生物学n或P/FDR，也不把标准化评分差称表达倍数。
同研究Spp1全身敲除小鼠中，肿瘤抗原TRP2阳性CD8的CD39阳性比例由WT均值53.27%降至34.45%，下降18.82个百分点；TIM3由36.11%降至23.08%，下降13.03个百分点。KO每组6、WT每组7，仅报告效应；原作者P另列。TRP2阴性对照和全部其他面板完整保存。

三、外部患者空间源数据把基质环境与CD8联系起来
同研究的人HCC CODEX源表提供15位患者αSMA富集区与非富集区的配对比较。本轮按患者复算，未把细胞或区域当样本量。
巨噬至总CD8距离：14/15位患者在αSMA富集区更近；患者内距离比中位数{g.median_ratio_positive_over_negative:.3f}，即约缩短{1-g.median_ratio_positive_over_negative:.1%}；配对检验P={g.P:.8f}，FDR={g.FDR:.8f}。
巨噬至CD39阳性CD8距离：15/15位患者更近；距离比中位数{k.median_ratio_positive_over_negative:.3f}，约缩短{1-k.median_ratio_positive_over_negative:.1%}；P={k.P:.8f}，FDR={k.FDR:.8f}。这是校正后显著的外部空间证据。
巨噬密度本身P=0.095668；Treg相关两项距离未显示名义显著。9项空间距离构成一个校正家族；NK相关两项源表14对而图注15，保留描述结果，缺项仅在保守BH计算中占据计划家族位置，不补造P或患者。全部原值和阴性结果随附。

四、本项目定向单细胞复核未复制预期的免疫抑制关联
依据ccRCC阻断实验，在计算前固定3个RNA接收端：PDCD1/ENTPD1共同检出比例（对应PD-1/CD39）、TOX检出比例、IFNG/GZMB平均logCP10K。沿用Phase305巨噬/CD8身份、原始计数和患者，不改变分群及来源门槛。主要模型控制作者批次和pT分期，敏感性再控制CD8测序深度。
Obradovic仍为11位未治疗患者。SPP1与PD-1/CD39共同检出呈反向名义关联：rho={rev.rho:.3f}，P={rev.P:.4f}，FDR={rev.FDR:.4f}；加深度控制后rho={depth.rho:.3f}，P={depth.P:.4f}。TYMP在深度模型中也呈反向名义信号：rho={tr.rho:.3f}，P={tr.P:.4f}，FDR={tr.FDR:.4f}。这些方向不支持“SPP1/TYMP越高，终末障碍表型越多”，不能包装成正面机制证据。
IFNG/GZMB部分效应方向与受抑假设一致，但均未达到原始P<0.05；TOX也无支持性名义信号。HCRN基线只有1位、进展只有5位具备巨噬/CD8配对，方向单列、不作P/FDR，不合并时点或重复患者。
这是一轮既有队列上的文献定向探索，不是新的独立队列验证，也不能抵消Phase305宽泛CD8端点的阴性结果。

五、原始数据通路已核查，当前瓶颈是样本映射
ccRCC干预计数OMIX010399-02约4.9MB，当前受控且页面显示不可下载。公开SRA/ENA的PRJNA1223824有60个测序run；其中B1–B7及SPP1_1–7这14个候选run压缩FASTQ共39.02GB。但按编号对应的7组年龄全部不一致，处理组没有明确字段，记录采集年2024也不同于原文干预队列2023。编号不能替代供者配对。
已保存60个BioSample元数据、14个实验元数据和逐组冲突表，未下载39GB后强行配对。作者补充表S1–S9是队列特征和试剂信息，未补齐个体干预映射。小鼠公开源表只有NanoString程序评分，作者仓库没有完整785基因计数矩阵，不能据此构造“扰动响应基因签名”再投射到肾癌单细胞。

六、方法与统计
外部表格保持原始文件及源单元格坐标；空间按存储行匹配αSMA阳/阴区域，配对t检验，距离家族9项保守BH，密度单一终点单列。小鼠源评分仅作组均值差，敲除小组仅作效应描述，原文P单独标注。
单细胞每患者至少100个既有合格CD8，所有端点基因须在贡献文库中被测量。检出比例以原始计数>0定义，实测零值保留；IFNG/GZMB继承原程序检出门槛。巨噬TYMP/SPP1资格沿用Phase305。两基因程序用患者原始计数求和后log1pCP10K均值。患者秩残差相关，n>=10且剩余自由度>=5时用双侧t近似P；每队列/时点/巨噬来源3端点×2模型构成6项BH家族。原P<0.05为名义信号，FDR<0.05才称校正后显著。

七、讨论与接续重点
本轮最有价值的推进，是把“免疫屏障”落实为可分开验证的两个问题：基质环境中的局部细胞空间关系，以及CD8是否被功能性抑制。外部干预支持SPP1可影响CD8，外部空间支持基质富集区域的细胞接近；空间接近本身不等于T细胞进入减少，也未证明SPP1驱动胶原排列。
小鼠全身Spp1敲除不等于巨噬特异敲除；外部HCC/肝转移结果不能替代ccRCC验证，CD39或TOX转录检出也不等于蛋白、抗原特异性或实测杀伤功能。所有新数据都不构成TYMP干预证据。
本项目上游TYMP—SPP1患者内正关联继续保留。下一步最有信息量的输入是可核实供者配对的ccRCC SPP1干预计数/流式，或含SPP1巨噬、CD8功能标记和基质形态的同组织空间蛋白数据。当前数据下不继续换阈值或拆亚群寻找显著性。未向作者发送消息或访问申请。

八、验收与完整结果
独立逐文库重建端点、正交投影复算相关、手工BH/配对t公式、源表坐标核验共{v['checks']}项通过。所有36个单细胞模型组合、5项刺激/救援评分、8项敲除表型和10项空间端点均保留。报告及统计表位于delivery；权威进度仅保存在sgs10项目目录。无后台计算，旧对象、C4、Phase302 Word和正式文章未覆盖，无新增生物学图件。

Obradovic全部模型
{summary}

来源
ccRCC原文：https://pmc.ncbi.nlm.nih.gov/articles/PMC12485372/
计数访问状态：https://ngdc.cncb.ac.cn/omix/release/OMIX010399
公开原始测序：https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1223824
小鼠扰动/人空间原文：https://www.nature.com/articles/s41467-025-59529-0
源表下载与SHA256详见resources/SOURCE_ACCESS_AND_UNIT_AUDIT.json；作者源表CC BY4.0，本轮重算及衍生统计均已与原始数据区分。
'''
stem='26210_Phase306_下游功能与空间证据';(D/(stem+'.txt')).write_text(text,encoding='utf-8-sig')
body=''.join('<p>'+html.escape(s).replace('\n','<br>')+'</p>' for s in text.split('\n\n'));(D/(stem+'.html')).write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>26210 下游证据补强</title><style>body{max-width:1000px;margin:40px auto;padding:0 25px;font:16px/1.9 "Microsoft YaHei",sans-serif;color:#24313c}p{margin:1.5em 0}table{border-collapse:collapse;font-size:12px}th,td{border:1px solid #bbc7d2;padding:6px}@media print{body{font-size:11pt;margin:0}}</style><body>'+body+'</body></html>',encoding='utf-8')
for folder in ['tables','resources','qc','code']:
 (D/folder).mkdir(exist_ok=True)
 for f in (P/folder).iterdir():
  if f.is_file() and not f.name.startswith('DOWNSTREAM_PROTOCOL_before'):shutil.copy2(f,D/folder/f.name)
(D/'source_data').mkdir(exist_ok=True);shutil.copy2(P/'downloads/Trehan_SourceData.xlsx',D/'source_data/Trehan_SourceData.xlsx')
manifest={str(f.relative_to(D)):sha(f) for f in D.rglob('*') if f.is_file() and f.name!='SHA256.json'};(D/'SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2));archive=P/(stem+'.zip')
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in D.rglob('*'):
  if f.is_file():z.write(f,f.relative_to(D))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
(P/'qc/DELIVERY_VALIDATED.json').write_text(json.dumps({'status':'PASS','zip_sha256':sha(archive),'zip_bytes':archive.stat().st_size,'files':len(manifest),'checks':v['checks']},indent=2))
message=f'''## Phase306终态（2026-10-02）：下游外部扰动/空间补强及定向单细胞复核完成
外部源数据新增：Trehan官方source xlsx通过SHA/CRC。纯CD8 Spp1刺激耗竭评分均值-1.842→1.275，加antiCD44降至0.567，抵消增量22.7%（每组仅2源值、供者不明，不补P/n）；KO6/WT7的TRP2+CD8 CD39和TIM3分别下降18.82/13.03个百分点，保留全部对照及作者P。外部HCC空间15患者：巨噬—CD8距离14/15更近，比例中位0.531，P{g.P:.8f}/FDR{g.FDR:.8f}；CD39+CD8距离15/15更近，比0.519，P{k.P:.8f}/FDR{k.FDR:.8f}。FDR距离9项，NK两项源14与图注15冲突描述单列；密度/Treg阴性保留。外部结果不替代RCC或TYMP因果，不称物理排斥。
ccRCC CIR-24-1146阻断文献支持CD8功能恢复（原文流式各组n16/转录组7供者14库），本轮未取得个体流式源数据。OMIX010399受控/不可下载；SRA60run可读，候选14run39.02GB的B1-7/SPP1_1-7全部同号年龄冲突，处理/配对未明确，采集年份与论文冲突，未硬配对/未下载FASTQ。60份样本+14份实验元数据、作者S1-S9和逐组冲突表保存；未联系作者。
本项目文献固定接收端（PDCD1/ENTPD1双检出、TOX检出、IFNG/GZMB）完成全部36模型。Obradovic11患者SPP1—PD1/CD39反向名义rho-0.727451/P0.026336/FDR0.158014，深度调后rho-0.327832/P0.427937；TYMP深度模型反向rho-0.710583/P0.048212/FDR0.289270。无预期方向原P<.05；不把反向信号包装成支持。HCRN基线1/进展5仅描述，不混时点。旧Phase305四项CD8阴性与上游阳性均保持。
权威目录phase306_downstream_perturbation_evidence_20261002_v1。冻结方案/完整表/源坐标/独立逐文库验收1161项PASS；analysis退出0。报告delivery/{stem}.txt和html，轻量包{stem}.zip。没有新增生物学图，旧C4/对象/正式文章/Phase302 Word保持。本轮补强分析完成，完整屏障机制尚未成立；后续需有效配对干预或同组织空间蛋白，停止同队列换阈值。无计算后台。

'''
backup=P/'resources'/('root_before_delivery_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'));backup.mkdir()
for name in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=R/name;old=f.read_bytes();(backup/name).write_bytes(old);new=message.encode()+old;assert f.read_bytes()==old;tmp=f.with_name(f.name+'.phase306.tmp');tmp.write_bytes(new);os.replace(tmp,f)
f=R/'RESULTS_MANIFEST.tsv';old=f.read_bytes();(backup/f.name).write_bytes(old);content=old.decode('utf-8-sig');fields=next(csv.reader(io.StringIO(content),delimiter='\t'));buf=io.StringIO();writer=csv.DictWriter(buf,fieldnames=fields,delimiter='\t',lineterminator='\n')
for ident,file,typ in [('phase306_report',D/(stem+'.txt'),'report'),('phase306_statistics',P/'tables/RECIPIENT_ASSOCIATIONS_ALL.tsv','patient_statistics'),('phase306_validation',P/'qc/FINAL_VALIDATION.json','validation'),('phase306_delivery',archive,'archive')]:
 writer.writerow({'result_id':ident,'stage':'Phase306','artifact_type':typ,'absolute_path':str(file),'generated_at':stamp,'validation_status':'PASS','validation_evidence':str(P/'qc/FINAL_VALIDATION.json'),'notes':'External functional/spatial support strengthened; no supportive new RCC recipient association; complete negatives retained'})
assert f.read_bytes()==old;tmp=f.with_name(f.name+'.phase306.tmp');tmp.write_bytes(old+(b'' if old.endswith(b'\n') else b'\n')+buf.getvalue().encode());os.replace(tmp,f)
(P/'STATUS.json').write_text(json.dumps({'status':'DOWNSTREAM_FIRST_ROUND_COMPLETE_VALIDATED','updated':stamp,'validation_checks':v['checks'],'active_jobs':[],'complete_mechanism_proven':False,'zip_sha256':sha(archive)},indent=2));print('DELIVERY',archive.name,archive.stat().st_size,sha(archive),len(manifest),flush=True)
