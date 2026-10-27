from pathlib import Path
import json,hashlib,datetime,csv,io,os,shutil,subprocess
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
sources=[
R/'phase280_TYMP_source_niche_function_20260928_v1/delivery/TYMP三方向深化研究报告.md',
R/'phase291_spatial_null_TYMP_matched_regions_20260930_v1/delivery/补充方法检验与完整结果.md',
R/'phase293_independent_barrier_followthrough_20260930_v1/delivery/Phase293_独立证据补强终态报告.md',
R/'phase295_literature_guided_sequential_trials_20260930_v1/delivery/Phase295_文献方法顺序尝试终态报告.md',
R/'phase296_TYMP_SPP1_barrier_segmented_chain_20260930_v1/delivery/Phase296_主线分段与参与效应检验.md']
with (P/'qc/memory_gate.log').open('a') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,check=True)
before={str(p):sha(p) for p in sources}
report="""# 当前线索整合：TYMP—SPP1相关巨噬—间质屏障
本轮仅综合既有证据，不新增统计、筛样、评分、图件或正文修改。服务器为唯一权威。用户要求：先把线索串起来，看看能组成什么。

## 目前可以形成的核心判断
现有生信最适合支撑“TYMP与SPP1相关髓系状态、间质基质重塑存在互补关联”的机制框架。研究主问题继续是肿瘤TYMP参与SPP1相关屏障；数据中最强的关联来自髓系内在及髓系—间质层面，不应将其偷换为肿瘤来源效应。
可以组织为一个中心明确、两端待连接的模型。不是没有故事，也不是已经证明串联因果。证据的互补性来自不同测量层级，不能把分散结果当作同一患者中完整链条的直接证明。

## 主线的层级结构
肿瘤TYMP与核苷/甲硫氨酸-SAM-多胺代谢背景：研究起点。
    待连接：肿瘤来源变化怎样影响SPP1相关巨噬状态。
SPP1相关巨噬状态：TYMP—SPP1巨噬胞内关联为支撑；C4为冻结研究载体，不要求C4专属。
    过程联系：髓系TYMP—间质胶原成熟关联；SPP1状态—间质活化的有限空间对应。
基质重塑相关生态位：胶原成熟、文献ECM基因方向对应；不是已测得的胶原排列。
    待连接：该生态位是否限制T/CD8进入。
SPP1相关免疫屏障形成/维持：保持为整体机制假说及研究目标。

上图式文字中“过程联系”也是关联，未赋予因果方向。客户肿瘤干预代谢方向只作为已提供实验背景，不对客户原始实验重新分析，也不补造未提供的SPP1、CAF或T细胞实验结论。

## 证据对照
|环节|已得到的线索|可使用的定位|
|---|---|---|
|TYMP与SPP1髓系状态|整体巨噬29患者24正，中位患者内rho0.049133，P0.000546/FDR0.002184；C4 15患者13正，rho0.134872，P0.007385/FDR0.029541|患者间方向一致、效应较小的胞内表达关联；不是肿瘤分泌TYMP刺激SPP1的证据|
|髓系TYMP与间质胶原成熟|原Phase292 13患者rho0.809604，P0.00895/FDR0.02685；独立GSE264586基线4患者rho1.0|过程端最强关联；13例原接收池为广义壁细胞，4例严格成纤维候选来自转移部位，仅描述|
|SPP1与间质过程|共同13患者SPP1—胶原rho0.666619，P0.13583/FDR0.339575，留一患者均正；独立4例rho0.8|正向描述线索，未达主队列名义显著，不能写成显著中介|
|SPP1状态与空间间质活化|GSE175540 400μm血管调整：15患者12正，中位rho0.224690，P0.0353394/FDR0.106018|名义显著、尺度/模型限定的空间协作线索；200μm主模型不支持，不能写成稳定跨尺度结果|
|间质反应与ECM参考|Tomko原文27基因21个方向对应；排除定义基因15/20；原文P<0.05并排除定义基因13/18|同一7患者的外部文献特征对应，非新增独立人群，不等于实测纤维排列|
|候选输出|NicheNet TGFB1全排名1/797，SPP1巨噬表达配体排名162/190|协作过程的候选解释；未证明SPP1状态特异TGFB1输出，不能将TGFB1改为替代主线|
|肿瘤TYMP到巨噬SPP1|共同13患者rho0.013862，P0.98529/FDR0.98529|当前数据不支持患者级上游联系；髓系阳性不能替代|
|肿瘤侧TYMP到空间联合生态位|TJRCC 200μm 9例8正，中位rho0.123034；严格空间零模型0/9原P<0.05|探索性描述；400μm反向/不稳，不能补上上游因果缺口|
|生态位到T/CD8进入|多队列/严格空间模型未获得稳定进入减少；GLMF两例总T方向可描述，CD8不同向|功能屏障环节尚未接上；静态内外分布也不是实测迁移|

各行n来自不同测量和资格人群，不能相加为总验证人数。Phase296重新统一检验家族/随机种子所得P/FDR不覆盖Phase292原统计；同一数据重复分析不是复制。

## 可以怎样安排结果叙事
第一部分：肿瘤TYMP及其代谢背景提出微环境参与问题。使用既有正式前两部分及用户已给实验方向，不在本次综合中重写或追加未经核查机制。
第二部分：TYMP与SPP1相关巨噬状态相联系。优先报告患者内小效应但方向一致的关联，再报告C4作为研究载体；不把来源从肿瘤悄然改成巨噬。
第三部分：TYMP相关髓系变化与间质基质处理相对应。13患者胶原成熟关联为主体，独立4患者方向为补充，说明接受细胞身份范围。
第四部分：SPP1相关状态与间质空间协作提供生态位线索。保留400μm名义信号与GLMF案例，集中说明200μm、严格背景模型和T/CD8终点未稳定支持。
讨论合并解释：这些结果支持在TYMP—SPP1相关髓系—基质重塑框架内继续理解屏障问题。正向关联与文献程序对应可支持“参与线索”，不能升级为TYMP经SPP1诱导CAF屏障或阻止T细胞进入已经证实。

建议生信结果总标题：TYMP与SPP1相关髓系状态及间质基质重塑的联系。
建议整体机制假说标题：肿瘤TYMP参与SPP1相关髓系—间质屏障形成的候选机制。
“假说标题”和“数据结果标题”必须分开。主线保持，证据强度不被标题放大。

## 可用于综合讨论的段落草稿
我们的分析显示，巨噬细胞内TYMP与SPP1表达在多数患者中呈同向联系，髓系TYMP同时与间质胶原成熟程序相关，且在独立小样本中观察到一致方向。部分空间分析进一步提示SPP1相关巨噬状态与间质活化相联系，接收端表达变化亦与已发表ECM重塑特征具有方向对应。这些结果共同提供了TYMP与SPP1相关髓系状态及基质重塑相联系的生信依据，并形成肿瘤TYMP参与相关屏障过程的候选机制框架。现有跨细胞关联尚未连接肿瘤TYMP与巨噬SPP1，T细胞进入减少亦未获得稳定支持，因此上述证据应解释为分层的参与线索，而非完整串联因果链。

## 改变解释的结果集中保留
控制SPP1后髓系TYMP—胶原成熟仍名义显著（rho0.742859/P0.02367/FDR0.09468）；这不证明SPP1介导，也不证明SPP1无作用。
TYMP×SPP1交互：肿瘤beta0.476258/P0.35851；髓系beta−0.402799/P0.05563，均不支持正向协同。不能将P0.05563称名义显著。
公共外源TYMP刺激GSE171542为非肿瘤破骨分化背景，SPP1/C4程序反向，不能当直接正向验证；情境不同也不应外推为肿瘤机制不存在。
27个差异/趋势ECM蛋白与24个排列相关蛋白不是同一集合。表达/方向对应不等于蛋白变化、纤维几何或功能屏障。
纯CAF身份、空间靶点缺失、分类器不相容及受限资料分别为测量资格问题，不写成生物学阴性。
文献先验网络排名不是实验处理效应。TGFB1候选结果不另立主线，也不填补TYMP→SPP1或T进入的缺口。

## 下一增量的明确位置
最值得补的是主线两端：肿瘤来源TYMP如何接到SPP1相关巨噬—间质过程，以及这个过程是否对应免疫进入受限。现有中段材料已经足以整理参与线索，无需继续在同13患者中轮换阈值和端点凑阳性。新增资料必须与缺口匹配，不以更强的髓系TYMP关联替代肿瘤来源证据。本轮为综合，不自行启动新增数据下载或统计。
"""
path=P/'delivery/当前线索整合与主线叙事.md';path.write_text(report)
matrix=[
('胞内TYMP-SPP1','患者内表达关联','29例24正;rho0.049133;P0.000546;FDR0.002184','支持小效应胞内联系','Phase280'),
('髓系TYMP-胶原成熟','跨细胞患者关联','13例rho0.809604;P0.00895;FDR0.02685','原广义壁细胞接收池','Phase292'),
('独立病例髓系TYMP-胶原','外部描述','4例rho1.0','转移部位严格成纤维候选;无P','Phase293'),
('SPP1状态-空间间质','模型限定患者关联','15例12正;P0.0353394;FDR0.106018','400μm血管模型;200μm不支持','Phase291'),
('肿瘤TYMP-SPP1','同患者关联','13例rho0.013862;P0.98529','当前数据不支持','Phase296'),
('SPP1-胶原','同患者关联','13例rho0.666619;P0.13583','仅正向描述','Phase296'),
('屏障-T进入','空间功能代理','未获稳定支持','不能写成功能屏障证实','Phase291/293')]
with (P/'delivery/EVIDENCE_MATRIX.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['link','measurement','evidence','boundary','source']);w.writerows(matrix)
assert all(sha(Path(p))==h for p,h in before.items())
now=datetime.datetime.now().astimezone().isoformat()
entry=f"""## Phase297（2026-09-30）：现有主线证据整合完成
用户要求把当前线索串起来。本轮只综合，无新统计/图/筛样/正文修改。
综合稿：{P.name}/delivery/当前线索整合与主线叙事.md；证据矩阵同目录EVIDENCE_MATRIX.tsv。
现有核心为TYMP与SPP1相关髓系状态、间质基质重塑的互补关联。整体研究目标仍为肿瘤TYMP参与SPP1相关屏障；髓系内在阳性不替代肿瘤来源。
中段：胞内TYMP-SPP1方向一致；髓系TYMP-胶原成熟13例关联及独立4例描述；SPP1状态-间质空间协作仅400μm血管模型名义支持。两端：肿瘤TYMP到SPP1及T/CD8进入受限未接上。模型是候选机制框架，非完整串联因果。
保留SPP1直接配体低排名、肿瘤关联/交互阴性、空间主模型阴性及外源TYMP反向。TGFB1仅候选输出，不另起主线。已提供讨论段落与结果组织建议，未改正式文章。

"""
names=['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']
old={n:(R/n).read_bytes() for n in names}
for n,b in old.items():
 with (P/'root_record_backups'/n).open('xb') as f:f.write(b)
new={}
for n,b in old.items():
 if n.endswith('.md'):new[n]=entry.encode()+b
 else:
  fields=b.decode('utf-8-sig').splitlines()[0].split('\t');out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
  for rid,typ,p in [('phase297_synthesis','report',path),('phase297_matrix','table',P/'delivery/EVIDENCE_MATRIX.tsv')]:
   w.writerow(dict(result_id=rid,stage='Phase297',artifact_type=typ,absolute_path=str(p),generated_at=now,validation_status='SYNTHESIS_SOURCE_CHECKED',validation_evidence=str(P/'qc/COMPLETION.json'),notes='No new statistics;SHA256='+sha(p)))
  new[n]=b+(b'' if b.endswith(b'\n') else b'\n')+out.getvalue().encode()
for n in names:assert (R/n).read_bytes()==old[n]
for n,b in new.items():
 temp=R/(n+'.phase297.tmp')
 with temp.open('xb') as f:f.write(b)
 shutil.copymode(R/n,temp);os.replace(temp,R/n)
assert all((R/n).read_bytes()==new[n] for n in names)
(P/'qc/COMPLETION.json').write_text(json.dumps(dict(status='COMPLETE',at=now,new_statistics=False,sources=before,source_files_unchanged=True,root_records_updated=names,outputs={str(path):sha(path),str(P/'delivery/EVIDENCE_MATRIX.tsv'):sha(P/'delivery/EVIDENCE_MATRIX.tsv')}),ensure_ascii=False,indent=2))
print('PHASE297_SYNTHESIS_COMPLETE')
