from pathlib import Path
import pandas as pd,numpy as np,json,datetime,hashlib,csv,io,os,shutil,subprocess
P=Path(__file__).resolve().parents[1];R=P.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
with (P/'logs/memory_gate.log').open('a') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,check=True)
a=pd.read_csv(P/'tables/SEGMENTED_ASSOCIATIONS.tsv',sep='\t');i=pd.read_csv(P/'tables/INTERACTION_EXPLORATORY.tsv',sep='\t');loo=pd.read_csv(P/'tables/LEAVE_OUT_EFFECTS.tsv',sep='\t');iloo=pd.read_csv(P/'tables/INTERACTION_LEAVE_OUT.tsv',sep='\t')
def row(e):
 z=a[a.edge_id.eq(e)].iloc[0];return f"rho={z.rho:.4f}，P={z.P:.5f}，FDR={z.FDR:.5f}"
lr=[]
for edge in a.edge_id:
 z=loo[loo.edge_id.eq(edge)&loo.omission_type.eq('patient')].dropna(subset=['rho'])
 lr.append(dict(edge=edge,evaluable=len(z),positive=int(z.rho.gt(0).sum()),negative=int(z.rho.lt(0).sum()),min_rho=z.rho.min(),max_rho=z.rho.max()))
pd.DataFrame(lr).to_csv(P/'tables/PATIENT_OMISSION_SUMMARY.tsv',sep='\t',index=False)
tbl='|检验|rho|原始P|FDR|\n|---|---:|---:|---:|\n'
for e,label in [('E01','肿瘤TYMP—巨噬SPP1'),('E03','巨噬SPP1—胶原成熟'),('E04','肿瘤TYMP—胶原成熟'),('E05','髓系TYMP—胶原成熟'),('E08','髓系TYMP—胶原，控制SPP1'),('E09','SPP1—胶原，控制髓系TYMP')]:
 z=a[a.edge_id.eq(e)].iloc[0];tbl+=f'|{label}|{z.rho:.4f}|{z.P:.5f}|{z.FDR:.5f}|\n'
report=f"""# Phase296：TYMP参与SPP1相关屏障的主线检验
用户要求不偏离主线继续推进。本轮坚持肿瘤TYMP为主要来源，SPP1相关巨噬—间质过程为第二段；TGFB1未成为替代主线。

## 同患者分段检验尚未连接肿瘤TYMP与SPP1—间质过程
在冻结的13位同患者同文库对象中，研究调整后肿瘤TYMP—巨噬SPP1接近零，SPP1—胶原成熟方向为正但原始P未达0.05。当前数据不支持这两段在患者层面形成统计支持的串联。
{tbl}
胶原成熟为既定间质转录过程端点，不是T细胞进入或实测免疫屏障。以上既有同定义边的效应与Phase292逐一一致；本轮统一到同13患者及新检验家族，不当新增独立复制。Monte Carlo种子和BH家族变化导致P/FDR不同，历史统计不覆盖。

## 髓系TYMP—间质关联控制SPP1后仍保留
控制巨噬SPP1后，髓系TYMP—胶原成熟为{row('E08')}，属于名义显著，未通过该家族FDR。再控制接收深度为{row('E18')}，仍为名义显著。不能将残余关联解释为“SPP1无作用”，也不能将它改写成SPP1介导已经成立。
SPP1—胶原控制肿瘤TYMP为{row('E07')}；控制髓系TYMP为{row('E09')}，均不达原P<0.05。固定C4比例—胶原为{row('E12')}；C4只是补充状态代理，不以该阴性否定全部SPP1相关状态。

## 参与效应的交互检验未支持正向协同
用户主张的参与关系不要求TYMP先升高SPP1 RNA，因此在看到线性分段结果后单独冻结探索性“TYMP×SPP1”分析；该序贯设计明确记录，未冒充最初预设。
肿瘤来源交互系数{i.iloc[0].interaction_beta:.4f}，P={i.iloc[0].P:.5f}，FDR={i.iloc[0].FDR:.5f}；髓系来源{i.iloc[1].interaction_beta:.4f}，P={i.iloc[1].P:.5f}，FDR={i.iloc[1].FDR:.5f}。当前数据不支持正向协同。髓系交互方向为负，P≥0.05，不称名义显著，也不靠“接近显著”纳入支持主线。
条件斜率、留一患者/研究全部结果保留在INTERACTION_EXPLORATORY.tsv与INTERACTION_LEAVE_OUT.tsv；交互是秩尺度模型效应，不是生物处理效应或因果协同。

## 方法、资格与完整性
以Phase292已冻结eligible13患者为共同分析集；无新筛样。研究组成见STUDY_COMPOSITION.tsv，9/13来自GSE207493。来源表达、C4比例和胶原端点沿用已核验表，合并严格按patient/study/library一对一。
分段模型为研究调整秩偏相关；条件项和接收深度秩转换，99999次研究内残差置换，种子296。有连续控制项的检验标为近似残差置换。原P<0.05为候选硬标准；BH家族core5、conditional4、state3、depth8各自保留，完整原P/FDR见SEGMENTED_ASSOCIATIONS.tsv。
交互使用标准化中心秩TYMP、SPP1及乘积，研究固定效应；零模型为研究和两个主效应，Freedman–Lane研究内近似残差置换99999次、种子29601，两个来源为家族2。目标是交互项的P，未使用协变量最小P。
n<10或残差df<5不计算推断P。留一患者/研究只报告效应，不追加检验。研究集中、小样本和横断面测量限制在解释中保留；不跨队列拼成同患者中介，不计算因果间接效应比例。
20个偏相关用独立最小二乘残差复算；3条同定义旧边精确核对；家族大小检查；2个交互系数用Frisch–Waugh残差公式复算。运行与交互脚本退出0。来源对象、旧C4、空间边界和历史结果只读，新图0。

## 主线保留与下一证据缺口
已经确立的巨噬胞内TYMP—SPP1关联，以及髓系TYMP—间质胶原成熟关联继续保留。它们不能替代客户肿瘤TYMP干预背景。本轮说明：继续在这13例上换线性串联、C4比例或交互写法，尚不能得到肿瘤TYMP参与SPP1屏障的新支持。
后续主问题仍是肿瘤TYMP与SPP1相关巨噬—间质屏障的对应证据，优先来自相容的独立患者/空间或公共干预资料；需同时包含来源和接收测量。TGFB1只作已获得的候选输出，不另立机制主线。公开数据不足不要求客户新增实验作为前提，也不将阴性条件分析写成机制不存在。

文件：PROTOCOL_FROZEN.md、INTERACTION_ADDENDUM_BEFORE_TESTING.md；tables/全部结果和患者表；qc/输入哈希及独立核验；logs/内存门禁、完整日志及退出状态。
"""
(P/'delivery/Phase296_主线分段与参与效应检验.md').write_text(report)
now=datetime.datetime.now().astimezone().isoformat()
entry=f"""## Phase296（2026-09-30）：按用户要求回归TYMP—SPP1屏障主线
报告：{P.name}/delivery/Phase296_主线分段与参与效应检验.md。
共同13患者：肿瘤TYMP—巨噬SPP1 rho0.0139/P0.98529；SPP1—胶原rho0.6666/P0.13583，当前数据不支持同患者串联。髓系TYMP—胶原控制SPP1后rho0.7429/P0.02367/FDR0.09468，名义显著；再控制深度P0.04162。保留胞内关联与过程关联，但不替代肿瘤来源主线，不声称SPP1中介。
另在结果后独立冻结探索性交互：肿瘤beta0.4763/P0.35851；髓系beta-0.4028/P0.05563，均无名义支持，不称正向协同。完整原P/FDR、20边、2交互、留一效应保存；24个主核验及2交互复算通过。
不将TGFB1提为新主线，不以C4特异性/比例作必要条件。旧对象/C4/边界/文章未改，无新图；下一增量需要相容独立来源—接收测量，不再重复同13例换端点/阈值。当前本轮计算完成，无后台。

"""
files=['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv'];old={n:(R/n).read_bytes() for n in files}
for n,b in old.items():
 with (P/'root_record_backups'/n).open('xb') as f:f.write(b)
arts=[('phase296_report','report',P/'delivery/Phase296_主线分段与参与效应检验.md'),('phase296_segmented','table',P/'tables/SEGMENTED_ASSOCIATIONS.tsv'),('phase296_interaction','table',P/'tables/INTERACTION_EXPLORATORY.tsv'),('phase296_protocol','protocol',P/'PROTOCOL_FROZEN.md'),('phase296_input','input',P/'tables/FROZEN_PATIENT_INPUTS.tsv')]
changes={}
for n,b in old.items():
 if n.endswith('.md'):changes[n]=entry.encode()+b
 else:
  fields=b.decode('utf-8-sig').splitlines()[0].split('\t');out=io.StringIO();w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
  for rid,typ,path in arts:w.writerow(dict(result_id=rid,stage='Phase296',artifact_type=typ,absolute_path=str(path),generated_at=now,validation_status='VALIDATED_EXPLORATORY',validation_evidence=str(P/'qc/FINAL_VALIDATION.json'),notes='SHA256='+sha(path)))
  changes[n]=b+(b'' if b.endswith(b'\n') else b'\n')+out.getvalue().encode()
for n in files:assert (R/n).read_bytes()==old[n]
for n,new in changes.items():
 tmp=R/(n+'.phase296.tmp')
 with tmp.open('xb') as f:f.write(new)
 shutil.copymode(R/n,tmp);os.replace(tmp,R/n)
checks=dict(run_exit0=(P/'logs/run.exit').read_text().strip()=='0',interaction_exit0=(P/'logs/interaction.exit').read_text().strip()=='0',edges20=len(a)==20,interactions2=len(i)==2,primary_checks24=len(pd.read_csv(P/'qc/INDEPENDENT_CHECKS.tsv',sep='\t'))==24,interaction_checks2=len(pd.read_csv(P/'qc/INTERACTION_INDEPENDENT_CHECKS.tsv',sep='\t'))==2,roots_written=all((R/n).read_bytes()==changes[n] for n in files),protocol_unchanged=sha(P/'PROTOCOL_FROZEN.md')==json.loads((P/'qc/INPUT_FREEZE.json').read_text())['protocol_sha256'])
assert all(checks.values())
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(dict(status='PASS',checks={k:bool(v) for k,v in checks.items()},at=now,artifacts=[dict(path=str(p),sha256=sha(p)) for _,_,p in arts]),indent=2))
(P/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE',scope='mainline segmented and exploratory interaction tests',patient_n=13,background_jobs=0,report=str(P/'delivery/Phase296_主线分段与参与效应检验.md')),indent=2))
print('PHASE296_COMPLETE final checks',len(checks))
print(pd.DataFrame(lr).query("edge in ['E01','E03','E08']").to_string(index=False))
