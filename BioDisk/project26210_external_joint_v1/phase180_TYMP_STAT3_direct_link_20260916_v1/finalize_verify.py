from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,shutil,csv
from scipy.stats import rankdata,wilcoxon
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase180_TYMP_STAT3_direct_link_20260916_v1";v=d/"run_v2"
checks=[]
def ck(n,ok,e=0):checks.append(dict(check=n,passed=bool(ok),error=float(e)))
pd0=pd.read_csv(v/"MALIGNANT_PATIENT_EFFECTS.tsv",sep="\t");s=pd.read_csv(v/"MALIGNANT_SUMMARY.tsv",sep="\t")
cache=np.load(v/"MALIGNANT_VALIDATION_CACHE.npz")
models={"raw":[],"depth":[2],"depth_genes":[2,3],"depth_cycle_hypoxia":[2,4,5],"depth_genes_cycle_hypoxia":[2,3,4,5]}
for key in cache.files:
 res,pat=key.split("__",1);a=cache[key]
 for adj,cols in models.items():
  C=np.column_stack([np.ones(len(a))]+[rankdata(a[:,i]) for i in cols]);Q,_=np.linalg.qr(C)
  x=rankdata(a[:,0]);y=rankdata(a[:,1]);x-=Q@(Q.T@x);y-=Q@(Q.T@y);rho=np.dot(x,y)/np.sqrt(np.dot(x,x)*np.dot(y,y))
  rr=pd0[(pd0.resource==res)&(pd0.scoring=="zmean")&(pd0.patient==pat)&(pd0.adjustment==adj)].rho.iloc[0]
  ck("cached_QR_"+key+adj,abs(rho-rr)<1e-10,abs(rho-rr))
for i,row in s.iterrows():
 a=pd0[(pd0.resource==row.resource)&(pd0.scoring==row.scoring)&(pd0.adjustment==row.adjustment)].rho.to_numpy()
 ck("summary_"+str(i),len(a)==row.n and np.isclose(np.median(a),row.median_rho) and sum(a>0)==row.positive and np.isclose(wilcoxon(a,method="approx",correction=True).pvalue,row.P))
b=pd.read_csv(v/"SAME_LIBRARY_PATIENT_VALUES.tsv",sep="\t");br=pd.read_csv(v/"SAME_LIBRARY_TYMP_MACRO_STAT3.tsv",sep="\t");group=[np.asarray(idx) for idx in b.groupby("study").groups.values()]
for _,row in br[br.adjustment=="study"].iterrows():
 x=rankdata(b.TYMP_pseudobulk);y=rankdata(b[row.resource+"_"+row.scoring])
 for idx in group:x[idx]-=x[idx].mean();y[idx]-=y[idx].mean()
 obs=np.dot(x,y)/np.sqrt(np.dot(x,x)*np.dot(y,y));rng=np.random.default_rng(262180);hits=0
 for j in range(9999):
  xp=x.copy()
  for idx in group:xp[idx]=rng.permutation(xp[idx])
  rr=np.dot(xp,y)/np.sqrt(np.dot(xp,xp)*np.dot(y,y));hits+=abs(rr)>=abs(obs)-1e-14
 pv=(hits+1)/10000;ck("study_perm_"+row.resource+row.scoring,abs(obs-row.rho)<1e-10 and pv==row.P)
for item in json.loads((v/"INPUT_MANIFEST.json").read_text()):
 h=hashlib.sha256()
 with Path(item["path"]).open("rb") as f:
  for block in iter(lambda:f.read(8*1024*1024),b""):h.update(block)
 ck("hash_"+item["path"],h.hexdigest()==item["sha256"])
pd.DataFrame(checks).to_csv(v/"INDEPENDENT_VALIDATION.tsv",sep="\t",index=False);assert all(x["passed"] for x in checks)
orig=pd.read_csv(v/"VALIDATION.tsv",sep="\t");assert orig["pass"].all()
total=len(checks)+len(orig)
bm=pd.read_csv(v/"BULK_TYMP_STAT3_META.tsv",sep="\t");bl=pd.read_csv(v/"BULK_TYMP_STAT3_LEAVE_COHORT.tsv",sep="\t");ml=pd.read_csv(v/"MALIGNANT_LEAVE_STUDY.tsv",sep="\t")
def tab(f):
 def fmt(x):return f"{x:.6g}" if isinstance(x,(float,np.floating)) else str(x)
 return "| "+" | ".join(f.columns)+" |\n|"+"|".join(["---"]*len(f.columns))+"|\n"+"\n".join("| "+" | ".join(fmt(x) for x in z)+" |" for z in f.itertuples(index=False,name=None))
pri=s[(s.resource=="collectri")&(s.scoring=="zmean")]
report="""# TYMP—STAT3活性代理补充验证
最终统计：run_v2；此前MTA—STAT3、STAT3—SPP1及完整短路径保留于RETAINED_PHASE179，内容与Phase179最终表逐字节一致。

## 结论
TYMP与STAT3相关转录状态的正向联系在恶性细胞内及公共bulk队列得到支持；同文库“肿瘤TYMP—巨噬STAT3”的跨区室联系仍弱。因此新增的是TYMP—STAT3关联证据，不是绕开MTA后已经打通了肿瘤到巨噬SPP1的整链。

## 一、恶性细胞内16/16患者方向一致
冻结11950恶性细胞基础，16患者符合预设TYMP检出和细胞数资格。主CollecTRI评分与归一化TYMP表达在深度调整后16/16正向，中位rho=0.14616，P=0.0004824，患者bootstrap中位数区间[0.08517,0.16311]。
"""+tab(pri[["adjustment","n","median_rho","positive","negative","P","CI_low","CI_high"]])+"""
增加检出基因数、周期及缺氧后的联合敏感性仍16/16正向，中位rho=0.11749，P=0.0004824。该结果支持肿瘤细胞内的状态联系，不是TYMP蛋白直接调控STAT3的证据。

## 二、四个公共队列TYMP—STAT3关联均正向
共255患者（97/52/50/56），主模型控制MTAP及髓系组成，合并rho=0.31856，95%区间[0.18305,0.44220]，P=8.0419e-06，四队列方向均正，I²=20.24%。
额外控制MTA后rho=0.34635、P=5.5302e-07；增加缺氧后rho=0.21483、P=0.0009105，均四队列正向。这说明当前测得的TYMP—代理关联在加入这些协变量后保留，不能推导为MTA在生物学上不参与。

主评分全部模型：
"""+tab(bm[(bm.resource=="collectri")&(bm.scoring=="zmean")])+"""
DoRothEA及秩均值对主模型方向提供支持；DoRothEA加缺氧后的P=0.0812（标准化均值）或0.1323（秩均值），不写成所有资源/调整都显著。

## 三、同文库肿瘤TYMP—巨噬STAT3尚未得到稳定支持
沿冻结同16患者、同肿瘤文库，将成熟巨噬原始计数聚合后评分。主CollecTRI评分raw rho=0.26176、P=0.32741；研究调整rho=0.11186、研究内置换P=0.7003。增加C4比例的敏感性rho=0.28872、P=0.3400。
"""+tab(br)+"""
DoRothEA未调整结果名义显著（rho=0.56765、P=0.02181），但研究调整后P=0.146；完整保留这一探索性信号，不替换预设主评分和主模型。

## 四、此前两段与短路径继续保留
|结果|原有结论|本次处理|
|---|---|---|
|整体巨噬STAT3代理—SPP1|31患者29正，rho中位0.08087，P=4.3194e-06|保留，不重写为恶性细胞结果|
|实测MTA—STAT3代理|主合并rho=-0.12332，P=0.10520，原FDR=0.15780|保留，仍不支持稳定正向|
|MTA—STAT3代理—SPP1|路径乘积区间跨零，P=0.76715|保留，不以新增TYMP关联替代|

当前逻辑可以写为：TYMP相关肿瘤状态与STAT3相关转录程序存在联系，STAT3相关程序也与整体巨噬SPP1表达相伴；MTA继续作为代谢候选与受体功能依赖假说研究。必须同时标明两种STAT3证据来自不同细胞区室，尚缺稳定的跨区室关联，不能将两段直接拼成已确认的TYMP→巨噬STAT3→SPP1中介。

## Methods
延用Phase179冻结CollecTRI205和DoRothEA62个有符号靶基因，不包含TYMP/SPP1/STAT3自身。评分为标准化log表达的有符号均值，秩均值为敏感性，固定12位小数后排名。STAT3指转录活性代理，不是STAT3 RNA或磷酸化实测。
恶性细胞TYMP采用log1p(CP10K)，患者≥50细胞且TYMP≥10细胞检出，患者内评分和秩相关，深度主调整。周期/缺氧使用预存Hallmark资源并剔除TYMP/SPP1/STAT3；完整检出基因数敏感性保留。患者层双侧Wilcoxon近似加连续性校正，5000次患者bootstrap中位效应区间。
跨区室分析在同16文库分别聚合巨噬原始计数，转换log2(CPM+1)，在患者间逐基因标准化后评分；不使用患者内z评分的均值作跨患者活性。按研究固定效应残差计算相关，主分析9999次研究内置换；C4协变量敏感性采用研究内残差置换后重新残差化。
公共队列直接复用Phase179最终run_v2评分和患者，TYMP与评分秩相关，调整MTAP/髓系代理；再加MTA或缺氧作敏感性。Fisher-z随机效应合并及留一队列描述。每个数据层仅一个预设主要关联，三个不同统计层面不混池FDR；资源、评分和协变量敏感性不升级为主检验。
本轮未进行新中介分析、重聚类或绘图。客户实验未参与统计。

## 数值验证与版本口径
最初诊断版本在恶性细胞侧使用原始TYMP计数；最终run_v2统一为前期采用的log1p(CP10K)表达，原计数结果保留但不作正文主结果。极少细胞患者的未定义评分在资格阶段排除，不构成推断样本。
原始计数聚合复核16文库细胞数和SPP1计数，相关残差以QR独立复算；另从缓存复算恶性细胞效应、患者摘要，并用组内中心化算法复算研究内置换P，核查所有输入哈希与保留表。
"""
report+=f"共{total}项技术检查通过（计算阶段{len(orig)}项，独立复核{len(checks)}项）。无新增图件。\n"
report+="\n## Discussion\n评分是外部靶程序代理，可能同时反映其他调控或共同肿瘤状态，不提供方向性因果。bulk混合细胞来源，同文库跨区室样本量16且研究结构影响关联。控制MTA后TYMP联系保留仅为观测条件关联，不能证明MTA旁路、必要性或非必要性。保留新正向事实和跨区室阴性，后续不强迫三个不同区室的相关拼成完整通路。\n"
(d/"TYMP_STAT3_RESULTS_CN.md").write_text(report,encoding="utf8")
(d/"TYMP_STAT3_MANUSCRIPT_MODULE_CN.md").write_text("""# TYMP相关肿瘤状态伴随STAT3相关转录程序
基于不包含TYMP、SPP1及STAT3自身的外部STAT3有符号靶基因评分，16名有效患者的恶性细胞内TYMP与STAT3活性代理均呈正向关联，深度调整后rho中位数为0.146（P=0.000482）；增加检出基因数、周期及缺氧调整后方向保留。四个公共配对队列共255例中，控制MTAP和髓系组成后的合并rho为0.319（95%CI 0.183–0.442，P=8.04×10⁻⁶），四队列方向一致。

同文库肿瘤TYMP与巨噬STAT3代理的研究调整关联未获支持（n=16，rho=0.112，置换P=0.700）。因此，该结果支持TYMP相关肿瘤转录状态与STAT3程序的联系，尚未证明跨细胞调控。此前整体巨噬STAT3—SPP1正向关联及MTA—STAT3端阴性/短路径未获支持结果均保留，不以新关联替代。
""",encoding="utf8")
stats="# TYMP—STAT3完整汇总统计（run_v2）\n\n"
for name in ["MALIGNANT_SUMMARY.tsv","MALIGNANT_LEAVE_STUDY.tsv","SAME_LIBRARY_TYMP_MACRO_STAT3.tsv","BULK_TYMP_STAT3_COHORT.tsv","BULK_TYMP_STAT3_META.tsv","BULK_TYMP_STAT3_LEAVE_COHORT.tsv"]:stats+="## "+name+"\n\n"+tab(pd.read_csv(v/name,sep="\t"))+"\n\n"
(d/"ALL_STATISTICS_CN.md").write_text(stats,encoding="utf8")
stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
status={"status":"COMPLETE_VALIDATED","statistics_version":"run_v2","checks":total,"previous_edges_retained":True,"malignant_primary_n":16,"same_library_n":16,"bulk_n":255,"new_mediation":False,"new_figures":0,"updated_at":stamp}
(d/"STATUS.json").write_text(json.dumps(status,indent=2))
bk=d/"root_record_backups";bk.mkdir(exist_ok=True)
entry=f"""
## Phase180 TYMP—STAT3关联补充完成（{stamp}）
用户要求测TYMP到STAT3，另外两段保留。沿Phase179外部靶集与最终评分，不含TYMP/SPP1/STAT3；最终run_v2采用归一化TYMP，原计数诊断保留。
恶性细胞16/16患者正，深度rho中位0.14616/P0.0004824；深度+基因数+周期+缺氧16/16正/rho0.11749。公共255例4/4正，MTAP+髓系合并rho0.31856/P8.04e-6，控制MTA后rho0.34635/P5.53e-7，主评分加缺氧仍支持。
同16文库肿瘤TYMP—巨噬STAT3主研究调整rho0.11186/置换P0.7003，不支持稳定跨区室联系；DoRothEA原始名义信号(P0.0218)研究调整后P0.146，保留但不替代主分析。
原MTA—STAT3、整体巨噬STAT3—SPP1和短路径统计逐字节保留。不跨区室拼中介。{total}检查通过，无新图/新中介。
报告：{d}/TYMP_STAT3_RESULTS_CN.md
正文：{d}/TYMP_STAT3_MANUSCRIPT_MODULE_CN.md
新增TYMP—STAT3肿瘤状态证据可纳入；MTA及原两段均继续保留，跨区室与MTA受体功能依赖仍为待解决环节。
"""
for n in ["PROJECT_STATUS.md","DECISIONS.md","CODEX_HANDOFF.md"]:
 p=r/n
 if not (bk/n).exists():shutil.copy2(p,bk/n)
 old=p.read_text()
 if "## Phase180 TYMP—STAT3关联补充完成" not in old:p.write_text(entry+"\n"+old,encoding="utf8")
p=r/"RESULTS_MANIFEST.tsv"
if not (bk/p.name).exists():shutil.copy2(p,bk/p.name)
if "phase180_TYMP_STAT3" not in p.read_text():
 with p.open("a",newline="") as f:csv.writer(f,delimiter="\t").writerow(["phase180_TYMP_STAT3","L2_supplement","report",str(d/"TYMP_STAT3_RESULTS_CN.md"),stamp,"PASS",str(total)+" checks","Malignant/bulk association supported;cross-compartment unsupported;prior edges retained"])
de=d/"delivery";de.mkdir(exist_ok=True)
for n in ["TYMP_STAT3_RESULTS_CN.md","TYMP_STAT3_MANUSCRIPT_MODULE_CN.md","ALL_STATISTICS_CN.md","STATUS.json"]:shutil.copy2(d/n,de/n)
for n in ["MALIGNANT_SUMMARY.tsv","SAME_LIBRARY_TYMP_MACRO_STAT3.tsv","BULK_TYMP_STAT3_META.tsv"]:shutil.copy2(v/n,de/n)
hh={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in de.iterdir()};(d/"DELIVERY_SHA256.json").write_text(json.dumps(hh,indent=2));print(json.dumps({"status":status,"hashes":hh},ensure_ascii=False))
