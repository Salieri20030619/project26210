from pathlib import Path
import json,hashlib,zipfile,shutil,sys,platform
from datetime import datetime
from zoneinfo import ZoneInfo
import pandas as pd,numpy as np,scipy,anndata
P=Path(__file__).resolve().parents[1];R=P.parent
q=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());assert q['status']=='PASS'
assert (P/'logs/validate.exit').read_text().strip()=='0'
assert q['round12_patient_sign_changes']==0 and q['round12_summary_P_changes']==0 and q['round12_summary_FDR_changes']==0
pr=json.loads((P/'resources/ANALYSIS_PROTOCOL.json').read_text())
s=pd.read_csv(P/'tables/ALL_EFFECT_SUMMARIES.tsv',sep='\t')
cx=pd.read_csv(P/'tables/TUMOUR_TYMP_MACRO_PROCESS.tsv',sep='\t')
def fmt(v):return 'NA' if pd.isna(v) else f'{v:.6g}'
def effect(branch,condition,source,target,model):
 return s[s.branch.eq(branch)&s.condition.eq(condition)&s.source.eq(source)&s.target.eq(target)&s.model.eq(model)].iloc[0]
def line(row):
 return f"{row.source}—{row.target}：n={row.patients}，{row.positive}正/{row.negative}负，中位rho={fmt(row.median_rho)}，P={fmt(row.P)}，FDR={fmt(row.FDR)}"
lines=[
'26210：TYMP参与过程优先验证及无预设方向的免疫特征探索',
'Phase314｜2026-10-04',
'',
'主要发现',
'在既有Obradovic未治疗肾癌11患者中，新增分析支持TYMP与CTSB、TGFB1的细胞内表达联系；结合既有TYMP—SPP1结果，可将蛋白处理及基质协作相关过程列为参与线索。免疫探索显示，TYMP和SPP1对应的转录特征有共同部分，也有明确的方向差异：TYMP与干扰素、MHC-I及吞噬程序正向，SPP1与干扰素和MHC-II程序负向。这些方向在严格成员复评中保留。',
'本轮按用户要求以TYMP参与为主，不再将重做末端空间排斥或建立完整串联因果闭环设为前提。参与过程和免疫探索采用分开的统计家族，完整保留反向、未达名义显著和不可推断结果。',
'',
'一　TYMP与CTSB、TGFB1的关联在宽松和严格成员中均得到支持',
'首先使用已审成熟巨噬簇的全部CD45+成员，共15056细胞/11患者；然后复评Phase305原严格成员13065细胞/同11患者。没有重新分群、挑患者或把两组成员当独立复制。',
'宽松主分析：']
for src,tgt in [('TYMP','SPP1'),('TYMP','CTSB'),('TYMP','TGFB1'),('SPP1','CTSB')]:
 for model in ['core','context']:
  lines.append(('基础调整；' if model=='core' else '加干扰素/缺氧调整；')+line(effect('participation','L2_broad',src,tgt,model)))
lines+=['',
'严格复评：']
for src,tgt in [('TYMP','SPP1'),('TYMP','CTSB'),('TYMP','TGFB1'),('SPP1','CTSB')]:
 for model in ['core','context']:
  lines.append(('基础调整；' if model=='core' else '加环境调整；')+line(effect('participation','L0_strict',src,tgt,model)))
lines+=['',
'宽松TYMP—SPP1基础模型9/11正，P=0.06543，未达到原P<0.05；加入预先固定的环境程序后11/11正并通过FDR。严格成员的基础模型准确复现Phase305全部11名患者效应，10/11正。该结果有成员与模型依赖，不表述为所有条件都显著。',
'TYMP—CTSB、TYMP—TGFB1在两套成员与两个模型中均正向并通过FDR，构成本轮优先过程线索。SPP1—CTSB在宽松成员为9/11正、P=0.06543；严格成员11/11正并通过FDR，必须保留这个差别。',
'RAB37宽松仅5患者、严格仅2患者达到单基因检出要求。TYMP—RAB37方向混合，SPP1—RAB37分别5/5及2/2负，仅描述，不计算P/FDR，不列为支持性分泌机制。',
'本轮是对Phase280/284过程候选在既有独立队列中的扩展；该队列已用于Phase305，不是再增加11例新患者。旧Phase305 P/FDR保留；本轮FDR随新的预先定义家族计算，不能互相覆盖。',
'',
'二　TYMP偏向干扰素、抗原呈递与吞噬相关转录特征，SPP1呈不同免疫轮廓',
'不预设阳性方向或免疫抑制结论，固定全部50个Hallmark程序，加6个GO免疫过程；每个程序移除TYMP、SPP1。仅在冻结方案后计算结果。下面为宽松主分析的固定免疫视图，括号内方向为独立患者正/负数：',
'程序\tTYMP rho（正/负）\tTYMP P / FDR\tSPP1 rho（正/负）\tSPP1 P / FDR']
for name in pr['immune_display_fixed']:
 a=effect('immune','L2_broad','TYMP',name,'core');b=effect('immune','L2_broad','SPP1',name,'core')
 lines.append(f"{name}\t{fmt(a.median_rho)}（{a.positive}/{a.negative}）\t{fmt(a.P)} / {fmt(a.FDR)}\t{fmt(b.median_rho)}（{b.positive}/{b.negative}）\t{fmt(b.P)} / {fmt(b.FDR)}")
lines+=['',
'以上固定免疫视图各模型均n=11。TYMP与干扰素α/γ、MHC-I、吞噬、补体程序正向，与TNFα/NF-κB程序负向；其中主要方向在严格成员保持。宽松炎症反应、IL-2/STAT5关联在严格成员未达到P<0.05，列为成员依赖。',
'SPP1与干扰素α/γ、MHC-II、促进T细胞活化的相关程序负向，与IL-2/STAT5程序弱正向；上述方向在严格成员保持。其MHC-I负向、补体正向仅严格成员达到P<0.05；宽松分别P=0.06543，不称主分析显著。',
'TYMP与T细胞负调节、正调节两个GO基因集均正向。这两个集合分别146和267基因，重叠53基因，不能据单个名称给TYMP贴统一免疫抑制标签；也不能把巨噬中的IL-2/STAT5评分称为实际T细胞增殖。SPP1的负调节GO集合本轮没有名义显著正向关联。',
'当前数据不支持两个目标在宽松主分析中均普遍增强IL-6/JAK/STAT3、TGF-β反应或淋巴细胞趋化程序。单基因TGFB1相关与整套TGF-β反应评分是不同终点，不能互换。',
'TYMP—SPP1同向关联与两者干扰素方向不同可以同时成立，说明其转录联系不等于两个基因代表完全相同的免疫状态。完整50个Hallmark及6个GO、全部宽松和严格结果见ALL_EFFECT_SUMMARIES.tsv，不只保留上述条目。',
'',
'三　肿瘤来源TYMP与髓系过程仍需保持来源区分',
'按原患者配对CD45−肿瘤样细胞TYMP与CD45+巨噬读出，以原计数伪汇总表达比较。SPP1/CTSB/TGFB1有8患者，RAB37宽松7例/严格4例，全部仅描述。以下为批次/pT调整结果：']
for _,row in cx[cx.model.eq('batch_pT_adjusted')].iterrows():
 lines.append(f"{row.condition}，{row.target}：n={row.patients}，rho={fmt(row.rho)}，P/FDR不作推断。")
lines+=['',
'这一来源配对没有形成一致正向关系。它不影响已成立的巨噬胞内联系，但不能被髓系阳性替代为肿瘤TYMP处理效应。客户肿瘤TYMP实验背景继续为研究起点；本轮未重分析客户湿实验。',
'',
'四　方法',
'数据为Phase305已质控的Obradovic原始计数及已审身份。宽松采用成熟巨噬簇2/3/5/6/8/9/10/15/17内所有CD45+细胞；严格成员完全沿用原lineage_Macrophage标志。患者内技术要求分别至少50/100细胞，两端单基因各至少10细胞检出。',
'按全部测量基因的细胞总计数归一化为CP10K后log1p。程序为可测基因的等权均值；至少80%集合基因可测，至少20%集合基因在max(3,ceil(1%细胞数))个细胞中检出。缺测不补零；实测零保留。各患者/文库资格见FEATURE_COVERAGE.tsv。',
'基础模型对两变量排序，控制深度、原始特征复杂度、文库和冻结髓系亚群的设计矩阵，再计算残差相关。参与过程的环境模型加入IFNα、IFNγ、缺氧程序，事先移除全部6个基因对涉及的基因。免疫探索以这些程序为结局，不再把它们作为协变量调掉。残差维度至少8。',
'患者内效应不产生细胞级P；按独立患者的正/负方向作双侧精确二项检验，有效非零效应至少10患者才推断。原P<0.05为名义筛选。参与家族6对×2模型×2成员=24项计划，其中RAB37相关8项仅描述，16项具有效P并在家族内BH；探索家族56程序×2目标×2成员=224项统一BH。报告全部原P、FDR、计划/实际家族大小，不以FDR作为进入候选的唯一门槛。',
'来源配对先按患者、细胞池加总原计数并计算log1pCP10K，单基因至少5细胞检出；未调整及作者批次/pT调整均保留。细胞数和文库数不替代患者数。',
'',
'五　讨论与下一步重点',
'本轮确立的是表达与转录程序关联，其中多数效应较小但患者方向一致。它支持TYMP参与SPP1相关状态及蛋白处理/基质协作过程的候选解释；不将其升级为SPP1分泌、CTSB酶活、活性TGFβ、实际吞噬/抗原呈递或T细胞迁移已被测得。该探索只覆盖巨噬细胞内的免疫转录特征，不能推广到全部免疫细胞或患者整体免疫结局。',
'研究推进继续优先TYMP参与，CTSB/TGFB1作为已有跨队列支持的过程候选；无方向免疫探索提示TYMP具有更广的免疫响应背景，解释SPP1联系时应区分共同与不同程序。下一有效增量应优先衔接来源明确的TYMP处理/过程资料，不再把末端空间重做设为主线门槛。',
'既有肿瘤来源未支持、外源TYMP刺激反向及末端免疫/空间结果均保留，各自按来源与测量层级解释；本轮不以未新增末端验证否定参与线索。',
'',
'六　验收与交付',
f"独立复核{q['checks']}项通过：原始选定计数逐项一致、全部评分重建、2728患者模型独立正交投影复算、患者符号检验组合数复核、BH独立函数核对；最大相关误差{q['max_effect_difference']:.3g}。主汇总248行，肿瘤来源描述16行，全部可复核。",
'验收核清两处定义/数值问题：209细胞原始特征复杂度与合并同名基因后的非零基因数相差1–2，已从全部受影响文库原始矩阵确认；模型继续使用原QC定义。独立重建与原评分在浮点精度内一致，但等价运算顺序可改变极少量并列秩；精确保存输入的独立SVD复算通过，另将全部评分取12位小数统一并列秩复评：患者方向、全部P/FDR均未改变，最大患者rho差0.00002784。以上均未改分析效应、成员或阈值，失败记录留档。验收启动器的Python兼容故障亦已修复并保存。',
'完整包包含报告、全部统计/资格表、代码、冻结协议、模型输入缓存、独立验收与来源记录；大原始矩阵及细胞计数对象保留服务器。未修改旧对象、Word或图件。',
'',
'资源出处',
'Obradovic原数据：https://data.mendeley.com/datasets/nc9bc8dn4m/1',
'作者样本/代码：https://github.com/Aleksobrad/single-cell-rcc-pipeline',
'50个Hallmark：沿用项目已下载的Enrichr MSigDB_Hallmark_2020固定版本，来源/哈希在资源登记中；没有假称使用最新版本。',
'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020',
'标准集合说明：https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp',
'6个GO条目从项目GO_BP_local.gmt按固定GO编号提取，原资源为GO.db/org.Hs.eg.db 3.22.0；软件与输入SHA均保存。'
]
report=P/'delivery/26210_Phase314_TYMP参与与免疫特征探索.txt'
report.write_text('\n'.join(lines)+'\n',encoding='utf-8-sig')
environment={'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'anndata':anndata.__version__}
(P/'qc/SOFTWARE_VERSIONS.json').write_text(json.dumps(environment,indent=2))
# Full comparison including non-significant and opposite directions.
primary=s[s.condition.eq('L2_broad')];sensitivity=s[s.condition.eq('L0_strict')]
comp=primary.merge(sensitivity,on=['branch','source','target','model'],suffixes=('_broad','_strict'),validate='one_to_one')
comp['median_direction_same']=np.sign(comp.median_rho_broad)==np.sign(comp.median_rho_strict)
comp.to_csv(P/'tables/BROAD_STRICT_COMPLETE_COMPARISON.tsv',sep='\t',index=False)
assets=[]
for directory in ['code','resources','tables','qc']:
 assets += [f for f in (P/directory).iterdir() if f.is_file()]
assets += list((P/'inputs').glob('*.npz'))+[P/'inputs/CACHE_MANIFEST.json',report]
assets += [f for f in (P/'logs').iterdir() if f.is_file() and f.suffix in ['.log','.exit']]
hashes={str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest() for f in assets}
manifest=P/'delivery/PACKAGE_SHA256.json';manifest.write_text(json.dumps(hashes,ensure_ascii=False,indent=2))
zipfile_path=P/'delivery/26210_Phase314_TYMP_participation_immune_results.zip'
with zipfile.ZipFile(zipfile_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in assets+[manifest]:z.write(f,str(f.relative_to(P)))
with zipfile.ZipFile(zipfile_path) as z:
 assert z.testzip() is None
 for name,h in hashes.items():assert hashlib.sha256(z.read(name)).hexdigest()==h
out={'status':'PASS','time':datetime.now(ZoneInfo('Asia/Taipei')).isoformat(),'zip':str(zipfile_path),'zip_bytes':zipfile_path.stat().st_size,'zip_sha256':hashlib.sha256(zipfile_path.read_bytes()).hexdigest(),'report':str(report),'report_sha256':hashlib.sha256(report.read_bytes()).hexdigest(),'files':len(hashes),'all_member_hashes_and_CRC':'PASS'}
(P/'qc/DELIVERY_VALIDATION.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
print(json.dumps(out,ensure_ascii=False),flush=True)
