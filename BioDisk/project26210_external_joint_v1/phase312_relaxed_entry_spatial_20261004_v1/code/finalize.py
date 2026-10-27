from pathlib import Path
import json,datetime,hashlib,zipfile,csv,shutil
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];root=P.parent;now=datetime.datetime.now().astimezone().isoformat()
for n in ['analyze_relaxed','validate_relaxed']:assert (P/'logs'/(n+'.exit')).read_text().strip()=='0'
v=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());assert v['status']=='PASS'
base=pd.read_csv(P/'tables/COMPLETE_DIRECTIONAL_MATCH.tsv',sep='\t');adj=pd.read_csv(P/'tables/ADJACENT_PATH_AND_DIRECT_EDGE.tsv',sep='\t')
pd.testing.assert_frame_equal(base,adj[base.columns]);assert np.array_equal(adj.adjacent_two_positive,(adj.SPP1_vs_fibroblast_fraction>0)&(adj.fibroblast_fraction_vs_nearest_CD8_distance>0))
check={'status':'PASS','original_direction_rows_unchanged':len(base),'adjacent_sign_flags_verified':len(adj),'new_models_or_thresholds':False}
(P/'qc/INTERPRETATION_TABLE_VALIDATION.json').write_text(json.dumps(check,indent=2))
for name,sha in json.loads((P/'qc/INPUT_SHA256.json').read_text()).items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha
status={'status':'COMPLETE_VALIDATED','time':now,'scope':'user-authorized relaxed-entry exploration before strict re-evaluation','primary_result':'broader coverage and weak positive adjacent unadjusted directions; adjusted downstream direction not retained','QC_cells_broad':{'myeloid':1456,'SPP1_detected_myeloid':59,'fibroblast_enriched':3055,'CD8_candidates':346},'primary50um_centers':1351,'primary50um_SPP1_detected':53,'primary_unadjusted_adjacent_rho':[0.053901,0.115273],'primary_adjusted_adjacent_rho':[0.070474,-0.025749],'models':{'primary_and_fixed_sensitivities':216,'matched_center_reference':216,'duplicated_strict_reference_rows':72},'validation_named_checks':v['named_checks'],'sample_unit':'one named region; no inferential P/FDR','old_objects_results':'unchanged','running_jobs':'none'}
(P/'STATUS.json').write_text(json.dumps(status,indent=2))
files=[P/'STATUS.json']
for d in ['code','tables','qc','logs','resources']:
 for f in (P/d).iterdir():
  if f.is_file() and not f.name.endswith('.pid') and not f.name.startswith('finalize') and f.name not in ['DELIVERY_VALIDATED.json','LOCAL_DELIVERY_VERIFIED.json']:files.append(f)
report=P/'delivery/26210_Phase312_宽松纳入探索结果.txt';files.append(report)
# Include finalizer as source, exclude only its changing logs.
files.append(P/'code/finalize.py');files=sorted(set(files))
manifest={str(f.relative_to(P)):{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files}
zpath=P/'delivery/26210_Phase312_relaxed_entry_results.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for f in files:z.write(f,str(f.relative_to(P)))
 z.writestr('SHA256_MANIFEST.json',json.dumps(manifest,indent=2,ensure_ascii=False))
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 for n,m in manifest.items():assert hashlib.sha256(z.read(n)).hexdigest()==m['sha256']
delivery={'status':'PASS','archive':str(zpath),'bytes':zpath.stat().st_size,'sha256':hashlib.sha256(zpath.read_bytes()).hexdigest(),'files':len(files),'CRC':'PASS','each_member_SHA256':'PASS','report_sha256':hashlib.sha256(report.read_bytes()).hexdigest()}
(P/'qc/DELIVERY_VALIDATED.json').write_text(json.dumps(delivery,indent=2));status['delivery']=delivery;(P/'STATUS.json').write_text(json.dumps(status,indent=2))
backup=P/'resources/root_before_final';backup.mkdir(exist_ok=False)
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','AGENTS.md']:shutil.copy2(root/n,backup/n)
note='''## Phase312终态（2026-10-04）：用户授权宽松优先，已出现相邻两段方向性线索
用户“先不从严格纳入开始，先看更宽松的条件下是否能符合要求”。已另开phase312_relaxed_entry_spatial_20261004_v1，优先L2整群髓系/成纤维富集群+T群CD8检出，再比较L1单标记和L0严格成员仅放低数量门槛；不重新聚类。技术覆盖中心30、成纤维/CD8各20、SPP1检出/未检出各10，全部条件与50/100μm等已先冻结，为用户要求的后验探索。
L2 QC髓系1456、SPP1检出59、成纤维3055、CD8 346；共同几何1449/3037/346，50μm窗口1351中心/53SPP1检出。L1髓系1384/54SPP1检出；L0 788/28。覆盖足够，全部432模型行可计算（216主与固定敏感性+216同中心，其中72严格参考重复列示），原严格16项NA被新门槛释放，历史56数值不改。
宽松Proseg主50μm：SPP1—成纤维原始rho0.053901、调整0.070474；成纤维—最近CD8距离原始0.115273、调整-0.025749。未调整相邻两段在两几何×两半径4组合均正，可保留为弱探索性方向对应；调整后末端未保持。直接SPP1—CD8距离-0.039349/-0.036479另列诊断，不额外设三条关系全部同时正向为保留相邻分段事实的硬门槛。
宽松原分割/核内的SPP1—成纤维调整-0.101023/-0.044900，SPP1—CD8距离+0.072397/+0.047143；方法差异完整保留。L0原计数/核内未调整三边均正的8组合保留，调整未保持。同734中心宽松参考群的Proseg相邻原始0.035533/0.150163、调整0.069645/-0.049574。不能跨计数方法挑各自正向边拼成链。
独立BallTree/直接距离/SVD残差/精确成员/62841输入行1147项PASS；解释表新增相邻两段标志与原值一致，源对象SHA未变。全部正式退出0，1区域无推断性P/FDR。用户要求“先宽松、后严格复评”已写AGENTS；旧限制不构成本轮授权放宽的阻碍，仍完整保留正负/方法差异，不把富集群叫纯CAF，不把上皮当髓系。
交付delivery/26210_Phase312_宽松纳入探索结果.txt与26210_Phase312_relaxed_entry_results.zip，CRC/全部SHA通过。旧结果/对象/Word/文章不改，无新图，无本轮计算后台。当前可保留两段弱方向性线索，尚不能称稳定屏障验证；后续遵照宽松优先再复评的工作顺序。

'''
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=root/n;f.write_text(note+f.read_text())
with (root/'RESULTS_MANIFEST.tsv').open('a',newline='') as f:
 w=csv.writer(f,delimiter='\t')
 for key,kind,path in [('report','relaxed_exploratory_report',report),('package','validated_package',zpath),('effects','descriptive_spatial_effects',P/'tables/ALL_SPATIAL_EFFECTS.tsv'),('members','relaxed_membership',P/'tables/ALL_CONDITION_MEMBERS.tsv.gz'),('adjacent','adjacent_direction_interpretation',P/'tables/ADJACENT_PATH_AND_DIRECT_EDGE.tsv'),('validation','independent_validation',P/'qc/FINAL_VALIDATION.json')]:
  w.writerow(['phase312_'+key,'Phase312',kind,str(path),now,'PASS',str(P/'qc/FINAL_VALIDATION.json'),'User-authorized relaxed exploratory entry; weak adjacent unadjusted correspondence; adjusted downstream not retained; all conditions reported'])
(P/'logs/finalize.exit').write_text('0\n');print(json.dumps(delivery,indent=2))
