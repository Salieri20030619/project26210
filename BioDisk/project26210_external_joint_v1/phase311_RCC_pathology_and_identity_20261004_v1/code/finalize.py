from pathlib import Path
import json,datetime,hashlib,zipfile,csv,shutil
P=Path(__file__).resolve().parents[1];root=P.parent;now=datetime.datetime.now().astimezone().isoformat()
stages=['install_proseg','proseg','cluster_identity','identity_controller','validate_segmentation','run_regional','validate_regional','source_and_score_audit','run_boundary','validate_boundary']
for n in stages:assert (P/'logs'/(n+'.exit')).read_text().strip()=='0',n
vals={}
for n in ['PROSEG_IDENTITY_VALIDATION','REGIONAL_INDEPENDENT_VALIDATION','BOUNDARY_REGIONAL_INDEPENDENT_VALIDATION']:
 a=json.loads((P/'qc'/(n+'.json')).read_text());assert a['status']=='PASS';vals[n]=a['named_checks']
assert sum(vals.values())==440
data=json.loads((P/'qc/REGIONAL_ANALYSIS_COMPLETE.json').read_text())
summary={'status':'COMPLETE_VALIDATED','time':now,'stage':'Phase311','scope':'single-region segmentation, identity, cellular source, descriptive spatial endpoints and numerical-boundary sensitivity','biological_effect':'no consistent new barrier support; primary SPP1-myeloid endpoints not assessable by fixed gate','sample_unit':'one named region; donor mapping unresolved; no inferential P/FDR','identity_counts':data['identity_counts'],'planned_models':72,'assessable_models':56,'numerical_sensitivity_models':72,'validation_checks':vals,'validation_total':440,'running_jobs':'none','old_objects_reports':'unchanged','no_new_biological_figures':True,'continuous_authorization':'persists for next qualified inputs; no routine reconfirmation required'}
(P/'STATUS.json').write_text(json.dumps(summary,indent=2))
assets=[]
for f in [root/'phase310_RCC_mapped_spatial_inputs_20261004_v1/downloads/pilot_tumor1_region6.csv.gz',*sorted((P/'resources/proseg_region6_v1').iterdir()),P/'resources/proseg_region6_identity.h5ad',P/'resources/proseg_region6_identity_score_boundary.h5ad',P/'resources/proseg206/bin/proseg']:
 if f.is_file():assets.append({'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(P/'qc/SERVER_ASSET_MANIFEST.json').write_text(json.dumps(assets,indent=2))
files=[P/'STATUS.json']
for d in ['code','tables','qc','metadata']:
 for f in (P/d).iterdir():
  if f.is_file() and f.name not in ['analysis_summary.html','DELIVERY_VALIDATED.json','LOCAL_DELIVERY_VERIFIED.json'] and f.suffix!='.png':files.append(f)
for f in (P/'logs').iterdir():
 if f.is_file() and not f.name.endswith('.pid') and not f.name.startswith('finalize'):files.append(f)
for f in (P/'resources').iterdir():
 if f.is_file():files.append(f)
files.extend(f for f in (P/'resources/validation_attempt1').iterdir() if f.is_file())
for n in ['maxpost-counts.csv.gz','expected-counts.csv.gz','cell-metadata.csv.gz','gene-metadata.csv.gz']:files.append(P/'resources/proseg_region6_v1'/n)
report=P/'delivery/26210_Phase311_重分割与SPP1空间来源及关联结果.txt';files.append(report)
files=sorted(set(files));manifest={str(f.relative_to(P)):{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files}
zpath=P/'delivery/26210_Phase311_resegmentation_and_spatial_results.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for f in files:z.write(f,str(f.relative_to(P)))
 z.writestr('SHA256_MANIFEST.json',json.dumps(manifest,indent=2,ensure_ascii=False))
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 for n,m in manifest.items():assert hashlib.sha256(z.read(n)).hexdigest()==m['sha256']
delivery={'status':'PASS','archive':str(zpath),'bytes':zpath.stat().st_size,'sha256':hashlib.sha256(zpath.read_bytes()).hexdigest(),'manifest_files':len(files),'CRC':'PASS','each_member_SHA256':'PASS','report':str(report),'report_sha256':hashlib.sha256(report.read_bytes()).hexdigest()}
(P/'qc/DELIVERY_VALIDATED.json').write_text(json.dumps(delivery,indent=2));summary['delivery']=delivery;(P/'STATUS.json').write_text(json.dumps(summary,indent=2))
backup=P/'resources/root_before_final';backup.mkdir(exist_ok=False)
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:shutil.copy2(root/n,backup/n)
note='''## Phase311终态（2026-10-04）：完成重分割、身份、区域空间端点及数值复核
用户“持续进行下一步”。已连续实施Proseg2.0.6（Bioconda SHA、隔离、8线程默认600迭代）→排除TYMP/SPP1的身份聚类→严格成员→72空间模型→问题修复与完整验收。权威目录phase311_RCC_pathology_and_identity_20261004_v1，旧对象/报告/C4只读。
同一预选tumor1_region6：2695925转录本、16697细胞、380基因，QC16238、11群；严格成熟髓系788/成纤维1911/CD8 206，共同几何785/1899/206。上皮具有AQP1/FXYD2/EPCAM，不直接叫癌细胞；群9未知，壁细胞不合并CAF。SPP1计数以上皮为主：原分配87.20%、核内88.62%、Proseg96.69%。严格髓系检出168/69/28；50μm窗口734中心中Proseg仅26检出，低于预设30，不降低阈值，SPP1主端点NA。原分配/核内SPP1—CD8距离调整rho0.099245/0.073940保留；SPP1—成纤维调整-0.021901/-0.003130。成纤维—CD8最近距离50μm原始0.068493、调整-0.088873；100μm0.073596/-0.089581。没有一致屏障方向，不能把该局部窗口阴性移植成整个研究无信号。1区域，无生物学P/FDR，无病理肿瘤进入推断。
独立源转录本汇总6344860矩阵条目逐项一致。验收发现Proseg2.0.6末次时间累计+1使分配权重最高1.00001，源码确认、原值/失败记录保留；923条在0.1边界（0.0382%）单列保守去除敏感性，严格髓系SPP1不变，72模型可评估状态不变、最大rho变化0.002600。76输入+182主空间+182敏感性=440项PASS，各正式阶段退出0；原始软件前景0.9参数未实现，不声称执行此门槛。
来源：已取2024预印本补充PDF，图示四肿瘤但无文件映射，不能据图补患者ID。另2025 10x RCC465534细胞与Phase309新690322不同输入，已取得experiment与H&E配准，病理GeoJSON尚未取得；S-BIAD2470为2025同数据重新打包，不加样本。10x与Figshare未解决的供者/病理映射完整保留，不批量下载伪独立新样本。
交付delivery/26210_Phase311_重分割与SPP1空间来源及关联结果.txt及26210_Phase311_resegmentation_and_spatial_results.zip；包含原/核内/重分割全部72模型、来源、冻结成员、计数、数值敏感性72模型及验收。ZIP CRC和所有SHA通过。无新增生物学图，未覆盖Word/正式文章。当前无本轮计算后台；持续推进授权保留，下一有效增量需可核验肾癌病理/供者/髓系身份或直接SPP1扰动配对输入，不从region6换阈值造阳性。Phase280/305上游阳性、Phase308及此前全部正负保持。

'''
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=root/n;f.write_text(note+f.read_text())
with (root/'RESULTS_MANIFEST.tsv').open('a',newline='') as f:
 w=csv.writer(f,delimiter='\t')
 entries=[('report','complete_technical_report',report),('package','validated_package',zpath),('identity','frozen_cell_identity',P/'tables/CELL_IDENTITY_FROZEN.tsv.gz'),('source','SPP1_cellular_source',P/'tables/SPP1_CELLULAR_SOURCE_AUDIT.tsv'),('effects','descriptive_spatial_effects',P/'tables/REGIONAL_ALL_SPATIAL_EFFECTS.tsv'),('boundary_effects','numerical_sensitivity',P/'tables/BOUNDARY_REGIONAL_ALL_SPATIAL_EFFECTS.tsv'),('object','new_annotated_counts',P/'resources/proseg_region6_identity.h5ad'),('raw_proseg','new_resegmentation',P/'resources/proseg_region6_v1')]
 for key,kind,fpath in entries:w.writerow(['phase311_'+key,'Phase311',kind,str(fpath),now,'PASS',str(P/'qc/PROSEG_IDENTITY_VALIDATION.json')+';'+str(P/'qc/REGIONAL_INDEPENDENT_VALIDATION.json'),'single-region descriptive; no inferential P/FDR; primary SPP1 gate unmet; all previous evidence retained'])
(P/'logs/finalize.exit').write_text('0\n');print(json.dumps(delivery,indent=2),flush=True)
