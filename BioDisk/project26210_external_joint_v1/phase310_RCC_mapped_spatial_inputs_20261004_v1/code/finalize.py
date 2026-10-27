from pathlib import Path
import json,csv,hashlib,zipfile,datetime,shutil
P=Path(__file__).resolve().parents[1];root=P.parent
v=json.loads((P/'qc/FINAL_VALIDATION.json').read_text());assert v['status']=='PASS'
for n in ['acquire_pilot_v2','reconstruct_pilot','validate_pilot','controller']:
 assert (P/'logs'/(n+'.exit')).read_text().strip()=='0',n
report=P/'delivery/26210_Phase310_肾癌空间输入恢复与试点重建.txt'
assert report.exists()
r=json.loads((P/'qc/RECONSTRUCTION_COMPLETE.json').read_text());r['status']='TECHNICAL_INPUT_VALIDATED';r['validation']='qc/FINAL_VALIDATION.json';(P/'qc/RECONSTRUCTION_COMPLETE.json').write_text(json.dumps(r,indent=2))
src=P/'downloads/pilot_tumor1_region6.csv.gz'
large={'source':str(src),'bytes':src.stat().st_size,'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'md5':hashlib.md5(src.read_bytes()).hexdigest(),'in_delivery':False}
(P/'qc/RAW_ASSET_MANIFEST.json').write_text(json.dumps(large,indent=2))
files=[]
for d in ['code','tables','qc','logs']:
 files += [f for f in (P/d).iterdir() if f.is_file() and not f.name.endswith('.pid') and f.name not in ['DELIVERY_VALIDATED.json','LOCAL_DELIVERY_VERIFIED.json']]
files += [f for f in (P/'resources').iterdir() if f.is_file()]
for n in ['figshare_current.json','figshare_file.json','proseg_readme.txt','paper_final.html','paper_final_text.txt','nicheverse_tree.json','nicheverse_reproduction.py']:
 files.append(P/'metadata'/n)
files.append(report)
archive=P/'delivery/26210_Phase310_spatial_input_qualification.zip'
manifest={str(f.relative_to(P)):{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files}
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for f in files:z.write(f,str(f.relative_to(P)))
 z.writestr('SHA256_MANIFEST.json',json.dumps(manifest,indent=2,ensure_ascii=False))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for name,m in manifest.items():assert hashlib.sha256(z.read(name)).hexdigest()==m['sha256']
now=datetime.datetime.now().astimezone().isoformat()
delivery={'status':'PASS','files':len(manifest),'archive':str(archive),'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'CRC':'PASS','all_file_sha256':'PASS','report_sha256':hashlib.sha256(report.read_bytes()).hexdigest()}
(P/'qc/DELIVERY_VALIDATED.json').write_text(json.dumps(delivery,indent=2))
status={'status':'TECHNICAL_INPUT_QUALIFICATION_COMPLETE','time':now,'biological_effect_analysis':'not_started_missing_donor_region_identity_and_cell_annotation','current_jobs':'none','effective_analyzed_donors':None,'downloaded_regions':1,'mapped_file_regions':10,'observed_tumor_labels':3,'author_claimed_patients':4,'validation':{'status':'PASS','named_checks':144},'delivery':delivery,'next':'Resolve public tumor-label to patient and tumor/adjacent region mapping; obtain qualified cell identity reference before effects.'}
(P/'STATUS.json').write_text(json.dumps(status,indent=2))
backup=P/'resources/root_before_final';backup.mkdir(exist_ok=False)
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv']:shutil.copy2(root/n,backup/n)
note='''## Phase310终态（2026-10-04）：肾癌空间入口打通并完成技术试点
用户“继续下一步”。已提取Figshare25685961官方folder_structure，10同名文件映射tumor1/3/4的6/2/2区域；作者四患者说明与三标签差异、肿瘤/邻近区域定义未核清，不把区域当患者。
按最小文件预选tumor1_region6技术试点，官方API代理下载103790480字节，MD5/gzip通过。3194464转录本，380非控制基因，16825合格计数细胞、技术QC16569。原分配/核内整数矩阵和核内中心代理已重建，不能称Proseg重分割、真实膜接触或病理边界。SPP1检出7112/核内4602细胞，尚非巨噬标签；TYMP完整文件未观察到，官方面板资格未补齐，不当零。CA9仅74/核内55，不能将预选小区域当患者肿瘤代表。
Nicheverse生产接口仍TLS失败，作者复现代码需仓库外注释对象；core示例不代替全队列。补充PDF未取到有效文件；官方正文及全部访问日志保留。初始摘要未分配标记修正为-1，qc/PILOT_INPUT_VALIDATION_v3.json为当前，重建本身从开始即正确排除。
独立全原始CSV/全部细胞和基因总量/核内坐标/64完整矩阵行/10文件映射144项PASS，获取v2、重建、验证及controller退出0。当前完成的是技术输入资格与重建，无新生物学效应/P/FDR或图；缺少可靠患者/组织区域映射和细胞身份参考，尚未启动关联分析。旧上游阳性、Phase308/309和全部阴性、旧C4/对象/Word/文章保持。
权威目录phase310_RCC_mapped_spatial_inputs_20261004_v1；报告delivery/26210_Phase310_肾癌空间输入恢复与试点重建.txt，轻量包26210_Phase310_spatial_input_qualification.zip，含重建计数/坐标/代码/全部技术表，原始文件仅服务器。包CRC/全部SHA通过。当前无本轮计算后台；下一输入门槛为公开患者与区域定义、可核查谱系身份，不直接批量下载其余区域或用此局部试点写屏障阳性。

'''
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=root/n;f.write_text(note+f.read_text())
with (root/'RESULTS_MANIFEST.tsv').open('a',newline='') as f:
 w=csv.writer(f,delimiter='\t')
 for key,kind,path in [('report','technical_report',report),('package','validated_package',archive),('mapping','source_mapping',P/'tables/FILE_REGION_MAP.tsv'),('raw','raw_transcripts',src),('assigned','sparse_counts',P/'resources/pilot_assigned_counts_qv20.npz'),('nuclear','sparse_counts',P/'resources/pilot_nuclear_counts_qv20.npz'),('validation','integrity_validation',P/'qc/FINAL_VALIDATION.json')]:
  w.writerow(['phase310_'+key,'Phase310',kind,str(path),now,'PASS',str(P/'qc/FINAL_VALIDATION.json'),'Technical qualification only; donor-region and identity unresolved; no biological effect tests'])
(P/'logs/finalize.exit').write_text('0\n')
print(json.dumps(delivery,indent=2),flush=True)
