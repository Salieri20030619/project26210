from pathlib import Path
import json,hashlib,datetime,shutil,csv,fcntl,os
p=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase258_spatial_GSEA_figure_20260926_v3');root=p.parent
qc=json.loads((p/'qc/FINAL_VALIDATION.json').read_text());assert qc['status']=='PASS'
assert (p/'build_EXIT_STATUS').read_text().strip()=='0' and (p/'finalize_EXIT_STATUS').read_text().strip()=='0'
now=datetime.datetime.now().astimezone().isoformat()
entry=f"""## {now} Phase258 空间大图加入GSEA，按用户建议移出独立H&E并合并面板
用户要求将类似参考Fig2G生成的功能图加入空间大图，并明确允许去掉单独原始H&E或合并。当前新版为{p}/delivery/26210_Spatial_A4_GSEA_v3.pdf，210×297mm单页A–H八面板、22张组织地图+4幅GSEA曲线。
A将原组织组成与TYMP表达合并；B代谢背景、C旧SPP1-state丰度、D代谢/状态高值区域；E为Lysosome与Cholesterol metabolism两通路×GSE175540(n23)/TJ-RCC(n12)；F状态/内皮区域；G六患者独立切片；H两例配对SPP1表达。独立原H&E两图移出，空间图的组织底图和比例尺保留。原四张免疫/侵袭功能地图由GSEA替换。Phase255旧版完整保留。
E采用Phase256主比较（SPP1-state q05前20%对其余，去七定义/目标基因），不是组成调整版本；由完整患者排序重建running ES和命中条带，沿用原NES/P/完整KEGG主家族FDR。四曲线与原结果独立核对通过，无新统计/模型。图内q显示FDR；GLMF两例不被拿来做推断性GSEA。
PDF、SVG、300dpi PNG、600dpi LZW TIFF、独立逐panel中文图注及绘图数据ZIP已生成。最小6pt，整页Poppler渲染及GSEA局部均已实际视觉检查，无文字图例挤压；FINAL_VALIDATION PASS。十张组织层表字节一致，几何/阈值/选片与旧图完全一致，旧图哈希、源输入及既有报告不变。构建/收尾退出0；元数据序列化修复留内部，图像像素和统计未改。
本地交付目标C:/Users/31937/Downloads/26210/Phase258_空间GSEA大图_v3/。图件独立交付，未插入报告；3G热点分析尚未启动。后续空间图审核以Phase258 v3为准；机制小节仍以Phase257 editorial_v2为准。全部进度只在服务器记录，无运行的本轮绘图任务。
"""
lock=(root/'phase258_record_update.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX)
backup=p/'qc/record_backup';backup.mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 f=root/name;old=f.read_text();assert 'Phase258 空间大图加入GSEA' not in old
 shutil.copyfile(f,backup/name);tmp=f.with_name(name+'.phase258.tmp');tmp.write_text(entry+'\n'+old,encoding='utf-8');os.replace(tmp,f)
f=root/'RESULTS_MANIFEST.tsv';shutil.copyfile(f,backup/f.name)
with f.open('a',newline='') as h:
 w=csv.writer(h,delimiter='\t',lineterminator='\n')
 for i,m in enumerate(json.loads((p/'DELIVERY_MANIFEST.json').read_text()),1):
  w.writerow([f'phase258_spatial_GSEA_v3_{i}',p.name,'figure_package',m['path'],now,'PASS',str(p/'qc/FINAL_VALIDATION.json'),'A4 eight panels,22 tissue maps,4 validated GSEA curves; SHA256='+m['sha256']])
fcntl.flock(lock,fcntl.LOCK_UN)
(p/'qc/RECORD_UPDATE_VALIDATION.json').write_text(json.dumps({'status':'PASS','time':now,'updated':['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv'],'entries':6},indent=2))
print('Phase258 recorded as current spatial figure.')
