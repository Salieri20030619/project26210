from pathlib import Path
import json,datetime,hashlib,shutil,csv,fcntl,os
p=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase259_spatial_GSEA_reference_style_20260926_v4');root=p.parent
q=json.loads((p/'qc/STYLE_ONLY_CHECK.json').read_text());assert q['status']=='PASS' and q['outside_E_changed_pixels']==0
assert (p/'finalize_EXIT_STATUS').read_text().strip()=='0'
now=datetime.datetime.now().astimezone().isoformat()
entry=f"""## {now} Phase259 E面板改为参考2G贯穿竖线渐变样式，v4完成
用户指出参考2G的竖线贯穿曲线区域并有颜色区分。本轮高分辨率核对原论文后，仅修改E：289条真实通路基因命中竖线从ES=0连接到该位置running ES；移除底部独立条带，曲线与竖线共用蓝—灰—紫渐变及0–0.7共同ES色标，增加小型色标与矩形框。明确本图颜色编码实际ES，不猜测原论文未说明的颜色定义。
A4八面板、22张组织地图和4曲线保留。与Phase258核对所有排序、基因集合、ES/NES/P/FDR完全一致；E以外PNG像素变化为0。整页PDF与局部实际视觉检查通过，STYLE_ONLY_CHECK与FINAL_VALIDATION PASS，构建及收尾退出0。
当前空间大图：{p}/delivery/26210_Spatial_A4_GSEA_v4.pdf。PDF/SVG/PNG300dpi/TIFF600dpi/中文图注/绘图数据ZIP六文件已同步C:/Users/31937/Downloads/26210/Phase259_空间GSEA参考样式_v4/，SHA256全部匹配。原v3及更早图保留，原报告未改，无新统计或模型，无运行绘图任务。

用户另问参考Fig5D是否适合加入空间。本轮仅核对：Fig5D是基于mIF的三元相对距离构成图，三轴CD8/APOE/POSTN，比较TIB与non-TIB；当前空间图尚无此类面板。项目可考虑用恶性表达区域作为参照，评估到SPP1-state/内皮/T细胞相关区域的距离构成；须明确区域而非单细胞测量、距离方向/归一化和目标资格，不能直接借用CAF/CD8身份或TIB分组。既有Phase246有表达区域图距离资料；三元图尚未计算、未生成、未加入。3G热点亦未启动。
"""
lock=(root/'phase259_record_update.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX)
backup=p/'qc/record_backup';backup.mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md']:
 f=root/name;old=f.read_text();assert 'Phase259 E面板改为参考2G贯穿竖线渐变样式' not in old
 shutil.copyfile(f,backup/name);tmp=f.with_name(name+'.phase259.tmp');tmp.write_text(entry+'\n'+old,encoding='utf-8');os.replace(tmp,f)
f=root/'RESULTS_MANIFEST.tsv';shutil.copyfile(f,backup/f.name)
with f.open('a',newline='') as h:
 w=csv.writer(h,delimiter='\t',lineterminator='\n')
 for i,m in enumerate(json.loads((p/'DELIVERY_MANIFEST.json').read_text()),1):
  w.writerow([f'phase259_spatial_GSEA_v4_{i}',p.name,'figure_package',m['path'],now,'PASS',str(p/'qc/FINAL_VALIDATION.json'),'Reference2G colored stems; unchanged data and non-E pixels; SHA256='+m['sha256']])
fcntl.flock(lock,fcntl.LOCK_UN)
(p/'DELIVERY_CONFIRMATION.json').write_text(json.dumps({'status':'PASS','time':now,'local_root':'C:/Users/31937/Downloads/26210/Phase259_空间GSEA参考样式_v4','files_verified':6,'verification':'Windows SHA256 matched server manifest before rename','current_figure':'26210_Spatial_A4_GSEA_v4.pdf'},ensure_ascii=False,indent=2))
(p/'qc/RECORD_UPDATE_VALIDATION.json').write_text(json.dumps({'status':'PASS','time':now,'updated':['PROJECT_STATUS.md','CODEX_HANDOFF.md','DECISIONS.md','RESULTS_MANIFEST.tsv'],'entries':6},indent=2))
print('Phase259 v4 delivery and current project records complete. Fig5D assessed only.')
