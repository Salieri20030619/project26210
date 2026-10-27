from pathlib import Path
import json,hashlib,shutil,datetime
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase192_sections28_212_figures_20260917_v1'
now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
(d/'record_backup').mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','DECISIONS.md','RESULTS_MANIFEST.tsv','CODEX_HANDOFF.md']:
 p=r/name
 if p.exists() and not (d/'record_backup'/name).exists():shutil.copy2(p,d/'record_backup'/name)
note=f"""
## phase192 | 2.8—2.12代表图 | {now}

已按用户新授权完成8张代表图（1/2/1/2/2），各有PDF、SVG及600 dpi PNG，并有中文图注。无散点图或箱线图；使用代谢反应图、密度状态图、表达/功能热图与原始分区轨迹。图源冻结，不重聚类、不拟合新轨迹、不新增推断性检验。所有PDF完成实际渲染与视觉检查；5研究原始图连通分量均未跨partition。保留组织MTA-SPP1名义显著、SPP1非C4特异、主轨迹覆盖不足等关键事实。
产物：{d}/figures；图注：{d}/FIGURE_LEGENDS_CN.md。
本地交付：Downloads/26210/22_2.8-2.12代表图_v1。当前正文仍为phase191 v3，尚未插图。
"""
for name in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 with (r/name).open('a',encoding='utf8') as f:f.write(note)
with (r/'RESULTS_MANIFEST.tsv').open('a',encoding='utf8') as f:
 f.write(f'\nphase192_representative_figures\tphase192\tfigures\t{d}/figures\t{now}\tPASS\t8 PDF rendered; SVG; 600dpi PNG; original graph checks\tindependent delivery; no new tests\n')
status={'status':'complete','time':now,'figures':8,'visual_review':'all eight; fixes rechecked; final PDF contact sheet reviewed','formats':['PDF','SVG','PNG 600dpi'],'report_inserted':False,'source_root':str(d)}
(d/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
delivery=d/'delivery';delivery.mkdir(exist_ok=True)
for f in (d/'figures').iterdir():shutil.copy2(f,delivery/f.name)
for name in ['FIGURE_LEGENDS_CN.md','README_CN.md']:shutil.copy2(d/name,delivery/name)
shutil.copy2(d/'qc/PDF_CONTACT_SHEET.jpg',delivery/'00_PREVIEW.jpg')
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in delivery.iterdir() if p.is_file() and p.name!='SHA256.json'}
(delivery/'SHA256.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2))
print('READY',len(hashes),'files')
