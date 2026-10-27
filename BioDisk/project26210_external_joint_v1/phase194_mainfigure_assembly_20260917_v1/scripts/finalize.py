from pathlib import Path
import json,hashlib,shutil,datetime
from PIL import Image
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase194_mainfigure_assembly_20260917_v1'
Image.MAX_IMAGE_PIXELS=None
qc=[]
for p in sorted((d/'figures').glob('*.png')):
 im=Image.open(p);assert min(im.info.get('dpi',(0,0)))>599
 qc.append({'file':p.name,'size':im.size,'dpi':im.info['dpi']})
for src,h in json.loads((d/'tables/SOURCE_SHA256.json').read_text()).items():
 assert hashlib.sha256(Path(src).read_bytes()).hexdigest()==h
(d/'qc/FINAL_VALIDATION.json').write_text(json.dumps({'pngs':qc,'source_files_unchanged':True,'visual_inspection':'all three final PDF renders inspected','new_statistics':False},indent=2))
delivery=d/'delivery';delivery.mkdir(exist_ok=True)
for f in (d/'figures').iterdir():shutil.copy2(f,delivery/f.name)
for f in (d/'qc').glob('Main*.jpg'):shutil.copy2(f,delivery/f.name)
for name in ['MAIN_FIGURE_GUIDE_CN.md','PANEL_INDEX_CN.md']:shutil.copy2(d/name,delivery/name)
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in delivery.iterdir() if p.name!='SHA256.json'}
(delivery/'SHA256.json').write_text(json.dumps(hashes,indent=2))
now=datetime.datetime.now().isoformat();(d/'record_backup').mkdir(exist_ok=True)
for name in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md','RESULTS_MANIFEST.tsv']:
 if (r/name).exists():shutil.copy2(r/name,d/'record_backup'/name)
note=f'\n## phase194 | 三张主图外部拼版 | {now}\n从2.1—2.12既有图件选取17组代表面板拼成3张主图（6/7/4），保留原图矢量内容；主图3纳入phase193等高线版。PDF/SVG/600dpi PNG及预览、中文拼图说明和面板来源表均已完成。三张PDF实际渲染已审阅。原数据和phase191正文未改动，源文件SHA256复核通过。未选面板及重要边界记录在MAIN_FIGURE_GUIDE_CN.md。交付路径：{d}/delivery；本地Downloads/26210/24_三张主图外部拼版_v1。\n'
for name in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 with (r/name).open('a',encoding='utf8') as f:f.write(note)
with (r/'RESULTS_MANIFEST.tsv').open('a') as f:f.write(f'\nphase194_mainfigures\tphase194\tfigures\t{d}/delivery\t{now}\tPASS\t3 vector PDFs; 17 panel groups; render inspection; SHA256\toutside report; no new analyses\n')
(d/'STATUS.json').write_text(json.dumps({'status':'complete','main_figures':3,'panel_groups':17,'report_modified':False,'formats':['PDF','SVG','PNG600dpi','JPGpreview']},indent=2))
print('READY',len(hashes),'files',sum(p.stat().st_size for p in delivery.iterdir())//1048576,'MB')
