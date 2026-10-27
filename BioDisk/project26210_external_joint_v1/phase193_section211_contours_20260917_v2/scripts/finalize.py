from pathlib import Path
import json,hashlib,shutil,datetime
from PIL import Image
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase193_section211_contours_20260917_v2'
im=Image.open(d/'figures/06_C4_expression_contours_v2.png');assert min(im.info['dpi'])>599
now=datetime.datetime.now().isoformat()
(d/'record_backup').mkdir(exist_ok=True)
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md','RESULTS_MANIFEST.tsv']:
 if (r/n).exists():shutil.copy2(r/n,d/'record_backup'/n)
note=f'\n## phase193 | 2.11等高线修订 | {now}\n用户指出原06图紫色连续填色区分不足。已改为7层等距表达等高线、冷暖填色、紫色虚线C4边界；显示主区间为各面板可见网格P2-P98，两端延伸。原始数据、局部均值和统计不变。PDF实际渲染已目视检查。替换phase192的06图，05热图不变。产物：{d}/figures；本地：Downloads/26210/23_2.11等高线修订_v2。\n'
for n in ['PROJECT_STATUS.md','DECISIONS.md','CODEX_HANDOFF.md']:
 with (r/n).open('a',encoding='utf8') as f:f.write(note)
with (r/'RESULTS_MANIFEST.tsv').open('a') as f:f.write(f'\nphase193_contours\tphase193\tfigure\t{d}/figures\t{now}\tPASS\tPDF render reviewed; PNG600dpi\tunchanged local means; 7 contour bands\n')
(d/'STATUS.json').write_text(json.dumps({'status':'complete','visual_review':'PDF rendered and inspected','input_sha256':hashlib.sha256((r/'phase192_sections28_212_figures_20260917_v1/inputs/MACRO_DISPLAY_METADATA.tsv.gz').read_bytes()).hexdigest()},indent=2))
delivery=d/'delivery';delivery.mkdir(exist_ok=True)
for p in (d/'figures').iterdir():shutil.copy2(p,delivery/p.name)
shutil.copy2(d/'FIGURE_LEGEND_CN.md',delivery/'FIGURE_LEGEND_CN.md')
(delivery/'SHA256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in delivery.iterdir() if p.name!='SHA256.json'},indent=2))
print('Delivery ready')
