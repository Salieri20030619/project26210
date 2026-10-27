from pathlib import Path
import json,datetime,hashlib,zipfile,subprocess
import numpy as np,pandas as pd
from PIL import Image
P=Path(__file__).resolve().parents[1]
assert (P/'build_EXIT_STATUS').read_text().strip()=='0'
gate=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
with (P/'logs/memory_gate.log').open('a') as f:f.write('finalize\t'+gate.stdout+gate.stderr)
assert gate.returncode==0
manifest=json.loads((P/'qc/INPUT_MANIFEST.json').read_text())
for m in manifest:
 p=Path(m['path']);s=p.stat()
 assert s.st_size==m['bytes'] and s.st_mtime_ns==m['mtime_ns'] and hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256']
summary=pd.read_csv(P/'tables/CASE_DISTANCE_SUMMARY.tsv',sep='\t')
for r in summary.itertuples():
 a=pd.read_csv(P/'plot_data'/f'{r.sample}_ALL_SPOT_DISTANCE_AUDIT.tsv.gz',sep='\t')
 d=pd.read_csv(P/'plot_data'/f'{r.sample}_TERNARY_POINTS.tsv',sep='\t')
 assert len(a)==r.QC_spots and a.spot_id.is_unique
 assert len(d)==r.plotted_spots and d.spot_id.is_unique
 assert r.source_spots==r.plotted_spots+r.unreachable_sources+r.allzero_sources
 assert set(d.spot_id)==set(a.loc[a.status=='PLOTTED','spot_id'])
 cols=['SPP1_state','Endothelium','T_cell']
 distance=d[['distance_'+c+'_um' for c in cols]].to_numpy()
 weights=d[['share_'+c for c in cols]].to_numpy()
 np.testing.assert_allclose(weights,distance/distance.sum(1)[:,None],atol=1e-12)
 np.testing.assert_allclose(weights.sum(1),1,atol=1e-12)
 np.testing.assert_allclose(d.ternary_x,weights[:,2]+weights[:,0]/2,atol=1e-12)
 np.testing.assert_allclose(d.ternary_y,weights[:,0]*np.sqrt(3)/2,atol=1e-12)
info=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/pdfinfo',str(P/'delivery/26210_Ternary_spatial_preview_v1.pdf')],capture_output=True,text=True,check=True).stdout
assert 'Pages:           1' in info
(P/'qc/PDF_INFO.txt').write_text(info)
pixels={}
for ext,dpi in [('png',300),('tiff',600)]:
 with Image.open(P/f'delivery/26210_Ternary_spatial_preview_v1.{ext}') as im:
  assert abs(im.size[0]-round(160/25.4*dpi))<=1 and abs(im.size[1]-round(170/25.4*dpi))<=1
  pixels[ext]=dict(size=list(im.size),dpi=[float(v) for v in im.info.get('dpi',[])])
rows=[]
for r in summary.itertuples():
 rows.append(f"| {r.sample} | {r.source_spots} | {r.plotted_spots} | {r.unreachable_sources} | {r.allzero_sources} |")
caption="""# 独立三元空间距离图：GLMF1与GLMF2

## 图题
恶性表达区域到SPP1-state、内皮及T细胞相关高值区域的相对距离构成。

## 布局与读图
上方叠加GLMF1（蓝色）与GLMF2（珊瑚色），下方按病例分别显示。同一病例所有合格观测点均保留；上方只以固定随机顺序交错绘制以减小覆盖顺序影响，没有抽样、抖动或平滑。图中透明度用于显示重叠，许多相同比例位置可有多个观测点重合。
每个点代表一个恶性表达区域中的spot，并非一个分割细胞。三轴为距离的相对构成，不是细胞比例：右斜轴SPP1-state、左斜轴Endoth.（Endothelium，内皮）和底轴T cells。沿坐标值增大的方向移动，表示到相应目标区域的距离占总距离比例更大，即相对更远；不能将靠近某顶点读成更靠近该类细胞。

令三个距离为dS、dE、dT，则坐标为dS/(dS+dE+dT)、dE/(dS+dE+dT)、dT/(dS+dE+dT)，每个点的三项之和为100%。图中几何顶点依次为上方SPP1-state、左下内皮、右下T细胞。相同比例可以对应不同绝对距离；每点的三个原始距离均保留于数据包。

## 输入与距离定义
GLMF1/GLMF2分别对应GSE210041的Patient1（GSM6415705）和Patient7（GSM6415706）。复用Phase234的冻结spot、坐标和cell2location映射。参照区域采用恶性成分比例>=0.20；三类目标分别为本片SPP1-state q05、内皮q05和T细胞q05位于第80百分位及以上的位置。
沿原first-ring组织邻接图，以实际物理坐标换算边长（中位100 μm），计算每个参照spot到各目标区域的最短图距离，不跨组织分量或空洞。参照点本身属于目标区域时距离记0。旧SPP1-state与单细胞C4保持各自定义，T细胞不改称CD8，内皮不改称成纤维细胞。

## 可绘图覆盖
| 病例 | 参照spot | 绘图spot | 至少一目标不可达 | 三项距离均为零 |
| --- | --- | --- | --- | --- |
"""+'\n'.join(rows)+"""

不可达目标的点没有有限完整三元坐标；三个距离均为零的点没有可定义分母，因此两类点单独保留在审计表，未添加任意小距离。来源、目标阈值及每点状态均可复核。全部病例与阈值在看图前固定。

## 统计定位
这是两名独立患者的描述性预览。图内数量是绘图spot数，不是生物学样本量；本轮没有进行spot层或患者层显著性检验。它展示的是表达/映射区域之间的距离构成，没有使用TIB/非TIB分组，也不承担单细胞接触或免疫排斥的推断。

## 验收与来源
scipy多源最短路与独立heapq实现逐位置一致；距离占比、三元坐标反算、输入哈希和绘图数量核对通过。整页PDF已经实际渲染并检查，坐标标签和病例图例清楚，无遮挡。
参考论文Ge等Science Advances 12, eaed0733（2026）Fig.5D的三元展示布局，数据与距离定义使用本项目冻结输入。图为160×170 mm独立文件，暂未加入当前空间大图；当前主图仍为Phase259 v4。
"""
(P/'delivery/26210_Ternary_spatial_preview_v1_caption.md').write_text(caption,encoding='utf-8')
validation=dict(status='PASS',time=datetime.datetime.now().astimezone().isoformat(),independent_patients=2,source_spots=int(summary.source_spots.sum()),plotted_spots=int(summary.plotted_spots.sum()),unreachable=int(summary.unreachable_sources.sum()),allzero=int(summary.allzero_sources.sum()),independent_shortest_path_checks=6,normalization='distance/sum_three_distances; no pseudocount',vertices={'SPP1_state':[.5,float(np.sqrt(3)/2)],'Endothelium':[0,0],'T_cell':[1,0]},axis_tick_meanings={'left':'Endothelium,100 at lower-left to0 at top','right':'SPP1-state,100 at top to0 at lower-right','bottom':'T-cell,0 at lower-left to100 at lower-right'},visual_check='Full-page PDF render reviewed; overlay and individual panels, grid percentages, rotated axis names, case labels and distance-share interpretation readable without overlap.',figure_mm=[160,170],export_pixels=pixels,smoothing=False,jitter=False,inferential_tests=0,input_files_unchanged=len(manifest),added_to_main_figure=False)
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2))
out=P/'delivery/26210_Ternary_spatial_preview_v1_data.zip'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for folder in ['plot_data','tables','code','qc']:
  for p in sorted((P/folder).rglob('*')):
   if p.is_file() and p.suffix!='.png' and '__pycache__' not in str(p):z.write(p,p.relative_to(P))
 z.write(P/'PLAN_CN.md','PLAN_CN.md')
with zipfile.ZipFile(out) as z:assert z.testzip() is None
files=[dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted((P/'delivery').iterdir()) if f.is_file()]
(P/'DELIVERY_MANIFEST.json').write_text(json.dumps(files,indent=2))
(P/'STATUS.json').write_text(json.dumps(dict(status='COMPLETE_VALIDATED',time=validation['time'],added_to_main_figure=False),indent=2))
print(json.dumps(files,indent=2))
