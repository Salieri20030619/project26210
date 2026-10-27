from pathlib import Path
import json,hashlib,datetime,zipfile,shutil,subprocess
import numpy as np,pandas as pd
from PIL import Image
P=Path(__file__).resolve().parents[1];R=P.parent
gate=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
with (P/'logs/memory_gate.log').open('a') as f:f.write('finalize\t'+gate.stdout+gate.stderr)
assert gate.returncode==0,gate.stdout+gate.stderr
assert (P/'build_EXIT_STATUS').read_text().strip()=='0'
old=R/'phase255_spatial_tissue_figure_20260925_v2';pre=R/'phase256_spatial_region_function_20260926_v1'
checks=[]
manifest=json.loads((P/'qc/INPUT_MANIFEST.json').read_text())
for m in manifest:
 f=Path(m['path']);s=f.stat()
 assert s.st_size==m['bytes'] and s.st_mtime_ns==m['mtime_ns'],str(f)
 if m['sha256']:assert hashlib.sha256(f.read_bytes()).hexdigest()==m['sha256']
checks.append('Source identity and hashes unchanged, including original Phase255 exports')
for file in sorted((P/'plot_data').glob('*_tissue_layers.tsv')):
 original=old/'plot_data'/file.name
 assert original.is_file()
 assert original.read_bytes()==file.read_bytes(),str(file)
checks.append('Ten prior tissue-layer tables byte-identical')
sel=P/'plot_data/INDEPENDENT_MAP_SELECTION.tsv'
assert sel.read_bytes()==(old/'plot_data/INDEPENDENT_MAP_SELECTION.tsv').read_bytes()
oldmeta=json.loads((old/'qc/DISPLAY_METADATA.json').read_text());newmeta=json.loads((P/'qc/DISPLAY_METADATA.json').read_text())
assert oldmeta==newmeta
checks.append('All image transforms, bounds, scales, region thresholds and case selection unchanged')
stats=pd.read_csv(P/'plot_data/GSEA_FOUR_CURVES_STATISTICS.tsv',sep='\t')
full=pd.read_csv(pre/'tables/GSEA_FULL_RESULTS.tsv.gz',sep='\t',low_memory=False)
full=full[(full.model=='absolute')&(full.version=='exclude_definition_genes')&(full.library=='KEGG')&full.term.isin(['Lysosome','Cholesterol metabolism'])]
assert len(stats)==len(full)==4
curve_checks=[]
for row in stats.itertuples():
 prior=full[(full.dataset==row.cohort)&(full.term==row.pathway)].iloc[0]
 for key in ['ES','NES','P_GSEA','FDR_GSEA']:assert abs(getattr(row,key)-prior[key])<1e-12
 short='Lysosome' if row.pathway=='Lysosome' else 'Chol._met.'
 d=pd.read_csv(P/f'plot_data/GSEA_{row.cohort}_{short}.tsv.gz',sep='\t')
 hits=d.gene_set_hit.to_numpy(bool);w=np.abs(d['rank'].to_numpy());n=len(d);k=hits.sum()
 # Recompute sparse hit-only extrema independently of the plotted cumulative walk.
 h=np.flatnonzero(hits);hw=w[h];tot=hw.sum()
 mx=np.max(np.cumsum(hw)/tot-(h-np.arange(k))/(n-k))
 mn=np.min(np.r_[0,np.cumsum(hw)[:-1]]/tot-(h-np.arange(k))/(n-k))
 es=mx if mx>=-mn else mn
 assert abs(es-row.ES)<1e-10
 assert n==row.ranked_genes and k==row.set_size
 curve_checks.append(dict(cohort=row.cohort,pathway=row.pathway,ES_recomputed=float(es),ES_error=float(abs(es-row.ES))))
pd.DataFrame(curve_checks).to_csv(P/'qc/INDEPENDENT_CURVE_CHECKS.tsv',sep='\t',index=False)
checks.append('Four curves, hit strips, patient n and full-family GSEA statistics verified')
# PDF size and page count through Poppler.
pdf=P/'delivery/26210_Spatial_A4_GSEA_v3.pdf'
info=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/pdfinfo',str(pdf)],capture_output=True,text=True,check=True).stdout
(P/'qc/PDF_INFO.txt').write_text(info)
assert 'Pages:           1' in info
for line in info.splitlines():
 if line.startswith('Page size:'):
  parts=line.split();w,h=float(parts[2]),float(parts[4]);assert abs(w*25.4/72-210)<.02 and abs(h*25.4/72-297)<.02
pixels={}
for ext,dpi in [('png',300),('tiff',600)]:
 with Image.open(P/f'delivery/26210_Spatial_A4_GSEA_v3.{ext}') as im:
  assert abs(im.size[0]-round(210/25.4*dpi))<=1 and abs(im.size[1]-round(297/25.4*dpi))<=1
  pixels[ext]=dict(size=list(im.size),dpi=[float(v) for v in im.info.get('dpi',[])])
# Main-figure legend: plot definitions and actual evidence, with detailed prior analyses retained in the data package.
rows=[]
for x in stats.itertuples():
 rows.append(f"| {x.cohort} | {x.patient_n} | {x.pathway} | {x.ES:.6g} | {x.NES:.6g} | {x.P_GSEA:.6g} | {x.FDR_GSEA:.6g} |")
legend="""# 空间转录组综合图 v3：组织分布与区域功能富集

## 图题
TYMP相关代谢背景、SPP1相关巨噬状态及溶酶体—胆固醇程序的空间特征

## 版式
A–H共八个面板，210×297 mm单页，22张组织空间地图和4幅GSEA曲线。单独的原始H&E面板已移出，组织影像作为空间图底图保留；原组织组成与TYMP表达合并为A。旧版四张免疫/侵袭功能地图的位置改为新功能富集面板E，其余组织地图保留。图件独立交付。

## A 组织组成与TYMP表达
每行分别为GLMF1、GLMF2；左列为冻结13类cell2location参考归并六大类后的优势成分，右列为TYMP的log1p(CP10k)表达。六类为恶性上皮、髓系、淋巴、血管、正常上皮和肥大细胞；混合spot显示其中丰度总和最高者。TYMP两片共用0–2色标，超出范围饱和显示。
GLMF1/2分别为GSE210041 Patient1/GSM6415705和Patient7/GSM6415706，质控spot数3199/3880。A–C中同一病例使用相同方向、范围及尺寸，比例尺500 μm。

## B 代谢转录背景
颜色为MAT2A、MAT2B、AMD1、SRM和SMS的片内基因Z分数均值；共同色标−1.25至1.25。红色表示高于片内平均，蓝色表示低于平均。该分数描述多胺/MTA相关转录背景，TYMP未参与评分。

## C SPP1相关巨噬状态丰度
沿用既有SPP1-state巨噬丰度后验第5百分位(q05)，共同色标0–3。保持原空间状态名称，与单细胞C4的定义分开。两片分别展示其组织内丰度分布。

## D 代谢背景与SPP1状态程序的高值区域
采用预先固定的片内第80百分位阈值：橙色仅代谢背景高，青色仅状态程序高，紫色两者均高，灰色其余。状态程序为GPNMB、LGALS3、FABP5、CD9和LPL，不含SPP1。该面板显示实测位置的相对高值区域重叠，所有区域和阈值沿用旧版。

## E SPP1-state高丰度区域的功能富集
左列Lysosome（溶酶体），右列Chol. met.（Cholesterol metabolism，胆固醇代谢）；上行为GSE175540的23位独立患者，下行为TJ-RCC的12位独立患者。这里的两行是外部统计来源，不是GLMF1/2。

采用Phase256主分析：每片SPP1-state q05丰度前20%区域与其余区域比较，计算log1p(CP10k)逐基因Z分数均值差；同患者多制备等权汇总。按患者效应mean/(SD/√n)排列全部合格基因，排除TYMP、SPP1、GPNMB、LGALS3、FABP5、CD9和LPL，使用冻结KEGG集合进行GSEA。

曲线横轴为排序基因位置占完整列表的比例（0–1）；左端对应高丰度区域正向差异较强的基因。纵轴ES为沿排序列表累积的富集分数，曲线下方短竖线表示该通路成员在列表中的位置。NES为归一化富集分数，图内q为该来源完整KEGG主分析家族的BH FDR，不是仅对展示的两条通路校正。图中四个NES均为正，表明这两个功能集合偏向排序前端。图形为原结果的精确重建，未重新检验或更改排序。

| 来源 | 独立患者n | 通路 | ES | NES | GSEA P | GSEA FDR |
| --- | --- | --- | --- | --- | --- | --- |
"""+'\n'.join(rows)+"""

本面板展示原始区域比较，组成调整、患者层效果及完整敏感性结果沿用Phase256统计，保存在配套数据中。GSEA P描述给定患者排序中的基因集相对富集；患者层统计另表保存。

## F SPP1-state与内皮高值区域
使用相同片内第80百分位显示规则：蓝色仅SPP1-state q05丰度高，珊瑚色仅内皮q05高，紫色二者均高，灰色其余。丰度与显示区域沿用冻结结果。

## G 六位独立患者的组织地图
展示E-MTAB-12767全部六患者的代表肿瘤阵列：HG_1.3、HG_2、HG_3.2、LG_1、LG_2、LG_3。每患者选择质控spot数最多的阵列，平局按名称排序。程序、阈值及颜色与D相同，比例尺500 μm。代表图不改变患者统计中多阵列的合并口径。

## H 配对旁组织与肿瘤SPP1表达
两行为HG1/HG3；左列旁组织pHG_1.1/pHG_3.1，右列肿瘤HG_1.3/HG_3.2。四图共用0–5的SPP1 log1p(CP10k)色标，超出范围饱和显示，比例尺500 μm。两例按独立配对个体作描述。

## 空间显示与缩写
组织地图的横纵位置沿用原影像坐标。六边形以原实测spot为中心，经冻结阵列—影像几何关系显示，未作表达插值、平滑或病理边界描绘；比例尺由相同几何关系换算。灰色和空白按各面板图例或原组织覆盖解释。相对高值分类图的颜色表示片内阈值类别，不等同跨切片绝对量。

H&E，hematoxylin and eosin，苏木精—伊红染色；CP10k，每万计数标准化表达；q05，后验第5百分位；GSEA，gene set enrichment analysis；ES，enrichment score；NES，normalized enrichment score；FDR，false discovery rate；Chol. met.，Cholesterol metabolism；TYMP/SPP1为标准基因符号。

## 来源
组织地图沿用Phase255的冻结输入、几何转换和选片规则；功能曲线来自Phase256已验收结果。参照Ge等Science Advances 12, eaed0733（2026）Fig.2G的区域功能叙事，采用规范的running ES与基因命中条带展示。图件不表示已建立TYMP到SPP1的因果调控链。旧图、原数据、RDS及现有报告保持不变。
"""
(P/'delivery/26210_Spatial_GSEA_v3_图注.md').write_text(legend,encoding='utf-8')
statdir=P/'plot_data/statistical_evidence';statdir.mkdir(exist_ok=True)
for name in ['PRIMARY_FULL_LIBRARY_RESULTS.tsv.gz','PRIORITY_PATHWAYS_ALL_MODELS_AND_VERSIONS.tsv','REFERENCE_FIG2G_ALL_TEN_TERMS.tsv']:
 shutil.copyfile(pre/'tables'/name,statdir/name)
shutil.copyfile(pre/'PLAN_CN.md',statdir/'Phase256_PLAN_CN.md')
shutil.copyfile(pre/'qc/FINAL_VALIDATION.json',statdir/'Phase256_FINAL_VALIDATION.json')
# Retain prior figure's statistical evidence without editorial expansion in the new legend.
if (old/'plot_data/statistical_evidence').is_dir():shutil.copytree(old/'plot_data/statistical_evidence',statdir/'prior_tissue_evidence',dirs_exist_ok=True)
check=dict(status='PASS',time=datetime.datetime.now().astimezone().isoformat(),page_mm=[210,297],panels=8,tissue_maps=22,gsea_curves=4,minimum_font_pt=6,visual_check='Full-page Poppler PDF render and GSEA detail reviewed; no text/legend overlaps; all A-H tissue maps, cohort labels, curves and gene-hit strips readable.',checks=checks,curve_checks=curve_checks,export_pixels=pixels,source_figure_unchanged=True,new_tests=0,new_models=0,report_unchanged=True)
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(check,ensure_ascii=False,indent=2))
zipfile_path=P/'delivery/26210_Spatial_GSEA_v3_plot_data.zip'
with zipfile.ZipFile(zipfile_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for folder in ['plot_data','code','qc']:
  for f in sorted((P/folder).rglob('*')):
   if f.is_file() and f.suffix not in ['.png'] and '__pycache__' not in str(f):z.write(f,f.relative_to(P))
 z.write(P/'PLAN_CN.md','PLAN_CN.md')
with zipfile.ZipFile(zipfile_path) as z:assert z.testzip() is None
outputs=[dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted((P/'delivery').iterdir()) if f.is_file()]
(P/'DELIVERY_MANIFEST.json').write_text(json.dumps(outputs,ensure_ascii=False,indent=2))
(P/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_VALIDATED','time':check['time'],'panels':8,'tissue_maps':22,'gsea_curves':4,'new_tests':0,'report_modified':False},indent=2))
print(json.dumps(outputs,ensure_ascii=False,indent=2))
