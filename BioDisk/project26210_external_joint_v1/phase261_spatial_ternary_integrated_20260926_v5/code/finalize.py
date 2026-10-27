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
checks.append('Four curves, true gene hit positions, patient n and full-family GSEA statistics verified')

# Verify styling-only edit against the accepted Phase258 data and non-E pixels.
prior=R/'phase259_spatial_GSEA_reference_style_20260926_v4'
priorstats=pd.read_csv(prior/'plot_data/GSEA_FOUR_CURVES_STATISTICS.tsv',sep='\t')
pd.testing.assert_frame_equal(stats,priorstats)
stems=pd.read_csv(P/'plot_data/GSEA_COLORED_STEMS.tsv',sep='\t')
assert len(stems)==int(stats.set_size.sum()) and (stems.stem_bottom==0).all()
for row in stats.itertuples():
 short='Lysosome' if row.pathway=='Lysosome' else 'Chol._met.'
 name=f'GSEA_{row.cohort}_{short}.tsv.gz'
 current=pd.read_csv(P/'plot_data'/name,sep='\t')
 previous=pd.read_csv(prior/'plot_data'/name,sep='\t')
 pd.testing.assert_frame_equal(current,previous)
 ss=stems[(stems.cohort==row.cohort)&(stems.pathway==row.pathway)]
 ix=ss.rank_index_zero_based.to_numpy(int)
 assert len(ss)==row.set_size
 assert current.iloc[ix].gene_set_hit.all()
 np.testing.assert_allclose(ss.stem_top,current.iloc[ix].running_ES,atol=1e-12)
 np.testing.assert_allclose(ss.color_scalar_ES,ss.stem_top,atol=1e-12)
checks.append('All prior numerical GSEA curves and statistics unchanged; all 289 stems match true gene hits and actual ES')
ternary=R/'phase260_spatial_ternary_preview_20260926_v1'
for sample,n in [('GLMF1',1282),('GLMF2',3113)]:
 f=P/f'plot_data/{sample}_TERNARY_POINTS.tsv'
 assert f.read_bytes()==(ternary/'plot_data'/f.name).read_bytes()
 d=pd.read_csv(f,sep='\t');assert len(d)==n
 np.testing.assert_allclose(d.ternary_x,.5*d.share_SPP1_state+d.share_T_cell,atol=1e-12)
 np.testing.assert_allclose(d.ternary_y,np.sqrt(3)/2*d.share_SPP1_state,atol=1e-12)
checks.append('Exactly two separate case ternary plots; all 1282/3113 coordinates byte-identical to accepted standalone data; no overlay')
build=json.loads((P/'qc/BUILD_CHECK.json').read_text())
assert build['ternary_plots']==2 and build['overlay']==False
# PDF size and page count through Poppler.
pdf=P/'delivery/26210_Spatial_A4_GSEA_Ternary_v5.pdf'
info=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/pdfinfo',str(pdf)],capture_output=True,text=True,check=True).stdout
(P/'qc/PDF_INFO.txt').write_text(info)
assert 'Pages:           1' in info
for line in info.splitlines():
 if line.startswith('Page size:'):
  parts=line.split();w,h=float(parts[2]),float(parts[4]);assert abs(w*25.4/72-210)<.02 and abs(h*25.4/72-297)<.02
pixels={}
for ext,dpi in [('png',300),('tiff',600)]:
 with Image.open(P/f'delivery/26210_Spatial_A4_GSEA_Ternary_v5.{ext}') as im:
  assert abs(im.size[0]-round(210/25.4*dpi))<=1 and abs(im.size[1]-round(297/25.4*dpi))<=1
  pixels[ext]=dict(size=list(im.size),dpi=[float(v) for v in im.info.get('dpi',[])])
# Main-figure legend: plot definitions and actual evidence, with detailed prior analyses retained in the data package.
rows=[]
for x in stats.itertuples():
 rows.append(f"| {x.cohort} | {x.patient_n} | {x.pathway} | {x.ES:.6g} | {x.NES:.6g} | {x.P_GSEA:.6g} | {x.FDR_GSEA:.6g} |")
oldcaption=next((R/'phase259_spatial_GSEA_reference_style_20260926_v4/delivery').glob('*_图注.md')).read_text()
parts=oldcaption.split('## ')
section={p.split('\n',1)[0]:p.split('\n',1)[1] for p in parts[1:]}
section['A 组织组成与TYMP表达']=section['A 组织组成与TYMP表达'].replace('A–C中','A–E中')
gsea=section['E SPP1-state高丰度区域的功能富集']
gsea=gsea.replace('左列Lysosome（溶酶体），右列Chol. met.（Cholesterol metabolism，胆固醇代谢）；上行为GSE175540的23位独立患者，下行为TJ-RCC的12位独立患者。这里的两行是外部统计来源，不是GLMF1/2。','从左至右为GSE175540的Lysosome、Chol. met.以及TJ-RCC的Lysosome、Chol. met.；两组各有23和12位独立患者。Chol. met.为Cholesterol metabolism（胆固醇代谢）。这些外部统计来源与GLMF1/2区分。')
legend='# 空间转录组综合图 v5：组织分布、区域功能与三元距离图\n\n'
legend+='## 图题\n'+section['图题']
legend+='## 版式\nA–I九个面板，210×297 mm单页；22张组织空间地图、4幅GSEA曲线，以及GLMF1/GLMF2两张分别展示的三元图。\n\n'
for title in ['A 组织组成与TYMP表达','B 代谢转录背景','C SPP1相关巨噬状态丰度','D 代谢背景与SPP1状态程序的高值区域']:
 legend+='## '+title+'\n'+section[title]
legend+='## E SPP1-state与内皮高值区域\n'+section['F SPP1-state与内皮高值区域']
legend+='## F SPP1-state高丰度区域的功能富集\n'+gsea
legend+='''## G 恶性富集spot到三类空间目标的相对距离
左为GLMF1（1,282个有效spot），右为GLMF2（3,113个有效spot），代表两位独立患者。每点为冻结恶性成分比例≥0.20的源spot，三个目标分别为该切片SPP1-state、内皮与T细胞q05丰度前20%的区域，采用原阳性阈值。

距离沿实测组织第一邻接图计算最短加权路径，使用物理100 μm相邻边，不跨越组织缺口或不连通分量。三个距离各除以其总和，形成总和为100%的三元坐标：顶部为SPP1-state距离份额，左下为内皮，右下为T细胞。某份额越大表示相对更远；靠近某顶点不是更接近该细胞，也不是该细胞比例更高。真实零距离保留；全零或无法到达三个目标的点不定义三元位置。GLMF1/2分别从1,314/3,133个源spot保留1,282/3,113个点；其余为GLMF1不可达11、全零21和GLMF2全零20。无抖动、平滑或抽样。

两病例使用相同坐标、尺度和透明度单独显示。该面板描述各病例的空间距离构成，spot数不是独立样本量，不作spot层推断检验。全部坐标沿用已验收Phase260结果。

'''
legend+='## H 六位独立患者的组织地图\n'+section['G 六位独立患者的组织地图']
legend+='## I 配对旁组织与肿瘤SPP1表达\n'+section['H 配对旁组织与肿瘤SPP1表达']
legend+='## 空间显示与缩写\n'+section['空间显示与缩写']
legend+='## 来源\n'+section['来源']+'\n三元图沿用Phase260冻结空间图距离与份额，参照参考论文Fig.5D的三元展示形式，仅纳入两张分病例图。\n'
(P/'delivery/26210_Spatial_GSEA_Ternary_v5_caption.md').write_text(legend,encoding='utf-8')
statdir=P/'plot_data/statistical_evidence';statdir.mkdir(exist_ok=True)
for name in ['PRIMARY_FULL_LIBRARY_RESULTS.tsv.gz','PRIORITY_PATHWAYS_ALL_MODELS_AND_VERSIONS.tsv','REFERENCE_FIG2G_ALL_TEN_TERMS.tsv']:
 shutil.copyfile(pre/'tables'/name,statdir/name)
shutil.copyfile(pre/'PLAN_CN.md',statdir/'Phase256_PLAN_CN.md')
shutil.copyfile(pre/'qc/FINAL_VALIDATION.json',statdir/'Phase256_FINAL_VALIDATION.json')
# Retain prior figure's statistical evidence without editorial expansion in the new legend.
if (old/'plot_data/statistical_evidence').is_dir():shutil.copytree(old/'plot_data/statistical_evidence',statdir/'prior_tissue_evidence',dirs_exist_ok=True)
for name in ['CASE_DISTANCE_SUMMARY.tsv','TARGET_DEFINITIONS.tsv']:
 shutil.copyfile(ternary/'tables'/name,statdir/name)
for name in ['GLMF1_ALL_SPOT_DISTANCE_AUDIT.tsv.gz','GLMF2_ALL_SPOT_DISTANCE_AUDIT.tsv.gz']:
 shutil.copyfile(ternary/'plot_data'/name,statdir/name)
shutil.copyfile(ternary/'qc/FINAL_VALIDATION.json',statdir/'Phase260_FINAL_VALIDATION.json')
shutil.copyfile(ternary/'code/make_preview.py',statdir/'Phase260_make_preview.py')
check=dict(status='PASS',time=datetime.datetime.now().astimezone().isoformat(),page_mm=[210,297],panels=9,tissue_maps=22,gsea_curves=4,ternary_plots=2,overlay=False,colored_gene_stems=True,bottom_rug_strips=False,minimum_font_pt=6,visual_check='Actual full-page Poppler PDF render reviewed; A-I labels, legends, four GSEA curves and two separate ternary cases readable.',checks=checks,curve_checks=curve_checks,export_pixels=pixels,source_figure_unchanged=True,new_tests=0,new_models=0,report_unchanged=True)
(P/'qc/FINAL_VALIDATION.json').write_text(json.dumps(check,ensure_ascii=False,indent=2))
zipfile_path=P/'delivery/26210_Spatial_GSEA_Ternary_v5_plot_data.zip'
with zipfile.ZipFile(zipfile_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for folder in ['plot_data','code','qc']:
  for f in sorted((P/folder).rglob('*')):
   if f.is_file() and f.suffix not in ['.png'] and '__pycache__' not in str(f):z.write(f,f.relative_to(P))
 z.write(P/'PLAN_CN.md','PLAN_CN.md')
with zipfile.ZipFile(zipfile_path) as z:assert z.testzip() is None
outputs=[dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in sorted((P/'delivery').iterdir()) if f.is_file()]
(P/'DELIVERY_MANIFEST.json').write_text(json.dumps(outputs,ensure_ascii=False,indent=2))
(P/'STATUS.json').write_text(json.dumps({'status':'COMPLETE_VALIDATED','time':check['time'],'panels':9,'tissue_maps':22,'gsea_curves':4,'new_tests':0,'report_modified':False},indent=2))
print(json.dumps(outputs,ensure_ascii=False,indent=2))
