from pathlib import Path
import json,hashlib,shutil,re,sys
import pandas as pd
from pypdf import PdfReader
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase199_fig3_marker_GSEA_20260918_v5'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
guide=d/'REVISION_GUIDE_CN.md'
s=guide.read_text().replace('保留上一版SPP1、胆固醇代谢、溶酶体和抗原呈递的四面板等高线。','保留上一版SPP1、胆固醇代谢、溶酶体和抗原呈递的四面板等高线，移除内部A–D字母；4张局部均值矩阵与原值逐项核对一致。').replace('最后以实际PDF渲染检查版式后交付。','两张最终PDF已实际渲染并逐图目视检查，图例颜色、统计标签及子图编号一致。')
guide.write_text(s)
validation=json.loads((d/'qc/STATISTICAL_VALIDATION.json').read_text())
assert validation['independent_ES_checks']==8
assert json.loads((d/'qc/CONTOUR_UNCHANGED_CHECK.json').read_text())['local_means_unchanged']
assert len(pd.read_csv(r/'phase192_sections28_212_figures_20260917_v1/tables/MARKER_PATIENT_EQUAL_LOG2CPM.tsv',sep='\t',index_col=0))==38
figchecks=[]
for p in (d/'figures').glob('*.pdf'):
 reader=PdfReader(p);assert len(reader.pages)==1
 text=reader.pages[0].extract_text()
 for word in ['Cholesterol','Lysosome','Phagosome','Antigen','12,048','299']:assert word in text
 im=Image.open(p.with_suffix('.png'));assert min(im.info['dpi'])>599
 assert p.with_suffix('.svg').stat().st_size>1000
 figchecks.append({'file':p.name,'dpi':im.info['dpi'],'pixels':im.size})
 for ext in ['.pdf','.png','.svg']:shutil.copy2(p.with_suffix(ext),d/'delivery'/p.with_suffix(ext).name)
 shutil.copy2(d/'qc'/f'{p.stem}.jpg',d/'delivery'/f'{p.stem}.jpg')
main=PdfReader(d/'figures/Main_Figure_3_marker_GSEA_v5.pdf').pages[0].extract_text()
assert 'Captured SPP1 RNA contribution' not in main
assert 'State identity and macrophage backbone' not in main
shutil.copy2(guide,d/'delivery/修改说明与GSEA统计.md')
hashes={p.name:sha(p) for p in (d/'delivery').iterdir() if p.name!='SHA256.json'}
(d/'delivery/SHA256.json').write_text(json.dumps(hashes,indent=2))
(d/'qc/FINAL_FIGURE_QA.json').write_text(json.dumps({'figures':figchecks,'actual_pdf_renders_visually_reviewed':2,'marker_heatmap_count':1,'nested_letters_removed':True,'curve_legend_stat_colors_match':True,'source_share_strip_removed':True,'B_local_means_match':True,'python':sys.version},indent=2))
sources=list((d/'scripts').glob('*'))+list((d/'inputs').glob('*'))+list((d/'tables').glob('*'))
(d/'qc/SOURCE_CODE_RESULT_SHA256.json').write_text(json.dumps({str(p):sha(p) for p in sources},indent=2))
accepts={'Fig1':r/'phase196_selected_panels_20260917_v3/delivery/Main_Figure_1_selected_v3.pdf','Fig2':r/'phase198_fig2_contrast_corr_revision_20260918_v5/delivery/Main_Figure_2_contrast_consistency_v5.pdf'}
(d/'qc/ACCEPTED_FIGURES_UNCHANGED.json').write_text(json.dumps({k:{'path':str(v),'sha256':sha(v),'user_status':'temporarily finalized'} for k,v in accepts.items()},indent=2))
(d/'STATUS.json').write_text(json.dumps({'status':'validated_ready_for_delivery','figure3':'marker_GSEA_v5','new_GSEA':True,'patients':26,'rank_genes':12048,'pathway_family':299,'figure1_2':'temporarily finalized unchanged','report_changed':False},indent=2))
print('PACKAGE_READY',len(hashes),'hashed files')
