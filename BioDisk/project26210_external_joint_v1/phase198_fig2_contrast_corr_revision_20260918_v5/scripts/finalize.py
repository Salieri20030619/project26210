from pathlib import Path
import json,hashlib,shutil,datetime,re,sys,platform
import numpy as np,pandas as pd
from pypdf import PdfReader
from PIL import Image
import matplotlib,pypdf,reportlab
Image.MAX_IMAGE_PIXELS=None
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase198_fig2_contrast_corr_revision_20260918_v5';prev=r/'phase197_top20_network_dynamics_20260917_v4'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
guide="""# Fig2 色阶、相关图及编号修订 v5

Fig1 按用户意见暂时定稿，沿用 Main_Figure_1_selected_v3；Fig3 沿用 Phase197。本次只修改 Fig2，尚未插入正文。

## A：增强同一通路内的相对差异
横轴为 F1–F7 的 cNMF 程序，分 TYMP−、TYMP+ 两组；纵轴及顺序仍为 Phase197 的 Top20。
颜色改为“每条通路内 mean AUCell 的 Z 分数”：将同一通路的两组、全部有合格配对患者的程序均值共同中心化，除以这些展示均值的样本标准差。两侧共同计算；没有分别拉伸 TYMP− 和 TYMP+。
原始患者等权 mean AUCell 和点面积均未改变。红色表示该通路内相对较高，蓝色表示相对较低；颜色不再表示跨通路绝对 AUCell 高低，也不表示 TYMP+ 相对 TYMP− 的效应量或显著性。低绝对活性通路也可呈红色。
点面积仍为 AUCell>0 的细胞比例。F3 无合格患者配对，用灰叉。程序患者数依次为6、9、0、2、3、2、2；Top20排序使用16名独立患者，原85通路全排名保留。本次不新增检验。

## B：恶性上皮轨迹
继续使用原11,950恶性上皮细胞、原UMAP及224段Monocle3主图；颜色为原全局伪时间顺序，无重新拟合。

## C：F1–通路相关及患者一致性
加入用户指定相关图的修订版，替换移出的生化示意图。
横纵轴为F1及嘧啶、半胱氨酸/甲硫氨酸、精氨酸/脯氨酸和甘油磷脂代谢。GPL为甘油磷脂缩写。
左下：20名患者内Spearman相关系数的中位数；圆面积为|rho|，蓝红为负正方向，统一色标−0.5至0.5，全部非对角原数值均在范围内。移除恒等于1的对角圆点。
右上：对应同一变量对中rho>0的患者数/20，青色比例带为正相关患者比例，灰色为其余患者；不是重复的相关系数，也不是显著患者比例。
F1对应四通路分别为：嘧啶15/20、Cys/Met15/20、Arg/Pro18/20、GPL12/20；原P/FDR直接保留于图底。
通路–通路的10对效应/一致性均由原逐患者表复核，未新增P值或重算FDR。F1检验继续沿用原四通路假设家族；GPL原P/FDR=0.49，当前数据不支持稳定的F1–GPL联系。
独立相关图同时交付，便于单独放大阅读。

## D：去除嵌套 A/B
保留Frozen macrophage states、C4 localization描述标题，删除内部A/B；复用原密度函数、坐标、冻结成员和图例。

## E：组织实测MTA及表达
继续采用原255例四队列、原患者顺序与原值。横轴为研究内按MTA排序的患者，纵轴为MTA和指定基因；颜色为队列内行Z分数。该图显示组织层面关系，不标记为C4特异关系。MTA–SPP1原P=0.02926、FDR=0.09755，仍为名义显著联系。

## 生化图独立交付
Standalone_MTA_biochemical_context 为原生化图的逐字节副本，移出完整Fig2；实线为生化反应，虚线为多步补救途径。它是生化背景，不代表本项目已经验证完整因果链。

## 交付
Main_Figure_2_contrast_consistency_v5 为新版主图。三张PDF各配SVG、600dpi PNG及JPG预览。
原始数值、排序、逐患者相关、原图与正文均保留。所有本次图件实际PDF渲染已审阅。
"""
(d/'REVISION_GUIDE_CN.md').write_text(guide)
# Validate unchanged raw values and explicit row scaling.
raw=pd.read_csv(r/'phase131_sections2_5_2_6_reference_correction_remote_v1/data/PATHWAY_DOTPLOT_SUMMARY.tsv',sep='\t')
got=pd.read_csv(d/'tables/TOP20_DISPLAY_VALUES.tsv',sep='\t')
keys=['pathway','factor','group']
a=got.set_index(keys).sort_index();b=raw.set_index(keys).loc[a.index]
assert np.allclose(a.mean_AUC,b.mean_AUC) and np.allclose(a.fraction,b.fraction)
for path,z in got.groupby('pathway'):
 assert np.allclose(z.row_z,(z.mean_AUC-z.mean_AUC.mean())/z.mean_AUC.std(ddof=1))
 assert abs(z.row_z.mean())<1e-12
assert sha(d/'tables/TOP20_RANKING_ALL85.tsv')==sha(prev/'tables/TOP20_RANKING_ALL85.tsv')
corr=pd.read_csv(d/'tables/CORRELATION_AND_CONSISTENCY.tsv',sep='\t')
assert len(corr)==10 and corr.n.eq(20).all()
eff=pd.read_csv(r/'phase131_sections2_5_2_6_reference_correction_remote_v1/data/CORRPLOT_PATIENT_EFFECTS.tsv',sep='\t')
for row in corr.itertuples():
 v=eff.loc[(eff.x==row.x)&(eff.y==row.y),'rho']
 assert row.positive==int((v>0).sum()) and abs(row.median_rho-v.median())<1e-12
main=d/'figures/Main_Figure_2_contrast_consistency_v5.pdf'
txt=PdfReader(main).pages[0].extract_text()
assert 'biochemical' not in txt.lower() and 'MAT2A' not in txt
ct=PdfReader(d/'panels/C4_no_nested_letters.pdf').pages[0].extract_text()
assert 'A   Frozen' not in ct and 'B   C4' not in ct
assert sha(d/'figures/Standalone_MTA_biochemical_context.pdf')==sha(r/'phase192_sections28_212_figures_20260917_v1/figures/01_MTA_metabolic_map.pdf')
norm=lambda x:re.sub(r'\s+','',x).lower()
sel=json.loads((d/'tables/TOP20_SELECTION.json').read_text())
assert all(norm(v) in norm(txt) for v in sel['paths'])
for p in (d/'figures').glob('*.pdf'):
 assert len(PdfReader(p).pages)==1
 im=Image.open(p.with_suffix('.png'));assert min(im.info['dpi'])>=599
 assert p.with_suffix('.svg').stat().st_size>1000
 for ext in ['.pdf','.png','.svg']:shutil.copy2(p.with_suffix(ext),d/'delivery'/p.with_suffix(ext).name)
 shutil.copy2(d/'qc'/f'{p.stem}.jpg',d/'delivery'/f'{p.stem}.jpg')
shutil.copy2(d/'REVISION_GUIDE_CN.md',d/'delivery/修改说明与图注.md')
status={'raw_AUC_and_fraction_unchanged':True,'row_z_verified':True,'top20_exact_ranking_preserved':True,'corr_pairs_verified':10,'patients_per_corr_pair':20,'new_inferential_tests':False,'C4_nested_letters_removed':True,'biochemical_map_removed_from_composite':True,'standalone_map_exact_copy':True,'PDF_visual_review':'3/3 actual PDF renders inspected','python':sys.version,'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__,'pypdf':pypdf.__version__,'reportlab':reportlab.Version}
(d/'qc/VALIDATION.json').write_text(json.dumps(status,ensure_ascii=False,indent=2))
accepted=r/'phase196_selected_panels_20260917_v3/delivery/Main_Figure_1_selected_v3.pdf'
(d/'qc/FIG1_PROVISIONAL_ACCEPTANCE.json').write_text(json.dumps({'user_decision':'Fig1 temporarily finalized; preserve existing version','path':str(accepted),'sha256':sha(accepted)},indent=2))
inputs=[r/'phase131_sections2_5_2_6_reference_correction_remote_v1/data'/n for n in ['PATHWAY_DOTPLOT_SUMMARY.tsv','CORRPLOT_SUMMARY.tsv','CORRPLOT_PATIENT_EFFECTS.tsv']]
inputs += [r/'phase192_sections28_212_figures_20260917_v1/inputs'/n for n in ['MACRO_DISPLAY_METADATA.tsv.gz','FEATURE_REGISTRY.json']]
inputs += list((d/'scripts').glob('*.py'))+list((d/'tables').glob('*'))
(d/'qc/INPUT_CODE_TABLE_SHA256.json').write_text(json.dumps({str(p):sha(p) for p in inputs},indent=2))
(d/'delivery/SHA256.json').write_text(json.dumps({p.name:sha(p) for p in (d/'delivery').iterdir()},indent=2))
(d/'STATUS.json').write_text(json.dumps({'status':'validated_ready_for_local_delivery','Fig1':'temporarily finalized unchanged','Fig2':'v5 ready','Fig3':'unchanged','report_changed':False},indent=2))
print('VALIDATED',len(list((d/'delivery').iterdir())),'files',flush=True)
