from pathlib import Path
import ast,json,subprocess
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase196_selected_panels_20260917_v3'
for s in ['panels','figures','tables','qc']:(d/s).mkdir(exist_ok=True)
o=r/'phase131_sections2_5_2_6_reference_correction_remote_v1'
ds=pd.read_csv(o/'data/PATHWAY_DOTPLOT_SUMMARY.tsv',sep='\t')
valid=json.loads((o/'qc/SECTION26_VALIDATION.json').read_text())
factors=['F'+str(i) for i in range(1,8)];npats=[valid['paired_patient_counts_per_program'][s] for s in factors]
b=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/09_metabolism_trajectory')
paths=pd.read_csv(b/'metabolism_v1/source_data/KEGG_AUCell_all_pathways_v1.tsv.gz',sep='\t',nrows=0).columns.tolist()[1:]
assert len(paths)==85
s=(o/'scripts/draw_section26_reference.py').read_text()
node=next(n for n in ast.parse(s).body if isinstance(n,ast.FunctionDef) and n.name=='dot_facets')
fn=ast.get_source_segment(s,node)
fn=fn.replace("(16,22)","(8.3,14.5)").replace("left=.405","left=.44").replace("pw=.22","pw=.19").replace("bottom=.075","bottom=.05").replace("top=.93","top=.956").replace("maxsize=118","maxsize=35").replace("fontsize=8 if full","fontsize=6.5 if full").replace("fontsize=7 if full","fontsize=6.7 if full").replace("fontsize=11)","fontsize=8)").replace("fontsize=9)","fontsize=7)").replace("fontsize=8)","fontsize=6.5)").replace("fontsize=8,title_fontsize=9","fontsize=6.5,title_fontsize=7")
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'pdf.fonttype':42})
def save(fig,name,title,foot,reference):
 fig.savefig(d/'panels/85_pathways_column.pdf',bbox_inches='tight',pad_inches=.06)
 plt.close(fig)
exec(compile(fn,'original_dotplot_function_adapted','exec'))
dot_facets(paths,'85_pathways_column','')
ds.to_csv(d/'tables/DOTPLOT_ORIGINAL_SUMMARY.tsv',sep='\t',index=False)
(d/'tables/DOTPLOT_VALIDATION.json').write_text(json.dumps({'pathways':paths,'n_pathways':85,'paired_patient_counts':dict(zip(factors,npats)),'row_order_unchanged':True,'data_changed':False,'size_semantics':'patient-equal fraction AUC>0','color_semantics':'patient-equal mean AUC','new_tests':False},indent=2))
print('COLUMN_READY')
