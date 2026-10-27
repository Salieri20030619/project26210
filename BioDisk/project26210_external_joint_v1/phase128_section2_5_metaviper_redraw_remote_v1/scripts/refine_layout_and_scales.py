from pathlib import Path
import json,subprocess,shutil
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
w=Path(__file__).resolve().parents[1];q=w/'qc';out=w/'figures';root=w.parent
base=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1')
genes=['CASP9','CHMP1A','PRSS8','PXN'];methods=['Monocle3','Slingshot','scTour']
a=pd.read_csv(w/'data/SELECTED_METAVIPER_ACTIVITY.tsv',sep='\t').set_index('exact_cell_id')
full=pd.read_csv(base/'07_pyscenic_v1/input/cell_metadata_v1.tsv.gz',sep='\t').set_index('exact_cell_id');m=full.loc[a.index]
points=pd.read_csv(q/'m11_mr_effect_activity_points_v1.tsv',sep='\t').set_index('protein').loc[genes]
mrpc=pd.read_csv(q/'m11_patient_cross_protein_correlations_v1.tsv',sep='\t');mrsummary=pd.read_csv(q/'m11_patient_correlation_summary_v1.tsv',sep='\t')
rglobal=json.loads((q/'INDEPENDENT_VALIDATION.json').read_text())['global_four_protein_r']
summ=pd.read_csv(q/'m12_f1_activity_patient_summary_v1.tsv',sep='\t')
groupmed=pd.read_csv(q/'PAIRED_PATIENT_MEDIANS.tsv',sep='\t');group_ids=groupmed.patient.unique()
endpoints=['continuous_F1_spearman','F1_dominant_median_difference','F1_Q4_minus_Q1_median_difference']
fit_table=pd.read_csv(w/'data/FITTED_PATIENT_BALANCED_CURVES.tsv.gz',sep='\t')
fd=pd.read_csv(w/'data/FITTED_CELL_PROJECTIONS.tsv.gz',sep='\t',index_col=0).loc[m.index]
fits={tuple(c.split('__',1)):fd[c].values for c in fd}
bs=pd.read_csv(q/'m13_patient_balanced_curve_summary_v1.tsv',sep='\t')
z=np.load(root/'phase123_section2_2_reference_redraw_remote_v1/data/plot_values_server_only.npz')
xyall=pd.DataFrame(z['xy'],index=z['cells']);xy=xyall.loc[m.index].values
s=(w/'scripts/draw_section25.py').read_text()
code=s[s.index('plt.rcParams.update'):s.index('# 1:')]+s[s.index('# 1:'):s.index('# 3:')]+s[s.index('# 4:'):s.index('# 9:')]+s[s.index('# 10:'):s.index("v={'status'")]
code=code.replace('The measured MR-effect / inferred-activity relationship is not positive','MR effect and inferred activity show a negative descriptive association')
code=code.replace("ax.set_title(g,loc='left',color=gc[g]);clean(ax)","ax.set_title(g,loc='left',color=gc[g],pad=24);clean(ax)")
code=code.replace('cmap=div,vmin=lo,vmax=hi','cmap=div,norm=TwoSlopeNorm(vmin=lo,vcenter=0,vmax=hi)')
code=code.replace('squeeze=False','squeeze=False,sharey="row"')
code=code.replace('top=.90,bottom=.11','top=.88 if rna_mode else .90,bottom=.20 if rna_mode else .11')
draft=w/'drafts_before_final_layout';draft.mkdir(exist_ok=True)
for f in out.iterdir():
 if f.is_file() and not (draft/f.name).exists():shutil.copy2(f,draft/f.name)
exec(code)
v=json.loads((q/'INDEPENDENT_VALIDATION.json').read_text());v.update(UMAP_color_center_NES=0,curve_y_axes_shared_within_protein=True)
(q/'INDEPENDENT_VALIDATION.json').write_text(json.dumps(v,indent=2))
for f in sorted(out.glob('*.pdf')):
 if f.name.startswith(('03_','09_')):continue
 subprocess.run(['pdftoppm','-scale-to','1700','-singlefile','-png',str(f),str(q/(f.stem+'_pdf'))],capture_output=True,check=True)
print('FINAL_ZERO_CENTERED_COLORS_AND_LAYOUTS_COMPLETE')
