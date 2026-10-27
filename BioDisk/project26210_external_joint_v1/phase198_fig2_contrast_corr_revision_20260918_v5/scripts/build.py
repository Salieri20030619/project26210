from pathlib import Path
import json,shutil
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase198_fig2_contrast_corr_revision_20260918_v5'
old=r/'phase197_top20_network_dynamics_20260917_v4'
# Reuse original Top20 drawing logic; only display normalization and height change.
s=(old/'scripts/redraw.py').read_text().split('# Membership network:')[0]
s=s.replace("d=r/'phase197_top20_network_dynamics_20260917_v4'","d=r/'phase198_fig2_contrast_corr_revision_20260918_v5'")
s=s.replace('LinearSegmentedColormap,Normalize','LinearSegmentedColormap,Normalize,PowerNorm')
s=s.replace('figsize=(8.7,8.8)','figsize=(8.7,11.0)')
s=s.replace("cmap='viridis',vmin=0,vmax=vmax","cmap='turbo',norm=PowerNorm(gamma=.45,vmin=0,vmax=vmax)")
s=s.replace("cb.ax.tick_params(labelsize=6)","cb.set_ticks([0,.01,.03,.06,.10,.20,.30,.40]);cb.ax.tick_params(labelsize=6)")
s=s.replace("'Patient-equal mean AUC'","'Patient-equal mean AUC (power scale)'")
s=s.replace("panel(fig,'top20_pathways')","fig.text(.43,.005,'Shared power scale: gamma = 0.45 | Original AUC values',fontsize=7)\npanel(fig,'top20_pathways')")
for name in ['TOP20_SELECTION.json','TOP20_RANKING_ALL85.tsv']:
 shutil.copy2(old/'tables'/name,d/'tables'/name)
(d/'scripts/draw_top20.py').write_text(s)
exec(compile(s,'draw_top20.py','exec'),{})
# Reuse exact density rendering and frozen metadata; remove only nested letters.
src=(r/'phase192_sections28_212_figures_20260917_v1/scripts/render_figures.py').read_text()
s=src.split('# 01 actual biochemical relations only')[0]+src.split('# 02\n')[1].split('# heatmap helper')[0]
s=s.replace("'A   Frozen macrophage states','B   C4 localization'","'Frozen macrophage states','C4 localization'")
s=s.replace("for ext in ['pdf','svg','png']:fig.savefig(d/'figures'/f'{name}.{ext}',dpi=600)","fig.savefig(r/'phase198_fig2_contrast_corr_revision_20260918_v5/panels/C4_no_nested_letters.pdf',dpi=400)")
s=s.replace("fig.savefig(d/'qc'/f'{name}.jpg',dpi=170)","pass")
(d/'scripts/draw_c4.py').write_text(s)
exec(compile(s,'draw_c4.py','exec'),{})
# Retain standalone biochemical schematic byte-for-byte.
shutil.copy2(r/'phase192_sections28_212_figures_20260917_v1/figures/01_MTA_metabolic_map.pdf',d/'figures/Standalone_MTA_biochemical_context.pdf')
print('PANELS_READY',flush=True)
