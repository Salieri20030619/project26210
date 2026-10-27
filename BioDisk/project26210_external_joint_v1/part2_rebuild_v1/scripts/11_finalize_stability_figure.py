from pathlib import Path
import json,shutil
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1')
C=R/'03_macrophage_clustering/run_v4_exploratory'
F=R/'10_figures/run_v3';F.mkdir(exist_ok=False)
for p in (R/'10_figures/run_v2').glob('*'):
 if p.suffix in ['.png','.svg']: shutil.copy2(p,F/p.name)
m=pd.read_csv(C/'RESOLUTION_STABILITY.tsv',sep='\t');s=json.loads((C/'STATUS.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
ax.plot(m.resolution,m.minimum_between_seed_ARI,'o-',label='Minimum seed-pair ARI',color='#355f8d')
ax.scatter([s['selected_resolution']],[s['subsample_regraph_ARI']],marker='s',s=55,color='#3a9d8f',label='80% regraph ARI (selected resolution)',zorder=5)
ax.axhline(.85,color='#355f8d',ls=':',lw=1,label='Seed ARI gate = 0.85')
ax.axhline(.8,color='#3a9d8f',ls=':',lw=1,label='Regraph ARI gate = 0.80')
ax.axvline(s['selected_resolution'],color='#a6514c',ls='--',label='Exploratory resolution = 0.8')
ax.set(xlabel='Leiden resolution',ylabel='Adjusted Rand index',ylim=(0,1.02),title='Purified candidate grid: global stability gates unmet')
ax.legend(frameon=False,fontsize=8,loc='lower center')
for ext in ['png','svg']: fig.savefig(F/('05_clustering_stability.'+ext),dpi=240,bbox_inches='tight')
(F/'FIGURE_FILES.json').write_text(json.dumps({'files':[p.name for p in sorted(F.iterdir())],'correction':'05 now displays measured selected-resolution subsample ARI 0.6328556; other six figures unchanged from run_v2'},indent=2))
print('FIGURES_RUN_V3_READY')
