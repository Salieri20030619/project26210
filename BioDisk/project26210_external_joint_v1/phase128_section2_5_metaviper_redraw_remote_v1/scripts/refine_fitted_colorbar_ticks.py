from pathlib import Path
import subprocess,json
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
w=Path(__file__).resolve().parents[1];root=w.parent;q=w/'qc';out=w/'figures'
genes=['CASP9','CHMP1A','PRSS8','PXN'];methods=['Monocle3','Slingshot','scTour']
a=pd.read_csv(w/'data/SELECTED_METAVIPER_ACTIVITY.tsv',sep='\t').set_index('exact_cell_id');m=a
z=np.load(root/'phase123_section2_2_reference_redraw_remote_v1/data/plot_values_server_only.npz')
xyall=pd.DataFrame(z['xy'],index=z['cells']);xy=xyall.loc[m.index].values
fd=pd.read_csv(w/'data/FITTED_CELL_PROJECTIONS.tsv.gz',sep='\t',index_col=0).loc[m.index]
fits={tuple(c.split('__',1)):fd[c].values for c in fd}
s=(w/'scripts/draw_section25.py').read_text()
exec(s[s.index('plt.rcParams.update'):s.index('# 1:')])
code=s[s.index('def umaps'):s.index('# 5:')].rsplit('\numaps(',1)[0]
code=code.replace('cmap=div,vmin=lo,vmax=hi','cmap=div,norm=TwoSlopeNorm(vmin=lo,vcenter=0,vmax=hi)')
code=code.replace('  cb.set_label','  if limits is not None:\n   ticks=[lo,lo/2,0,hi/2,hi];cb.set_ticks(ticks);cb.set_ticklabels([format(t,".2f") for t in ticks]);cb.ax.minorticks_off()\n  cb.set_label')
exec(code)
fitlimits=json.loads((q/'INDEPENDENT_VALIDATION.json').read_text())['fitted_maps_scales']
for i,method in enumerate(methods,6):
 stem=f'0{i}_F_{method}_fitted_activity_UMAP_v1'
 umaps({g:fits[(method,g)] for g in genes},method+': patient-local fitted activity projected on the frozen UMAP',stem,fitlimits)
 subprocess.run(['pdftoppm','-scale-to','1700','-singlefile','-png',str(out/(stem+'.pdf')),str(q/(stem+'_pdf'))],capture_output=True,check=True)
print('FIT_MAP_TICKS_REFINED')
