from pathlib import Path
import ast,json
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase195_mainfigures_compact_20260917_v2';o=r/'phase130_section2_4_BE_correction_remote_v1'
raw=(o/'scripts/redraw_E_gene_network.py').read_text();fn=next(n for n in ast.parse(raw).body if isinstance(n,ast.FunctionDef) and n.name=='network')
tfs=['ILF2','HOXD3','HNF4A'];colors={'Both methods':'#806080','GENIE3 only':'#B57756','GRNBoost2 only':'#397C9C'}
nodecolor={'TF':'#6A263E','Target':'#8898A3','Shared target':'#C5913D'};shared={'PTER','PTN','TRPV4'}
edges=pd.read_csv(o/'qc/GENE_NETWORK_DISPLAY_EDGES.tsv',sep='\t')
nodes=pd.read_csv(o/'qc/GENE_NETWORK_NODE_POSITIONS.tsv',sep='\t')
assert len(edges)==71 and len(nodes)==71
exec(compile(ast.Module(body=[fn],type_ignores=[]),'original_network_function','exec'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(7,7));fig.subplots_adjust(left=.01,right=.99,top=.99,bottom=.11)
pos=network(ax,edges,nodes.gene.tolist(),'',size=9)
assert all(np.allclose(pos[z.gene],[z.x,z.y]) for z in nodes.itertuples())
ax.set_xlim(-1.31,1.31);ax.set_ylim(-1.29,1.29)
handles=[Line2D([],[],color=c,lw=1,label=k) for k,c in colors.items()]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.044),ncol=3,frameon=False,fontsize=7)
handles=[Line2D([],[],ls='',marker='D' if k=='TF' else 'o',color=c,label=k,ms=5) for k,c in nodecolor.items()]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.008),ncol=3,frameon=False,fontsize=7)
fig.savefig(d/'panels/network_compact.pdf',bbox_inches='tight',pad_inches=.06)
(d/'tables/NETWORK_VALIDATION.json').write_text(json.dumps({'nodes':71,'edges':71,'node_positions_exactly_preserved':True,'original_function_reused':True},indent=2))
# Eliminate rounding a nonzero contribution to an apparent exact zero.
p=d/'scripts/redraw_panels.py';s=p.read_text().replace("f'{vv[i,j]*100:.0f}'","f'{vv[i,j]*100:.1f}'")
p.write_text(s)
start=s.index('# Contribution:');end=s.index('# Regenerate dynamics')
exec(compile(s[start:end],str(p),'exec'),dict(globals(),save=lambda fig,n:(fig.savefig(d/'panels'/f'{n}.pdf',bbox_inches='tight',pad_inches=.045),plt.close(fig))))
p=d/'scripts/assemble_compact.py';s=p.read_text().replace("'network':(d/'inputs/user_selected_regulon_network.pdf',[0,0,1,1])","'network':(d/'panels/network_compact.pdf',[0,0,1,1])")
p.write_text(s)
print('Network and contribution revised')
