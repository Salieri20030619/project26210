from pathlib import Path
import json,subprocess
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MPath
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
w=Path(__file__).resolve().parents[1];q=w/'qc';out=w/'figures'
methods=['GENIE3','GRNBoost2'];tfs=['ILF2','HOXD3','HNF4A'];mc={'GENIE3':'#AB7764','GRNBoost2':'#5E91AB'}
source=Path('/mnt/BioDisk/wangmingyuan/project26210_part1_cnmf_pyscenic_phase26_20260814_v1/analysis_v3')
gd=pd.read_csv(source/'patient_cnmf_group_regulon_effects_v1.tsv.gz',sep='\t');group_summary=pd.read_csv(source/'cnmf_group_regulon_patient_summary_v1.tsv',sep='\t')
targets=pd.read_csv(w/'data/CANDIDATE_REGULON_TARGETS.tsv',sep='\t')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':11,'axes.linewidth':.65,'pdf.fonttype':42,'savefig.facecolor':'white'})
def clean(ax):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3)
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=360,bbox_inches='tight',pad_inches=.15)
 plt.close(fig)
s=(w/'scripts/draw_section24.py').read_text()
part=s[s.index('# D:'):s.index('# F:')]
part=part.replace("top=.80,bottom=.23","top=.80,bottom=.28")
part=part.replace("xlabel='Patient median AUCell: F1-dominant minus other cells'","xlabel='Patient median AUCell difference\\n(F1-dominant minus other cells)'")
part=part.replace("ax.ticklabel_format(axis='x',style='sci',scilimits=(-2,2))","ax.xaxis.set_major_locator(MaxNLocator(5));ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))")
part=part.replace("figsize=(14.8,6.1)","figsize=(15.6,7.4)").replace("top=.82,bottom=.17","top=.86,bottom=.14")
part=part.replace("ylim=(-1.38,1.5)","ylim=(-1.87,1.5)")
part=part.replace("ax.text(.5,-.035,'All non-self targets:","ax.text(0,-1.77,'All non-self targets:")
part=part.replace(",transform=ax.transAxes,ha='center',fontsize=8)",",ha='center',fontsize=8)")
part=part.replace("bbox_to_anchor=(.5,.075)","bbox_to_anchor=(.5,.065)")
exec(part)
for stem in ['05_D_candidate_F1_group_exploratory_v1','06_E_candidate_circular_regulon_networks_v1']:
 subprocess.run(['pdftoppm','-scale-to','1800','-singlefile','-png',str(out/(stem+'.pdf')),str(q/(stem+'_pdf'))],check=True,capture_output=True)
print('TWO_LAYOUT_REFINEMENTS_COMPLETE')
