from geometry import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
allr=json.loads((P/'transforms/ALL_CANDIDATE_TRANSFORMS.json').read_text());fig,axes=plt.subplots(4,2,figsize=(10,18))
for row,label in enumerate(['R114-T-neg','R114-T-pos','Y12-T-neg','Y12-T-pos']):
 sample=label.rsplit('-',1)[0].replace('-','_');a=np.load(P/'inputs'/(sample+'_ST.npz'));b=np.load(P/'inputs'/(label+'_SM.npz'));img=warp(a['rgb'],a['T'],order=1);rr=sorted([x for x in allr if x['label']==label and x['model']=='similarity' and x['variant']=='HE'],key=lambda x:(-x['dice'],x['surface_mean']))
 for rank,r in enumerate(rr[:2]):
  ax=axes[row,rank];ax.imshow(img);ax.contour(a['HE'],levels=[.5],colors=['cyan'],linewidths=.9);ax.contour(warp(b['source'],np.array(r['matrix']))>.5,levels=[.5],colors=['magenta' if rank==0 else 'orange'],linewidths=1);angle=((r['angle_degrees']+180)%360)-180;ax.set_title(f"{label} | {'preferred' if rank==0 else 'runner-up'}\nDice {r['dice']:.3f}; rotation {angle:.1f} deg; mirror {r['flip']}",fontsize=11);ax.axis('off')
fig.legend(handles=[Line2D([0],[0],color='cyan',label='H&E tissue boundary'),Line2D([0],[0],color='magenta',label='Preferred MSI contour'),Line2D([0],[0],color='orange',label='Alternative MSI contour')],loc='lower center',ncol=3,frameon=False);fig.suptitle('Morphology-only registration; anatomical correspondence unverified',fontsize=13);fig.tight_layout(rect=[0,.035,1,.975]);fig.savefig(P/'diagnostics/registration_overlays.png',dpi=160);plt.close(fig)
