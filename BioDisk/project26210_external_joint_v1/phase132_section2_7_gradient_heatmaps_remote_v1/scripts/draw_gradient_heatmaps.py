from pathlib import Path
import json,subprocess,sys
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,ListedColormap,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram
w=Path(__file__).resolve().parents[1];q=w/'qc';out=w/'figures'
v=json.loads((q/'VALIDATION.json').read_text());assert v['status']=='PASS'
z=np.load(w/'data/GRADIENT_HEATMAP_VALUES.npz');genes=z['genes'];grid=z['grid'];leaves=z['leaves'];cluster=z['clusters'];tree=z['tree'];n=len(genes)
states=pd.read_csv(w/'data/PATIENT_STATE_COMPOSITION.tsv.gz',sep='\t')
state_colors={'Base':'#8CA8B6','Hypoxia':'#CB8B4B','Cycle':'#8766AB','TYMP_state':'#C94C64'}
clcolors={1:'#9B7890',2:'#47879A',3:'#87A86F',4:'#C3A987'}
cmap=LinearSegmentedColormap.from_list('reference_gradient',[(0,'#130078'),(.23,'#5400A8'),(.48,'#8900AC'),(.70,'#EF7949'),(1,'#F5F520')])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.65,'pdf.fonttype':42,'savefig.facecolor':'white'})
selected=['TYMP','TK1','SLC29A1','SLC29A2','ODC1','SRM','SMS','MTAP','CHKA','PCYT1A','PCYT2','CEPT1']
registry=[]
for num,method in enumerate(['Monocle3','Slingshot','scTour'],11):
 fig=plt.figure(figsize=(14,11));left,bottom,width,height=.19,.135,.47,.66;top=bottom+height
 ax=fig.add_axes([left,bottom,width,height]);ordered=z[method][leaves]
 # Horizontal interpolation only; no smoothing or interpolation between gene rows.
 xfine=np.linspace(grid[0],grid[-1],1400);smooth=np.stack([np.interp(xfine,grid,row) for row in ordered])
 im=ax.imshow(smooth,cmap=cmap,vmin=-2.5,vmax=2.5,aspect='auto',interpolation='nearest',extent=[grid[0],grid[-1],n-.5,-.5],rasterized=True)
 ax.set_xticks([.025,.25,.5,.75,.975],['0.025','0.25','0.50','0.75','0.975']);ax.set_xlabel('Within-patient pseudotime rank (early to late)',labelpad=8);ax.set_yticks([])
 for s in ax.spines.values():s.set_visible(False)
 boundaries=np.flatnonzero(np.diff(cluster)!=0)+.5
 for val in boundaries:ax.axhline(val,color='white',lw=1.2)
 ta=fig.add_axes([.081,bottom,.083,height]);dendrogram(tree,orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#727A80',ax=ta)
 for coll in ta.collections:coll.set_linewidth(.25)
 ta.invert_yaxis();ta.axis('off')
 strip=fig.add_axes([.171,bottom,.012,height]);strip.imshow((cluster-1)[:,None],cmap=ListedColormap(list(clcolors.values())),vmin=0,vmax=3,aspect='auto',interpolation='nearest');strip.axis('off')
 for c in clcolors:
  ix=np.flatnonzero(cluster==c);y=bottom+height*(1-(ix.mean()+.5)/n)
  fig.text(.018,y,f'C{c}\n({len(ix)})',ha='left',va='center',fontweight='bold',fontsize=9,color=clcolors[c],linespacing=1.4)
 fig.text(.018,top+.012,'Gene\ncluster',va='bottom',fontsize=9,fontweight='bold')
 # Top state strip represents actual patient-equal composition, without inventing one state per bin.
 composition=states[states.method==method].groupby(['rank','state']).fraction.mean().unstack().reindex(columns=list(state_colors)).sort_index()
 assert composition.shape==(len(grid),4) and np.allclose(composition.index,grid,atol=1e-12,rtol=0)
 assert np.allclose(composition.sum(axis=1),1)
 stax=fig.add_axes([left,.856,width,.052]);stax.stackplot(grid,composition.values.T,colors=list(state_colors.values()),linewidth=0)
 stax.set_xlim(grid[0],grid[-1]);stax.set_ylim(0,1);stax.set_xticks([]);stax.set_yticks([0,1],['0','1'],fontsize=7);stax.tick_params(length=2)
 for s in stax.spines.values():s.set_visible(False)
 fig.text(left-.027,.882,'Cell-state\ncomposition',ha='right',va='center',fontsize=9)
 pax=fig.add_axes([left,.825,width,.018]);pax.imshow(grid[None,:],cmap='RdBu_r',vmin=0,vmax=1,aspect='auto',extent=[grid[0],grid[-1],0,1]);pax.axis('off')
 fig.text(left-.008,.834,'Pseudotime',ha='right',va='center',fontsize=9)
 # Locate each named gene in the real row ordering; repel text labels only.
 ordergenes=genes[leaves].tolist();labels=[]
 for gene in selected:
  row=ordergenes.index(gene);y=bottom+height*(1-(row+.5)/n);labels.append([gene,row,y,y])
 labels.sort(key=lambda r:r[2]);gap=.023
 for i in range(1,len(labels)):labels[i][3]=max(labels[i][3],labels[i-1][3]+gap)
 labels[-1][3]=min(labels[-1][3],top-.008)
 for i in range(len(labels)-2,-1,-1):labels[i][3]=min(labels[i][3],labels[i+1][3]-gap)
 assert labels[0][3]>=bottom
 for gene,row,y,yp in labels:
  fig.lines.append(plt.Line2D([left+width+.002,.692,.724],[y,y,yp],transform=fig.transFigure,color='#66727D',lw=.5))
  fig.text(.730,yp,gene,fontsize=8.5,va='center',fontweight='bold' if gene=='TYMP' else 'normal')
 fig.text(.692,top+.018,'Selected enzymes',fontsize=9,fontweight='bold')
 lx=.84
 ca=fig.add_axes([lx,.751,.022,.145]);cb=fig.colorbar(im,cax=ca,extend='both',ticks=[-2,-1,0,1,2]);cb.ax.tick_params(labelsize=8,length=3);cb.set_label('Row z-score',fontsize=10)
 fig.text(lx,.704,'Cell state',fontsize=10,fontweight='bold')
 for i,(state,col) in enumerate(state_colors.items()):
  yy=.676-i*.027;fig.add_artist(Rectangle((lx,yy-.006),.013,.015,facecolor=col,transform=fig.transFigure,edgecolor='none'))
  fig.text(lx+.02,yy+.001,state.replace('_',' '),fontsize=9,va='center')
 fig.text(lx,.548,'Pseudotime rank',fontsize=10,fontweight='bold')
 pca=fig.add_axes([lx,.505,.12,.016]);pcb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap='RdBu_r'),cax=pca,orientation='horizontal',ticks=[0,1]);pcb.ax.set_xticklabels(['Early','Late'],fontsize=8);pcb.ax.tick_params(length=2)
 fig.text(lx,.433,'Gene cluster',fontsize=10,fontweight='bold')
 for i,(c,col) in enumerate(clcolors.items()):
  yy=.405-i*.028;fig.add_artist(Rectangle((lx,yy-.006),.013,.015,facecolor=col,transform=fig.transFigure,edgecolor='none'))
  fig.text(lx+.02,yy+.001,f'C{c} ({int(np.sum(cluster==c))})',fontsize=9,va='center')
 fig.text(lx,.264,'Rows and clusters\nfixed from Monocle3\nfor all three methods.',fontsize=9,va='top',linespacing=1.6,color='#53616D')
 fig.text(.025,.966,f'{method} | smooth gene-expression dynamics',fontsize=17,fontweight='bold')
 fig.text(.025,.935,f'{n:,} genes | 20 independent patients | 11,921 malignant cells | fixed gene set across methods',fontsize=10,color='#53616D')
 fig.text(.19,.048,'Each row: log1p RNA trend, standardized across rank within the method. Equal weight per patient; Gaussian bandwidth = 0.075.',fontsize=8)
 fig.text(.19,.029,'Columns: 200 fitted rank positions, interpolated horizontally for display. Colors saturate at +/-2.5; no new statistical tests.',fontsize=8)
 name=f'{num:02d}_{method}_gradient_heatmap_v1'
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=360)
 plt.close(fig)
 subprocess.run(['pdftoppm','-r','125','-singlefile','-png',str(out/(name+'.pdf')),str(q/(name+'_pdf'))],check=True)
 registry.append({'method':method,'name':name,'genes':n,'selected_gene_label_rows':{gene:int(row) for gene,row,_,_ in labels}})
 print('SAVED',name,flush=True)
(q/'FIGURE_REGISTRY.json').write_text(json.dumps(registry,indent=2))
(q/'SOFTWARE_VERSIONS.json').write_text(json.dumps({'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},indent=2))
