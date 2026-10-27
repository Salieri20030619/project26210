from pathlib import Path
import pandas as pd,numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
r=Path(__file__).resolve().parents[2];work=Path(__file__).resolve().parents[1];d=Path('C:/Users/31937/Downloads/26210/02_2.1恶性上皮亚克隆_v1');d.mkdir(parents=True,exist_ok=True)
s=r/'phase32_part1_m03_cnv_clone_trees_local_v1';src=r/'phase43_part1_client_figures_local_v1/source_remote_v1'
c=pd.read_csv(s/'input_snapshot_local_v1/formal_cell_cnv_clone_assignments_v1.tsv.gz',sep='\t');b=pd.read_csv(src/'tumor_epithelial_malignancy_calls_v1.tsv.gz',sep='\t');assert len(c)==11950 and len(b)==38670
c=c.merge(b[['exact_cell_id','infercnv_cnv_gene_burden','infercnv_high_cnv_threshold']],on='exact_cell_id',validate='one_to_one');c['relative_burden']=np.log2((c.infercnv_cnv_gene_burden+1)/(c.infercnv_high_cnv_threshold+1));b['relative_burden']=np.log2((b.infercnv_cnv_gene_burden+1)/(b.infercnv_high_cnv_threshold+1))
e=pd.read_csv(s/'results_local_v1/m03_tree_edges_v1.tsv',sep='\t');nodes=pd.read_csv(s/'results_local_v1/m03_tree_nodes_cnmf_annotation_v1.tsv',sep='\t');trans=pd.read_csv(s/'results_local_v1/m03_tree_edge_bin_transitions_v1.tsv.gz',sep='\t')
patients=['GSE207493__RCC104','GSE207493__RCC86','GSE156632__RCC4','GSE156632__RCC5'];colors=['#6BA5BA','#D97A69','#E4AD55','#8AAF88','#AF7EA4','#8B93BF','#B7B66C'];rng=np.random.default_rng(26210)
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'pdf.fonttype':42,'axes.linewidth':.7})
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(d/(name+'.'+ext),dpi=300,bbox_inches='tight',pad_inches=.12)
 plt.close(fig)
def clean(ax):ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3)
def distribution(ax,arr,col):
 for j,(v,co) in enumerate(zip(arr,col)):
  v=np.array(v);v=v[np.isfinite(v)]
  if len(v)>2 and np.ptp(v)>1e-10:
   part=ax.violinplot([v],positions=[j],widths=.78,showextrema=False)
   part['bodies'][0].set_facecolor(co);part['bodies'][0].set_edgecolor('#555555');part['bodies'][0].set_alpha(.7);part['bodies'][0].set_linewidth(.7)
  pts=v if len(v)<70 else rng.choice(v,70,replace=False);ax.scatter(j+rng.uniform(-.13,.13,len(pts)),pts,s=4,color=co,edgecolors='none',alpha=.5,rasterized=True)
  q=np.quantile(v,[.25,.5,.75]);ax.plot([j,j],[q[0],q[2]],lw=3,color='#555555');ax.plot([j-.075,j+.075],[q[1],q[1]],lw=1.5,color='white')
 clean(ax)
# Patient-level separation panels; complete cases retained, discordant/indeterminate shown separately.
fig,axs=plt.subplots(1,2,figsize=(10,4.8),gridspec_kw={'width_ratios':[1.05,1.6]})
st=['non_malignant_consensus','malignant_consensus'];arr=[b.loc[b.malignancy_status_v1==t,'relative_burden'].dropna() for t in st];distribution(axs[0],arr,['#78AE9A','#BA78A5']);axs[0].axhline(0,ls='--',c='#888888',lw=.6);axs[0].set_xticks([0,1],[f'Non-malignant\nconsensus\n{len(arr[0]):,} cells',f'Malignant\nconsensus\n{len(arr[1]):,} cells']);axs[0].set_ylabel('Relative CNV burden\nlog2((burden + 1)/(threshold + 1))');axs[0].set_title('CNV-defined epithelial groups',loc='left',fontweight='bold')
status=b.malignancy_status_v1.value_counts();names={'malignant_consensus':'Malignant consensus','malignant_infercnv_only':'Malignant inferCNV only','non_malignant_consensus':'Non-malignant consensus','method_discordant':'Method-discordant','indeterminate':'Indeterminate'}
cols={'malignant_consensus':'#BA78A5','malignant_infercnv_only':'#DBB56C','non_malignant_consensus':'#78AE9A','method_discordant':'#BCBCBC','indeterminate':'#E1E1E1'}
axs[1].barh(range(len(status)),status.values,color=[cols[t] for t in status.index]);axs[1].set_yticks(range(len(status)),[names[t] for t in status.index],fontsize=8);axs[1].invert_yaxis();axs[1].set_xlim(0,27000)
for j,v in enumerate(status.values):axs[1].text(v+300,j,f'{v:,}',va='center',fontsize=8)
axs[1].set_xlabel('Tumour-tissue epithelial cells');axs[1].set_title('All CNV classifications retained',loc='left',fontweight='bold');clean(axs[1]);fig.tight_layout(w_pad=2);save(fig,'04_CNV_groups_violin_and_classification_v1')
fig,axs=plt.subplots(4,2,figsize=(11.5,13));summary=[]
for row,pt in enumerate(patients):
 pc=c[(c.infercnv_unit==pt)&(c.clone_cells>=20)].copy();clones=sorted(pc.formal_inferred_cnv_clone_v1.unique());clcols=colors[:len(clones)]
 for col,field in enumerate(['relative_burden','TYMP_log1p_CP10K_v1']):
  ax=axs[row,col];arrays=[pc.loc[pc.formal_inferred_cnv_clone_v1==cl,field].dropna() for cl in clones];distribution(ax,arrays,clcols);ax.set_xticks(range(len(clones)),[cl.rsplit('_',1)[-1]+'\n'+str(len(v))+' cells' for cl,v in zip(clones,arrays)],fontsize=7)
  ax.set_ylabel('Relative CNV burden' if col==0 else 'TYMP log1p CP10K');ax.set_title(pt.replace('__','  ')+' | '+('CNV' if col==0 else 'TYMP'),loc='left',fontweight='bold')
  if col==0:ax.axhline(0,ls='--',lw=.5,c='#999999')
  for cl,v in zip(clones,arrays):summary.append({'patient':pt,'clone':cl,'feature':field,'cells':len(v),'median':float(np.median(v))})
fig.tight_layout(h_pad=2,w_pad=2);save(fig,'05_subclone_CNV_TYMP_distributions_v1')
# Neutral-root-anchored topology rendered as descending branches, preserving every edge.
def tree(ax,pt):
 pe=e[e.infercnv_unit==pt];pn=nodes[nodes.infercnv_unit==pt].set_index('node_id');children={}
 for row in pe.itertuples():children.setdefault(row.parent_node,[]).append(row.child_node);children.setdefault(row.child_node,[])
 for nd in children:children[nd].sort()
 leaves=[]
 def find(n):
  if not children[n]:leaves.append(n)
  else:
   for k in children[n]:find(k)
 find('NEUTRAL_ROOT');pos={};lx={n:i*2.6 for i,n in enumerate(leaves)}
 def assign(n,depth):
  for k in children[n]:assign(k,depth+1)
  xx=lx[n] if not children[n] else np.mean([pos[k][0] for k in children[n]]);pos[n]=(xx,-depth*1.6)
 assign('NEUTRAL_ROOT',0);ndorder=sorted(n for n in pos if n!='NEUTRAL_ROOT');co={n:colors[i%len(colors)] for i,n in enumerate(ndorder)}
 for edge in pe.itertuples():
  pa,ch=edge.parent_node,edge.child_node;x1,y1=pos[pa];x2,y2=pos[ch];ax.plot([x1,x2],[y1,y2],lw=8,c='#888888',solid_capstyle='round',zorder=1);ax.plot([x1,x2],[y1,y2],lw=6,c=co[ch],solid_capstyle='round',zorder=2)
 for n in ndorder:
  xx,yy=pos[n];row=pn.loc[n];ax.scatter(xx,yy,s=75,c=co[n],edgecolor='white',lw=.6,zorder=3);ax.text(xx+.16,yy+.12,n.rsplit('_',1)[-1],fontsize=8,fontweight='bold',zorder=6,bbox={'facecolor':'white','edgecolor':'none','alpha':.8,'pad':.3})
  ax.text(xx+.16,yy-.20,f'{row.dominant_cnmf_factor}; {int(row.clone_cells)} cells',fontsize=7,zorder=6,bbox={'facecolor':'white','edgecolor':'none','alpha':.8,'pad':.3})
 rootx,rooty=pos['NEUTRAL_ROOT'];ax.plot([rootx,rootx],[rooty,rooty+.7],lw=8,c='#B7B7B7',solid_capstyle='round');ax.text(rootx,rooty+.85,'Neutral reference',ha='center',fontsize=8)
 ax.set_xlim(-1,max(lx.values())+1.9);ax.set_ylim(min(v[1] for v in pos.values())-.7,1.45);ax.axis('off');ax.set_title(pt.replace('__','  '),fontweight='bold',loc='left',pad=8)
 # Event table beside the branch schematic avoids obscuring the true topology.
 changed=trans[(trans.infercnv_unit==pt)&(trans.state_delta_from_parent!=0)].copy();changed['direction']=np.where(changed.state_delta_from_parent>0,'gain','loss');event=[]
 for n in ndorder:
  ss=changed[changed.child_node==n].groupby(['chr_number','direction']).size().sort_values(ascending=False)
  if len(ss):key=ss.index[0];event.append(f"{n.rsplit('_',1)[-1]}: chr{int(key[0])} {key[1]} ({ss.iloc[0]} bins)")
 return event
fig,axs=plt.subplots(2,2,figsize=(12,12))
for ax,pt in zip(axs.ravel(),patients):tree(ax,pt)
fig.subplots_adjust(wspace=.20,hspace=.25,bottom=.07);fig.text(.1,.025,'Colours identify inferred clones; labels show dominant cNMF factor and cell count. Branch depth is not clinical time.',fontsize=9,color='#555555');save(fig,'06_four_CNV_clone_trees_v1')
for pt in patients:
 fig=plt.figure(figsize=(7.4,8.5));ax=fig.add_axes([.06,.30,.88,.64]);ev=tree(ax,pt);fig.text(.09,.24,'Largest chromosome-level shift per branch',fontweight='bold',fontsize=10);fig.text(.09,.215,'\n'.join(ev),va='top',fontsize=9,linespacing=1.5);fig.text(.09,.025,'Relative to the parent inferred state; no mutation or LOH inference.',fontsize=8,color='#555555');save(fig,'Tree_'+pt+'_v1')
(work/'qc/local_plot_summary.json').write_text(json.dumps({'cells':11950,'all_tumour_epithelial':38670,'matched':len(c),'examples':patients,'clone_distributions':summary},indent=2),encoding='utf8')
print('Local figures generated')
