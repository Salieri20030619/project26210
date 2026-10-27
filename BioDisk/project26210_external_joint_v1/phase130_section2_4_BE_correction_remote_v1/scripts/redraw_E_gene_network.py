from pathlib import Path
import json,hashlib,subprocess,sys
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
w=Path(__file__).resolve().parents[1];q=w/'qc';old=w.parent/'phase127_section2_4_pyscenic_redraw_remote_v1'
source=old/'data/CANDIDATE_REGULON_TARGETS.tsv';previous=old/'qc/CIRCULAR_DISPLAY_EDGES.tsv'
raw=pd.read_csv(source,sep='\t');before=pd.read_csv(previous,sep='\t');tfs=['ILF2','HOXD3','HNF4A'];methods=['GENIE3','GRNBoost2']
assert not raw.duplicated(['tf','target','method']).any()
targets=list(dict.fromkeys(before.target));nodes=set(targets)|set(tfs)
# Gene-only induced subgraph on the prior displayed genes: include every real edge among these nodes.
induced=raw[raw.tf.isin(nodes)&raw.target.isin(nodes)&(raw.tf!=raw.target)].copy()
pivot=induced.pivot(index=['tf','target'],columns='method',values='weight').reset_index()
pivot['support']=np.where(pivot[methods].notna().all(axis=1),'Both methods',np.where(pivot.GENIE3.notna(),'GENIE3 only','GRNBoost2 only'))
assert len(nodes)==71 and len(targets)==68 and len(pivot)==71 and len(induced)==89
previous_edges=set(map(tuple,before[['tf','target']].drop_duplicates().values));current_edges=set(map(tuple,pivot[['tf','target']].values))
assert previous_edges<=current_edges and len(current_edges-previous_edges)==3
assert not ({'GENIE3','GRNBoost2'}&nodes)
shared=set(pivot.groupby('target').tf.nunique().loc[lambda s:s>1].index)
assert shared=={'PTER','PTN','TRPV4'}
for r in before.itertuples():
 observed=raw[(raw.tf==r.tf)&(raw.target==r.target)&(raw.method==r.method)].weight.iloc[0]
 assert np.isclose(observed,r.weight,atol=1e-12)
induced.to_csv(q/'GENE_NETWORK_SOURCE_METHOD_EDGES.tsv',sep='\t',index=False)
pivot.to_csv(q/'GENE_NETWORK_DISPLAY_EDGES.tsv',sep='\t',index=False)
pd.DataFrame(sorted(current_edges-previous_edges),columns=['tf','target']).to_csv(q/'ADDED_TRUE_CROSS_REGULON_EDGES.tsv',sep='\t',index=False)
colors={'Both methods':'#806080','GENIE3 only':'#B57756','GRNBoost2 only':'#397C9C'}
nodecolor={'TF':'#6A263E','Target':'#8898A3','Shared target':'#C5913D'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'savefig.facecolor':'white'})
registry=[]
def network(ax,sub,order,title,size=8.8):
 assert set(order)==set(sub.tf)|set(sub.target) and len(order)==len(set(order))
 angles=np.pi/2-2*np.pi*np.arange(len(order))/len(order)
 pos={g:np.array([np.cos(t),np.sin(t)]) for g,t in zip(order,angles)}
 # Actual directed gene-to-gene connections; methods only determine edge appearance.
 for support in ['GRNBoost2 only','GENIE3 only','Both methods']:
  for r in sub[sub.support==support].itertuples():
   a=pos[r.tf];b=pos[r.target];path=MPath([a*.975,a*.13,b*.13,b*.975],[MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4])
   arrow=FancyArrowPatch(path=path,arrowstyle='-|>',mutation_scale=6.0 if len(order)>40 else 7.5,lw=1.15 if support=='Both methods' else .8,color=colors[support],alpha=.83,zorder=1)
   ax.add_patch(arrow)
 for gene,t in zip(order,angles):
  kind='TF' if gene in tfs else 'Shared target' if gene in shared else 'Target'
  xy=pos[gene];is_tf=kind=='TF';ax.scatter(*xy,s=115 if is_tf else 32,marker='D' if is_tf else 'o',color=nodecolor[kind],edgecolors='white',lw=.7,zorder=4)
  deg=np.rad2deg(t)%360;left=90<deg<270;rotation=deg+180 if left else deg
  ax.text(xy[0]*1.055,xy[1]*1.055,gene,rotation=rotation,rotation_mode='anchor',ha='right' if left else 'left',va='center',fontsize=size+2 if is_tf else size,fontweight='bold' if is_tf else 'normal',color=nodecolor[kind],zorder=5)
 ax.set_aspect('equal');ax.set_xlim(-1.54,1.54);ax.set_ylim(-1.45,1.45);ax.axis('off');ax.set_title(title,loc='left',fontsize=12,fontweight='bold',pad=5)
 return pos
edge_handles=[Line2D([],[],color=c,lw=1.8,label=k) for k,c in colors.items()]
node_handles=[Line2D([],[],ls='',marker='D' if k=='TF' else 'o',color=c,label={'TF':'TF gene','Target':'Target gene','Shared target':'Target of multiple displayed TFs'}[k],ms=7) for k,c in nodecolor.items()]
def save(fig,name,title,subtitle):
 fig.text(.045,.969,title,fontsize=17,fontweight='bold',color='#273844',ha='left')
 fig.text(.045,.940,subtitle,fontsize=10,color='#56636D',ha='left')
 l1=fig.legend(handles=edge_handles,loc='lower center',bbox_to_anchor=(.5,.078),ncol=3,frameon=False,title='Edge support',fontsize=9,title_fontsize=10)
 fig.add_artist(l1)
 fig.legend(handles=node_handles,loc='lower center',bbox_to_anchor=(.5,.029),ncol=3,frameon=False,fontsize=9)
 fig.text(.5,.013,'Arrows: TF gene -> target gene | predicted regulon edges; methods are edge evidence, not network nodes.',ha='center',fontsize=8)
 for ext in ['png','pdf']:fig.savefig(w/'figures'/(name+'.'+ext),dpi=360)
 plt.close(fig)
 pdf=w/'figures'/(name+'.pdf');subprocess.run(['pdftoppm','-r','125','-singlefile','-png',str(pdf),str(q/(name+'_pdf'))],check=True)
 registry.append({'name':name,'title':title,'subtitle':subtitle});print('SAVED',name,flush=True)
# Preserve the original target membership/order to avoid selecting a new gene story.
order=[]
for tf in tfs:order += [tf]+list(dict.fromkeys(before.loc[before.tf==tf,'target']))
fig,ax=plt.subplots(figsize=(12,12));fig.subplots_adjust(left=.03,right=.97,bottom=.16,top=.9)
pos=network(ax,pivot,order,'',size=8.5)
pd.DataFrame([{'gene':g,'x':float(v[0]),'y':float(v[1]),'kind':'TF' if g in tfs else 'Shared target' if g in shared else 'Target'} for g,v in pos.items()]).to_csv(q/'GENE_NETWORK_NODE_POSITIONS.tsv',sep='\t',index=False)
save(fig,'06_E_gene_to_gene_regulon_network_v2','E  Gene-to-gene connections among candidate regulons','3 TF genes + 68 target genes | 71 unique directed edges | all supported edges among the previously displayed genes')
# Individual networks expose each source gene and retain the same induced-network edges.
for i,tf in enumerate(tfs,1):
 sub=pivot[pivot.tf==tf].copy();orig=list(dict.fromkeys(before.loc[before.tf==tf,'target']));extra=sorted(set(sub.target)-set(orig));node_order=[tf]+orig+extra
 fig,ax=plt.subplots(figsize=(10,10));fig.subplots_adjust(left=.04,right=.96,bottom=.17,top=.88)
 network(ax,sub,node_order,'',size=10)
 counts=sub.support.value_counts().to_dict()
 save(fig,f'06_E{i}_{tf}_gene_network_v2',f'E{i}  {tf} gene-to-gene regulon network',f'{len(sub)} target genes | '+', '.join(f'{k}: {counts.get(k,0)}' for k in colors))
validation={'status':'PASS','gene_nodes':71,'TF_nodes':tfs,'target_nodes':68,'unique_directed_gene_edges':71,'original_method_edge_records':89,'previous_gene_edges_preserved':68,'additional_source_verified_edges':[{'tf':tf,'target':target} for tf,target in sorted(current_edges-previous_edges)],'shared_target_genes':sorted(shared),'algorithm_nodes':0,'edge_support_counts':pivot.support.value_counts().to_dict(),'individual_network_targets':pivot.groupby('tf').target.nunique().to_dict(),'edge_weights_changed':False,'models_rerun':False,'new_statistical_tests':False}
(q/'E_NETWORK_VALIDATION.json').write_text(json.dumps(validation,indent=2));(q/'E_FIGURE_REGISTRY.json').write_text(json.dumps(registry,indent=2))
manifest=[{'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [source,previous,Path(__file__)]]
(q/'E_SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(validation),flush=True)
