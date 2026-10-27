from pathlib import Path
import json,datetime,hashlib,subprocess,heapq,zipfile
import numpy as np,pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import dijkstra,connected_components
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parents[1];R=P.parent;A=R/'phase234_part4_replan_execution_20260921_v1'
TARGETS=[('SPP1_state','V1_SPP1_q05'),('Endothelium','endothelial_q05'),('T_cell','Tcell_q05')]
COLORS={'GLMF1':'#4B94B5','GLMF2':'#D47C6C'}
def gate(stage):
 x=subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
 with (P/'logs/memory_gate.log').open('a') as f:f.write(stage+'\t'+x.stdout+x.stderr)
 if x.returncode:raise RuntimeError(x.stdout+x.stderr)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def meta(f):
 f=Path(f);s=f.stat();return dict(path=str(f),bytes=s.st_size,mtime_ns=s.st_mtime_ns,sha256=sha(f))
def dump(f,o):(P/f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)),encoding='utf-8')
def manual_distance(G,targets):
 n=G.shape[0];dist=np.full(n,np.inf);dist[targets]=0
 heap=[(0.,int(j)) for j in targets];heapq.heapify(heap)
 while heap:
  v,u=heapq.heappop(heap)
  if v>dist[u]:continue
  for at in range(G.indptr[u],G.indptr[u+1]):
   j=G.indices[at];cand=v+G.data[at]
   if cand<dist[j]:
    dist[j]=cand;heapq.heappush(heap,(cand,int(j)))
 return dist
def main():
 gate('freeze_and_distance')
 reg=pd.read_csv(A/'tables/INPUT_SECTION_REGISTER.tsv',sep='\t')
 arrays={};summ=[];manifest=[];target_rows=[];qc=[]
 for sid in ['GLMF1','GLMF2']:
  gate('distance_'+sid)
  fp=A/'inputs'/f'{sid}_INPUTS.tsv.gz';gp=A/'graphs'/f'{sid}_ring.npz'
  before=[meta(fp),meta(gp)];manifest+=before
  d=pd.read_csv(fp,sep='\t');G=sparse.load_npz(gp).tocsr()
  assert G.shape==(len(d),len(d)) and d.spot_id.is_unique and (G-G.T).nnz==0
  xy=d[['physical_x_um','physical_y_um']].to_numpy()
  rr,cc=G.nonzero();w=np.linalg.norm(xy[rr]-xy[cc],axis=1)
  assert np.isfinite(w).all() and (w>0).all() and 90<=np.median(w)<=110
  WG=sparse.csr_matrix((w,(rr,cc)),shape=G.shape)
  assert np.max(np.abs((WG-WG.T).data),initial=0)<1e-9
  nc,component=connected_components(WG,directed=False)
  source=d.malignant_fraction.to_numpy()>=.2
  assert source.sum()>=50
  distances=[];target_masks={}
  for label,col in TARGETS:
   v=d[col].to_numpy(float);assert np.isfinite(v).all() and len(np.unique(v))>=3
   threshold=float(np.quantile(v,.8));mask=v>=threshold
   assert threshold>0 and mask.sum()>=20
   targets=np.flatnonzero(mask)
   dist=dijkstra(WG,directed=False,indices=targets,min_only=True)
   reference=manual_distance(WG,targets)
   assert np.array_equal(np.isfinite(dist),np.isfinite(reference))
   err=float(np.max(np.abs(dist[np.isfinite(dist)]-reference[np.isfinite(reference)])))
   assert err<1e-8
   distances.append(dist);target_masks[label]=mask
   target_rows.append(dict(sample=sid,target=label,input_column=col,quantile=.8,threshold=threshold,target_spots=int(mask.sum())))
   qc.append(dict(sample=sid,target=label,n_spots=len(d),independent_shortest_path_max_error=err))
  D=np.column_stack(distances)
  reachable=np.isfinite(D).all(1);total=D.sum(1)
  zero=reachable&(total<=1e-9)
  valid=source&reachable&~zero
  assert valid.sum()>=50
  W=D[valid]/total[valid,None]
  assert (W>=0).all() and (W<=1).all();np.testing.assert_allclose(W.sum(1),1,atol=1e-12)
  pos=np.column_stack([.5*W[:,0]+W[:,2],np.sqrt(3)/2*W[:,0]])
  # Recover all three barycentric coordinates from plotted XY independently.
  a=pos[:,1]/(np.sqrt(3)/2);c=pos[:,0]-.5*a;b=1-a-c
  np.testing.assert_allclose(np.column_stack([a,b,c]),W,atol=1e-12)
  status=np.full(len(d),'NOT_SOURCE',dtype=object);status[source&~reachable]='UNREACHABLE_TARGET'
  status[source&zero]='ALL_THREE_DISTANCES_ZERO';status[valid]='PLOTTED'
  obs=pd.DataFrame({'sample':sid,'spot_id':d.spot_id,'component':component,'malignant_fraction':d.malignant_fraction,'source_region':source,'status':status})
  for i,(label,_) in enumerate(TARGETS):
   obs['distance_'+label+'_um']=D[:,i]
   obs['target_'+label]=target_masks[label]
  obs.to_csv(P/'plot_data'/f'{sid}_ALL_SPOT_DISTANCE_AUDIT.tsv.gz',sep='\t',index=False)
  plot=pd.DataFrame({'sample':sid,'spot_id':d.loc[valid,'spot_id'].to_numpy(),'component':component[valid],'ternary_x':pos[:,0],'ternary_y':pos[:,1]})
  for i,(label,_) in enumerate(TARGETS):
   plot['distance_'+label+'_um']=D[valid,i];plot['share_'+label]=W[:,i]
  plot.to_csv(P/'plot_data'/f'{sid}_TERNARY_POINTS.tsv',sep='\t',index=False)
  arrays[sid]=plot
  record=dict(sample=sid,patient=reg.loc[reg['sample']==sid,'patient'].iloc[0],QC_spots=len(d),source_spots=int(source.sum()),plotted_spots=int(valid.sum()),unreachable_sources=int((source&~reachable).sum()),allzero_sources=int((source&zero).sum()),connected_components=int(nc),median_edge_um=float(np.median(w)))
  for i,(label,_) in enumerate(TARGETS):
   record['median_distance_'+label+'_um']=float(np.median(D[valid,i]))
   record['median_share_'+label]=float(np.median(W[:,i]))
   record['zero_distance_'+label]=int((D[valid,i]<=1e-9).sum())
  summ.append(record)
  assert before==[meta(fp),meta(gp)]
  print('DISTANCE_DONE',json.dumps(record),flush=True)
 pd.DataFrame(summ).to_csv(P/'tables/CASE_DISTANCE_SUMMARY.tsv',sep='\t',index=False)
 pd.DataFrame(target_rows).to_csv(P/'tables/TARGET_DEFINITIONS.tsv',sep='\t',index=False)
 pd.DataFrame(qc).to_csv(P/'qc/INDEPENDENT_DISTANCE_CHECKS.tsv',sep='\t',index=False)
 dump('qc/INPUT_MANIFEST.json',manifest)
 gate('render_ternary')
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.6,'savefig.facecolor':'white'})
 width,height=160,170
 fig=plt.figure(figsize=(width/25.4,height/25.4))
 def tx(x,y,s,**kw):return fig.text(x/width,y/height,s,**kw)
 tx(80,164,'Spatial distance composition',fontsize=11,fontweight='bold',ha='center')
 tx(80,158,'Tumor-associated regions | GLMF1 and GLMF2',fontsize=8,ha='center',color='#4A4A4A')
 def triangle(rect):
  ax=fig.add_axes([rect[0]/width,rect[1]/height,rect[2]/width,rect[3]/height])
  h=np.sqrt(3)/2
  vertices=np.array([[.5,h],[0,0],[1,0]])
  def xy(a,b,c):return np.array([.5*a+c,h*a])
  for t in [.2,.4,.6,.8]:
   for p1,p2 in [(xy(t,1-t,0),xy(t,0,1-t)),(xy(1-t,t,0),xy(0,t,1-t)),(xy(1-t,0,t),xy(0,1-t,t))]:
    ax.plot([p1[0],p2[0]],[p1[1],p2[1]],color='#DBDEE2',lw=.45,zorder=0)
  boundary=np.vstack([vertices,vertices[0]]);ax.plot(boundary[:,0],boundary[:,1],color='#44484E',lw=.85,zorder=2)
  ax.text(.91,.48,'SPP1-state',ha='center',va='center',rotation=-60,fontsize=8.2,fontweight='bold')
  ax.text(.09,.48,'Endoth.',ha='center',va='center',rotation=60,fontsize=8.2,fontweight='bold')
  ax.text(.5,-.10,'T cells',ha='center',va='top',fontsize=8.2,fontweight='bold')
  # Each edge shows the percentage of its adjacent labeled coordinate.
  for i in [0,20,40,60,80,100]:
   t=i/100
   ax.text(t,-.017,str(i),ha='center',va='top',fontsize=6.4,color='#50555A')
   ax.text(.5*t-.025,h*t+.008,str(100-i),ha='right',va='center',fontsize=6.4,color='#50555A')
   ax.text(.5+.5*t+.025,h*(1-t)+.008,str(100-i),ha='left',va='center',fontsize=6.4,color='#50555A')
  ax.set_xlim(-.16,1.16);ax.set_ylim(-.16,1.03);ax.set_aspect('equal');ax.axis('off')
  return ax
 ax=triangle((43,87,74,64))
 tx(80,151,'Overlay',ha='center',fontsize=8.5,fontweight='bold')
 joined=pd.concat([arrays['GLMF1'],arrays['GLMF2']],ignore_index=True)
 order=np.random.default_rng(262260).permutation(len(joined));joined=joined.iloc[order]
 color=np.array([COLORS[s] for s in joined['sample']])
 ax.scatter(joined.ternary_x,joined.ternary_y,s=5.5,c=color,alpha=.22,lw=0,rasterized=True,zorder=1)
 handles=[Line2D([0],[0],marker='o',linestyle='None',color=COLORS[s],markersize=4,label=s) for s in ['GLMF1','GLMF2']]
 fig.legend(handles=handles,loc='center left',bbox_to_anchor=(119/width,124/height),frameon=False,fontsize=7.7,handletextpad=.35,borderaxespad=0)
 for k,sid in enumerate(['GLMF1','GLMF2']):
  xx=3+k*80;ax=triangle((xx,16,74,64));d=arrays[sid]
  ax.scatter(d.ternary_x,d.ternary_y,s=6.5,c=COLORS[sid],alpha=.35,lw=0,rasterized=True,zorder=1)
  tx(xx+37,83,f'{sid}  |  {len(d):,} spots',ha='center',fontsize=8.2,fontweight='bold',color=COLORS[sid])
 tx(80,9,'Coordinates: distance shares (%)',ha='center',fontsize=7.5)
 tx(80,4,'Larger share = relatively farther; two patients, descriptive',ha='center',fontsize=7.1,color='#4A4A4A')
 fig.canvas.draw();renderer=fig.canvas.get_renderer();bad=[]
 for obj in fig.findobj(matplotlib.text.Text):
  if not obj.get_visible() or not obj.get_text():continue
  assert obj.get_fontsize()>=6
  bb=obj.get_window_extent(renderer)
  if bb.width and (bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5):bad.append(obj.get_text())
 assert not bad,bad
 base=P/'delivery/26210_Ternary_spatial_preview_v1'
 fig.savefig(str(base)+'.pdf',dpi=450);fig.savefig(str(base)+'.svg',dpi=450)
 fig.savefig(str(base)+'.png',dpi=300);fig.savefig(str(base)+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
 dump('qc/BUILD_CHECK.json',dict(status='PASS',figure_mm=[160,170],patients=2,panels=3,targets=[t[0] for t in TARGETS],normalization='d_i/(d_SPP1+d_Endothelium+d_Tcell)',true_points_only=True,jitter=False,smoothing=False,sampled_points=False,overlay_drawing_order_seed=262260,inferential_tests=0,new_models=0,minimum_font_pt=6.4,text_bounds_pass=True))
 dump('STATUS.json',dict(status='BUILT_AWAITING_VISUAL_QA',time=datetime.datetime.now().astimezone().isoformat()))
 print('BUILD_PASS',flush=True)
if __name__=='__main__':main()
