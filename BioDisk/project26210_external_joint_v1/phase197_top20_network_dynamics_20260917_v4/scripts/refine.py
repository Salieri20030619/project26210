from pathlib import Path
d=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase197_top20_network_dynamics_20260917_v4')
p=d/'scripts/redraw.py';s=p.read_text().replace("figsize=(8.7,6.4)","figsize=(8.7,8.8)")
s+="""
# Short, wide tissue panel: unchanged 255 samples and row values.
bulk=pd.read_csv(r/'phase192_sections28_212_figures_20260917_v1/tables/TISSUE_MTA_ORDERED_PATIENTS.tsv',sep='\\t',index_col=0)
features=['MET:5-methylthioadenosine (MTA)','SPP1','LST1','AIF1','CSF1R','CD68','FCER1G','SMS','MTAP']
fig=plt.figure(figsize=(12.5,3.2));start=0
for i,(co,dat) in enumerate(bulk.groupby('cohort',sort=False)):
 ax=fig.add_axes([.07+.86*start/255,.25,.86*len(dat)/255-.006,.62])
 vv=dat[features].values.T.copy();vv[0]=np.log2(vv[0]+1)
 zz=(vv-vv.mean(1,keepdims=True))/vv.std(1,keepdims=True)
 im=ax.imshow(zz,aspect='auto',cmap=cm,vmin=-2.5,vmax=2.5,interpolation='nearest')
 ax.set_title(co.replace('_',' ')+f' | n = {len(dat)}',fontsize=9)
 ax.set_yticks(range(9),['MTA']+features[1:] if i==0 else ['']*9,fontsize=7)
 ax.set_xticks([]);ax.tick_params(length=0)
 ax.add_patch(Rectangle((-.5,.5),len(dat),1,fill=False,ec=purple,lw=.9))
 start+=len(dat)
fig.text(.5,.14,'Patients ordered by measured MTA within each cohort',ha='center',fontsize=8)
ca=fig.add_axes([.951,.30,.011,.48]);fig.colorbar(im,cax=ca,label='Within-cohort row z-score',ticks=[-2,0,2]);ca.tick_params(labelsize=7)
panel(fig,'tissue_compact')
"""
p.write_text(s)
p=d/'scripts/assemble.py';s=p.read_text().replace("'expression_dynamics_curves']:","'expression_dynamics_curves','tissue_compact']:")
s=s.replace("'Measured MTA and tissue expression','2G',8,ye,324","'Measured MTA and tissue expression','tissue_compact',8,ye,324")
p.write_text(s)
