from geometry import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from skimage.transform import resize
import csv,traceback

def oldnorm(mask,size=128):
 y,x=np.where(mask);x0=x.min();y0=y.min();h=y.max()-y0+1;w=x.max()-x0+1;s=min((size-16)/h,(size-16)/w);nh=max(1,int(round(h*s)));nw=max(1,int(round(w*s)));sx=nw/w;sy=nh/h;ox=(size-nw)//2;oy=(size-nh)//2
 return np.array([[sx,0,ox+.5*sx-.5-sx*x0],[0,sy,oy+.5*sy-.5-sy*y0],[0,0,1.]])
def oldmatrix(sample,pol,rawmask,sf):
 pos=pd.read_csv(ST/sample/'outs/spatial/tissue_positions_list.csv',header=None,names=['barcode','in_tissue','ar','ac','y','x']);xy=pos[['x','y']].to_numpy();spacing=np.median(cKDTree(xy).query(xy,k=2)[0][:,1]);lo=xy.min(0)-spacing;hi=xy.max(0)+spacing;yy,xx=np.meshgrid(np.linspace(lo[1],hi[1],256),np.linspace(lo[0],hi[0],256),indexing='ij');dist,idx=cKDTree(xy).query(np.c_[xx.ravel(),yy.ravel()]);mask=((pos.in_tissue.to_numpy()[idx]==1)&(dist<.66*spacing)).reshape(256,256);mask=ndimage.binary_closing(mask,iterations=2)
 grid=np.array([[255/(hi[0]-lo[0]),0,-255*lo[0]/(hi[0]-lo[0])],[0,255/(hi[1]-lo[1]),-255*lo[1]/(hi[1]-lo[1])],[0,0,1.]])
 old=pd.read_csv(OLD/'tables/PAIRING_ORIENTATION_DIAGNOSTIC.tsv',sep='\t');r=old[(old['sample']==sample)&(old.polarity==pol)].sort_values('dice',ascending=False).iloc[0];O=np.eye(3)
 if r.flip:O=np.array([[-1.,0,127],[0,1,0],[0,0,1]])@O
 R=np.array([[0.,1,0],[-1,0,127],[0,0,1.]])
 for _ in range(int(r.rotation90)):O=R@O
 Tr=np.array([[1.,0,r.translation_x],[0,1,r.translation_y],[0,0,1.]])
 return np.diag([sf,sf,1.])@np.linalg.inv(oldnorm(mask)@grid)@Tr@O@oldnorm(rawmask)

def main():
 assert (P/'COMPUTE_EXIT_STATUS').read_text().strip()=='0'
 allr=json.loads((P/'transforms/ALL_CANDIDATE_TRANSFORMS.json').read_text());kept=json.loads((P/'transforms/RETAINED_CANDIDATES.json').read_text());held=pd.read_csv(P/'tables/HELDOUT_BOUNDARY_STABILITY.tsv',sep='\t');summary=[];mapping=[];exports=[];baseline=[]
 fig,axes=plt.subplots(4,2,figsize=(10,18));plt.rcParams.update({'font.size':10})
 labels=['R114-T-neg','R114-T-pos','Y12-T-neg','Y12-T-pos']
 for row,label in enumerate(labels):
  sample=label.rsplit('-',1)[0].replace('-','_');pol=label.rsplit('-',1)[1];a=np.load(P/'inputs'/(sample+'_ST.npz'));b=np.load(P/'inputs'/(label+'_SM.npz'));T=a['T'];S=b['S'];sf=float(a['hires_scale']);pitch=float(a['pitch']);rr=sorted([x for x in allr if x['label']==label and x['model']=='similarity' and x['variant']=='HE'],key=lambda x:(-x['dice'],x['surface_mean']));best=rr[0];ret=[x for x in kept if x['label']==label];variants=[max([x for x in allr if x['label']==label and x['model']=='similarity' and x['variant']==v],key=lambda x:x['dice']) for v in ['HE','HE_low','HE_high','provider','erode','dilate']]
  stable=sum(x['flip']==best['flip'] and x['rotation_seed']==best['rotation_seed'] for x in variants);aff=max([x for x in allr if x['label']==label and x['model']=='affine'],key=lambda x:x['dice']);hh=held[(held.label==label)&(held.flip==best['flip'])&(held.rotation_seed==best['rotation_seed'])]
  normimg=warp(a['rgb'],T,order=1);target=a['HE'];positive=b['tic']>0;rawxy=b['xy'][positive];coords=a['xy_hires'];intissue=a['in_tissue']==1;barcodes=a['barcode'];qc_ids=set(pd.read_csv(P.parent/'phase234_part4_replan_execution_20260921_v1/inputs'/(sample+'_INPUTS.tsv.gz'),sep='\t',usecols=['spot_id']).spot_id)
  tree=cKDTree(rawxy);tested=[];newM=np.array(best['matrix']);Hbest=np.array(best['SM_zero_based_to_ST_hires']);oldH=oldmatrix(sample,pol,b['raw_mask'],sf);oldM=T@oldH@np.linalg.inv(S);oldw=warp(b['source'],oldM)>.5
  def dice(w,t):return 2*np.count_nonzero(w&t)/(w.sum()+t.sum())
  olddice=dice(oldw,target);neww=warp(b['source'],newM)>.5
  baseline.append(dict(label=label,old_transform_Dice_on_current_HE=olddice,new_primary_Dice_on_current_HE=dice(neww,target),delta=dice(neww,target)-olddice,old_transform_Dice_on_provider=dice(oldw,a['provider']),new_primary_Dice_on_provider=dice(neww,a['provider']),comparison='same current target masks for old/new transforms'))
  for rank,r in enumerate(rr[:2]):
   ax=axes[row,rank];ax.imshow(normimg);ax.contour(target,levels=[.5],colors=['cyan'],linewidths=.9);w=warp(b['source'],np.array(r['matrix']))>.5;ax.contour(w,levels=[.5],colors=['magenta' if rank==0 else 'orange'],linewidths=1)
   angle=((r['angle_degrees']+180)%360)-180;ax.set_title(f"{label} | {'preferred' if rank==0 else 'runner-up'}\nDice {r['dice']:.3f}; rotation {angle:.1f} deg; mirror {r['flip']}",fontsize=11);ax.axis('off')
  for r in ret:
   H=np.array(r['SM_zero_based_to_ST_hires']);one=np.array([[1.,0,-1],[0,1,-1],[0,0,1.]]);oneH=H@one;full=np.diag([1/sf,1/sf,1.])@oneH
   exported=dict(candidate_id=r['candidate_id'],label=label,rank=r['rank'],status='GEOMETRY_CANDIDATE_NOT_ANATOMICALLY_VALIDATED',SM_zero_based_xy_to_ST_hires_xy=H.tolist(),SM_original_one_based_xy_to_ST_hires_xy=oneH.tolist(),SM_original_one_based_xy_to_ST_fullres_xy=full.tolist(),ST_hires_xy_to_SM_original_one_based_xy=np.linalg.inv(oneH).tolist(),ST_hires_scale=sf,source_coordinates='imzML x/y begin at 1; supplied xy matrix also has explicit zero-based version',target_coordinates='x right y down; ST hires image pixels or fullres as labelled',physical_SM_unit='unknown; no micrometre conversion claimed',dice=r['dice'],angle_clockwise_degrees=((r['angle_degrees']+180)%360)-180,mirror_x_before_rotation=r['flip'],model=r['model'],molecular_targets_used=False)
   (P/'transforms'/(r['candidate_id']+'.json')).write_text(json.dumps(exported,indent=2));exports.append(exported)
   q=apply(coords,np.linalg.inv(H));dist,ix=tree.query(q);native_indices=np.flatnonzero(positive)[ix];within=(dist<=1.5)&intissue;delta=np.linalg.norm(apply(rawxy[ix],H)-coords,axis=1);qc=np.array([x in qc_ids for x in barcodes]);tested.append(dict(candidate=r['candidate_id'],rank=r['rank'],QC_spots=int(qc.sum()),QC_spots_mapped=int((within&qc).sum()),QC_mapped_fraction=float((within&qc).sum()/qc.sum()),tissue_spots_mapped=int(within.sum()),median_distance_hires_pixels=float(np.median(delta[within&qc])),median_distance_SM_grid=float(np.median(dist[within&qc]))))
   frame=pd.DataFrame(dict(sample=sample,polarity=pol,candidate_id=r['candidate_id'],rank=r['rank'],spot_id=barcodes,in_tissue=intissue,QC_pass=qc,ST_x_hires=coords[:,0],ST_y_hires=coords[:,1],SM_continuous_x_one_based=q[:,0]+1,SM_continuous_y_one_based=q[:,1]+1,SM_pixel_index=native_indices,SM_nearest_x_one_based=rawxy[ix,0]+1,SM_nearest_y_one_based=rawxy[ix,1]+1,distance_SM_grid=dist,distance_ST_hires=delta,included=within,validated_registration=False));mapping.append(frame)
   probes=rawxy[np.linspace(0,len(rawxy)-1,min(1000,len(rawxy))).astype(int)];assert np.max(np.abs(apply(apply(probes,H),np.linalg.inv(H))-probes))<1e-8;assert np.allclose(apply(probes+1,oneH),apply(probes,H));assert np.allclose(apply(probes+1,full)*sf,apply(probes,H))
  pd.DataFrame(tested).to_csv(P/'tables'/(label+'_mapping_coverage.tsv'),sep='\t',index=False)
  probes=rawxy[np.linspace(0,len(rawxy)-1,min(2000,len(rawxy))).astype(int)];positions=np.array([apply(probes,np.array(r['SM_zero_based_to_ST_hires'])) for r in ret]);deviation=np.linalg.norm(positions-positions[0],axis=2)/pitch;bestcov=tested[0]
  summary.append(dict(label=label,sample=sample,polarity=pol,best_Dice=best['dice'],runner_up_Dice=rr[1]['dice'],orientation_Dice_margin=best['dice']-rr[1]['dice'],preferred_mirror=best['flip'],preferred_clockwise_degrees=((best['angle_degrees']+180)%360)-180,retained_candidates=len(ret),preferred_orientation_wins_mask_variants=stable,total_mask_variants=6,affine_best_Dice=aff['dice'],primary_surface_mean_spot_pitches=best['surface_mean_spot_pitches'],primary_surface_p95_spot_pitches=best['surface_p95_spot_pitches'],heldout_max_median_shift_spot_pitches=float(hh.mapping_displacement_median_spot_pitches.max()),heldout_max_p95_shift_spot_pitches=float(hh.mapping_displacement_p95_spot_pitches.max()),alternative_max_median_shift_spot_pitches=float(np.median(deviation,axis=1).max()),alternative_max_p95_shift_spot_pitches=float(np.percentile(deviation,95,axis=1).max()),QC_spots=bestcov['QC_spots'],QC_spots_mapped=bestcov['QC_spots_mapped'],QC_mapped_fraction=bestcov['QC_mapped_fraction'],SM_foreground_on_acquisition_border=int(b['raw_mask'][0].sum()+b['raw_mask'][-1].sum()+b['raw_mask'][:,0].sum()+b['raw_mask'][:,-1].sum()),validated_registration=False))
 fig.legend(handles=[Line2D([0],[0],color='cyan',label='H&E tissue boundary'),Line2D([0],[0],color='magenta',label='Preferred MSI contour'),Line2D([0],[0],color='orange',label='Alternative MSI contour')],loc='lower center',ncol=3,frameon=False);fig.suptitle('Morphology-only registration diagnostics; anatomical correspondence unverified',fontsize=13);fig.tight_layout(rect=[0,.035,1,.975]);fig.savefig(P/'diagnostics/registration_overlays.png',dpi=160);plt.close(fig)
 pd.DataFrame(summary).to_csv(P/'tables/REGISTRATION_SUMMARY.tsv',sep='\t',index=False);pd.DataFrame(baseline).to_csv(P/'tables/SAME_TARGET_BASELINE_COMPARISON.tsv',sep='\t',index=False);pd.concat(mapping,ignore_index=True).to_csv(P/'transforms/CANDIDATE_ST_TO_SM_MAPS.tsv.gz',sep='\t',index=False);(P/'transforms/COORDINATE_CONVENTIONS.json').write_text(json.dumps(dict(SM='original imzML one-based integer x/y; zero-based matrices clearly identified',ST='hires and fullres image pixel frames, x right y down',candidate_count=len(exports),validated_registration=False),indent=2));print(pd.DataFrame(summary).to_string(index=False));print('BASELINE',pd.DataFrame(baseline).to_string(index=False))
if __name__=='__main__':main()
