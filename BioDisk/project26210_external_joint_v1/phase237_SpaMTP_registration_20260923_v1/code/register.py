from geometry import *
import concurrent.futures,traceback,os
CACHE={}
def data(label):
 if label not in CACHE:
  sample=label.rsplit('-',1)[0].replace('-','_');a=np.load(P/'inputs'/(sample+'_ST.npz'));b=np.load(P/'inputs'/(label+'_SM.npz'));CACHE[label]=dict(st={k:a[k] for k in ['HE','HE_low','HE_high','provider','T','pitch']},sm={k:b[k] for k in ['source','erode','dilate','S','xy','raw_mask']})
 return CACHE[label]
def task(job):
 label,variant,flip,k,model=job;d=data(label);target=d['st'][variant] if variant in ['HE_low','HE_high','provider'] else d['st']['HE'];source=d['sm'][variant] if variant in ['erode','dilate'] else d['sm']['source'];g=geometry(source,target)
 res=fit(g,flip,k,affine=model=='affine');res.update(label=label,sample=label.rsplit('-',1)[0].replace('-','_'),polarity=label.rsplit('-',1)[1],variant=variant)
 H=np.linalg.inv(d['st']['T'])@np.array(res['matrix'])@d['sm']['S'];res['SM_zero_based_to_ST_hires']=H.tolist();res['surface_mean_hires_pixels']=res['surface_mean']/d['st']['T'][0,0];res['surface_mean_spot_pitches']=res['surface_mean_hires_pixels']/float(d['st']['pitch']);res['surface_p95_spot_pitches']=res['surface_p95']/d['st']['T'][0,0]/float(d['st']['pitch']);return res

def holdout(job):
 ref,fold=job;d=data(ref['label']);g=geometry(d['sm']['source'],d['st']['HE']);M=np.array(ref['matrix']);sa=apply(g['sp'],M)-g['tc'];ta=g['tp']-g['tc'];ss=(np.floor((np.arctan2(sa[:,1],sa[:,0])+np.pi)/(np.pi/4)).astype(int)%8)%4;ts=(np.floor((np.arctan2(ta[:,1],ta[:,0])+np.pi)/(np.pi/4)).astype(int)%8)%4
 fitr=fit(g,ref['flip'],ref['rotation_seed'],initial=np.array(ref['parameters']),train=(g['sp'][ss!=fold],g['tp'][ts!=fold]));A=np.array(fitr['matrix']);da,db=distance_values(A,g,g['sp'][ss==fold],g['tp'][ts==fold]);test=np.r_[da,db];dis=np.linalg.norm(apply(g['sp'],A)-apply(g['sp'],M),axis=1);conversion=d['st']['T'][0,0]*float(d['st']['pitch'])
 return dict(label=ref['label'],flip=ref['flip'],rotation_seed=ref['rotation_seed'],fold=fold,training_source_points=int((ss!=fold).sum()),training_target_points=int((ts!=fold).sum()),heldout_boundary_mean_spot_pitches=float(test.mean()/conversion),heldout_boundary_p95_spot_pitches=float(np.percentile(test,95)/conversion),mapping_displacement_median_spot_pitches=float(np.median(dis)/conversion),mapping_displacement_p95_spot_pitches=float(np.percentile(dis,95)/conversion),angle_degrees=fitr['angle_degrees'],dice=fitr['dice'],matrix=fitr['matrix'])

def emit(event,**kw):
 item=dict(time=datetime.datetime.now().astimezone().isoformat(),event=event,**kw);print(json.dumps(item),flush=True)
 with (P/'logs/events.jsonl').open('a') as f:f.write(json.dumps(item)+'\n')
def main():
 # Numerical check for raster x/y conventions and mirrored affine inversion.
 test=np.zeros((N,N),bool);test[55:195,55:150]=True;test[145:215,150:195]=True;test[75:105,150:180]=True
 angle=np.deg2rad(100);A=np.eye(3);A[:2,:2]=1.03*np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])@np.diag([-1,1]);center=np.array([127.5,127.5]);A[:2,2]=center+np.array([2,-3])-A[:2,:2]@center
 tar=warp(test,A)>.5;g=geometry(test,tar);r=fit(g,True,1);B=np.array(r['matrix']);pts=g['sp'];error=np.median(np.linalg.norm(apply(pts,A)-apply(pts,B),axis=1));assert r['dice']>.97 and error<2,(r,error)
 roundtrip=np.max(np.abs(apply(apply(pts,A),np.linalg.inv(A))-pts));assert roundtrip<1e-9
 (P/'metadata/TRANSFORM_NUMERICAL_QA.json').write_text(json.dumps(dict(status='PASS',synthetic_dice=r['dice'],known_transform_median_error_pixels=float(error),inverse_roundtrip_max_error=float(roundtrip)),indent=2))
 labels=['R114-T-neg','R114-T-pos','Y12-T-neg','Y12-T-pos'];jobs=[(l,v,f,k,'similarity') for l in labels for v in ['HE','HE_low','HE_high','provider','erode','dilate'] for f in [False,True] for k in range(4)];jobs += [(l,'HE',f,k,'affine') for l in labels for f in [False,True] for k in range(4)]
 results=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
  futures={pool.submit(task,j):j for j in jobs}
  for fut in concurrent.futures.as_completed(futures):
   r=fut.result();results.append(r)
   if len(results)%8==0:
    emit('REGISTRATION_PROGRESS',finished=len(results),total=len(jobs));(P/'STATUS.json').write_text(json.dumps(dict(status='REGISTRATION_RUNNING',finished=len(results),total=len(jobs)),indent=2))
 (P/'transforms/ALL_CANDIDATE_TRANSFORMS.json').write_text(json.dumps(results,indent=2));pd.DataFrame([{k:v for k,v in r.items() if k not in ['matrix','parameters','SM_zero_based_to_ST_hires']} for r in results]).to_csv(P/'tables/ALL_REGISTRATION_CANDIDATES.tsv',sep='\t',index=False)
 candidates=[];holdjobs=[]
 for label in labels:
  rr=sorted([r for r in results if r['label']==label and r['model']=='similarity' and r['variant']=='HE'],key=lambda r:(-r['dice'],r['surface_mean']));best=rr[0]
  kept=[r for r in rr if r['dice']>=best['dice']-.03]
  for rank,r in enumerate(kept,1):r=dict(r,rank=rank,candidate_id=label+'_f'+str(int(r['flip']))+'_r'+str(r['rotation_seed']));candidates.append(r)
  for r in rr[:2]:
   for fold in range(4):holdjobs.append((r,fold))
  emit('PRIMARY_CANDIDATES',label=label,best_dice=best['dice'],best_angle=best['angle_degrees'],best_flip=best['flip'],retained=len(kept),second_dice=rr[1]['dice'])
 (P/'transforms/RETAINED_CANDIDATES.json').write_text(json.dumps(candidates,indent=2))
 held=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
  for r in pool.map(holdout,holdjobs):held.append(r)
 (P/'transforms/HELDOUT_TRANSFORMS.json').write_text(json.dumps(held,indent=2));pd.DataFrame([{k:v for k,v in r.items() if k!='matrix'} for r in held]).to_csv(P/'tables/HELDOUT_BOUNDARY_STABILITY.tsv',sep='\t',index=False)
 (P/'REGISTRATION_COMPUTED').write_text(datetime.datetime.now().astimezone().isoformat()+'\n');emit('REGISTRATION_COMPUTED',candidates=len(results),retained=len(candidates),heldout_refits=len(held))
code=1
if __name__=='__main__':
 try:main();code=0
 except BaseException:(P/'logs/ERROR.txt').write_text(traceback.format_exc());print(traceback.format_exc(),flush=True)
 finally:(P/'COMPUTE_EXIT_STATUS').write_text(str(code)+'\n')
 raise SystemExit(code)
