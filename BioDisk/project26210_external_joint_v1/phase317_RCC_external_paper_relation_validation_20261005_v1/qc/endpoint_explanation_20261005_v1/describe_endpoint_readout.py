from pathlib import Path
import pandas as pd,numpy as np,json,subprocess,datetime
P=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
reg=pd.read_csv(P/'resources/EXTERNAL_SECTION_REGISTER.tsv',sep='\t');rows=[]
for a in reg.itertuples():
 subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],check=True,stdout=subprocess.DEVNULL)
 d=pd.read_csv(P/'inputs'/a.sample/'METADATA.tsv',sep='\t',float_precision='round_trip')
 valid=d.boundary_input_valid.astype(str).str.lower().eq('true');inside=valid&d.inside_morphology_supported.astype(str).str.lower().eq('true');outside=valid&~inside
 if inside.sum()<50 or outside.sum()<50:continue
 for gene in ['CD3D','CD8A']:
  x=d.loc[inside,'g_'+gene].to_numpy();y=d.loc[outside,'g_'+gene].to_numpy();fx=(x>0).mean();fy=(y>0).mean();mx=x.mean();my=y.mean();px=x[x>0].mean() if (x>0).any() else np.nan;py=y[y>0].mean() if (y>0).any() else np.nan
  assert np.isclose(mx,fx*px) and np.isclose(my,fy*py)
  rows.append({'sample':a.sample,'dataset':a.dataset,'patient':a.patient,'preparation':a.preparation,'gene':gene,'mean_delta':mx-my,'detection_delta':fx-fy,'positive_spot_mean_delta':px-py,'log2_median_library_ratio':np.log2(d.loc[inside,'library_total'].median()/d.loc[outside,'library_total'].median()),'log2_median_feature_ratio':np.log2(d.loc[inside,'detected_genes'].median()/d.loc[outside,'detected_genes'].median())})
z=pd.DataFrame(rows);z.to_csv(O/'SECTION_DIAGNOSTICS.tsv',sep='\t',index=False)
pa=z.groupby(['dataset','patient','gene'],as_index=False).agg({c:'mean' for c in ['mean_delta','detection_delta','positive_spot_mean_delta','log2_median_library_ratio','log2_median_feature_ratio']});pa.to_csv(O/'PATIENT_DESCRIPTIVE_DIAGNOSTICS.tsv',sep='\t',index=False)
out=[]
for key,d in pa.groupby(['dataset','gene']):
 out.append({'dataset':key[0],'gene':key[1],'patients':len(d),'inside_library_higher':int(d.log2_median_library_ratio.gt(0).sum()),'median_library_ratio':float(2**d.log2_median_library_ratio.median()),'inside_feature_higher':int(d.log2_median_feature_ratio.gt(0).sum()),'median_feature_ratio':float(2**d.log2_median_feature_ratio.median()),'mean_lower_but_detection_higher':int((d.mean_delta.lt(0)&d.detection_delta.gt(0)).sum()),'inside_positive_spot_intensity_lower':int(d.positive_spot_mean_delta.lt(0).sum())})
(O/'SUMMARY.json').write_text(json.dumps({'time':datetime.datetime.now().astimezone().isoformat(),'type':'descriptive explanation audit only; no inferential tests or model changes','summary':out},indent=2));print(json.dumps(out,indent=2))
