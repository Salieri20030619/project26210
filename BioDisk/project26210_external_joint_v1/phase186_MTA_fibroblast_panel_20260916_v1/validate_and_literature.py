from pathlib import Path
import subprocess,xml.etree.ElementTree as E,pandas as pd,numpy as np,json,statistics
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase186_MTA_fibroblast_panel_20260916_v1')
z=subprocess.run(['curl','-fLsS','--connect-timeout','10','--max-time','35','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12823564/fullTextXML','-o',str(r/'sources/PMC12823564.xml')],capture_output=True,text=True)
print('LITERATURE',z.returncode,z.stderr,flush=True)
if z.returncode==0:
 x=E.parse(r/'sources/PMC12823564.xml')
 for e in x.iter('p'):
  t=' '.join(e.itertext())
  if any(k in t for k in ['IL-6 treatment','recombinant IL-6','Data availability','GSE','All experiments','independent experiments']):print(t[:2600])
 for e in x.iter('ref'):
  if e.attrib.get('id') in ['CR31','CR6']:print('REFERENCE',e.attrib,' '.join(e.itertext()))
s=pd.read_csv(r/'PANEL_SAMPLE_VALUES.tsv',sep='\t');a=pd.read_csv(r/'PANEL_EFFECTS.tsv',sep='\t')
checks=[];loo=[]
for row in a.to_dict('records'):
 g=s[s.gene==row['gene']];col=row['method']+'_normalized'
 cc=np.log2(g[g.group=='Control'][col].to_numpy()+1);tt=np.log2(g[g.group=='MTA'][col].to_numpy()+1)
 d=statistics.median(tt)-statistics.median(cc);assert abs(d-row['delta_log2_median'])<1e-10;checks.append(True)
 ds=[np.median(np.delete(tt,i))-np.median(np.delete(cc,j)) for i in range(4) for j in range(4)]
 loo.append({'gene':row['gene'],'method':row['method'],'leave_one_library_each_group_delta_min':min(ds),'max':max(ds)})
pd.DataFrame(loo).to_csv(r/'LEAVE_ONE_LIBRARY_SENSITIVITY.tsv',sep='\t',index=False)
pivot=a.pivot(index='gene',columns='method',values='delta_log2_median')
print(pivot.to_string())
print('DIRECTION_CONCORDANCE',int((np.sign(pivot.median_ratio)==np.sign(pivot.upper_quartile)).sum()))
print(pd.read_csv(r/'LIBRARY_QC.tsv',sep='\t').to_string(index=False))
(r/'INDEPENDENT_RECALCULATION.json').write_text(json.dumps({'group_median_checks':len(checks),'all_pass':all(checks),'direction_concordance':int((np.sign(pivot.median_ratio)==np.sign(pivot.upper_quartile)).sum())},indent=2))
