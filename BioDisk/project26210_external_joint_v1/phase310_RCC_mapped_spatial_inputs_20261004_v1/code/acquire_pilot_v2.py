from pathlib import Path
import subprocess,json,hashlib,gzip,csv,datetime,collections,traceback,os
P=Path(__file__).resolve().parents[1]
a=json.loads((P/'metadata/figshare_current.json').read_text())
rows=[]
for f in a['files']:
 region=a['folder_structure'][str(f['id'])].split('/')[-1];group=region.split('_')[0]
 rows.append(dict(file_id=f['id'],region=region,tumor_label=group,size=f['size'],md5=f['computed_md5'],url=f['download_url'],donor_identity='unresolved_three_labels_vs_four_patients'))
with (P/'tables/FILE_REGION_MAP.tsv').open('w') as o:
 w=csv.DictWriter(o,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
# Smallest file is a technical schema/panel pilot, never a representative biological region.
f=min(rows,key=lambda x:x['size']);out=P/'downloads'/('pilot_'+f['region']+'.csv.gz')
status={'scope':'technical input qualification only','selection':'smallest bytes, chosen before any expression/geometry','record':f,'time':datetime.datetime.now().astimezone().isoformat()}
(P/'resources/PILOT_INPUT_PROTOCOL.json').write_text(json.dumps(status,indent=2))
rc=1
try:
 cmd=['curl','--proxy','http://127.0.0.1:7897','--connect-timeout','20','--max-time','600','--retry','2','--retry-all-errors','--retry-delay','3','-fL','-sS','-o',str(out)+'.partial','https://api.figshare.com/v2/file/download/'+str(f['file_id'])]
 print('DOWNLOAD_START',f,flush=True)
 subprocess.run(cmd,check=True)
 tmp=Path(str(out)+'.partial');assert tmp.stat().st_size==f['size']
 md5=hashlib.md5(tmp.read_bytes()).hexdigest();assert md5==f['md5'];os.replace(tmp,out)
 genes=collections.Counter();assigned=collections.Counter();n=0;heads=[]
 with gzip.open(out,'rt') as z:
  reader=csv.DictReader(z);columns=reader.fieldnames
  for row in reader:
   if n<5:heads.append(row)
   genes[row.get('feature_name',row.get('gene','UNKNOWN'))]+=1
   assigned[row.get('cell_id','MISSING')=='UNASSIGNED']+=1;n+=1
 result={'status':'PASS','file':str(out),'bytes':out.stat().st_size,'md5':md5,'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'gzip_full_read':'PASS','rows':n,'columns':columns,'head':heads,'features':dict(genes),'unassigned_flag_counts':dict(assigned)}
 (P/'qc/PILOT_INPUT_VALIDATION.json').write_text(json.dumps(result,indent=2))
 print('PILOT_COMPLETE',n,columns,dict(genes),flush=True);rc=0
except Exception:
 traceback.print_exc()
finally:
 (P/'logs/acquire_pilot_v2.exit').write_text(str(rc)+'\n')
raise SystemExit(rc)
