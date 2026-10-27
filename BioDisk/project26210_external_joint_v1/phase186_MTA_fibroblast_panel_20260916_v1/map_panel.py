from pathlib import Path
import json,struct,subprocess,concurrent.futures
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase186_MTA_fibroblast_panel_20260916_v1')
b=(r/'sources/zip_tail.bin').read_bytes();base=30392672;p=30504747-base;rows=[]
while b[p:p+4]==b'PK\x01\x02':
 a=struct.unpack_from('<4s6H3L5H2L',b,p);nl,el,cl=a[10:13];name=b[p+46:p+46+nl].decode('utf-8')
 rows.append({'name':name,'compressed':a[8],'uncompressed':a[9],'method':a[4],'crc32':a[7],'offset':a[16]});p+=46+nl+el+cl
assert len(rows)==170
(r/'ZIP_DIRECTORY.json').write_text(json.dumps(rows,indent=2))
print('DATA_FILES',json.dumps([x for x in rows if not x['name'].endswith('/') and not '__MACOSX' in x['name']],indent=2)[:12000],flush=True)
panel=['HGF','TGFB1','TGFB2','TGFB3','IL6','CSF1','WNT5A','PTGS2','PTGES','SPP1','MTAP']
(r/'PRESET_PANEL.json').write_text(json.dumps({'genes':panel,'fixed_from':'Phase185 section5','normalization':'median-ratio primary; upper-quartile sensitivity','statistics':'unpaired descriptive group medians/ranges; no P/FDR; no independent donor assertion'},indent=2))
def get(g):
 p=r/'sources'/f'ensembl_{g}.json'
 z=subprocess.run(['curl','-fLsS','--connect-timeout','10','--max-time','30',f'https://rest.ensembl.org/lookup/symbol/homo_sapiens/{g}?content-type=application/json','-o',str(p)],capture_output=True,text=True)
 if z.returncode==0:
  a=json.loads(p.read_text());return {'symbol':g,'ensembl_id':a['id'],'display_name':a['display_name']}
 return {'symbol':g,'error':z.stderr}
out=list(concurrent.futures.ThreadPoolExecutor(max_workers=6).map(get,panel))
(r/'GENE_MAPPING.json').write_text(json.dumps(out,indent=2));print('MAPPING',json.dumps(out))
