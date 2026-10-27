import pathlib,subprocess,json,hashlib,sys,datetime
C=pathlib.Path(__file__).resolve().parent
meta=json.loads((C/'OFFICIAL_ARCHIVE_METADATA.json').read_text())
q=C/'SRR32389763.sra.part'
cmd=['curl','-q','--fail','--location','--silent','--show-error','--connect-timeout','60','--speed-limit','1024','--speed-time','300','--retry','10','--retry-delay','15','--retry-all-errors','--continue-at','-','--proxy','http://127.0.0.1:7897','--noproxy','','--output',str(q),meta['locations'][0]['link']]
r=subprocess.run(cmd)
qa=dict(download_exit_code=r.returncode,expected_bytes=meta['size'],expected_md5=meta['md5'],time=datetime.datetime.now().astimezone().isoformat())
if q.exists():
 h=hashlib.md5()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 qa.update(observed_bytes=q.stat().st_size,observed_md5=h.hexdigest())
ok=r.returncode==0 and qa.get('observed_bytes')==meta['size'] and qa.get('observed_md5')==meta['md5']
qa['status']='PASS' if ok else 'FAIL'
if ok:q.rename(C/'SRR32389763.sra')
(C/'ARCHIVE_QA.json').write_text(json.dumps(qa,indent=2));(C/'DOWNLOAD_EXIT_STATUS').write_text('0\n' if ok else '1\n');print(json.dumps(qa),flush=True)
sys.exit(0 if ok else 1)
