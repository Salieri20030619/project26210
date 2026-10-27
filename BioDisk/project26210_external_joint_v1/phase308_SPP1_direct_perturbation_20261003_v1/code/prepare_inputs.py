from pathlib import Path
import requests,json,hashlib,datetime,os,re,tarfile,gzip,io,shutil,time
P=Path(__file__).resolve().parents[1];R=P.parent
def session():
 s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'};return s
proto={'created':datetime.datetime.now().astimezone().isoformat(),'dataset':'GSE274609','experiment':'WT6620c1PDAC cells orthotopically implanted in WT or globalhostSpp1KO mice; notmacrophage-specific andnotRCC. AuthorreportsapCAFs majorstromalSPP1source. TumorcellsnotgeneticallyKO.',
'unit':'3biologicalpools/group; eachpoolcontains2tumors,6tumors/group. Do not claim6individuallymeasuredmice/group,12independenttranscriptomes,orcell-levelP. Unpairedgroupcomparison; same-numberKO/WTnotpaired.',
'QC':'Perlibrary integercounts, unique barcodes. min200genes/min500UMIs; <=6000genes; mitochondrial<=20%; upperUMI99.5percentile. Scrublet expected_doublet_rate0.06 perlibrary; exclude predicted doublets. NoSpp1/Tymp/effectselection.',
'identity':'Log1pCP10K,HVG3000 excludingSpp1/Tymp,PCA40/neighbors20/Leiden0.7. Review marker tables before freezing lineages. Purefibro/mural/mesothelial-like andmacro/DC/monocyte separated; no pureCAFclaimwithoutmarkers. TargetSpp1doesnotdefinebroadlineage.',
'primary_endpoints':['CapturedCD8fractionofallQCcells','CapturedCD8fractionofimmune','Fibroblast_collagen_maturation','Fibroblast_focal_adhesion','Fibroblast_actomyosin'],
'secondary_endpoints':['CapturedTfraction','CD8cytotoxicity','CD8dysfunction','Macrophage_Cxcl9_Cxcl10','Fibroblast_Postn_Acta2_Lrrc15','Fibroblast_Pi16_Dpt','Spp1_and_Tymp_lineage_expression'],
'programs':'ReusePhase292frozenstromalpanels translatedviaMGI1to1 orthology; CD8cytoGzmbPrf1Nkg7Ctsw, dysfunctionPdcd1Ctla4Lag3Havcr2TigitToxCxcl13Entpd1. Mouse lacksdirect1to1GNLY; recordcoverage. myCAF3andresting2 fromsourceexperiment assecondaryreadouts, notnewclusterdefinitions.',
'aggregation':'RawRNApseudobulk perpool/lineage tolog1pCP10K, equalweightfixedgenes. Lineage>=100macro/CD8 and>=20fibro. >=80%panelmeasured; >=50%panelgenes detectedin>=5cells. Sourcegenes>=5detectedcells. Fractions useallQCandmarkedimmune denominators; capturedcomposition notabsolute infiltration.',
'statistics':'Descriptiveonly:3pools/arm below10. Reportallpoolvalues,medianKOminusmedianWT,ratioofmedianswheredefined,between-groupoverlap andall9cross-pooldirectioncomparisons(namedcomparisons,notn). No biologicalP,FDR,cellDEscreenorpooledcellP.',
'boundaries':'HostSPP1dependence externalexperiment; noTYMPperturbation,noRCCcausality,noMacrophage-specificattribution ordirectcollagen geometry. Allnegative andreverseendpointspreserved. Noextrafigures.'}
f=P/'resources/PERTURBATION_PROTOCOL.json'
if not f.exists():f.write_text(json.dumps(proto,indent=2))
samples=[]
for block in (P/'metadata/GSE274609.soft').read_text().split('^SAMPLE = ')[1:]:
 gsm=block.splitlines()[0].strip();title=re.search(r'!Sample_title = (.+)',block).group(1);group='KO' if 'knockout' in title else 'WT';rep=int(re.search(r'replicate (\d+)',title).group(1));samples.append({'GSM':gsm,'library':group+str(rep),'group':group,'replicate':rep,'pool_tumors':2,'paired':False})
assert len(samples)==6;import pandas as pd
pd.DataFrame(samples).to_csv(P/'tables/SAMPLE_PLAN.tsv',sep='\t',index=False)
stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S');backup=P/'resources'/('root_before_start_'+stamp);backup.mkdir()
message='## Phase308进行中（2026-10-03）：启动直接宿主SPP1敲除单细胞验证\n用户要求继续下一阶段。已找到GSE274609：WT6620c1肿瘤植入WT或全身宿主Spp1KO；每组6肿瘤，两肿瘤混一库，实际3池/组，不是12独立转录组。原文提示apCAF主要SPP1来源，因此不称巨噬特异敲除。优先评估纯成纤维基质程序及CD8组成/功能；仅描述方向效应，不计算细胞P。665MB归档按代理规则获取；QC/端点先冻结，身份复核后再计算效应。JITC巨噬特异Spp1研究仅按需提供，GSE316918为滑膜OPN刺激情境，留来源审计。新目录phase308_SPP1_direct_perturbation_20261003_v1；旧C4/对象/所有阴性与阳性保持。\n\n'
for n in ['CODEX_HANDOFF.md','PROJECT_STATUS.md','DECISIONS.md']:
 f=R/n;old=f.read_bytes();(backup/n).write_bytes(old);assert f.read_bytes()==old;tmp=f.with_name(f.name+'.phase308.tmp');tmp.write_bytes(message.encode()+old);os.replace(tmp,f)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274609/suppl/GSE274609_RAW.tar';s=session();head=s.head(url,timeout=45);head.raise_for_status();size=int(head.headers.get('Content-Length',0));assert 0<size<2_000_000_000;dest=P/'downloads/GSE274609_RAW.tar';partial=dest.with_suffix('.tar.partial')
print('DOWNLOAD_PLAN',size,flush=True)
if not dest.exists():
 for attempt in range(4):
  try:
   rr=s.get(url,stream=True,timeout=(25,90));rr.raise_for_status()
   with partial.open('wb') as out:
    for b in rr.iter_content(8*1024*1024):
     if b:out.write(b)
   assert partial.stat().st_size==size;os.replace(partial,dest);break
  except Exception as e:
   print('RETRY',attempt,type(e).__name__,flush=True)
   if attempt==3:raise RuntimeError('Matrix download failed') from None
assert dest.stat().st_size==size;h=hashlib.sha256()
with dest.open('rb') as inf:
 for b in iter(lambda:inf.read(8*1024*1024),b''):h.update(b)
inner=[]
with tarfile.open(dest) as t:
 members=[m for m in t if m.isfile()];assert len(members)==6
 for m in members:
  gsm=m.name.split('_')[0];sample=next(x for x in samples if x['GSM']==gsm);data=t.extractfile(m).read();name=Path(m.name).name;f=P/'downloads'/name
  if not f.exists():f.write_bytes(data)
  decompressed=gzip.decompress(data);lib=P/'inputs'/sample['library'];lib.mkdir(exist_ok=True)
  with tarfile.open(fileobj=io.BytesIO(decompressed)) as sub:
   files=[x for x in sub if x.isfile()];kept=[]
   for item in files:
    fn=Path(item.name).name
    if fn in ['matrix.mtx.gz','features.tsv.gz','genes.tsv.gz','barcodes.tsv.gz','matrix.mtx','features.tsv','genes.tsv','barcodes.tsv']:
     raw=sub.extractfile(item).read()
     if fn.endswith('.gz'):gzip.decompress(raw)
     (lib/fn).write_bytes(raw);kept.append(fn)
   assert len(kept)==3,kept
  inner.append({'GSM':gsm,'library':sample['library'],'archive':name,'sha256':hashlib.sha256(data).hexdigest(),'members':kept});print('INPUT_VERIFIED',sample['library'],kept,flush=True)
(P/'qc/INPUT_ARCHIVES_VALIDATED.json').write_text(json.dumps({'status':'PASS','url':url,'route':'proxy','bytes':size,'sha256':h.hexdigest(),'http_etag':head.headers.get('ETag'),'gzipCRC':'all6innerarchivesandcompressedmatricespassed','libraries':inner},indent=2));(P/'logs/prepare_inputs.exit').write_text('0\n');(P/'STATUS.json').write_text(json.dumps({'status':'INPUTS_READY_FOR_QC','time':datetime.datetime.now().astimezone().isoformat()},indent=2));print('INPUTS_COMPLETE',size,h.hexdigest(),flush=True)
