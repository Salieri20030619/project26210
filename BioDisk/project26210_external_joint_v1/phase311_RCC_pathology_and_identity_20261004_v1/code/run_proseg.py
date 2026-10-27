from pathlib import Path
import json,datetime,subprocess,hashlib,os,traceback
P=Path(__file__).resolve().parents[1];src=P.parent/'phase310_RCC_mapped_spatial_inputs_20261004_v1/downloads/pilot_tumor1_region6.csv.gz'
out=P/'resources/proseg_region6_v1';out.mkdir(exist_ok=False)
exe=P/'resources/proseg206/bin/proseg'
cmd=[str(exe),'--xenium','--cell-id-unassigned=-1','--min-qv','20','--excluded-genes','(?i)^(NegControl|NegativeControl|Blank|Unassigned|antisense)','--nthreads','8','--output-path',str(out),'--output-maxpost-counts','maxpost-counts.csv.gz','--output-gene-metadata','gene-metadata.csv.gz',str(src)]
protocol={'stage':'single-region independent segmentation and lineage qualification','frozen_at':datetime.datetime.now().astimezone().isoformat(),'raw_input':str(src),'raw_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'command':cmd,'software':'proseg 2.0.6 pinned official Bioconda SHA256','sampling':'author defaults 150/150/300 schedule, 100 recorded samples; stochastic, no exposed seed','replicate_unit':'one named region only; donor identity unresolved, no population inference','original_region_selection':'smallest file technical pilot before expression; not representative donor/tumor sampling','segmentation':'all observed noncontrol genes; target excluded from subsequent identity PCA and clustering; not independent measurement','identity_plan':{'matrix':'max posterior integer point estimate','QC':{'counts_min':20,'detected_genes_min':5},'normalization':'total-count 1e4 log1p','features':'all measured noncontrol genes excluding TYMP and SPP1, standardize clipped at10','PCA':30,'neighbors':15,'leiden_resolution':0.6,'random_seed':26210,'annotation':'inspect marker coherence, keep ambiguous/unknown; no malignant identity based solely on EPCAM; no CAF from collagen alone','markers':{'T':['CD3D','CD3E','CD3G','TRAC','CD247'],'CD8':['CD8A','CD8B'],'myeloid':['CD68','C1QA','C1QB','C1QC','CD163','AIF1','LYZ','FCGR3A'],'fibroblast':['COL1A1','COL1A2','COL3A1','LUM','PDGFRA','DPT','FBLN1'],'mural':['MYH11','ACTA2','ACTG2','TAGLN','MYL9','PDGFRB'],'endothelial':['PECAM1','VWF','KDR','EMCN','EGFL7','PLVAP','MMRN2'],'epithelial':['EPCAM','KRT8','KRT7','KRT18','KRT19','CDH1'],'RCC_candidate':['CA9','EGLN3','NDUFA4L2'],'renal_normal':['AQP1','FXYD2','CLIC6'],'B':['MS4A1','CD79A','CD79B','CD19'],'plasma':['MZB1','JCHAIN','IGKC','IGHG1'],'mast':['TPSAB1','CPA3','KIT']}},'effect_tests':'none until independent identity and spatial endpoint gates frozen; no biological P/FDR for this region','old_assets':'readonly'}
(P/'resources/SEGMENTATION_IDENTITY_PROTOCOL.json').write_text(json.dumps(protocol,indent=2))
rc=1
try:
 print('START',protocol['frozen_at'],cmd,flush=True)
 env=os.environ.copy();env.update({'OMP_NUM_THREADS':'8','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
 with (P/'logs/proseg_raw.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
 rc=r.returncode;print('PROSEG_EXIT',rc,datetime.datetime.now().astimezone().isoformat(),flush=True)
except Exception:traceback.print_exc()
finally:(P/'logs/proseg.exit').write_text(str(rc)+'\n')
raise SystemExit(rc)
