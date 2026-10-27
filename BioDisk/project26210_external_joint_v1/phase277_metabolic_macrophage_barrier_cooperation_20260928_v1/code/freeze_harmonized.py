from pathlib import Path
import pandas as pd,json,hashlib,datetime,h5py
P=Path('phase277_metabolic_macrophage_barrier_cooperation_20260928_v1');H=P/'harmonized_platform_v1'
assert not H.exists()
for d in ['code','tables','inputs','qc','logs']:(H/d).mkdir(parents=True,exist_ok=True)
programs=json.loads((P/'PROGRAMS_FROZEN.json').read_text())['programs'];q=pd.read_csv(P/'tables/PROGRAM_ELIGIBILITY.tsv',sep='\t');common=set(programs['C4_19'])
for r in q[q.program=='C4_19'].itertuples():common &= set(r.present_genes.split(';'))
genes=[g for g in programs['C4_19'] if g in common];removed=sorted(set(programs['C4_19'])-set(genes));assert len(genes)==17 and removed==['CES1','FABP3']
files=[Path('/mnt/BioDisk/ludingyu/project26210/external/GSE210041/unpacked')/f'{gsm}_{sample}_filtered_feature_bc_matrix.h5' for gsm,sample in [('GSM6415705','GLMF1'),('GSM6415706','GLMF2')]]
files+=sorted(Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase54_gse175540_spatial_v1/input_extract_v1').glob('*ffpe*filtered_feature_bc_matrix.h5'));assert len(files)==14;raw=[]
for f in files:
 with h5py.File(f) as h:
  names=[v.decode() if isinstance(v,bytes) else str(v) for v in h['matrix/features/name'][:]];assert all(g not in names for g in removed)
  raw.append(dict(path=str(f),features=len(names),CES1_present=False,FABP3_present=False,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
pd.DataFrame(raw).to_csv(H/'tables/RAW_PLATFORM_MISSING_GENES.tsv',sep='\t',index=False)
plan=dict(at=datetime.datetime.now().astimezone().isoformat(),reason='CES1 and FABP3 absent from original filtered feature lists of both GLMF and all12 GSE175540 FFPE sections; not zero biological expression or QC deletion',previous_results_known=True,original19gene_results_preserved=True,genes_selected_only_by_assay_intersection=True,C4_platform17=genes,C4_platform15_noAPOEshared=[g for g in genes if g not in ['MARCO','S100A9']],removed=removed,models='same Phase277 geometry,qualification,core/state adjustment,outcomes,interactions,sign/Wilcoxon/bootstrap; no relaxation of50segments/quadrant/condition thresholds',families='technical extension separate C4 interaction8,disjoint sensitivity8,metabolic context6 percohort; all66 endpoints retained; not replacement original primary or independent replication',identity='descriptive platform17 paired C4-other scores and reference lineage specificity; not new heldout validation',no_new_pathways_or_regions=True,figures=0)
(H/'ADDENDUM_BEFORE_HARMONIZED_MODELS.json').write_text(json.dumps(plan,indent=2));print('ADDENDUM_FROZEN',len(raw),genes)
