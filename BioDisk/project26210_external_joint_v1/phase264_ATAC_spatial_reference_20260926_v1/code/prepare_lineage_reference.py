from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
P=Path(__file__).resolve().parents[1]
ref=Path('/mnt/BioDisk/wangmingyuan/project26210_part4_tymp_reset_phase17_20260812_v1/reference_tymp_continuous_v1/reference_tymp_continuous_v1_signatures_study_balanced_mean_counts.tsv.gz')
groups={'Malignant':['Malignant'],'Normal_epithelial':['Epithelial-Normal'],'Myeloid':['C1QC-state Macrophage','SPP1-state Macrophage','Other Macrophage','Monocyte','Dendritic cell'],'Lymphoid':['T cell','B cell','NK cell'],'Endothelial':['Endothelial'],'SMC_Pericyte':['SMC/Pericyte'],'Mast':['Mast cell']}
canonical={'Malignant':['CA9','EPCAM','KRT8','KRT18','KRT19','PAX8','NDUFA4L2'],
'Normal_epithelial':['LRP2','CUBN','SLC34A1','ALDOB','UMOD','AQP2'],
'Myeloid':['LST1','FCER1G','TYROBP','CSF1R','CTSS','C1QA','CD68'],
'Lymphoid':['CD3D','CD3E','TRAC','IL7R','LCK','CD79A','MS4A1','NKG7','GNLY','KLRD1'],
'Endothelial':['PECAM1','VWF','KDR','EMCN','ESAM','RAMP2'],
'SMC_Pericyte':['RGS5','ACTA2','PDGFRB','MCAM','CSPG4','NOTCH3'],
'Mast':['KIT','TPSAB1','TPSB2','CPA3','MS4A2','GATA2']}
regions=pd.read_csv(P/'reference/protein_coding_gene_activity_regions.tsv',sep='\t')
a=pd.read_csv(ref,sep='\t',index_col=0);a=a.loc[a.index.isin(regions.gene)]
norm=a/a.sum(axis=0)*1e4
means=pd.DataFrame({k:norm[v].mean(axis=1) for k,v in groups.items()})
all_holdout=set(sum(canonical.values(),[]));features=[]
for k in groups:
 others=means.drop(columns=k).max(axis=1)
 ratio=(means[k]+.1)/(others+.1)
 select=ratio[(ratio>=2)&(means[k]>=.5)&(~ratio.index.isin(all_holdout))].sort_values(ascending=False).head(40)
 features.extend(dict(lineage=k,gene=g,reference_fold_specificity=float(v)) for g,v in select.items())
f=pd.DataFrame(features);assert min(f.groupby('lineage').size())>=15
means.to_csv(P/'reference/broad_lineage_RNA_centroids.tsv',sep='\t')
f.to_csv(P/'reference/broad_lineage_training_features.tsv',sep='\t',index=False)
plan=dict(source=str(ref),source_sha256=hashlib.sha256(ref.read_bytes()).hexdigest(),groups=groups,canonical_holdout=canonical,
selection='Top 40 protein-coding centroid markers with >=2-fold over largest other lineage and >=0.5 mean CP10K; canonical validation genes excluded.',
classification='Centered log-CP10K cosine to broad RNA centroids; cell call requires >=0.05 margin, same heldout canonical top lineage and >=2 nonzero canonical genes.',
qc='Use called ATAC barcodes, >=1000 peak fragments and >=500 detected gene-activity features; classify otherwise unassigned.',
pseudobulk='At least 100 confidently assigned cells and at least 20000 chromosome-1 unique fragments per lineage before footprinting.',
projection='Sum cell2location q05 abundances for matching broad lineages, normalize by total mapped abundance; calculate coverage; only project spots with >=80% represented abundance. Report omitted lineages and do not fill them as inactive.',
statistics='One ATAC donor: no biological inferential P-values; any HINT internal test values are model/region-level output only, not patient-level replication.')
(P/'metadata/LINEAGE_TYPING_PLAN.json').write_text(json.dumps(plan,indent=2))
print('LINEAGE_FEATURES_FROZEN',f.groupby('lineage').size().to_dict(),flush=True)

