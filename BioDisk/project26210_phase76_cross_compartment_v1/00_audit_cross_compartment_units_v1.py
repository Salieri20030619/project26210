#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT = Path('/mnt/BioDisk/wangmingyuan')
OUT = ROOT / 'project26210_phase76_cross_compartment_v1' / 'preflight_v1'
OUT.mkdir(parents=True, exist_ok=True)

sources = {
    'GSE207493': ROOT / 'project26210_external_joint_v1_phase64_large_immune_validation_v1/results_v1/GSE207493/gse207493_cell_metadata_scores_v1.tsv.gz',
    'GSE264586': ROOT / 'project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2/GSE264586_phase67_metadata_v1.tsv.gz',
    'GSE294109': ROOT / 'project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2/GSE294109_phase67_metadata_v1.tsv.gz',
}

rows = []
annotation_rows = []
for dataset, path in sources.items():
    d = pd.read_csv(path, sep='\t')
    if dataset == 'GSE207493':
        d['patient_key'] = dataset + '::' + d['sample'].astype(str)
        d['timepoint'] = 'Single'
        d['unit_id'] = d['patient_key'] + '::Single'
        d['is_tumor'] = d['cell_type'].eq('Tumor')
        d['is_myeloid'] = d['cell_type'].isin(['Myeloid', 'Dendritic'])
        annotation_rows.extend({'dataset': dataset, 'annotation_level': 'cell_type', 'label': str(k), 'n_cells': int(v)} for k, v in d['cell_type'].value_counts().items())
    else:
        d['unit_id'] = d['patient_key'].astype(str) + '::' + d['timepoint'].astype(str)
        d['is_tumor'] = d['phase7_fine_annotation'].eq('Epithelial_tumor_like')
        d['is_myeloid'] = d['phase7_broad_annotation'].eq('Myeloid')
        annotation_rows.extend({'dataset': dataset, 'annotation_level': 'phase7_fine_annotation', 'label': str(k), 'n_cells': int(v)} for k, v in d.loc[d['phase7_broad_annotation'].eq('Epithelial'), 'phase7_fine_annotation'].value_counts().items())
    for unit_id, g in d.groupby('unit_id', sort=True):
        rows.append({
            'dataset': dataset,
            'unit_id': unit_id,
            'patient_key': str(g['patient_key'].iloc[0]),
            'timepoint': str(g['timepoint'].iloc[0]),
            'n_all_cells': len(g),
            'n_tumor_cells': int(g['is_tumor'].sum()),
            'n_myeloid_cells': int(g['is_myeloid'].sum()),
        })

units = pd.DataFrame(rows)
for threshold in [20, 30, 50, 100]:
    units[f'eligible_ge{threshold}'] = (units['n_tumor_cells'] >= threshold) & (units['n_myeloid_cells'] >= threshold)
units.to_csv(OUT / 'cross_compartment_unit_counts_v1.tsv', sep='\t', index=False)
pd.DataFrame(annotation_rows).to_csv(OUT / 'annotation_inventory_v1.tsv', sep='\t', index=False)

summary = []
for dataset, g in units.groupby('dataset', sort=True):
    summary.append({
        'dataset': dataset,
        'n_units': len(g),
        'total_tumor_cells': int(g['n_tumor_cells'].sum()),
        'total_myeloid_cells': int(g['n_myeloid_cells'].sum()),
        **{f'eligible_units_ge{t}': int(g[f'eligible_ge{t}'].sum()) for t in [20, 30, 50, 100]},
    })
summary = pd.DataFrame(summary)
summary.to_csv(OUT / 'cross_compartment_preflight_summary_v1.tsv', sep='\t', index=False)
print(summary.to_csv(sep='\t', index=False))
