import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r'D:\project26210_external_joint_v1\phase114_unweighted_seven_feature_analysis_local_v1')
DOWNLOADS = ROOT / 'downloads_v1'
INPUTS = ROOT / 'inputs_v1'
QC = ROOT / 'qc_v1'

DATASETS = {
    'RC18': ('matched_Harmonized_RawData_RC18.csv', 'matched_tpm_RC18.csv', 144, 97, 47),
    'RC20': ('matched_Harmonized_RawData_RC20.csv', 'matched_tpm_RC20.csv', 76, 52, 24),
    'CPTAC': ('matched_Harmonized_Met_CPTAC.csv', 'matched_tpm_CPTAC.csv', 50, 50, 0),
    'CPTAC_val': ('matched_Harmonized_Met_CPTAC_val.csv', 'matched_tpm_CPTAC_val.csv', 71, 56, 15),
}

UNITEDMET_FEATURES = {
    '5-MTA': '5-methylthioadenosine (MTA)',
    'Thymidine': 'thymidine',
    "2'-Deoxyuridine": "2'-deoxyuridine",
    'Methionine': 'methionine',
    'SAM': 'S-adenosylmethionine (SAM)',
    'Thymine': 'thymine',
    'Spermidine': 'spermidine',
}

NCI60_FEATURES = {
    '5-MTA': "5'-s-methyl-5'-thioadenosine",
    "2'-Deoxyuridine": "2'-deoxyuridine",
    'Methionine': 'methionine',
    'SAM': "S-(5'-adenosyl)-L-methionine",
    'Spermidine': 'Spermidine',
}

EXPECTED_COVERAGE = {
    'RC18': set(UNITEDMET_FEATURES),
    'RC20': set(UNITEDMET_FEATURES),
    'CPTAC': {'5-MTA', 'Methionine', 'Thymine'},
    'CPTAC_val': {'5-MTA', 'Methionine', 'Thymine'},
}


def digest(path, algorithm):
    h = hashlib.new(algorithm)
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def normalize_cell(value):
    text = str(value).upper().strip()
    if text.endswith('/ATCC'):
        text = text[:-5]
    return ''.join(char for char in text if char.isalnum())


def write_tsv(path, rows):
    if not rows:
        raise RuntimeError(f'no rows for {path}')
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False, quoting=csv.QUOTE_MINIMAL)


def main():
    QC.mkdir(parents=True, exist_ok=True)
    record = json.loads((INPUTS / 'zenodo_record_11286535_v1.json').read_text(encoding='utf-8'))
    expected_files = {
        item['key']: {
            'size': int(item['size']),
            'md5': item['checksum'].split(':', 1)[1],
        }
        for item in record['files']
    }
    inventory = []
    file_pass = True
    for name, expected in sorted(expected_files.items()):
        path = DOWNLOADS / name
        exists = path.is_file()
        size = path.stat().st_size if exists else -1
        md5 = digest(path, 'md5') if exists else ''
        sha256 = digest(path, 'sha256') if exists else ''
        passed = exists and size == expected['size'] and md5 == expected['md5']
        file_pass &= passed
        inventory.append({
            'file': name,
            'absolute_path': str(path),
            'expected_bytes': expected['size'],
            'observed_bytes': size,
            'expected_md5': expected['md5'],
            'observed_md5': md5,
            'sha256': sha256,
            'status': 'PASS' if passed else 'FAIL',
        })
    write_tsv(QC / 'public_input_file_inventory_v1.tsv', inventory)

    mapping = pd.read_csv(INPUTS / 'MasterMapping_ccRCC.csv', dtype=str).fillna('')
    mapping_rows = []
    feature_rows = []
    unitedmet_pass = True

    for dataset, (met_name, rna_name, expected_total, expected_tumor, expected_normal) in DATASETS.items():
        mapd = mapping.loc[mapping['Dataset'].eq(dataset)].copy()
        met_path = DOWNLOADS / met_name
        rna_path = DOWNLOADS / rna_name
        met_header = pd.read_csv(met_path, nrows=0).columns.tolist()
        rna_header = pd.read_csv(rna_path, nrows=0).columns.tolist()
        met_id_col = met_header[0]
        rna_id_col = rna_header[0]

        available_columns = [met_id_col] + [
            original for original in UNITEDMET_FEATURES.values() if original in met_header
        ]
        met = pd.read_csv(met_path, usecols=available_columns, dtype={met_id_col: str})
        rna = pd.read_csv(rna_path, usecols=[rna_id_col, 'TYMP', 'MTAP'], dtype={rna_id_col: str})

        met_ids = met[met_id_col].astype(str)
        rna_ids = rna[rna_id_col].astype(str)
        ids_match = (
            len(mapd) == len(met) == len(rna)
            and set(mapd['MetabID']) == set(met_ids)
            and set(mapd['RNAID']) == set(rna_ids)
            and not mapd['MetabID'].duplicated().any()
            and not mapd['RNAID'].duplicated().any()
            and not met_ids.duplicated().any()
            and not rna_ids.duplicated().any()
        )
        tumor_n = int(mapd['TN'].eq('Tumor').sum())
        normal_n = int(mapd['TN'].eq('Normal').sum())
        count_pass = (
            len(mapd) == expected_total
            and tumor_n == expected_tumor
            and normal_n == expected_normal
        )
        mapping_status = ids_match and count_pass
        unitedmet_pass &= mapping_status
        mapping_rows.append({
            'dataset': dataset,
            'mapping_n': len(mapd),
            'metabolomics_n': len(met),
            'rna_n': len(rna),
            'tumor_n': tumor_n,
            'normal_n': normal_n,
            'metabolomics_id_unique': not met_ids.duplicated().any(),
            'rna_id_unique': not rna_ids.duplicated().any(),
            'mapping_ids_exact': ids_match,
            'expected_counts_exact': count_pass,
            'status': 'PASS' if mapping_status else 'FAIL',
        })

        tumor_met_ids = set(mapd.loc[mapd['TN'].eq('Tumor'), 'MetabID'])
        tumor_rna_ids = set(mapd.loc[mapd['TN'].eq('Tumor'), 'RNAID'])
        met_tumor = met.loc[met_ids.isin(tumor_met_ids)].copy()
        rna_tumor = rna.loc[rna_ids.isin(tumor_rna_ids)].copy()
        rna_finite = all(
            np.isfinite(pd.to_numeric(rna_tumor[gene], errors='coerce')).all()
            for gene in ['TYMP', 'MTAP']
        )
        expected_originals = [
            UNITEDMET_FEATURES[canonical] for canonical in EXPECTED_COVERAGE[dataset]
        ]
        expected_numeric = met_tumor[expected_originals].apply(pd.to_numeric, errors='coerce')
        complete_feature_case_n = int(np.isfinite(expected_numeric).all(axis=1).sum())
        complete_case_gate = complete_feature_case_n >= 10
        mapping_rows[-1]['complete_feature_case_n'] = complete_feature_case_n
        mapping_rows[-1]['complete_case_gate_n_at_least_10'] = complete_case_gate
        unitedmet_pass &= (
            rna_finite
            and len(met_tumor) == expected_tumor
            and len(rna_tumor) == expected_tumor
            and complete_case_gate
        )

        observed_coverage = set()
        for canonical, original in UNITEDMET_FEATURES.items():
            exact_count = sum(column == original for column in met_header)
            present = exact_count == 1
            finite_n = 0
            unique_n = 0
            variance_nonzero = False
            if present:
                values = pd.to_numeric(met_tumor[original], errors='coerce')
                finite = values[np.isfinite(values)]
                finite_n = len(finite)
                unique_n = finite.nunique(dropna=True)
                variance_nonzero = finite_n >= 10 and unique_n > 1 and float(finite.var(ddof=0)) > 0
                observed_coverage.add(canonical)
            expected_present = canonical in EXPECTED_COVERAGE[dataset]
            passed = present == expected_present and (not expected_present or variance_nonzero)
            unitedmet_pass &= passed
            feature_rows.append({
                'dataset': dataset,
                'canonical_feature': canonical,
                'platform_original_name': original,
                'exact_column_count': exact_count,
                'expected_present': expected_present,
                'tumor_finite_n': finite_n,
                'tumor_unique_n': unique_n,
                'tumor_variance_nonzero': variance_nonzero,
                'status': 'PASS' if passed else 'FAIL',
            })
        if observed_coverage != EXPECTED_COVERAGE[dataset]:
            unitedmet_pass = False

    write_tsv(QC / 'unitedmet_sample_mapping_audit_v1.tsv', mapping_rows)
    write_tsv(QC / 'metabolite_identity_coverage_v1.tsv', feature_rows)

    nci_raw = pd.read_csv(
        INPUTS / 'NCI60_WEB_DATA_METABOLON_ALL.TXT',
        usecols=['SAMPLENAME', 'TITLE', 'pname', 'cellname', 'Value'],
    )
    nci_existing = pd.read_csv(INPUTS / 'nci60_all_lineages_preprocessed_v1.tsv', sep='\t')
    selected = nci_raw.loc[nci_raw['TITLE'].isin(NCI60_FEATURES.values())].copy()
    selected['cell_key'] = selected['cellname'].map(normalize_cell)
    selected['Value'] = pd.to_numeric(selected['Value'], errors='coerce')
    selected['canonical_feature'] = selected['TITLE'].map({value: key for key, value in NCI60_FEATURES.items()})
    agg = selected.groupby(['cell_key', 'canonical_feature'], as_index=False).agg(
        value=('Value', 'mean'),
        n_technical_repeats=('SAMPLENAME', 'nunique'),
        tissue=('pname', 'first'),
    )
    wide = agg.pivot(index='cell_key', columns='canonical_feature', values='value')
    existing_keys = set(nci_existing['cell_key'].astype(str))
    nci_rows = []
    nci_pass = True
    for canonical, original in NCI60_FEATURES.items():
        subset = agg.loc[agg['canonical_feature'].eq(canonical)]
        present_keys = set(subset['cell_key'])
        finite_n = int(np.isfinite(subset['value']).sum())
        status = (
            present_keys == existing_keys
            and finite_n == 58
            and subset['value'].nunique(dropna=True) > 1
            and subset['n_technical_repeats'].min() >= 2
        )
        nci_pass &= status
        nci_rows.append({
            'canonical_feature': canonical,
            'platform_original_name': original,
            'raw_exact_title_count': int(nci_raw['TITLE'].eq(original).sum()),
            'independent_cell_lines': len(present_keys),
            'finite_cell_line_means': finite_n,
            'min_technical_repeats': int(subset['n_technical_repeats'].min()),
            'unique_cell_line_means': int(subset['value'].nunique(dropna=True)),
            'status': 'PASS' if status else 'FAIL',
        })
    mta_compare = wide['5-MTA'].reindex(nci_existing['cell_key']).to_numpy() - nci_existing['mta_mean_raw'].to_numpy()
    mta_delta = float(np.nanmax(np.abs(mta_compare)))
    renal_keys = set(nci_existing.loc[nci_existing['is_renal'].astype(str).str.upper().eq('TRUE'), 'cell_key'])
    nci_global_checks = {
        'five_features_exact': bool(set(wide.columns) == set(NCI60_FEATURES)),
        'all_58_complete': bool(wide.reindex(sorted(existing_keys)).notna().all(axis=None)),
        'mta_reproduces_phase102': bool(mta_delta <= 1e-12),
        'renal_n_8': bool(len(renal_keys) == 8),
        'cell_7860_present_once': bool(int((nci_existing['cell_key'] == '7860').sum()) == 1),
    }
    nci_pass &= all(nci_global_checks.values())
    write_tsv(QC / 'nci60_five_feature_input_audit_v1.tsv', nci_rows)

    status = {
        'status': 'PASS' if file_pass and unitedmet_pass and nci_pass else 'FAIL',
        'zenodo_record': 11286535,
        'official_files_expected': len(expected_files),
        'official_file_identity_pass': file_pass,
        'unitedmet_mapping_and_coverage_pass': unitedmet_pass,
        'nci60_five_feature_input_pass': nci_pass,
        'nci60_global_checks': nci_global_checks,
        'nci60_mta_max_abs_delta_phase102': mta_delta,
        'complete_case_counts': {
            row['dataset']: int(row['complete_feature_case_n']) for row in mapping_rows
        },
        'total_tumor_counts': {'RC18': 97, 'RC20': 52, 'CPTAC': 50, 'CPTAC_val': 56},
        'nci60_counts': {'NCI60': 58, 'NCI60_renal': 8},
        'analysis_authorized_if_pass': 'Phase112 steps 3-7 only; conditional steps 8-11 remain gated',
    }
    (QC / 'INPUT_AUDIT_STATUS_v1.json').write_text(
        json.dumps(status, ensure_ascii=False, indent=2) + '\n', encoding='utf-8'
    )
    print(json.dumps(status, ensure_ascii=False, indent=2))
    raise SystemExit(0 if status['status'] == 'PASS' else 1)


if __name__ == '__main__':
    main()
