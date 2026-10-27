[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$phaseRoot = $PSScriptRoot
$syncRoot = Join-Path $phaseRoot 'server_results_sync_v4'
$checksPath = Join-Path $phaseRoot 'gse242299_h5ad_recovery_validation_checks_v1.tsv'
$statusPath = Join-Path $phaseRoot 'GSE242299_H5AD_RECOVERY_LOCAL_VALIDATION_PASS_v1.md'
$checks = [System.Collections.Generic.List[object]]::new()

function Add-Check {
    param([string]$Check, [bool]$Passed, [string]$Observed, [string]$Expected)
    $checks.Add([pscustomobject]@{
        check = $Check
        status = if ($Passed) { 'PASS' } else { 'FAIL' }
        observed = $Observed
        expected = $Expected
    })
}

$required = @(
    'GSE242299_H5AD_RECOVERY_PASS_v3',
    'H5AD_SOURCE_MANIFEST_v3.tsv',
    'h5ad_dataset_inventory_v3.tsv',
    'h5ad_library_cell_counts_v3.tsv',
    'h5ad_matrix_candidate_audit_v3.tsv',
    'h5ad_obs_candidate_columns_v3.tsv',
    'h5ad_recovery_metadata_v3.json',
    '11_finalize_gse242299_h5ad_recovery_v3.log'
)
foreach ($name in $required) {
    Add-Check "required_file:$name" (Test-Path -LiteralPath (Join-Path $syncRoot $name) -PathType Leaf) 'present_or_absent' 'present'
}

$marker = Get-Content -LiteralPath (Join-Path $syncRoot 'GSE242299_H5AD_RECOVERY_PASS_v3') -Raw
$metadata = Get-Content -LiteralPath (Join-Path $syncRoot 'h5ad_recovery_metadata_v3.json') -Raw | ConvertFrom-Json
$matrix = Import-Csv -LiteralPath (Join-Path $syncRoot 'h5ad_matrix_candidate_audit_v3.tsv') -Delimiter "`t"
$libraries = Import-Csv -LiteralPath (Join-Path $syncRoot 'h5ad_library_cell_counts_v3.tsv') -Delimiter "`t"
$source = Import-Csv -LiteralPath (Join-Path $syncRoot 'H5AD_SOURCE_MANIFEST_v3.tsv') -Delimiter "`t"
$log = Get-Content -LiteralPath (Join-Path $syncRoot '11_finalize_gse242299_h5ad_recovery_v3.log') -Raw

Add-Check 'marker_status' ($marker -match '(?m)^status=PASS$') (($marker -split "`n")[0]) 'status=PASS'
Add-Check 'metadata_status' ($metadata.status -eq 'PASS') ([string]$metadata.status) 'PASS'
Add-Check 'dimensions' (($metadata.n_obs -eq 50236) -and ($metadata.n_vars -eq 33538)) "$($metadata.n_obs)x$($metadata.n_vars)" '50236x33538'
Add-Check 'library_column' (($metadata.selected_library_column -eq 'library_name') -and (@($metadata.exact_library_columns).Count -eq 1)) ([string]$metadata.selected_library_column) 'library_name as the sole exact 31-library field'
Add-Check 'selected_count_matrix' (($metadata.selected_count_matrix_path -eq '/X') -and ($metadata.count_matrix_rule -eq 'accept_unique_integer_X_when_raw_and_layers_absent')) "$($metadata.selected_count_matrix_path);$($metadata.count_matrix_rule)" '/X; accept unique integer X when raw/layers absent'
Add-Check 'author_code_commit' ($metadata.author_code_commit -eq '1c9daebb33449f32a06978227aa5fcf0db726c3a') ([string]$metadata.author_code_commit) '1c9daebb33449f32a06978227aa5fcf0db726c3a'

Add-Check 'matrix_row_count' ($matrix.Count -eq 1) ([string]$matrix.Count) '1'
if ($matrix.Count -eq 1) {
    $m = $matrix[0]
    $matrixPass = (
        $m.matrix_path -eq '/X' -and $m.storage -eq 'csr_matrix' -and
        [int]$m.shape_rows -eq 50236 -and [int]$m.shape_columns -eq 33538 -and
        [long]$m.stored_values -eq 26714153 -and $m.finite_all -eq 'true' -and
        $m.nonnegative_all -eq 'true' -and $m.integer_all -eq 'true' -and
        $m.structural_checks -eq 'true' -and $m.raw_count_candidate -eq 'true' -and
        [double]$m.minimum -eq 1 -and [double]$m.maximum -eq 1251
    )
    Add-Check 'matrix_safety_and_structure' $matrixPass "$($m.matrix_path);$($m.shape_rows)x$($m.shape_columns);nnz=$($m.stored_values);finite=$($m.finite_all);nonnegative=$($m.nonnegative_all);integer=$($m.integer_all);structural=$($m.structural_checks);range=$($m.minimum)-$($m.maximum)" '/X CSR; 50236x33538; nnz=26714153; all safety checks true; range 1-1251'
}

$uniqueLibraries = @($libraries.library_id | Sort-Object -Unique)
$cellSum = [long](($libraries | Measure-Object -Property cells -Sum).Sum)
$t2 = $libraries | Where-Object library_id -eq 'T2_1'
$n3 = $libraries | Where-Object library_id -eq 'N3_1'
Add-Check 'library_count_and_uniqueness' (($libraries.Count -eq 31) -and ($uniqueLibraries.Count -eq 31)) "$($libraries.Count) rows;$($uniqueLibraries.Count) unique" '31 rows;31 unique'
Add-Check 'library_cell_conservation' (($cellSum -eq 50236) -and -not ($libraries.status -contains 'FAIL')) "sum=$cellSum;failures=$(@($libraries | Where-Object status -ne 'PASS').Count)" 'sum=50236;failures=0'
Add-Check 'affected_library_recovery' (([int]$t2.cells -eq 4906) -and ([int]$n3.cells -eq 240)) "T2_1=$($t2.cells);N3_1=$($n3.cells)" 'T2_1=4906;N3_1=240'

Add-Check 'source_manifest_row_count' ($source.Count -eq 1) ([string]$source.Count) '1'
if ($source.Count -eq 1) {
    $s = $source[0]
    $sourcePass = (
        $s.status -eq 'PASS' -and [long]$s.decompressed_size_bytes -eq 7345987312 -and
        $s.decompressed_sha256 -match '^[0-9a-f]{64}$' -and
        [long]$s.compressed_size_bytes -eq 5350598000 -and
        $s.compressed_sha256 -eq '8539723bf40f3d8647c9773c386db56256b305f69e0cdf80848cec5bae428f94'
    )
    Add-Check 'source_manifest_integrity' $sourcePass "decompressed=$($s.decompressed_size_bytes):$($s.decompressed_sha256);compressed=$($s.compressed_size_bytes):$($s.compressed_sha256);status=$($s.status)" 'validated sizes, 64-char decompressed SHA256, exact official compressed SHA256, PASS'
}
Add-Check 'atomic_finalize_log' (($log -match '(?m)^finished_at=') -and ($log -notmatch 'Traceback|RuntimeError|BLOCKED')) 'finish marker and error scan evaluated' 'finished_at present; no Traceback/RuntimeError/BLOCKED'

$checks | Export-Csv -LiteralPath $checksPath -Delimiter "`t" -NoTypeInformation -Encoding UTF8
$failed = @($checks | Where-Object status -ne 'PASS')
if ($failed.Count -gt 0) {
    throw "GSE242299 H5AD independent validation failed: $($failed.check -join ', ')"
}

$generatedAt = (Get-Date).ToString('yyyy-MM-ddTHH:mm:ssK')
$status = @"
# GSE242299 official H5AD recovery — independent local validation

- Status: **PASS**
- Generated: `$generatedAt`
- Checks: $($checks.Count)/$($checks.Count) PASS
- Matrix: `/X`, CSR, 50,236 cells × 33,538 genes, 26,714,153 stored nonzero values
- Count safety: finite, nonnegative, exactly integer-valued, structurally valid
- Library field: `library_name`, exact 31-library match; cell total 50,236
- Recovered defective-MTX libraries: T2_1 = 4,906 cells; N3_1 = 240 cells
- Official compressed source SHA256: `8539723bf40f3d8647c9773c386db56256b305f69e0cdf80848cec5bae428f94`
- Decompressed H5AD SHA256: `$($source[0].decompressed_sha256)`
- Author-code snapshot commit: `1c9daebb33449f32a06978227aa5fcf0db726c3a`
- Large H5AD remains on `sgs10`; only small audit evidence was synchronized locally.
"@
Set-Content -LiteralPath $statusPath -Value $status -Encoding UTF8
Write-Output "PASS checks=$($checks.Count) status=$statusPath"
