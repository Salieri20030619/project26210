[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$phaseRoot = $PSScriptRoot
$sourceDir = Join-Path $phaseRoot 'scope_freeze_local_v2'
$sourceScope = Join-Path $sourceDir 'full_cohort_library_scope_v2.tsv'
$recoveryPass = Join-Path $phaseRoot 'GSE242299_H5AD_RECOVERY_LOCAL_VALIDATION_PASS_v2.md'
$finalDir = Join-Path $phaseRoot 'scope_freeze_local_v3'
$partDir = "$finalDir.part"

if (-not (Test-Path -LiteralPath $sourceScope -PathType Leaf)) { throw "Missing v2 scope: $sourceScope" }
if (-not (Test-Path -LiteralPath $recoveryPass -PathType Leaf)) { throw "Missing independent H5AD recovery PASS: $recoveryPass" }
if (Test-Path -LiteralPath $finalDir) { throw "Refusing existing final directory: $finalDir" }
if (Test-Path -LiteralPath $partDir) { throw "Refusing existing part directory: $partDir" }
New-Item -ItemType Directory -Path $partDir | Out-Null

$rows = @(Import-Csv -LiteralPath $sourceScope -Delimiter "`t")
foreach ($row in $rows) {
    $row.scope_version = 'full_cohort_library_scope_v3'
    if ($row.study -eq 'GSE242299') {
        $row.input_source = 'NCBI_GEO_official_series_H5AD_validated_raw_X'
        $row.input_path_or_archive = 'public_geo_alternatives_v1/GSE242299_all_cells_50236_33538.h5ad.gz'
        $row.archive_members = 'GSE242299_all_cells_50236_33538.h5ad.gz::/X'
        $row.inclusion_reason = 'complete_31_of_31_recovered_from_validated_series_H5AD'
    }
}

$scopeOut = Join-Path $partDir 'full_cohort_library_scope_v3.tsv'
$rows | Export-Csv -LiteralPath $scopeOut -Delimiter "`t" -NoTypeInformation -Encoding UTF8

$expectedIncluded = @{
    GSE156632 = 12
    GSE159115 = 11
    GSE178481 = 26
    GSE207493 = 19
    GSE242299 = 31
}
$summary = foreach ($study in @('GSE156632','GSE159115','GSE178481','GSE207493','GSE242299')) {
    $studyRows = @($rows | Where-Object study -eq $study)
    $included = @($studyRows | Where-Object inclusion_status -eq 'include')
    [pscustomobject]@{
        study = $study
        public_libraries = $studyRows.Count
        proposed_included_libraries = $included.Count
        proposed_excluded_libraries = @($studyRows | Where-Object inclusion_status -eq 'exclude').Count
        proposed_patients = @($included.source_patient_id | Sort-Object -Unique).Count
        expected_included_libraries = $expectedIncluded[$study]
        status = if ($included.Count -eq $expectedIncluded[$study]) { 'PASS' } else { 'FAIL' }
    }
}
$summaryOut = Join-Path $partDir 'full_cohort_scope_summary_v3.tsv'
$summary | Export-Csv -LiteralPath $summaryOut -Delimiter "`t" -NoTypeInformation -Encoding UTF8

$checks = [System.Collections.Generic.List[object]]::new()
function Add-Check {
    param([string]$Check, [bool]$Passed, [string]$Observed, [string]$Expected)
    $checks.Add([pscustomobject]@{check=$Check;status=if($Passed){'PASS'}else{'FAIL'};observed=$Observed;expected=$Expected})
}
$includedRows = @($rows | Where-Object inclusion_status -eq 'include')
$excludedRows = @($rows | Where-Object inclusion_status -eq 'exclude')
$g242 = @($rows | Where-Object study -eq 'GSE242299')
$other = @($rows | Where-Object study -ne 'GSE242299')
Add-Check 'total_scope_rows' ($rows.Count -eq 102) ([string]$rows.Count) '102'
Add-Check 'included_and_excluded_rows' (($includedRows.Count -eq 99) -and ($excludedRows.Count -eq 3)) "include=$($includedRows.Count);exclude=$($excludedRows.Count)" 'include=99;exclude=3'
Add-Check 'study_summary' (-not ($summary.status -contains 'FAIL')) (($summary | ForEach-Object { "$($_.study)=$($_.proposed_included_libraries)" }) -join ';') '12;11;26;19;31 and all PASS'
Add-Check 'scope_version' (@($rows | Where-Object scope_version -ne 'full_cohort_library_scope_v3').Count -eq 0) 'evaluated all rows' 'all v3'
Add-Check 'gse242299_consistent_h5ad_source' (
    $g242.Count -eq 31 -and
    @($g242 | Where-Object input_source -ne 'NCBI_GEO_official_series_H5AD_validated_raw_X').Count -eq 0 -and
    @($g242 | Where-Object input_path_or_archive -ne 'public_geo_alternatives_v1/GSE242299_all_cells_50236_33538.h5ad.gz').Count -eq 0 -and
    @($g242 | Where-Object archive_members -ne 'GSE242299_all_cells_50236_33538.h5ad.gz::/X').Count -eq 0
) '31 rows evaluated' '31/31 use one validated official H5AD /X source'

# Compare all non-version fields for non-GSE242299 rows to the v2 source.
$v2Other = @((Import-Csv -LiteralPath $sourceScope -Delimiter "`t") | Where-Object study -ne 'GSE242299')
$compareFields = @('study','geo_accession','source_patient_id','library_id','tissue','public_member_count','input_source','input_path_or_archive','archive_members','inclusion_status','inclusion_reason')
$conservationFailures = 0
for ($i = 0; $i -lt $v2Other.Count; $i++) {
    foreach ($field in $compareFields) {
        if ([string]$v2Other[$i].$field -cne [string]$other[$i].$field) { $conservationFailures++ }
    }
}
Add-Check 'other_study_exact_field_conservation' ($conservationFailures -eq 0) ([string]$conservationFailures) '0 differing compared fields'

$checksOut = Join-Path $partDir 'full_cohort_scope_validation_v3.tsv'
$checks | Export-Csv -LiteralPath $checksOut -Delimiter "`t" -NoTypeInformation -Encoding UTF8
$failed = @($checks | Where-Object status -ne 'PASS')
if ($failed.Count -gt 0) { throw "Scope v3 validation failed: $($failed.check -join ', ')" }

$scopeHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $scopeOut).Hash.ToLowerInvariant()
$summaryHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $summaryOut).Hash.ToLowerInvariant()
$marker = @"
status=PASS
scope_version=full_cohort_library_scope_v3
rows=102
included_libraries=99
excluded_libraries=3
study_specific_patients=53
gse242299_libraries=31
gse242299_count_source=NCBI_GEO_official_series_H5AD_validated_raw_X::/X
scope_sha256=$scopeHash
summary_sha256=$summaryHash
"@
Set-Content -LiteralPath (Join-Path $partDir 'FULL_COHORT_SCOPE_PASS_v3') -Value $marker -Encoding UTF8
Move-Item -LiteralPath $partDir -Destination $finalDir
Write-Output "PASS scope=$scopeHash summary=$summaryHash final=$finalDir"
