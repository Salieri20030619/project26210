$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$phaseDir = $PSScriptRoot
$projectRoot = Split-Path -Parent $phaseDir
$syncDir = Join-Path $phaseDir 'server_results_sync_v2'
$outDir = Join-Path $phaseDir 'scope_freeze_local_v2'
$scopeFile = Join-Path $outDir 'full_cohort_library_scope_v2.tsv'
$summaryFile = Join-Path $outDir 'full_cohort_scope_summary_v2.tsv'
$statusFile = Join-Path $outDir 'FULL_COHORT_SCOPE_PASS_v2'

if (Test-Path -LiteralPath $outDir) {
    throw "Refusing existing output directory: $outDir"
}
New-Item -ItemType Directory -Path $outDir -ErrorAction Stop | Out-Null

function New-ScopeRow {
    param(
        [string]$Study,
        [string]$GeoAccession,
        [string]$Patient,
        [string]$LibraryId,
        [string]$Tissue,
        [int]$PublicMemberCount,
        [string]$InputSource,
        [string]$InputPathOrArchive,
        [string]$ArchiveMembers,
        [string]$InclusionStatus,
        [string]$InclusionReason
    )
    [pscustomobject]@{
        study = $Study
        geo_accession = $GeoAccession
        source_patient_id = $Patient
        library_id = $LibraryId
        tissue = $Tissue
        public_member_count = $PublicMemberCount
        input_source = $InputSource
        input_path_or_archive = $InputPathOrArchive
        archive_members = $ArchiveMembers
        inclusion_status = $InclusionStatus
        inclusion_reason = $InclusionReason
        scope_version = 'full_cohort_library_scope_v2'
    }
}

$rows = [System.Collections.Generic.List[object]]::new()

# GSE156632: all 12 public count CSV files.
$g156Members = Get-Content -LiteralPath (Join-Path $syncDir 'GSE156632_tar_members_v1.txt')
foreach ($member in $g156Members) {
    if ($member -notmatch '^(GSM\d+)_(RCC\d+)([tn])\.csv\.gz$') {
        throw "Unexpected GSE156632 archive member: $member"
    }
    $gsm = $Matches[1]
    $patient = $Matches[2]
    $suffix = $Matches[3]
    $tissue = if ($suffix -eq 't') { 'Tumor' } else { 'Normal' }
    $rows.Add((New-ScopeRow -Study 'GSE156632' -GeoAccession $gsm -Patient $patient `
        -LibraryId "$patient$suffix" -Tissue $tissue -PublicMemberCount 1 `
        -InputSource 'NCBI_GEO_processed_raw_counts' `
        -InputPathOrArchive 'public_geo_archives_v1/GSE156632_RAW.tar' -ArchiveMembers $member `
        -InclusionStatus 'include' -InclusionReason 'complete_12_of_12_ccRCC_scope'))
}

# GSE159115: use the curated patient mapping for 11 ccRCC-compatible records and explicitly exclude 3 chRCC-related records.
$provenance = Import-Csv -Delimiter "`t" -LiteralPath (Join-Path $syncDir 'geo_sample_provenance_v1.tsv')
$g159Provenance = @{}
foreach ($row in ($provenance | Where-Object { $_.study -eq 'GSE159115' })) {
    $g159Provenance[$row.geo_accession] = $row
}
$g159Excluded = @{
    'GSM4819730' = [pscustomobject]@{ patient = 'SS_2014'; tissue = 'Normal'; reason = 'exclude_chRCC_patient_benign_cortex' }
    'GSM4819731' = [pscustomobject]@{ patient = 'SS_2014'; tissue = 'Normal'; reason = 'exclude_chRCC_patient_benign_medulla' }
    'GSM4819732' = [pscustomobject]@{ patient = 'SS_2016'; tissue = 'Tumor'; reason = 'exclude_chRCC_tumor' }
}
$g159Members = Get-Content -LiteralPath (Join-Path $syncDir 'GSE159115_tar_members_v1.txt')
foreach ($member in $g159Members) {
    if ($member -notmatch '^(GSM\d+)_(SI_\d+)_filtered_gene_bc_matrices_h5\.h5$') {
        throw "Unexpected GSE159115 archive member: $member"
    }
    $gsm = $Matches[1]
    $libraryId = $Matches[2]
    if ($g159Provenance.ContainsKey($gsm)) {
        $p = $g159Provenance[$gsm]
        $patient = $p.source_patient_id
        $tissue = $p.source_tissue
        $status = 'include'
        $reason = 'ccRCC_compatible_scope'
    } elseif ($g159Excluded.ContainsKey($gsm)) {
        $e = $g159Excluded[$gsm]
        $patient = $e.patient
        $tissue = $e.tissue
        $status = 'exclude'
        $reason = $e.reason
    } else {
        throw "GSE159115 record lacks include/exclude mapping: $gsm"
    }
    $rows.Add((New-ScopeRow -Study 'GSE159115' -GeoAccession $gsm -Patient $patient `
        -LibraryId $libraryId -Tissue $tissue -PublicMemberCount 1 `
        -InputSource 'NCBI_GEO_processed_10x_H5' `
        -InputPathOrArchive 'public_geo_archives_v1/GSE159115_RAW.tar' -ArchiveMembers $member `
        -InclusionStatus $status -InclusionReason $reason))
}

# GSE178481 and GSE207493: complete local RDS sets already matched to the frozen transfer manifest.
$transfer = Import-Csv -Delimiter "`t" -LiteralPath (Join-Path $projectRoot 'local_audit_v1\manifest_reconciliation.tsv')
foreach ($row in ($transfer | Where-Object { $_.dataset -eq 'GSE178481' } | Sort-Object object_name)) {
    if ($row.object_name -notmatch '^(GSM\d+)_RCC-((?:PR|BM)\d+)-(.+)$') {
        throw "Unexpected GSE178481 object name: $($row.object_name)"
    }
    $gsm = $Matches[1]
    $patient = $Matches[2]
    $sampleType = $Matches[3]
    $tissue = if ($sampleType -match '^Normal$') { 'Normal' } else { 'Tumor' }
    $rows.Add((New-ScopeRow -Study 'GSE178481' -GeoAccession $gsm -Patient $patient `
        -LibraryId $row.object_name -Tissue $tissue -PublicMemberCount 1 `
        -InputSource 'validated_local_raw_count_RDS' -InputPathOrArchive $row.manifest_source_file `
        -ArchiveMembers '' -InclusionStatus 'include' -InclusionReason 'complete_26_of_26_ccRCC_scope'))
}

$g207Gsm = @{
    'RCC81'='GSM6290518'; 'RCC84'='GSM6290520'; 'RCC86'='GSM6290522'; 'RCC87'='GSM6290524'
    'RCC94'='GSM6290526'; 'RCC96'='GSM6290528'; 'RCC99'='GSM6290530'; 'RCC100'='GSM6290532'
    'RCC101'='GSM6290534'; 'RCC103'='GSM6290536'; 'RCC104'='GSM6290538'; 'RCC106'='GSM6290540'
    'RCC112'='GSM6290542'; 'RCC113'='GSM6290544'; 'RCC114'='GSM6290546'; 'RCC115'='GSM6290548'
    'RCC116'='GSM6290550'; 'RCC119'='GSM6290552'; 'RCC120'='GSM6290554'
}
foreach ($row in ($transfer | Where-Object { $_.dataset -eq 'GSE207493' } | Sort-Object object_name)) {
    if ($row.object_name -notmatch '^GSE207493__(RCC\d+)$') {
        throw "Unexpected GSE207493 object name: $($row.object_name)"
    }
    $patient = $Matches[1]
    if (-not $g207Gsm.ContainsKey($patient)) { throw "Missing GSE207493 GSM mapping: $patient" }
    $rows.Add((New-ScopeRow -Study 'GSE207493' -GeoAccession $g207Gsm[$patient] -Patient $patient `
        -LibraryId $row.object_name -Tissue 'Tumor' -PublicMemberCount 1 `
        -InputSource 'validated_local_raw_count_RDS' -InputPathOrArchive $row.manifest_source_file `
        -ArchiveMembers '' -InclusionStatus 'include' -InclusionReason 'complete_19_of_19_scRNA_scope'))
}

# GSE242299: all 31 public emulsions, each represented by barcode/features/matrix files.
$g242Members = Get-Content -LiteralPath (Join-Path $syncDir 'GSE242299_tar_members_v1.txt')
$g242Groups = $g242Members | Group-Object { if ($_ -match '^(GSM\d+)_([TN]\d+_\d+)_') { "$($Matches[1])|$($Matches[2])" } else { throw "Unexpected GSE242299 member: $_" } }
foreach ($group in $g242Groups) {
    $parts = $group.Name -split '\|'
    $gsm = $parts[0]
    $libraryId = $parts[1]
    if ($group.Count -ne 3) { throw "GSE242299 library does not have 3 matrix members: $($group.Name)" }
    $expectedSuffixes = @('_barcodes.tsv.gz','_features.tsv.gz','_matrix.mtx.gz')
    foreach ($suffix in $expectedSuffixes) {
        if (@($group.Group | Where-Object { $_.EndsWith($suffix) }).Count -ne 1) {
            throw "Missing or duplicate GSE242299 member $suffix for $($group.Name)"
        }
    }
    if ($libraryId -notmatch '^([TN])(\d+)_\d+$') { throw "Unexpected GSE242299 library id: $libraryId" }
    $tissue = if ($Matches[1] -eq 'T') { 'Tumor' } else { 'Normal' }
    $patient = "P$($Matches[2])"
    $rows.Add((New-ScopeRow -Study 'GSE242299' -GeoAccession $gsm -Patient $patient `
        -LibraryId $libraryId -Tissue $tissue -PublicMemberCount 3 `
        -InputSource 'NCBI_GEO_processed_MTX' `
        -InputPathOrArchive 'public_geo_archives_v1/GSE242299_RAW.tar' `
        -ArchiveMembers (($group.Group | Sort-Object) -join '|') `
        -InclusionStatus 'include' -InclusionReason 'complete_31_of_31_emulsion_scope'))
}

$scope = @($rows | Sort-Object study, geo_accession)
if ($scope.Count -ne 102) { throw "Public library row count mismatch: $($scope.Count) != 102" }
if (@($scope | Where-Object { $_.inclusion_status -eq 'include' }).Count -ne 99) { throw 'Included row count mismatch' }
if (@($scope | Where-Object { $_.inclusion_status -eq 'exclude' }).Count -ne 3) { throw 'Excluded row count mismatch' }
if (@($scope | Group-Object { "$($_.study)|$($_.geo_accession)" } | Where-Object Count -ne 1).Count -ne 0) { throw 'Duplicate study/GSM key' }

$expectedIncluded = @{ GSE156632=12; GSE159115=11; GSE178481=26; GSE207493=19; GSE242299=31 }
$summary = foreach ($study in ($expectedIncluded.Keys | Sort-Object)) {
    $z = @($scope | Where-Object { $_.study -eq $study })
    $included = @($z | Where-Object { $_.inclusion_status -eq 'include' })
    $excluded = @($z | Where-Object { $_.inclusion_status -eq 'exclude' })
    $patients = @($included.source_patient_id | Sort-Object -Unique).Count
    $status = if ($included.Count -eq $expectedIncluded[$study]) { 'PASS' } else { 'FAIL' }
    [pscustomobject]@{
        study = $study
        public_libraries = $z.Count
        proposed_included_libraries = $included.Count
        proposed_excluded_libraries = $excluded.Count
        proposed_patients = $patients
        expected_included_libraries = $expectedIncluded[$study]
        status = $status
    }
}
if (@($summary | Where-Object status -ne 'PASS').Count -ne 0) { throw 'Per-study scope validation failed' }
if (($summary.proposed_patients | Measure-Object -Sum).Sum -ne 53) { throw 'Study-specific patient count sum mismatch' }

$scope | Export-Csv -Delimiter "`t" -NoTypeInformation -Encoding utf8 -LiteralPath $scopeFile
$summary | Export-Csv -Delimiter "`t" -NoTypeInformation -Encoding utf8 -LiteralPath $summaryFile
@(
    'status=PASS'
    'public_library_records=102'
    'included_libraries=99'
    'excluded_libraries=3'
    'study_specific_patients=53'
    'previous_libraries=75'
    'previous_patients=40'
    'downstream_started=false'
    ('generated_at=' + (Get-Date -Format 'yyyy-MM-ddTHH:mm:ssK'))
) | Set-Content -Encoding utf8 -LiteralPath $statusFile

Get-Content -LiteralPath $statusFile
$summary | Format-Table -AutoSize
