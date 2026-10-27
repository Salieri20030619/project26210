[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$phaseRoot = $PSScriptRoot
$projectRoot = Split-Path $phaseRoot -Parent
$source = Join-Path $projectRoot 'phase48_part1_isolated_rebuild_pilot_local_v1\scope_freeze_local_v3\full_cohort_library_scope_v3.tsv'
$finalDir = Join-Path $phaseRoot 'scope_v1'
$partDir = "$finalDir.part"
if (Test-Path $finalDir) { throw "Refusing existing final directory: $finalDir" }
if (Test-Path $partDir) { throw "Refusing existing part directory: $partDir" }
New-Item -ItemType Directory -Path $partDir | Out-Null

$rows = @(Import-Csv -LiteralPath $source -Delimiter "`t" | Where-Object inclusion_status -eq 'include')
$old178 = @('GSM5392398','GSM5392399','GSM5392400','GSM5392401','GSM5392402','GSM5392403','GSM5392404','GSM5392405','GSM5392406','GSM5392407','GSM5392408')
$old207 = @('RCC86','RCC87','RCC94','RCC96','RCC99','RCC100','RCC101','RCC103','RCC104','RCC106','RCC112','RCC113')
$old242Missing = @('N3_1','N4_1')

$out = foreach ($row in $rows) {
    $isOld = switch ($row.study) {
        'GSE156632' { $true }
        'GSE159115' { $true }
        'GSE178481' { $row.geo_accession -in $old178 }
        'GSE207493' { $row.source_patient_id -in $old207 }
        'GSE242299' { $row.library_id -notin $old242Missing }
        default { throw "Unexpected study: $($row.study)" }
    }
    [pscustomobject]@{
        study = $row.study
        geo_accession = $row.geo_accession
        source_patient_id = $row.source_patient_id
        library_id = $row.library_id
        tissue = $row.tissue
        comparison_group = if ($isOld) { 'old75' } else { 'added24' }
        incremental_action = if ($isOld) { 'reuse_old_part1_result' } elseif ($row.study -in @('GSE178481','GSE207493')) { 'reuse_existing_full_cohort_scores' } else { 'targeted_official_h5ad_composition_and_tymp_check' }
        input_source = $row.input_source
        input_path_or_archive = $row.input_path_or_archive
        archive_members = $row.archive_members
    }
}
$scopePath = Join-Path $partDir 'incremental_library_scope_v1.tsv'
$out | Export-Csv -LiteralPath $scopePath -Delimiter "`t" -NoTypeInformation -Encoding UTF8

$summary = foreach ($study in @('GSE156632','GSE159115','GSE178481','GSE207493','GSE242299')) {
    $z = @($out | Where-Object study -eq $study)
    $old = @($z | Where-Object comparison_group -eq 'old75')
    $added = @($z | Where-Object comparison_group -eq 'added24')
    [pscustomobject]@{
        study = $study
        old_libraries = $old.Count
        added_libraries = $added.Count
        expanded_libraries = $z.Count
        old_patients = @($old.source_patient_id | Sort-Object -Unique).Count
        added_library_patients = @($added.source_patient_id | Sort-Object -Unique).Count
        newly_introduced_patients = @($added.source_patient_id | Where-Object { $_ -notin $old.source_patient_id } | Sort-Object -Unique).Count
    }
}
$summaryPath = Join-Path $partDir 'incremental_scope_summary_v1.tsv'
$summary | Export-Csv -LiteralPath $summaryPath -Delimiter "`t" -NoTypeInformation -Encoding UTF8

$checks = @(
    [pscustomobject]@{check='expanded_total';status=if($out.Count -eq 99){'PASS'}else{'FAIL'};observed=$out.Count;expected=99},
    [pscustomobject]@{check='old_total';status=if(@($out|Where-Object comparison_group -eq 'old75').Count -eq 75){'PASS'}else{'FAIL'};observed=@($out|Where-Object comparison_group -eq 'old75').Count;expected=75},
    [pscustomobject]@{check='added_total';status=if(@($out|Where-Object comparison_group -eq 'added24').Count -eq 24){'PASS'}else{'FAIL'};observed=@($out|Where-Object comparison_group -eq 'added24').Count;expected=24},
    [pscustomobject]@{check='per_study_delta';status=if((($summary.added_libraries -join ',') -eq '0,0,15,7,2')){'PASS'}else{'FAIL'};observed=($summary.added_libraries -join ',');expected='0,0,15,7,2'},
    [pscustomobject]@{check='new_patients';status=if((($summary.newly_introduced_patients|Measure-Object -Sum).Sum -eq 13)){'PASS'}else{'FAIL'};observed=(($summary.newly_introduced_patients|Measure-Object -Sum).Sum);expected=13},
    [pscustomobject]@{check='gse242299_added_set';status=if(((@($out|Where-Object {$_.study -eq 'GSE242299' -and $_.comparison_group -eq 'added24'}).library_id|Sort-Object) -join ',') -eq 'N3_1,N4_1'){'PASS'}else{'FAIL'};observed=((@($out|Where-Object {$_.study -eq 'GSE242299' -and $_.comparison_group -eq 'added24'}).library_id|Sort-Object) -join ',');expected='N3_1,N4_1'}
)
$checksPath = Join-Path $partDir 'incremental_scope_validation_v1.tsv'
$checks | Export-Csv -LiteralPath $checksPath -Delimiter "`t" -NoTypeInformation -Encoding UTF8
if ($checks.status -contains 'FAIL') { throw 'Incremental scope validation failed' }

$scopeHash = (Get-FileHash -Algorithm SHA256 $scopePath).Hash.ToLowerInvariant()
Set-Content -LiteralPath (Join-Path $partDir 'INCREMENTAL_SCOPE_PASS_v1') -Encoding UTF8 -Value @"
status=PASS
old_libraries=75
added_libraries=24
expanded_libraries=99
newly_introduced_patients=13
added_by_study=GSE178481:15|GSE207493:7|GSE242299:2
gse242299_added=N3_1|N4_1
scope_sha256=$scopeHash
"@
Move-Item -LiteralPath $partDir -Destination $finalDir
Write-Output "PASS scope_sha256=$scopeHash"
