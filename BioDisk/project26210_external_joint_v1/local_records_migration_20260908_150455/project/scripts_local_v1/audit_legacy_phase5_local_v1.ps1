[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$Root = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$AuditDir = Join-Path $Root 'local_audit_v1'
$Rscript = 'D:\R-4.5.2\bin\Rscript.exe'
$ExpectedR = Join-Path $Root 'scripts_v1\phase5_merge_by_dataset_v1.R'
$RunnerPatterns = @('scripts_v1\*phase5*.ps1','scripts_v1\*phase5*.sh')

function Write-TsvAtomic {
    param($Rows,[string]$Path)
    $part = "$Path.part"
    $Rows | Export-Csv -LiteralPath $part -Delimiter "`t" -NoTypeInformation -Encoding UTF8
    Move-Item -LiteralPath $part -Destination $Path -Force
}

$rows = [System.Collections.Generic.List[object]]::new()
if (Test-Path -LiteralPath $ExpectedR -PathType Leaf) {
    $parseOut = Join-Path $AuditDir 'legacy_phase5_R_parse_output.txt'
    $output = & $Rscript --vanilla -e "f <- commandArgs(TRUE)[1]; tryCatch({parse(file=f, keep.source=TRUE); cat('PARSE_PASS\n'); q('no',0,FALSE)}, error=function(e){cat('PARSE_FAIL\n',conditionMessage(e),'\n',sep=''); q('no',2,FALSE)})" --args $ExpectedR 2>&1
    $exit = $LASTEXITCODE
    Set-Content -LiteralPath $parseOut -Value $output -Encoding UTF8
    $text = Get-Content -LiteralPath $ExpectedR -Raw
    $serverRootToken = '/mnt/' + 'BioDisk'
    $rows.Add([PSCustomObject]@{
        artifact=$ExpectedR; artifact_type='legacy_R'; exists=$true
        parse_status=if($exit -eq 0){'PASS'}else{'FAIL'}
        linux_absolute_path_hits=([regex]::Matches($text,[regex]::Escape($serverRootToken))).Count
        windows_incompatible_command_hits=([regex]::Matches($text,'(?im)\b(nohup|chmod|conda|find|awk)\b')).Count
        evidence=$parseOut
        note='Legacy file found and fully parsed.'
    })
} else {
    $rows.Add([PSCustomObject]@{
        artifact=$ExpectedR; artifact_type='legacy_R'; exists=$false
        parse_status='NOT_FOUND'
        linux_absolute_path_hits='NA'
        windows_incompatible_command_hits='NA'
        evidence='Recursive disk search and support-archive listing contained no phase5/merge_by_dataset artifact.'
        note='Cannot parse a file that was not transferred; it was not launched or reconstructed in place.'
    })
}

$runnerFiles = @()
foreach ($pattern in $RunnerPatterns) { $runnerFiles += @(Get-ChildItem -Path (Join-Path $Root $pattern) -File -ErrorAction SilentlyContinue) }
if ($runnerFiles.Count -eq 0) {
    $rows.Add([PSCustomObject]@{
        artifact=(Join-Path $Root 'scripts_v1\*phase5* runner'); artifact_type='legacy_runner'; exists=$false
        parse_status='NOT_FOUND'; linux_absolute_path_hits='NA'; windows_incompatible_command_hits='NA'
        evidence='Recursive disk search and support-archive listing contained no Phase 5 runner.'
        note='No legacy runner was launched.'
    })
} else {
    foreach ($runner in $runnerFiles) {
        $text = Get-Content -LiteralPath $runner.FullName -Raw
        $serverRootToken = '/mnt/' + 'BioDisk'
        $parseStatus = 'NOT_APPLICABLE'
        if ($runner.Extension -eq '.ps1') {
            $tokens=$null; $errors=$null
            [System.Management.Automation.Language.Parser]::ParseFile($runner.FullName,[ref]$tokens,[ref]$errors) | Out-Null
            $parseStatus = if($errors.Count -eq 0){'PASS'}else{'FAIL'}
        }
        $rows.Add([PSCustomObject]@{
            artifact=$runner.FullName; artifact_type='legacy_runner'; exists=$true; parse_status=$parseStatus
            linux_absolute_path_hits=([regex]::Matches($text,[regex]::Escape($serverRootToken))).Count
            windows_incompatible_command_hits=([regex]::Matches($text,'(?im)\b(nohup|chmod|conda|find|awk)\b')).Count
            evidence='local file scan'; note='Legacy runner was audited but not launched.'
        })
    }
}

$archive = Join-Path $Root 'transfer_control\project_support_v1.tar.gz'
$archivePhase5 = @(& tar -tzf $archive | Where-Object { $_ -match '(?i)phase5|merge_by_dataset' })
$archiveEvidence = Join-Path $AuditDir 'support_archive_phase5_listing.txt'
if ($archivePhase5.Count -eq 0) { Set-Content -LiteralPath $archiveEvidence -Value 'NO_PHASE5_ARTIFACTS' -Encoding UTF8 } else { Set-Content -LiteralPath $archiveEvidence -Value $archivePhase5 -Encoding UTF8 }

Write-TsvAtomic -Rows $rows -Path (Join-Path $AuditDir 'legacy_phase5_script_audit.tsv')
Write-Output "LEGACY_PHASE5_R=$(if(Test-Path -LiteralPath $ExpectedR){'FOUND'}else{'NOT_FOUND'})"
Write-Output "LEGACY_PHASE5_RUNNERS=$($runnerFiles.Count)"
Write-Output "ARCHIVE_PHASE5_ENTRIES=$($archivePhase5.Count)"
