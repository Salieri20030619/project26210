[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$Root = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$ExpectedRoot = 'D:\project26210_external_joint_v1'
$AuditDir = Join-Path $Root 'local_audit_v1'
$TransferManifest = Join-Path $Root 'transfer_control\RDS_FILES.tsv'
$ClinicalManifest = Join-Path $Root 'phase4e_clinical_metadata_harmonization_v1\tables\object_clinical_manifest.tsv'
$ObjectAudit = Join-Path $Root 'phase4d_95_object_compatibility_audit_v1\qc\object_level_audit.tsv'
$ExpectedDatasets = @('GSE178481','GSE202374','GSE207493','GSE264586','GSE294109')

if ($Root.TrimEnd('\') -ne $ExpectedRoot) {
    throw "Project root mismatch: actual=$Root expected=$ExpectedRoot"
}
foreach ($required in @($TransferManifest, $ClinicalManifest, $ObjectAudit)) {
    if (-not (Test-Path -LiteralPath $required -PathType Leaf)) { throw "Required input is missing: $required" }
}
if (-not (Test-Path -LiteralPath $AuditDir)) { New-Item -ItemType Directory -Path $AuditDir | Out-Null }

function Write-TsvAtomic {
    param([Parameter(Mandatory)]$Rows, [Parameter(Mandatory)][string]$Path)
    $Part = "$Path.part"
    $Rows | Export-Csv -LiteralPath $Part -Delimiter "`t" -NoTypeInformation -Encoding UTF8
    Move-Item -LiteralPath $Part -Destination $Path -Force
}

function Write-TextAtomic {
    param([Parameter(Mandatory)][string]$Text, [Parameter(Mandatory)][string]$Path)
    $Part = "$Path.part"
    Set-Content -LiteralPath $Part -Value $Text -Encoding UTF8
    Move-Item -LiteralPath $Part -Destination $Path -Force
}

$GeneratedRoots = @(
    (Join-Path $Root 'local_audit_v1'),
    (Join-Path $Root 'logs_local_v1'),
    (Join-Path $Root 'phase5_merge_by_dataset_local_v1')
)
$Inventory = Get-ChildItem -LiteralPath $Root -Recurse -File -Force | Where-Object {
    $full = $_.FullName
    -not ($GeneratedRoots | Where-Object { $full.StartsWith($_, [System.StringComparison]::OrdinalIgnoreCase) })
} | Sort-Object FullName | ForEach-Object {
    [PSCustomObject]@{
        absolute_path   = $_.FullName
        relative_path   = $_.FullName.Substring($Root.Length).TrimStart('\')
        filename        = $_.Name
        extension       = $_.Extension.ToLowerInvariant()
        size_bytes      = $_.Length
        last_write_time = $_.LastWriteTime.ToString('yyyy-MM-ddTHH:mm:sszzz')
    }
}
$InventoryPath = Join-Path $AuditDir 'file_inventory.tsv'
Write-TsvAtomic -Rows $Inventory -Path $InventoryPath

$Transfer = @(Import-Csv -LiteralPath $TransferManifest -Delimiter "`t")
$Clinical = @(Import-Csv -LiteralPath $ClinicalManifest -Delimiter "`t")
$PriorAudit = @(Import-Csv -LiteralPath $ObjectAudit -Delimiter "`t")

$ClinicalIndex = @{}
foreach ($row in $Clinical) { $ClinicalIndex["$($row.dataset)`t$($row.object_name)"] = $row }
$AuditIndex = @{}
foreach ($row in $PriorAudit) { $AuditIndex["$($row.dataset)`t$($row.object_name)"] = $row }
$ObjectDup = @{}; $Transfer | Group-Object object_name | ForEach-Object { $ObjectDup[$_.Name] = $_.Count }
$SourceDup = @{}; $Transfer | Group-Object source_file | ForEach-Object { $SourceDup[$_.Name] = $_.Count }
$RelDup = @{}; $Transfer | Group-Object local_relpath | ForEach-Object { $RelDup[$_.Name] = $_.Count }

$Reconciliation = foreach ($row in $Transfer) {
    $localRel = $row.local_relpath -replace '/', '\'
    $localPath = [System.IO.Path]::GetFullPath((Join-Path $Root $localRel))
    $exists = Test-Path -LiteralPath $localPath -PathType Leaf
    $actualSize = if ($exists) { (Get-Item -LiteralPath $localPath).Length } else { $null }
    $key = "$($row.dataset)`t$($row.object_name)"
    $clinicalRow = $ClinicalIndex[$key]
    $auditRow = $AuditIndex[$key]
    [PSCustomObject]@{
        dataset                    = $row.dataset
        object_name                = $row.object_name
        manifest_source_file       = $row.source_file
        local_relpath              = $row.local_relpath
        local_absolute_path        = $localPath
        expected_size_bytes        = [Int64]$row.size_bytes
        actual_size_bytes          = $actualSize
        exists                     = $exists
        zero_bytes                 = ($exists -and $actualSize -eq 0)
        size_matches_manifest      = ($exists -and $actualSize -eq [Int64]$row.size_bytes)
        manifest_cells             = if ($null -ne $clinicalRow) { [Int64]$clinicalRow.cells } else { $null }
        prior_audit_cells          = if ($null -ne $auditRow) { [Int64]$auditRow.cells } else { $null }
        cells_agree_between_records = ($null -ne $clinicalRow -and $null -ne $auditRow -and [Int64]$clinicalRow.cells -eq [Int64]$auditRow.cells)
        duplicate_object_name      = ($ObjectDup[$row.object_name] -gt 1)
        duplicate_source_file      = ($SourceDup[$row.source_file] -gt 1)
        duplicate_local_relpath    = ($RelDup[$row.local_relpath] -gt 1)
        expected_sha256            = $row.sha256
    }
}
$ReconciliationPath = Join-Path $AuditDir 'manifest_reconciliation.tsv'
Write-TsvAtomic -Rows $Reconciliation -Path $ReconciliationPath

$ManifestPaths = @{}; foreach ($row in $Reconciliation) { $ManifestPaths[$row.local_absolute_path.ToLowerInvariant()] = $true }
$CandidateExtensions = @('.rds','.rda','.rdata','.h5ad','.h5','.hdf5','.mtx')
$CandidateFiles = @(Get-ChildItem -LiteralPath $Root -Recurse -File -Force | Where-Object {
    $name = $_.Name.ToLowerInvariant()
    ($CandidateExtensions -contains $_.Extension.ToLowerInvariant()) -or $name.EndsWith('.mtx.gz')
})
$ExtraCandidates = @($CandidateFiles | Where-Object { -not $ManifestPaths.ContainsKey($_.FullName.ToLowerInvariant()) } | ForEach-Object {
    [PSCustomObject]@{ absolute_path=$_.FullName; relative_path=$_.FullName.Substring($Root.Length).TrimStart('\'); size_bytes=$_.Length; reason='candidate_data_file_not_in_transfer_manifest' }
})
if ($ExtraCandidates.Count -eq 0) { $ExtraCandidates = @([PSCustomObject]@{absolute_path='';relative_path='';size_bytes='';reason='NONE'}) }
Write-TsvAtomic -Rows $ExtraCandidates -Path (Join-Path $AuditDir 'extra_candidate_objects.tsv')

$TempFiles = @(Get-ChildItem -LiteralPath $Root -Recurse -File -Force | Where-Object {
    $_.Name -match '(?i)(\.part|\.partial|\.tmp|\.crdownload|\.download|\.aria2)$'
} | ForEach-Object {
    [PSCustomObject]@{ absolute_path=$_.FullName; relative_path=$_.FullName.Substring($Root.Length).TrimStart('\'); size_bytes=$_.Length; last_write_time=$_.LastWriteTime.ToString('yyyy-MM-ddTHH:mm:sszzz') }
})
if ($TempFiles.Count -eq 0) { $TempFiles = @([PSCustomObject]@{absolute_path='';relative_path='';size_bytes='';last_write_time=''}) }
Write-TsvAtomic -Rows $TempFiles -Path (Join-Path $AuditDir 'incomplete_download_candidates.tsv')

# Hash only a targeted subset because all 95 files already have checksums and sizes agree.
$HashTargets = @($Reconciliation | Group-Object dataset | ForEach-Object { $_.Group | Sort-Object object_name | Select-Object -First 1 })
$Largest = $Reconciliation | Sort-Object actual_size_bytes -Descending | Select-Object -First 1
if ($null -ne $Largest -and -not ($HashTargets.object_name -contains $Largest.object_name)) { $HashTargets += $Largest }
$HashAudit = foreach ($row in $HashTargets) {
    $actualHash = if ($row.exists) { (Get-FileHash -LiteralPath $row.local_absolute_path -Algorithm SHA256).Hash.ToLowerInvariant() } else { '' }
    [PSCustomObject]@{
        dataset=$row.dataset; object_name=$row.object_name; local_absolute_path=$row.local_absolute_path
        size_bytes=$row.actual_size_bytes; expected_sha256=$row.expected_sha256; actual_sha256=$actualHash
        sha256_matches=($row.exists -and $actualHash -eq $row.expected_sha256.ToLowerInvariant())
    }
}
Write-TsvAtomic -Rows $HashAudit -Path (Join-Path $AuditDir 'targeted_sha256_audit.tsv')

$ControlChecks = foreach ($name in @('RDS_FILES.tsv','project_support_v1.tar.gz')) {
    $path = Join-Path $Root "transfer_control\$name"
    $expectedLine = Get-Content -LiteralPath (Join-Path $Root 'transfer_control\CONTROL_SHA256SUMS.txt') | Where-Object { $_ -match "  $([regex]::Escape($name))$" }
    $expected = if ($expectedLine) { ($expectedLine -split '\s+')[0].ToLowerInvariant() } else { '' }
    $actual = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
    [PSCustomObject]@{file=$path;expected_sha256=$expected;actual_sha256=$actual;sha256_matches=($expected -eq $actual)}
}
Write-TsvAtomic -Rows $ControlChecks -Path (Join-Path $AuditDir 'control_sha256_audit.tsv')

$nativeCode = @'
using System;
using System.Runtime.InteropServices;
public static class LocalAuditNativeStatus {
  [StructLayout(LayoutKind.Sequential, CharSet=CharSet.Auto)]
  public class MEMORYSTATUSEX { public uint dwLength=64; public uint dwMemoryLoad; public ulong ullTotalPhys; public ulong ullAvailPhys; public ulong ullTotalPageFile; public ulong ullAvailPageFile; public ulong ullTotalVirtual; public ulong ullAvailVirtual; public ulong ullAvailExtendedVirtual; }
  [DllImport("kernel32.dll", CharSet=CharSet.Auto, SetLastError=true)] public static extern bool GlobalMemoryStatusEx([In, Out] MEMORYSTATUSEX lpBuffer);
  [DllImport("kernel32.dll", CharSet=CharSet.Auto, SetLastError=true)] public static extern bool GetDiskFreeSpaceEx(string lpDirectoryName, out ulong freeBytesAvailable, out ulong totalNumberOfBytes, out ulong totalNumberOfFreeBytes);
}
'@
Add-Type -TypeDefinition $nativeCode
$mem = New-Object LocalAuditNativeStatus+MEMORYSTATUSEX
if (-not [LocalAuditNativeStatus]::GlobalMemoryStatusEx($mem)) { throw 'GlobalMemoryStatusEx failed' }
[UInt64]$diskAvail=0; [UInt64]$diskTotal=0; [UInt64]$diskFree=0
if (-not [LocalAuditNativeStatus]::GetDiskFreeSpaceEx('D:\',[ref]$diskAvail,[ref]$diskTotal,[ref]$diskFree)) { throw 'GetDiskFreeSpaceEx failed' }
$reg = Get-ItemProperty -LiteralPath 'HKLM:\SOFTWARE\Microsoft\Windows NT\CurrentVersion'

$Rscript = 'D:\R-4.5.2\bin\Rscript.exe'
if (-not (Test-Path -LiteralPath $Rscript -PathType Leaf)) {
    $cmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue | Select-Object -First 1
    if ($null -eq $cmd) { $Rscript = '' } else { $Rscript = $cmd.Source }
}
$RValues = @{}
if ($Rscript) {
    $rLines = & $Rscript --vanilla -e "ip <- installed.packages(); cat('r_version=', R.version.string, '\n', sep=''); cat('r_home=', R.home(), '\n', sep=''); for (p in c('Seurat','SeuratObject','Matrix','data.table')) cat('package_',p,'=',if(p %in% rownames(ip)) ip[p,'Version'] else 'MISSING','\n',sep=''); q('no', status=0, runLast=FALSE)"
    foreach ($line in $rLines) { if ($line -match '^([^=]+)=(.*)$') { $RValues[$Matches[1]] = $Matches[2] } }
}
$Environment = @(
    [PSCustomObject]@{category='path';key='project_root';value=$Root;unit='path';status=if($Root -eq $ExpectedRoot){'PASS'}else{'FAIL'}},
    [PSCustomObject]@{category='os';key='product_name';value=$reg.ProductName;unit='text';status='INFO'},
    [PSCustomObject]@{category='os';key='display_version';value=$reg.DisplayVersion;unit='text';status='INFO'},
    [PSCustomObject]@{category='os';key='build';value=("$($reg.CurrentBuild).$($reg.UBR)");unit='version';status='INFO'},
    [PSCustomObject]@{category='os';key='architecture';value=[System.Runtime.InteropServices.RuntimeInformation]::OSArchitecture.ToString();unit='text';status='INFO'},
    [PSCustomObject]@{category='powershell';key='version';value=$PSVersionTable.PSVersion.ToString();unit='version';status='INFO'},
    [PSCustomObject]@{category='disk';key='D_total';value=$diskTotal;unit='bytes';status='INFO'},
    [PSCustomObject]@{category='disk';key='D_free';value=$diskFree;unit='bytes';status=if($diskFree -gt 50GB){'PASS'}else{'WARN'}},
    [PSCustomObject]@{category='memory';key='physical_total';value=$mem.ullTotalPhys;unit='bytes';status='INFO'},
    [PSCustomObject]@{category='memory';key='physical_available';value=$mem.ullAvailPhys;unit='bytes';status=if($mem.ullAvailPhys -gt 8GB){'PASS'}else{'WARN'}},
    [PSCustomObject]@{category='memory';key='load';value=$mem.dwMemoryLoad;unit='percent';status='INFO'},
    [PSCustomObject]@{category='R';key='Rscript';value=$Rscript;unit='path';status=if($Rscript){'PASS'}else{'MISSING'}},
    [PSCustomObject]@{category='R';key='version';value=$RValues['r_version'];unit='version';status=if($RValues['r_version']){'PASS'}else{'MISSING'}},
    [PSCustomObject]@{category='R';key='R_home';value=$RValues['r_home'];unit='path';status='INFO'},
    [PSCustomObject]@{category='R_package';key='Seurat';value=$RValues['package_Seurat'];unit='version';status=if($RValues['package_Seurat'] -and $RValues['package_Seurat'] -ne 'MISSING'){'PASS'}else{'MISSING'}},
    [PSCustomObject]@{category='R_package';key='SeuratObject';value=$RValues['package_SeuratObject'];unit='version';status=if($RValues['package_SeuratObject'] -and $RValues['package_SeuratObject'] -ne 'MISSING'){'PASS'}else{'MISSING'}},
    [PSCustomObject]@{category='R_package';key='Matrix';value=$RValues['package_Matrix'];unit='version';status=if($RValues['package_Matrix'] -and $RValues['package_Matrix'] -ne 'MISSING'){'PASS'}else{'MISSING'}},
    [PSCustomObject]@{category='R_package';key='data.table';value=$RValues['package_data.table'];unit='version';status=if($RValues['package_data.table'] -and $RValues['package_data.table'] -ne 'MISSING'){'PASS'}else{'MISSING'}}
)
Write-TsvAtomic -Rows $Environment -Path (Join-Path $AuditDir 'environment_audit.tsv')

$DatasetRows = $Reconciliation | Group-Object dataset | ForEach-Object {
    [PSCustomObject]@{dataset=$_.Name;source_objects=$_.Count;manifest_cells=($_.Group|Measure-Object manifest_cells -Sum).Sum;files_present=@($_.Group|Where-Object exists).Count;size_matches=@($_.Group|Where-Object size_matches_manifest).Count}
} | Sort-Object dataset
Write-TsvAtomic -Rows $DatasetRows -Path (Join-Path $AuditDir 'dataset_manifest_summary.tsv')

$Missing = @($Reconciliation | Where-Object { -not $_.exists })
$SizeMismatch = @($Reconciliation | Where-Object { -not $_.size_matches_manifest })
$Zero = @($Reconciliation | Where-Object zero_bytes)
$DupObject = @($Reconciliation | Where-Object duplicate_object_name)
$DupSource = @($Reconciliation | Where-Object duplicate_source_file)
$UnexpectedDatasets = @($DatasetRows.dataset | Where-Object { $ExpectedDatasets -notcontains $_ })
$MissingDatasets = @($ExpectedDatasets | Where-Object { $DatasetRows.dataset -notcontains $_ })
$RealExtraCount = @($ExtraCandidates | Where-Object {$_.reason -ne 'NONE'}).Count
$RealTempCount = @($TempFiles | Where-Object {$_.absolute_path}).Count
$InputComplete = ($Transfer.Count -eq 95 -and $DatasetRows.Count -eq 5 -and $MissingDatasets.Count -eq 0 -and $UnexpectedDatasets.Count -eq 0 -and $Missing.Count -eq 0 -and $SizeMismatch.Count -eq 0 -and $Zero.Count -eq 0 -and $DupObject.Count -eq 0 -and $DupSource.Count -eq 0 -and $RealExtraCount -eq 0 -and $RealTempCount -eq 0 -and @($HashAudit|Where-Object {-not $_.sha256_matches}).Count -eq 0)

$datasetLines = ($DatasetRows | ForEach-Object { "| $($_.dataset) | $($_.source_objects) | $($_.manifest_cells) | $($_.files_present) | $($_.size_matches) |" }) -join "`r`n"
$Summary = @"
# Local input audit summary

Generated: $([DateTimeOffset]::Now.ToString('yyyy-MM-ddTHH:mm:sszzz'))

## Manifest

- Transfer manifest: $TransferManifest
- Fields: $((($Transfer | Select-Object -First 1).PSObject.Properties.Name) -join ', ')
- Rows: $($Transfer.Count)
- Clinical manifest: $ClinicalManifest
- Clinical rows: $($Clinical.Count)
- Prior object audit: $ObjectAudit
- Prior audit rows: $($PriorAudit.Count)

| dataset | source objects | manifest cells | files present | size matches |
|---|---:|---:|---:|---:|
$datasetLines

## Local reconciliation

- Candidate single-cell data files found: $($CandidateFiles.Count)
- Manifest files present: $(@($Reconciliation | Where-Object exists).Count) / $($Transfer.Count)
- Exact size matches: $(@($Reconciliation | Where-Object size_matches_manifest).Count) / $($Transfer.Count)
- Missing files: $($Missing.Count)
- Zero-byte manifest files: $($Zero.Count)
- Size mismatches: $($SizeMismatch.Count)
- Duplicate object names: $($DupObject.Count)
- Duplicate source-file records: $($DupSource.Count)
- Candidate data files outside manifest: $RealExtraCount
- Incomplete-download filename candidates: $RealTempCount
- Targeted per-object SHA256 checks: $(@($HashAudit | Where-Object sha256_matches).Count) / $($HashAudit.Count) PASS
- Control-file SHA256 checks: $(@($ControlChecks | Where-Object sha256_matches).Count) / $($ControlChecks.Count) PASS
- Missing expected datasets: $(if($MissingDatasets.Count){$MissingDatasets -join ', '}else{'NONE'})
- Unexpected datasets: $(if($UnexpectedDatasets.Count){$UnexpectedDatasets -join ', '}else{'NONE'})

## Decision

- Lightweight input gate: $(if($InputComplete){'PASS'}else{'BLOCKED_INPUT_INCOMPLETE'})
- Source-object content is not inferred from filenames. Full readRDS() and raw-count validation will occur once per source object during Phase 5, with incremental evidence written to the Phase 5 audit tables.
- No source data file was modified by this audit.
"@
Write-TextAtomic -Text $Summary -Path (Join-Path $AuditDir 'input_audit_summary.md')

$Gate = @(
    "status=$(if($InputComplete){'PASS'}else{'BLOCKED_INPUT_INCOMPLETE'})",
    "manifest_rows=$($Transfer.Count)",
    "datasets=$($DatasetRows.Count)",
    "manifest_files_present=$(@($Reconciliation | Where-Object exists).Count)",
    "size_matches=$(@($Reconciliation | Where-Object size_matches_manifest).Count)",
    "targeted_sha256_pass=$(@($HashAudit | Where-Object sha256_matches).Count)",
    "targeted_sha256_total=$($HashAudit.Count)",
    "generated_at=$([DateTimeOffset]::Now.ToString('yyyy-MM-ddTHH:mm:sszzz'))"
) -join "`r`n"
Write-TextAtomic -Text $Gate -Path (Join-Path $AuditDir 'INPUT_GATE_STATUS')

Write-Output "INPUT_GATE=$(if($InputComplete){'PASS'}else{'BLOCKED_INPUT_INCOMPLETE'})"
Write-Output "MANIFEST_ROWS=$($Transfer.Count)"
Write-Output "DATASETS=$($DatasetRows.Count)"
Write-Output "CANDIDATE_DATA_FILES=$($CandidateFiles.Count)"
Write-Output "MANIFEST_FILES_PRESENT=$(@($Reconciliation | Where-Object exists).Count)"
Write-Output "DISK_FREE_BYTES=$diskFree"
Write-Output "MEMORY_AVAILABLE_BYTES=$($mem.ullAvailPhys)"
