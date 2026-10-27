[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$Root = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$ExpectedRoot = 'D:\project26210_external_joint_v1'
$Rscript = 'D:\R-4.5.2\bin\Rscript.exe'
$RFile = Join-Path $PSScriptRoot 'phase5_merge_by_dataset_local_v2.R'
$MonitorFile = Join-Path $PSScriptRoot 'monitor_phase5_merge_local_v2.ps1'
$LogDir = Join-Path $Root 'logs_local_v2'
$ResultDir = Join-Path $Root 'phase5_merge_by_dataset_local_v2'
$TableDir = Join-Path $ResultDir 'tables'
$ObjectDir = Join-Path $ResultDir 'objects_raw'
$InputAuditDir = Join-Path $Root 'local_audit_v1'
$AuditDir = Join-Path $Root 'local_audit_v2'

if ($Root.TrimEnd('\') -ne $ExpectedRoot) { throw "Project root mismatch: actual=$Root expected=$ExpectedRoot" }
foreach ($required in @($Rscript,$RFile,$MonitorFile,(Join-Path $InputAuditDir 'INPUT_GATE_STATUS'))) {
    if (-not (Test-Path -LiteralPath $required -PathType Leaf)) { throw "Required launch input is missing: $required" }
}
if ((Get-Content -LiteralPath (Join-Path $InputAuditDir 'INPUT_GATE_STATUS') | Select-String '^status=').Line -ne 'status=PASS') {
    throw 'Local input gate is not PASS; Phase 5 will not start.'
}
if (Test-Path -LiteralPath (Join-Path $ResultDir 'RUN_LOCK')) { throw 'RUN_LOCK exists; refusing a duplicate Phase 5 launch.' }
if (Test-Path -LiteralPath (Join-Path $ResultDir 'PHASE5_MERGE_STATUS')) { throw 'PHASE5_MERGE_STATUS already exists; use a higher version instead of overwriting.' }
$existingResults = @(Get-ChildItem -LiteralPath $ObjectDir -File -ErrorAction SilentlyContinue) + @(Get-ChildItem -LiteralPath $TableDir -File -ErrorAction SilentlyContinue)
if ($existingResults.Count -gt 0) { throw "Phase 5 result directory is not empty; refusing overwrite. Files=$($existingResults.Count)" }

# Parse this runner and monitor with the PowerShell AST parser.
foreach ($psFile in @($MyInvocation.MyCommand.Path,$MonitorFile)) {
    $tokens=$null; $errors=$null
    [System.Management.Automation.Language.Parser]::ParseFile($psFile,[ref]$tokens,[ref]$errors) | Out-Null
    if ($errors.Count -gt 0) {
        $context = ($errors | ForEach-Object { "$($_.Extent.File):$($_.Extent.StartLineNumber):$($_.Extent.StartColumnNumber) $($_.Message)" }) -join [Environment]::NewLine
        throw "PowerShell parse failed:`n$context"
    }
}

# Fully parse the R script and refuse Linux server paths or commands in the local implementation.
$parseOutput = & $Rscript --vanilla -e "f <- commandArgs(trailingOnly=TRUE)[1]; invisible(parse(file=f, keep.source=TRUE)); cat('R_PARSE=PASS\n'); q('no',0,FALSE)" $RFile 2>&1
if ($LASTEXITCODE -ne 0 -or ($parseOutput -notcontains 'R_PARSE=PASS')) { throw "R parse failed: $($parseOutput -join [Environment]::NewLine)" }
$rText = Get-Content -LiteralPath $RFile -Raw
$serverRootToken = '/mnt/' + 'BioDisk'
if ($rText.Contains($serverRootToken)) { throw 'Local R script contains a forbidden Linux server-root path.' }
if ($rText -match '(?im)\b(nohup|chmod|conda|find|awk)\b') { throw 'Local R script contains a Windows-incompatible shell command token.' }

# Detect an active prior run from its exact record and PID before launching.
$priorRecords = @(Get-ChildItem -LiteralPath $LogDir -Filter 'phase5_run_*.json' -File -ErrorAction SilentlyContinue)
foreach ($recordFile in $priorRecords) {
    try {
        $record = Get-Content -LiteralPath $recordFile.FullName -Raw | ConvertFrom-Json
        $process = Get-Process -Id ([int]$record.pid) -ErrorAction SilentlyContinue
        if ($null -ne $process -and $process.Path -eq $Rscript -and $record.r_script -eq $RFile -and $record.state -notmatch 'EXITED|COMPLETED|FAILED') {
            throw "An existing recorded Phase 5 command is active: PID=$($record.pid) record=$($recordFile.FullName)"
        }
    } catch [System.Management.Automation.RuntimeException] {
        throw
    } catch {
        # An unreadable stale record is evidence to stop, not a reason to guess and relaunch.
        throw "Cannot safely audit prior run record $($recordFile.FullName): $($_.Exception.Message)"
    }
}

$nativeCode = @'
using System;
using System.Runtime.InteropServices;
public static class Phase5LaunchNativeStatus {
  [StructLayout(LayoutKind.Sequential, CharSet=CharSet.Auto)]
  public class MEMORYSTATUSEX { public uint dwLength=64; public uint dwMemoryLoad; public ulong ullTotalPhys; public ulong ullAvailPhys; public ulong ullTotalPageFile; public ulong ullAvailPageFile; public ulong ullTotalVirtual; public ulong ullAvailVirtual; public ulong ullAvailExtendedVirtual; }
  [DllImport("kernel32.dll", CharSet=CharSet.Auto, SetLastError=true)] public static extern bool GlobalMemoryStatusEx([In, Out] MEMORYSTATUSEX lpBuffer);
  [DllImport("kernel32.dll", CharSet=CharSet.Auto, SetLastError=true)] public static extern bool GetDiskFreeSpaceEx(string lpDirectoryName, out ulong freeBytesAvailable, out ulong totalNumberOfBytes, out ulong totalNumberOfFreeBytes);
}
'@
Add-Type -TypeDefinition $nativeCode
$mem = New-Object Phase5LaunchNativeStatus+MEMORYSTATUSEX
if (-not [Phase5LaunchNativeStatus]::GlobalMemoryStatusEx($mem)) { throw 'GlobalMemoryStatusEx failed' }
[UInt64]$diskAvail=0; [UInt64]$diskTotal=0; [UInt64]$diskFree=0
if (-not [Phase5LaunchNativeStatus]::GetDiskFreeSpaceEx('D:\',[ref]$diskAvail,[ref]$diskTotal,[ref]$diskFree)) { throw 'GetDiskFreeSpaceEx failed' }

$priorAudit = Import-Csv -LiteralPath (Join-Path $Root 'phase4d_95_object_compatibility_audit_v1\qc\object_level_audit.tsv') -Delimiter "`t"
$resourceRows = $priorAudit | Group-Object dataset | ForEach-Object {
    $nnz = [double](($_.Group | Measure-Object nonzero_values -Sum).Sum)
    $cells = [double](($_.Group | Measure-Object cells -Sum).Sum)
    $matrixBytes = 12.0*$nnz + 4.0*($cells+1)
    $maxSource = ($_.Group | ForEach-Object {12.0*[double]$_.nonzero_values + 4.0*([double]$_.cells+1)} | Measure-Object -Maximum).Maximum
    $peak = [Math]::Ceiling(1.6*$matrixBytes + $maxSource + 1GB)
    [PSCustomObject]@{dataset=$_.Name;prior_nnz=[Int64]$nnz;estimated_merged_dgC_bytes=[Int64]$matrixBytes;max_source_dgC_bytes=[Int64]$maxSource;estimated_peak_memory_bytes=[Int64]$peak;required_with_25pct_margin_bytes=[Int64][Math]::Ceiling(1.25*$peak);available_memory_bytes=[Int64]$mem.ullAvailPhys;memory_gate=([double]$mem.ullAvailPhys -ge 1.25*$peak)}
}
$estimatedOutput = [Math]::Ceiling((($resourceRows | Measure-Object estimated_merged_dgC_bytes -Sum).Sum)*1.15)
$estimatedWork = [Math]::Ceiling((($resourceRows | Measure-Object estimated_merged_dgC_bytes -Maximum).Maximum)*1.25)
$diskRequired = $estimatedOutput + $estimatedWork + 10GB
if (@($resourceRows | Where-Object {-not $_.memory_gate}).Count -gt 0) {
    $bad = $resourceRows | Where-Object {-not $_.memory_gate} | Sort-Object required_with_25pct_margin_bytes -Descending | Select-Object -First 1
    throw "Memory preflight failed before launch: dataset=$($bad.dataset) available=$($mem.ullAvailPhys) required=$($bad.required_with_25pct_margin_bytes)"
}
if ([double]$diskFree -lt $diskRequired) { throw "Disk preflight failed before launch: free=$diskFree required=$diskRequired" }

$runId = [DateTimeOffset]::Now.ToString('yyyyMMdd_HHmmss_fff')
$preflightPath = Join-Path $AuditDir "phase5_launch_preflight_$runId.tsv"
$resourceRows | Select-Object *,@{n='free_disk_bytes';e={[Int64]$diskFree}},@{n='disk_required_bytes';e={[Int64]$diskRequired}},@{n='disk_gate';e={$true}} | Export-Csv -LiteralPath $preflightPath -Delimiter "`t" -NoTypeInformation -Encoding UTF8
$stdoutLog = Join-Path $LogDir "phase5_merge_$runId.stdout.log"
$stderrLog = Join-Path $LogDir "phase5_merge_$runId.stderr.log"
$pidFile = Join-Path $LogDir "phase5_pid_$runId.txt"
$runRecord = Join-Path $LogDir "phase5_run_$runId.json"
$progressFile = Join-Path $ResultDir "PHASE5_PROGRESS_$runId.tsv"

$env:PHASE5_RUN_ID = $runId
$env:PHASE5_AVAILABLE_MEMORY_BYTES = [string]$mem.ullAvailPhys
$env:PHASE5_FREE_DISK_BYTES = [string]$diskFree
# The managed host can expose both Path and PATH. Windows PowerShell 5.1 Start-Process
# builds a case-insensitive environment dictionary and fails on that duplicate pair.
$childPath = [Environment]::GetEnvironmentVariable('Path','Process')
[Environment]::SetEnvironmentVariable('PATH',$null,'Process')
[Environment]::SetEnvironmentVariable('Path',$childPath,'Process')
$arguments = @('--vanilla',$RFile)
$process = Start-Process -FilePath $Rscript -ArgumentList $arguments -RedirectStandardOutput $stdoutLog -RedirectStandardError $stderrLog -WindowStyle Hidden -PassThru
Set-Content -LiteralPath $pidFile -Value $process.Id -Encoding ASCII
$recordObject = [ordered]@{
    run_id=$runId; pid=$process.Id; state='STARTED'; started_at=[DateTimeOffset]::Now.ToString('o')
    rscript=$Rscript; r_script=$RFile; arguments=($arguments -join ' ')
    stdout_log=$stdoutLog; stderr_log=$stderrLog; pid_file=$pidFile; progress_file=$progressFile
    status_file=(Join-Path $ResultDir 'PHASE5_MERGE_STATUS'); preflight_file=$preflightPath
    available_memory_bytes=[Int64]$mem.ullAvailPhys; free_disk_bytes=[Int64]$diskFree
}
$recordObject | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $runRecord -Encoding UTF8

$deadline = (Get-Date).AddSeconds(45)
$enteredDataset = $false
do {
    Start-Sleep -Seconds 2
    $live = Get-Process -Id $process.Id -ErrorAction SilentlyContinue
    if (Test-Path -LiteralPath $stdoutLog) {
        $enteredDataset = [bool](Select-String -LiteralPath $stdoutLog -Pattern 'DATASET_START dataset=' -Quiet)
    }
    if ($null -eq $live -or $enteredDataset) { break }
} while ((Get-Date) -lt $deadline)

$live = Get-Process -Id $process.Id -ErrorAction SilentlyContinue
if ($null -eq $live) {
    $recordObject.state = 'EXITED_DURING_STARTUP_CHECK'
    $recordObject | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $runRecord -Encoding UTF8
    $stderrTail = if(Test-Path -LiteralPath $stderrLog){(Get-Content -LiteralPath $stderrLog -Tail 30) -join [Environment]::NewLine}else{''}
    $stdoutTail = if(Test-Path -LiteralPath $stdoutLog){(Get-Content -LiteralPath $stdoutLog -Tail 30) -join [Environment]::NewLine}else{''}
    throw "Phase 5 exited during startup validation. stdout=$stdoutTail stderr=$stderrTail"
}
$recordObject.state = if($enteredDataset){'RUNNING_DATASET'}else{'RUNNING_STARTUP'}
$recordObject.dataset_marker_seen = $enteredDataset
$recordObject.startup_checked_at = [DateTimeOffset]::Now.ToString('o')
$recordObject | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $runRecord -Encoding UTF8

Write-Output 'POWERSHELL_PARSE=PASS'
Write-Output 'R_PARSE=PASS'
Write-Output "RUN_ID=$runId"
Write-Output "PID=$($process.Id)"
Write-Output "STATE=$($recordObject.state)"
Write-Output "DATASET_MARKER_SEEN=$enteredDataset"
Write-Output "STDOUT_LOG=$stdoutLog"
Write-Output "STDERR_LOG=$stderrLog"
Write-Output "RUN_RECORD=$runRecord"
Write-Output "MONITOR_COMMAND=& '$MonitorFile' -RunRecord '$runRecord'"
