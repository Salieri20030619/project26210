[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)]
    [string]$RunRecord,
    [int]$Tail = 30
)

$ErrorActionPreference = 'Stop'
$RunRecord = [System.IO.Path]::GetFullPath($RunRecord)
if (-not (Test-Path -LiteralPath $RunRecord -PathType Leaf)) { throw "Run record does not exist: $RunRecord" }
$record = Get-Content -LiteralPath $RunRecord -Raw | ConvertFrom-Json
$process = Get-Process -Id ([int]$record.pid) -ErrorAction SilentlyContinue
$alive = $null -ne $process

Write-Output "RUN_ID=$($record.run_id)"
Write-Output "PID=$($record.pid)"
Write-Output "PROCESS_ALIVE=$alive"
Write-Output "RECORDED_STATE=$($record.state)"
Write-Output "RSCRIPT=$($record.rscript)"
Write-Output "R_SCRIPT=$($record.r_script)"
Write-Output "STDOUT_LOG=$($record.stdout_log)"
Write-Output "STDERR_LOG=$($record.stderr_log)"
Write-Output "PROGRESS_FILE=$($record.progress_file)"
Write-Output "STATUS_FILE=$($record.status_file)"

if (Test-Path -LiteralPath $record.progress_file) {
    $progress = @(Import-Csv -LiteralPath $record.progress_file -Delimiter "`t")
    if ($progress.Count -gt 0) {
        $last = $progress[-1]
        Write-Output "LATEST_EVENT=$($last.event)"
        Write-Output "CURRENT_DATASET=$($last.dataset)"
        Write-Output "CURRENT_OBJECT=$($last.object_name)"
        Write-Output "LATEST_DETAIL=$($last.detail)"
        $completed = @($progress | Where-Object {$_.event -eq 'DATASET_PASS'}).Count
        Write-Output "DATASETS_COMPLETED=$completed"
    }
}
if (Test-Path -LiteralPath $record.status_file) {
    Write-Output 'PHASE5_STATUS_BEGIN'
    Get-Content -LiteralPath $record.status_file
    Write-Output 'PHASE5_STATUS_END'
}
if (Test-Path -LiteralPath $record.stdout_log) {
    Write-Output 'STDOUT_TAIL_BEGIN'
    Get-Content -LiteralPath $record.stdout_log -Tail $Tail
    Write-Output 'STDOUT_TAIL_END'
}
if (Test-Path -LiteralPath $record.stderr_log) {
    Write-Output 'STDERR_TAIL_BEGIN'
    Get-Content -LiteralPath $record.stderr_log -Tail $Tail
    Write-Output 'STDERR_TAIL_END'
}
