[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$runner = Join-Path $PSScriptRoot 'monitor_cnv_then_finalize_v2.ps1'
try {
    & $runner `
        -RemoteChainVersion 'v3' `
        -RemoteValidationRun '20260812_final_v5' `
        -MonitorVersion 'v4' `
        -SyncVersion 'v5'
    exit $LASTEXITCODE
}
catch {
    $monitorRoot = Join-Path $PSScriptRoot 'monitor_cnv_finalize_v4'
    $status = Join-Path $monitorRoot 'MONITOR_STATUS.tsv'
    $log = Join-Path $monitorRoot 'monitor.log'
    $line = "$(Get-Date -Format o)`tlauncher`tFAIL`t$($_.Exception.Message -replace "`r?`n", ' ')"
    Add-Content -LiteralPath $status -Value $line -Encoding utf8
    Add-Content -LiteralPath $log -Value $line -Encoding utf8
    throw
}

