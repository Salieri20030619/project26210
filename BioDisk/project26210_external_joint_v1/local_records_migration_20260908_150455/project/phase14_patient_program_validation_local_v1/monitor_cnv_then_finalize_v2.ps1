[CmdletBinding()]
param(
    [string]$RemoteChainVersion = 'v2',
    [string]$RemoteValidationRun = '20260812_final_v4',
    [string]$MonitorVersion = 'v2',
    [string]$SyncVersion = 'v3'
)

$ErrorActionPreference = 'Stop'
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$MonitorRoot = Join-Path $PSScriptRoot "monitor_cnv_finalize_$MonitorVersion"
$Log = Join-Path $MonitorRoot 'monitor.log'
$Status = Join-Path $MonitorRoot 'MONITOR_STATUS.tsv'
$RemoteRoot = '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1'
$RemoteChain = "$RemoteRoot/phase11_external_cnv_consensus_v1/results/FINAL_BACKGROUND_CHAIN_STATUS_$RemoteChainVersion.tsv"
$RemoteConsensus = "$RemoteRoot/phase11_external_cnv_consensus_v1/results/runs/20260812_113000_final_v2"
$RemoteValidation = "$RemoteRoot/phase11_external_cnv_consensus_v1/results/rds_validation_runs/$RemoteValidationRun"
$SyncRoot = Join-Path $ProjectRoot "phase11_external_cnv_consensus_local_v1\final_server_sync_$SyncVersion"
$ProvisionalRun = Join-Path $PSScriptRoot 'runs\20260812_082000_provisional_v3'
$FinalRun = Join-Path $PSScriptRoot 'runs\20260812_final_cnv_gated_v4'
$AnalysisOut = Join-Path $PSScriptRoot 'patient_analysis_final_cnv_v4'
$Rscript = 'D:\R-4.5.2\bin\Rscript.exe'
$SshArgs = @('-F','C:\Users\31937\.ssh\config','-o','BatchMode=yes','-o','ConnectTimeout=20','-o','ConnectionAttempts=1','-o','ServerAliveInterval=10','-o','ServerAliveCountMax=2','sgs10')
$ScpBase = @('-F','C:\Users\31937\.ssh\config','-o','BatchMode=yes','-o','ConnectTimeout=20','-o','ConnectionAttempts=1','-o','ServerAliveInterval=10','-o','ServerAliveCountMax=2')

New-Item -ItemType Directory -Path $MonitorRoot -Force | Out-Null
if (Test-Path -LiteralPath $Status) { throw "Refusing existing monitor status: $Status" }
"time`tstage`tresult`tdetail" | Set-Content -LiteralPath $Status -Encoding utf8
function Record([string]$Stage,[string]$Result,[string]$Detail) {
    $line = "$(Get-Date -Format o)`t$Stage`t$Result`t$Detail"
    Add-Content -LiteralPath $Status -Value $line -Encoding utf8
    Add-Content -LiteralPath $Log -Value $line -Encoding utf8
}

Record 'monitor' 'START' "poll_interval_seconds=600;remote_chain=$RemoteChainVersion;validation_run=$RemoteValidationRun"
$ConsecutiveFailures = 0
while ($true) {
    $savedPreference = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try {
        $remote = & ssh @SshArgs "if test -f '$RemoteChain'; then cat '$RemoteChain'; else echo RUNNING; fi" 2>&1
        $exit = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $savedPreference
    }
    if ($exit -ne 0) {
        $ConsecutiveFailures++
        Record 'poll' 'NETWORK_RETRY' "exit=$exit consecutive=$ConsecutiveFailures"
        Start-Sleep -Seconds 600
        continue
    }
    $ConsecutiveFailures = 0
    $text = ($remote | Out-String)
    if ($text -match "(?m)^OVERALL\s+FAIL(?:\s|$)") {
        Record 'remote_chain' 'FAIL' (($text -replace "`r?`n", '; ').Trim())
        exit 3
    }
    if ($text -notmatch "(?m)^OVERALL\s+PASS(?:\s|$)") {
        Record 'poll' 'RUNNING' 'server isolated environment/validation chain not terminal'
        Start-Sleep -Seconds 600
        continue
    }
    Record 'remote_chain' 'PASS' 'isolated infercnv environment and 111 RDS validation complete'
    break
}

if (Test-Path -LiteralPath $SyncRoot) { throw "Refusing existing sync root: $SyncRoot" }
New-Item -ItemType Directory -Path $SyncRoot | Out-Null
$downloads = @(
    @("sgs10:${RemoteConsensus}/external_cnv_cell_consensus_incremental_v1.tsv.gz", (Join-Path $SyncRoot 'external_cnv_cell_consensus_final_v2.tsv.gz')),
    @("sgs10:${RemoteConsensus}/RUN_METADATA.json", (Join-Path $SyncRoot 'RUN_METADATA_final_v2.json')),
    @("sgs10:${RemoteConsensus}/FILE_MANIFEST.tsv", (Join-Path $SyncRoot 'FILE_MANIFEST_final_v2.tsv')),
    @("sgs10:${RemoteValidation}/TERMINAL_RDS_VALIDATION_STATUS.json", (Join-Path $SyncRoot "TERMINAL_RDS_VALIDATION_STATUS_$RemoteChainVersion.json")),
    @("sgs10:${RemoteValidation}/terminal_rds_validation_v1.tsv", (Join-Path $SyncRoot "terminal_rds_validation_$RemoteChainVersion.tsv")),
    @("sgs10:${RemoteRoot}/phase9_cnv_alignment_v1/results/cnv_gate_v1/phase9_all_cells_cnv_gate_v1.tsv.gz", (Join-Path $SyncRoot 'phase9_all_cells_cnv_gate_v1.tsv.gz'))
)
foreach ($pair in $downloads) {
    $part = "$($pair[1]).part"
    while ($true) {
        if (Test-Path -LiteralPath $part) { Remove-Item -LiteralPath $part -Force }
        $savedPreference = $ErrorActionPreference
        $ErrorActionPreference = 'Continue'
        try {
            & scp @ScpBase $pair[0] $part *>> $Log
            $scpExit = $LASTEXITCODE
        }
        finally {
            $ErrorActionPreference = $savedPreference
        }
        if ($scpExit -eq 0 -and (Test-Path -LiteralPath $part) -and (Get-Item -LiteralPath $part).Length -gt 0) {
            Move-Item -LiteralPath $part -Destination $pair[1]
            break
        }
        Record 'sync' 'NETWORK_RETRY' "source=$($pair[0]);exit=$scpExit"
        Start-Sleep -Seconds 600
    }
}
Record 'sync' 'PASS' "files=$($downloads.Count)"

& $Rscript (Join-Path $PSScriptRoot '05_finalize_cnv_gated_program_inputs_v4.R') $ProjectRoot $ProvisionalRun $SyncRoot $FinalRun *>> $Log
if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath (Join-Path $FinalRun 'FINAL_CNV_GATED_INPUT_STATUS_v4.tsv'))) {
    throw 'Final CNV-gated input builder failed'
}
Record 'cnv_gate' 'PASS' $FinalRun

$env:PROJECT26210_ANALYSIS_RUN = $FinalRun
$env:PROJECT26210_ANALYSIS_OUT = $AnalysisOut
$env:PROJECT26210_ANALYSIS_VERSION = 'v4'
$env:PROJECT26210_FINAL_MALIGNANCY = '1'
& python (Join-Path $PSScriptRoot '03_run_patient_level_analysis_v1.py') *>> $Log
if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath (Join-Path $AnalysisOut 'PATIENT_ANALYSIS_STATUS_v4.json'))) {
    throw 'Final patient analysis failed'
}
Record 'patient_analysis_final' 'PASS' $AnalysisOut
Record 'OVERALL' 'PASS' 'remote validation synced; final CNV gate and patient v4 complete'
