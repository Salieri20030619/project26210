[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$RunRoot = Join-Path $PSScriptRoot 'continuation_final_v6'
$Status = Join-Path $RunRoot 'CONTINUATION_STATUS.tsv'
$RStdout = Join-Path $RunRoot 'cnv_gate.stdout.log'
$RStderr = Join-Path $RunRoot 'cnv_gate.stderr.log'
$PyStdout = Join-Path $RunRoot 'patient_analysis.stdout.log'
$PyStderr = Join-Path $RunRoot 'patient_analysis.stderr.log'
$Rscript = 'D:\R-4.5.2\bin\Rscript.exe'
$Python = (Get-Command python -ErrorAction Stop).Source
$ProvisionalRun = Join-Path $PSScriptRoot 'runs\20260812_082000_provisional_v3'
$SyncRoot = Join-Path $ProjectRoot 'phase11_external_cnv_consensus_local_v1\final_server_sync_v5'
$FinalRun = Join-Path $PSScriptRoot 'runs\20260812_final_cnv_gated_v6'
$AnalysisOut = Join-Path $PSScriptRoot 'patient_analysis_final_cnv_v6'

if (Test-Path -LiteralPath $RunRoot) { throw "Refusing existing continuation root: $RunRoot" }
if (Test-Path -LiteralPath $FinalRun) { throw "Refusing existing final run: $FinalRun" }
if (Test-Path -LiteralPath $AnalysisOut) { throw "Refusing existing analysis output: $AnalysisOut" }
foreach ($required in @(
    'external_cnv_cell_consensus_final_v2.tsv.gz',
    'FILE_MANIFEST_final_v2.tsv',
    'phase9_all_cells_cnv_gate_v1.tsv.gz',
    'RUN_METADATA_final_v2.json',
    'TERMINAL_RDS_VALIDATION_STATUS_v3.json',
    'terminal_rds_validation_v3.tsv'
)) {
    $path = Join-Path $SyncRoot $required
    if (-not (Test-Path -LiteralPath $path) -or (Get-Item -LiteralPath $path).Length -le 0) {
        throw "Missing or empty synchronized input: $path"
    }
}

New-Item -ItemType Directory -Path $RunRoot | Out-Null
"time`tstage`tresult`tdetail" | Set-Content -LiteralPath $Status -Encoding utf8
function Record([string]$Stage,[string]$Result,[string]$Detail) {
    Add-Content -LiteralPath $Status -Encoding utf8 -Value "$(Get-Date -Format o)`t$Stage`t$Result`t$Detail"
}

Record 'preflight' 'PASS' 'six synchronized files nonempty; server RDS validation 111/111 PASS; explicit data.table logical filter'
Record 'cnv_gate' 'RUNNING' $FinalRun
$r = Start-Process -FilePath $Rscript -ArgumentList @(
    (Join-Path $PSScriptRoot '05_finalize_cnv_gated_program_inputs_v4.R'),
    $ProjectRoot,
    $ProvisionalRun,
    $SyncRoot,
    $FinalRun
) -WorkingDirectory $ProjectRoot -WindowStyle Hidden -RedirectStandardOutput $RStdout -RedirectStandardError $RStderr -PassThru -Wait
if ($r.ExitCode -ne 0 -or -not (Test-Path -LiteralPath (Join-Path $FinalRun 'FINAL_CNV_GATED_INPUT_STATUS_v4.tsv'))) {
    Record 'cnv_gate' 'FAIL' "exit=$($r.ExitCode)"
    exit 3
}
Record 'cnv_gate' 'PASS' $FinalRun

$env:PROJECT26210_ANALYSIS_RUN = $FinalRun
$env:PROJECT26210_ANALYSIS_OUT = $AnalysisOut
$env:PROJECT26210_ANALYSIS_VERSION = 'v6'
$env:PROJECT26210_FINAL_MALIGNANCY = '1'
Record 'patient_analysis' 'RUNNING' $AnalysisOut
$p = Start-Process -FilePath $Python -ArgumentList @(
    (Join-Path $PSScriptRoot '03_run_patient_level_analysis_v1.py')
) -WorkingDirectory $ProjectRoot -WindowStyle Hidden -RedirectStandardOutput $PyStdout -RedirectStandardError $PyStderr -PassThru -Wait
if ($p.ExitCode -ne 0 -or -not (Test-Path -LiteralPath (Join-Path $AnalysisOut 'PATIENT_ANALYSIS_STATUS_v6.json'))) {
    Record 'patient_analysis' 'FAIL' "exit=$($p.ExitCode)"
    exit 4
}
Record 'patient_analysis' 'PASS' $AnalysisOut
Record 'OVERALL' 'PASS' 'final CNV gate and patient-level v6 complete'

