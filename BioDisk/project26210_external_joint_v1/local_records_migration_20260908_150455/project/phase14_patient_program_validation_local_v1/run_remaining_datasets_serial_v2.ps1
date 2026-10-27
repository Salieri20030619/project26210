[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$Rscript = 'D:\R-4.5.2\bin\Rscript.exe'
$RunId = '20260812_082000_provisional_v3'
$RunRoot = Join-Path $PSScriptRoot "runs\$RunId"
$ScoreScript = Join-Path $PSScriptRoot '01_score_dataset_programs_v1.R'
$ValidateScript = Join-Path $PSScriptRoot '02_validate_dataset_programs_v1.R'
$Datasets = @('GSE202374', 'GSE178481', 'GSE264586', 'GSE294109', 'GSE207493')

if (-not (Test-Path -LiteralPath $Rscript -PathType Leaf)) { throw "Rscript absent: $Rscript" }
if (-not (Test-Path -LiteralPath $ScoreScript -PathType Leaf)) { throw "Score script absent" }
if (-not (Test-Path -LiteralPath $ValidateScript -PathType Leaf)) { throw "Validator absent" }

$LogRoot = Join-Path $RunRoot 'serial_logs'
New-Item -ItemType Directory -Path $LogRoot -Force | Out-Null
$Rows = [System.Collections.Generic.List[object]]::new()

foreach ($Dataset in $Datasets) {
    $OutputDir = Join-Path $RunRoot $Dataset
    $LogPath = Join-Path $LogRoot "${Dataset}_score_v3.log"
    $ValidationLog = Join-Path $LogRoot "${Dataset}_validate_v3.log"
    if (Test-Path -LiteralPath $OutputDir) { throw "Refusing existing dataset output: $OutputDir" }
    if (Test-Path -LiteralPath $LogPath) { throw "Refusing existing log: $LogPath" }
    if (Test-Path -LiteralPath $ValidationLog) { throw "Refusing existing log: $ValidationLog" }

    $Started = (Get-Date).ToString('o')
    $ErrorActionPreference = 'Continue'
    & $Rscript $ScoreScript $ProjectRoot $Dataset $OutputDir *>&1 | Tee-Object -FilePath $LogPath
    $ScoreExit = $LASTEXITCODE
    $ErrorActionPreference = 'Stop'
    $StatusPath = Join-Path $OutputDir "status\${Dataset}_STATUS_v1.tsv"
    $ScorePass = $false
    if (Test-Path -LiteralPath $StatusPath -PathType Leaf) {
        $Status = Import-Csv -Delimiter "`t" -LiteralPath $StatusPath
        $ScorePass = [bool]($Status | Where-Object { $_.stage -eq 'OVERALL' -and $_.result -eq 'PASS_PROVISIONAL_MALIGNANCY' })
    }
    if (-not $ScorePass) { throw "Scoring did not produce a valid PASS status for $Dataset (exit=$ScoreExit)" }

    $ErrorActionPreference = 'Continue'
    & $Rscript $ValidateScript $Dataset $OutputDir *>&1 | Tee-Object -FilePath $ValidationLog
    $ValidateExit = $LASTEXITCODE
    $ErrorActionPreference = 'Stop'
    $ValidationPath = Join-Path $OutputDir "status\${Dataset}_INDEPENDENT_VALIDATION_v1.tsv"
    $ValidatePass = $false
    if (Test-Path -LiteralPath $ValidationPath -PathType Leaf) {
        $Validation = Import-Csv -Delimiter "`t" -LiteralPath $ValidationPath
        $ValidatePass = [bool]($Validation | Where-Object { $_.check -eq 'OVERALL' -and $_.result -eq 'PASS' })
    }
    if (-not $ValidatePass) { throw "Independent validation failed for $Dataset (exit=$ValidateExit)" }

    $Rows.Add([pscustomobject]@{
        dataset = $Dataset
        started = $Started
        completed = (Get-Date).ToString('o')
        score_process_exit = $ScoreExit
        score_status = 'PASS_PROVISIONAL_MALIGNANCY'
        validator_process_exit = $ValidateExit
        validation_status = 'PASS'
        output_dir = $OutputDir
    })
}

$Summary = Join-Path $RunRoot 'SERIAL_ALL_DATASETS_STATUS_v3.tsv'
$Part = "$Summary.part"
if (Test-Path -LiteralPath $Summary) { throw "Refusing existing summary: $Summary" }
if (Test-Path -LiteralPath $Part) { throw "Refusing existing partial summary: $Part" }
$Rows | Export-Csv -Delimiter "`t" -NoTypeInformation -Encoding utf8 -LiteralPath $Part
Move-Item -LiteralPath $Part -Destination $Summary
