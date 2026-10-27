[CmdletBinding()]
param(
  [string]$Root = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'
$Root = (Resolve-Path -LiteralPath $Root).Path
$RFile = Join-Path $PSScriptRoot 'phase6a_validation_readiness_local_v1.R'
$LogDir = Join-Path $Root 'logs_validation_local_v1'
New-Item -ItemType Directory -Path $LogDir -Force | Out-Null

$RscriptCandidates = @(@(
  $env:R_SCRIPT,
  'D:\R-4.5.2\bin\Rscript.exe',
  (Get-Command Rscript.exe -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Source -First 1)
) | Where-Object { $_ -and (Test-Path -LiteralPath $_) } | Select-Object -Unique)
if (-not $RscriptCandidates) { throw 'Rscript.exe was not found' }
$Rscript = $RscriptCandidates[0]

$tokens = $null
$errors = $null
[System.Management.Automation.Language.Parser]::ParseFile($PSCommandPath, [ref]$tokens, [ref]$errors) | Out-Null
if ($errors.Count -gt 0) { throw ('PowerShell parse failed: ' + (($errors | ForEach-Object Message) -join '; ')) }

$parseExpression = "parse(file=tail(commandArgs(trailingOnly=TRUE), 1), keep.source=TRUE); cat('R_PARSE=PASS\n')"
$parseOutput = & $Rscript --vanilla -e $parseExpression --args $RFile 2>&1
if ($LASTEXITCODE -ne 0) { throw ('R parse failed: ' + ($parseOutput -join [Environment]::NewLine)) }

$runId = [DateTimeOffset]::Now.ToString('yyyyMMdd_HHmmss_fff')
$env:PHASE6A_RUN_ID = $runId
$logPath = Join-Path $LogDir "phase6a_validation_readiness_$runId.log"
$stdoutLog = Join-Path $LogDir "phase6a_validation_readiness_$runId.stdout.log"
$stderrLog = Join-Path $LogDir "phase6a_validation_readiness_$runId.stderr.log"
$runRecord = Join-Path $LogDir "phase6a_validation_readiness_$runId.json"
[PSCustomObject]@{
  run_id = $runId
  root = $Root
  rscript = $Rscript
  r_script = $RFile
  log = $logPath
  stdout_log = $stdoutLog
  stderr_log = $stderrLog
  started_at = [DateTimeOffset]::Now.ToString('o')
  r_parse = 'PASS'
  powershell_parse = 'PASS'
} | ConvertTo-Json -Depth 3 | Set-Content -LiteralPath $runRecord -Encoding UTF8

$argumentList = @('--vanilla', ('"' + $RFile + '"'), ('"' + $Root + '"'))
$process = Start-Process -FilePath $Rscript -ArgumentList $argumentList -PassThru -Wait -WindowStyle Hidden -RedirectStandardOutput $stdoutLog -RedirectStandardError $stderrLog
$exitCode = $process.ExitCode
@(
  '=== STDOUT ==='
  (Get-Content -LiteralPath $stdoutLog -Raw -ErrorAction SilentlyContinue)
  '=== STDERR ==='
  (Get-Content -LiteralPath $stderrLog -Raw -ErrorAction SilentlyContinue)
) | Set-Content -LiteralPath $logPath -Encoding UTF8
if ($exitCode -ne 0) { throw "Phase 6A failed with exit code $exitCode; log=$logPath" }

$statusPath = Join-Path $Root "singlecell_validation_local_v1\runs\$runId\PHASE6A_VALIDATION_READINESS_STATUS"
if (-not (Test-Path -LiteralPath $statusPath)) { throw "Phase 6A status file missing: $statusPath" }
Write-Output "RUN_ID=$runId"
Write-Output "LOG=$logPath"
Write-Output "STATUS=$statusPath"
