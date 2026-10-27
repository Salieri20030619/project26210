[CmdletBinding()]
param([string]$Root = (Split-Path -Parent $PSScriptRoot))
$ErrorActionPreference = 'Stop'
$Root = (Resolve-Path -LiteralPath $Root).Path
$Rscript = 'D:\R-4.5.2\bin\Rscript.exe'
$RFile = Join-Path $PSScriptRoot 'phase6_process_cluster_local_v1.R'
$BaseDir = Join-Path $Root 'phase6_processing_local_v1'
$LogDir = Join-Path $Root 'logs_validation_local_v1'
$Lock = Join-Path $BaseDir 'RUN_LOCK'
New-Item -ItemType Directory -Force -Path $BaseDir,$LogDir | Out-Null
if (Test-Path -LiteralPath $Lock) { throw "Active or stale lock exists: $Lock" }

$tok=$null; $err=$null
[System.Management.Automation.Language.Parser]::ParseFile($PSCommandPath,[ref]$tok,[ref]$err)|Out-Null
if($err.Count){throw (($err|ForEach-Object Message)-join '; ')}
$expr="invisible(parse(file=tail(commandArgs(trailingOnly=TRUE),1),keep.source=TRUE));cat('R_PARSE=PASS\n')"
$parse=& $Rscript --vanilla -e $expr --args $RFile 2>&1
if($LASTEXITCODE -ne 0){throw ($parse -join [Environment]::NewLine)}

$runId=[DateTimeOffset]::Now.ToString('yyyyMMdd_HHmmss_fff')
$env:PHASE6_RUN_ID=$runId
$runDir=Join-Path $BaseDir "runs\$runId"
$stdout=Join-Path $LogDir "phase6_process_cluster_$runId.stdout.log"
$stderr=Join-Path $LogDir "phase6_process_cluster_$runId.stderr.log"
$pidFile=Join-Path $LogDir "phase6_process_cluster_$runId.pid.txt"
$record=Join-Path $LogDir "phase6_process_cluster_$runId.json"
Set-Content -LiteralPath $Lock -Value $runId -Encoding UTF8
$arguments=@('--vanilla',('"'+$RFile+'"'),('"'+$Root+'"'))
$childPath=[Environment]::GetEnvironmentVariable('Path','Process')
[Environment]::SetEnvironmentVariable('PATH',$null,'Process')
[Environment]::SetEnvironmentVariable('Path',$childPath,'Process')
try {
  $proc=Start-Process -FilePath $Rscript -ArgumentList $arguments -PassThru -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr
} catch {
  if(Test-Path -LiteralPath $Lock){Remove-Item -LiteralPath $Lock -Force}
  throw
}
Set-Content -LiteralPath $pidFile -Value $proc.Id -Encoding UTF8
[PSCustomObject]@{run_id=$runId;pid=$proc.Id;root=$Root;rscript=$Rscript;r_script=$RFile;stdout=$stdout;stderr=$stderr;run_dir=$runDir;started_at=[DateTimeOffset]::Now.ToString('o');r_parse='PASS';powershell_parse='PASS'}|ConvertTo-Json -Depth 3|Set-Content -LiteralPath $record -Encoding UTF8

$deadline=(Get-Date).AddSeconds(45)
$entered=$false
while((Get-Date)-lt $deadline){Start-Sleep -Seconds 3;if(Test-Path -LiteralPath (Join-Path $runDir 'PHASE6_PROCESSING_PROGRESS.tsv')){$entered=$true;break};if($proc.HasExited){break}}
if(-not $entered){throw "Process did not enter first dataset; pid=$($proc.Id); stderr=$stderr"}
Write-Output "RUN_ID=$runId"
Write-Output "PID=$($proc.Id)"
Write-Output "STDOUT=$stdout"
Write-Output "STDERR=$stderr"
Write-Output "RUN_DIR=$runDir"
