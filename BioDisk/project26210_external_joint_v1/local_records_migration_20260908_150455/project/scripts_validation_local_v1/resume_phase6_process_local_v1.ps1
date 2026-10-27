[CmdletBinding()]
param([Parameter(Mandatory=$true)][string]$RunId,[string]$Root=(Split-Path -Parent $PSScriptRoot))
$ErrorActionPreference='Stop'
$Root=(Resolve-Path -LiteralPath $Root).Path
$Rscript='D:\R-4.5.2\bin\Rscript.exe'
$RFile=Join-Path $PSScriptRoot 'phase6_process_cluster_local_v1.R'
$BaseDir=Join-Path $Root 'phase6_processing_local_v1'
$RunDir=Join-Path $BaseDir "runs\$RunId"
$Lock=Join-Path $BaseDir 'RUN_LOCK'
$LogDir=Join-Path $Root 'logs_validation_local_v1'
if(-not(Test-Path -LiteralPath $RunDir)){throw "Run directory missing: $RunDir"}
if(Test-Path -LiteralPath $Lock){throw "Lock exists: $Lock"}
$expr="invisible(parse(file=tail(commandArgs(trailingOnly=TRUE),1),keep.source=TRUE));cat('R_PARSE=PASS\n')"
$parse=& $Rscript --vanilla -e $expr --args $RFile 2>&1
if($LASTEXITCODE -ne 0){throw ($parse -join [Environment]::NewLine)}
$launch=[DateTimeOffset]::Now.ToString('yyyyMMdd_HHmmss_fff')
$stdout=Join-Path $LogDir "phase6_process_cluster_${RunId}_resume_${launch}.stdout.log"
$stderr=Join-Path $LogDir "phase6_process_cluster_${RunId}_resume_${launch}.stderr.log"
$record=Join-Path $LogDir "phase6_process_cluster_${RunId}_resume_${launch}.json"
Set-Content -LiteralPath $Lock -Value $RunId -Encoding UTF8
$env:PHASE6_RUN_ID=$RunId
$childPath=[Environment]::GetEnvironmentVariable('Path','Process')
[Environment]::SetEnvironmentVariable('PATH',$null,'Process')
[Environment]::SetEnvironmentVariable('Path',$childPath,'Process')
$arguments=@('--vanilla',('"'+$RFile+'"'),('"'+$Root+'"'))
try{$proc=Start-Process -FilePath $Rscript -ArgumentList $arguments -PassThru -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr}catch{if(Test-Path $Lock){Remove-Item -LiteralPath $Lock -Force};throw}
[PSCustomObject]@{run_id=$RunId;launch_id=$launch;pid=$proc.Id;stdout=$stdout;stderr=$stderr;run_dir=$RunDir;resumed_at=[DateTimeOffset]::Now.ToString('o')}|ConvertTo-Json|Set-Content -LiteralPath $record -Encoding UTF8
Write-Output "RUN_ID=$RunId"
Write-Output "PID=$($proc.Id)"
Write-Output "STDOUT=$stdout"
Write-Output "STDERR=$stderr"
