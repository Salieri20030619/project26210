[CmdletBinding()]
param([Parameter(Mandatory=$true)][string]$RunId,[string]$Root=(Split-Path -Parent $PSScriptRoot))
$Root=(Resolve-Path -LiteralPath $Root).Path
$runDir=Join-Path $Root "phase6_processing_local_v1\runs\$RunId"
$record=Get-Content -LiteralPath (Join-Path $Root "logs_validation_local_v1\phase6_process_cluster_$RunId.json") -Raw -Encoding UTF8|ConvertFrom-Json
$proc=Get-Process -Id $record.pid -ErrorAction SilentlyContinue
[PSCustomObject]@{run_id=$RunId;pid=$record.pid;process_exists=[bool]$proc;run_dir=$runDir;status_exists=(Test-Path (Join-Path $runDir 'PHASE6_PROCESSING_STATUS'))}|Format-List
$progress=Join-Path $runDir 'PHASE6_PROCESSING_PROGRESS.tsv'
if(Test-Path $progress){Get-Content -LiteralPath $progress -Encoding UTF8}
if(Test-Path $record.stderr){Get-Content -LiteralPath $record.stderr -Tail 20 -Encoding UTF8}
