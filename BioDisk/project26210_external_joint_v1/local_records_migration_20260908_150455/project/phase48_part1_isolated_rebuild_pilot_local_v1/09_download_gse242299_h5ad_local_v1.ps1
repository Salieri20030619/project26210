$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$phaseDir = $PSScriptRoot
$outputDir = Join-Path $phaseDir 'gse242299_official_h5ad_local_v1'
$finalFile = Join-Path $outputDir 'GSE242299_all_cells_50236_33538.h5ad.gz'
$partFile = "$finalFile.part"
$manifestFile = Join-Path $outputDir 'GSE242299_official_h5ad_download_manifest_local_v1.tsv'
$url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242299/suppl/GSE242299_all_cells_50236_33538.h5ad.gz'
$expectedSize = 5350598000L

if (-not (Test-Path -LiteralPath $outputDir)) {
    New-Item -ItemType Directory -Path $outputDir -ErrorAction Stop | Out-Null
}
if (Test-Path -LiteralPath $finalFile) {
    throw "Refusing existing final file: $finalFile"
}

Write-Output "attempt_started=$((Get-Date).ToString('o'))"
Write-Output "url=$url"
Write-Output "expected_size=$expectedSize"
& curl.exe --fail --location --retry 5 --retry-delay 10 --continue-at - --output $partFile $url
if ($LASTEXITCODE -ne 0) {
    throw "curl failed with exit code $LASTEXITCODE"
}

$observedSize = (Get-Item -LiteralPath $partFile).Length
if ($observedSize -ne $expectedSize) {
    throw "Size mismatch: $observedSize != $expectedSize"
}
$sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $partFile).Hash.ToLowerInvariant()
Move-Item -LiteralPath $partFile -Destination $finalFile -ErrorAction Stop
@(
    "file`tsize_bytes`tsha256`turl`tstatus"
    "$(Split-Path -Leaf $finalFile)`t$observedSize`t$sha256`t$url`tPASS"
) | Set-Content -LiteralPath $manifestFile -Encoding utf8
Write-Output "sha256=$sha256"
Write-Output "final_file=$finalFile"
Write-Output "attempt_finished=$((Get-Date).ToString('o'))"
