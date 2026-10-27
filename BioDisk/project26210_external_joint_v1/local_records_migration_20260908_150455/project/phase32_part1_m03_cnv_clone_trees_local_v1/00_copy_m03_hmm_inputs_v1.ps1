$ErrorActionPreference = 'Stop'

$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$inputRoot = Join-Path $PSScriptRoot 'input_snapshot_local_v1\hmm_profiles'
$sshConfig = 'C:\Users\31937\.ssh\config'
$remoteRoot = '/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/infercnv_v2/patients'
$units = @(
    'GSE156632__RCC2',
    'GSE156632__RCC4',
    'GSE156632__RCC5',
    'GSE159115__SS_2017',
    'GSE159115__SS_2023',
    'GSE207493__RCC100',
    'GSE207493__RCC104',
    'GSE207493__RCC106',
    'GSE207493__RCC112',
    'GSE207493__RCC86',
    'GSE207493__RCC87',
    'GSE207493__RCC94',
    'GSE207493__RCC96',
    'GSE207493__RCC99'
)
$files = @(
    '17_HMM_predHMMi6.leiden.hmm_mode-subclusters.genes_used.dat',
    '17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat'
)

if (Test-Path -LiteralPath $inputRoot) {
    throw "Refusing to overwrite existing HMM snapshot: $inputRoot"
}
New-Item -ItemType Directory -Path $inputRoot | Out-Null

foreach ($unit in $units) {
    $unitDir = Join-Path $inputRoot $unit
    New-Item -ItemType Directory -Path $unitDir | Out-Null
    foreach ($name in $files) {
        $part = Join-Path $unitDir ($name + '.part')
        $final = Join-Path $unitDir $name
        $remote = "sgs10:${remoteRoot}/${unit}/infercnv_run/${name}"
        & scp -F $sshConfig $remote $part
        if ($LASTEXITCODE -ne 0) {
            throw "scp failed for $remote"
        }
        Move-Item -LiteralPath $part -Destination $final
    }
}

Write-Output "M03 HMM input snapshot copied: $($units.Count) patients, $($files.Count * $units.Count) files"
