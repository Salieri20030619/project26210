$ErrorActionPreference = 'Stop'
$DestinationRoot = Join-Path 'C:\Users\31937\Downloads' 'project26210_spatial_minimal_tables'
New-Item -ItemType Directory -Force -Path $DestinationRoot | Out-Null
$Files = @(
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/customer_style/tables/GLMF1_GSVA_scores.csv'
        Relative = 'results\part4_spatial\customer_style\tables\GLMF1_GSVA_scores.csv'
        Bytes = [int64]321168
        SHA256 = '5ff9255423742ee0b6dc3d31952e82029f314fcf350a62eeadec09ef7ea1b7ae'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/advanced_patch_misty_glmf1_v1/GLMF1_misty_contribution_stats_v1.tsv'
        Relative = 'results\part4_spatial\advanced_patch_misty_glmf1_v1\GLMF1_misty_contribution_stats_v1.tsv'
        Bytes = [int64]4657
        SHA256 = 'a44e30f67e0fbeb3c9ad9c28ad415ddcec6ec8851cb86968c163f277796e3e87'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/advanced_patch_misty_glmf1_v1/GLMF1_misty_improvement_stats_v1.tsv'
        Relative = 'results\part4_spatial\advanced_patch_misty_glmf1_v1\GLMF1_misty_improvement_stats_v1.tsv'
        Bytes = [int64]6099
        SHA256 = 'a7977ce423f1be7287f15d7eb4a7bd3db1f878ffb7c9cdc5758380be99f09699'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/advanced_patch_spatial_programs_glmf1_v1/GLMF1_hotspot_local_correlation_z.tsv'
        Relative = 'results\part4_spatial\advanced_patch_spatial_programs_glmf1_v1\GLMF1_hotspot_local_correlation_z.tsv'
        Bytes = [int64]4749849
        SHA256 = '7570dff4af9eadf75ce986596448ed655339449c85d60f68b71ec5f0a6e96971'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/redrawn_figures/tables/GLMF1_abundance_correlations.csv'
        Relative = 'results\part4_spatial\redrawn_figures\tables\GLMF1_abundance_correlations.csv'
        Bytes = [int64]2153
        SHA256 = '3747f7a397922a5670564066065a2f1842145caba74210f8614819e68df7d791'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/redrawn_figures/tables/GLMF1_colocalization_correlations.csv'
        Relative = 'results\part4_spatial\redrawn_figures\tables\GLMF1_colocalization_correlations.csv'
        Bytes = [int64]1059
        SHA256 = '9d8aacb699a9b37b1008523e46986368d22efceb2f65f177961c65ef21d2e20a'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/advanced_patch_bayesspace_glmf1_v1/GLMF1_boundary_distance_profiles.tsv'
        Relative = 'results\part4_spatial\advanced_patch_bayesspace_glmf1_v1\GLMF1_boundary_distance_profiles.tsv'
        Bytes = [int64]2153
        SHA256 = 'd827d537199522f6f74aec199b0c741a788f3956e04c1d8077a703c8d3ef916d'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/customer_style/tables/GLMF2_GSVA_scores.csv'
        Relative = 'results\part4_spatial\customer_style\tables\GLMF2_GSVA_scores.csv'
        Bytes = [int64]390637
        SHA256 = 'd867d4856410919659d1512af2d93a0955a1e7acc490c0367db5118d4e655bcf'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/advanced_patch_misty_glmf2_v1/GLMF2_misty_contribution_stats_v1.tsv'
        Relative = 'results\part4_spatial\advanced_patch_misty_glmf2_v1\GLMF2_misty_contribution_stats_v1.tsv'
        Bytes = [int64]4740
        SHA256 = '25de7b901407925199f7d226b1daff86e8c7d8d8dc3daa9bfddd5bd02af97364'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/advanced_patch_misty_glmf2_v1/GLMF2_misty_improvement_stats_v1.tsv'
        Relative = 'results\part4_spatial\advanced_patch_misty_glmf2_v1\GLMF2_misty_improvement_stats_v1.tsv'
        Bytes = [int64]6099
        SHA256 = 'f37b9106f779f37e01b1f0850b76fca61476bcc69dfedb821f695ec77ea1cbc2'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/advanced_patch_spatial_programs_glmf2_v2/GLMF2_hotspot_local_correlation_z.tsv'
        Relative = 'results\part4_spatial\advanced_patch_spatial_programs_glmf2_v2\GLMF2_hotspot_local_correlation_z.tsv'
        Bytes = [int64]4751263
        SHA256 = '0325ce61911d42d521c7310c33dc4d94e976a907852a09feaa607138b6170901'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/redrawn_figures/tables/GLMF2_abundance_correlations.csv'
        Relative = 'results\part4_spatial\redrawn_figures\tables\GLMF2_abundance_correlations.csv'
        Bytes = [int64]2173
        SHA256 = 'c2fd26d1359ed1394cacf01f8c205c6f76f4bc189e34fba474d26691382761a7'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/redrawn_figures/tables/GLMF2_colocalization_correlations.csv'
        Relative = 'results\part4_spatial\redrawn_figures\tables\GLMF2_colocalization_correlations.csv'
        Bytes = [int64]1066
        SHA256 = '9c89bdbbfcfec7b4372cb6069562163c7b0b051b1ad7a6c466ee3e6ed6f3db99'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/advanced_patch_bayesspace_glmf2_v1/GLMF2_boundary_distance_profiles.tsv'
        Relative = 'results\part4_spatial\advanced_patch_bayesspace_glmf2_v1\GLMF2_boundary_distance_profiles.tsv'
        Bytes = [int64]2188
        SHA256 = '1191b254283b2b288dcde9ec496dc0a5ed955bc39f274451ee4d45225e3e4030'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/advanced_patch_misty_qc_v1/verification_summary_v1.json'
        Relative = 'results\part4_spatial\advanced_patch_misty_qc_v1\verification_summary_v1.json'
        Bytes = [int64]491
        SHA256 = 'dc6435cad9c019dca127d9b65102b57dd7406a9409431ba58d5412fb23d8569a'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/pathology_annotation/annotation_assets/annotation_manifest.json'
        Relative = 'results\part4_spatial\pathology_annotation\annotation_assets\annotation_manifest.json'
        Bytes = [int64]1388
        SHA256 = '9a9100fa497013aad8296d9309c2df83b25b84dc60a2b8974549bd05bedf2c1d'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/customer_style/tables/spatial_lag_TYMPHigh_association.csv'
        Relative = 'results\part4_spatial\customer_style\tables\spatial_lag_TYMPHigh_association.csv'
        Bytes = [int64]387
        SHA256 = '1e2425f524b320279a47143048f715a07f1365f989df4766d81cb441de385fed'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/redrawn_figures/tables/TYMP_TAM_colocalization_correlations.csv'
        Relative = 'results\part4_spatial\redrawn_figures\tables\TYMP_TAM_colocalization_correlations.csv'
        Bytes = [int64]2056
        SHA256 = '5f8174e9a3f34115ecfcb52fc4989636945dacdb2c778bdc81c3d99f43be75f7'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/spatial_restricted_validation_v1/spatial_restricted_relation_tests.csv'
        Relative = 'results\part4_spatial\spatial_restricted_validation_v1\spatial_restricted_relation_tests.csv'
        Bytes = [int64]3145
        SHA256 = '32e7b818414c092a4b7fde758524637b08143283916f7ccb77cfe9bda2b06fd7'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/stats/TYMP_TAM_colocalization_correlations.csv'
        Relative = 'results\part4_spatial\stats\TYMP_TAM_colocalization_correlations.csv'
        Bytes = [int64]3041
        SHA256 = '48147e45a2c55db45bc15b7583cfae4657c0b84ef8801e3a23eccbc33d1a1b98'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/tymp_myeloid_state_ordering_v1/tables/sample_pseudotime_associations_v1.tsv'
        Relative = 'results\part4_spatial\tymp_myeloid_state_ordering_v1\tables\sample_pseudotime_associations_v1.tsv'
        Bytes = [int64]44349
        SHA256 = '6f58ee3650a42cc8f59d253dbe135b81c2bd1048b7630e599943b706c750717c'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/tymp_spp1_state_v1/tables/patient_state_associations_v1.tsv'
        Relative = 'results\part4_spatial\tymp_spp1_state_v1\tables\patient_state_associations_v1.tsv'
        Bytes = [int64]90100
        SHA256 = '582f5c70563c80d85816781d6de7391d9514ed3b7a114945b4030e844603d55b'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/customer_style/tables/tumor_NAT_boundary_distance_template.csv'
        Relative = 'results\part4_spatial\customer_style\tables\tumor_NAT_boundary_distance_template.csv'
        Bytes = [int64]184088
        SHA256 = '4c3c4d2bb931487b3a8dd704a60b5ce3babc53ea9ebd761169411df30a726603'
    }
    [pscustomobject]@{
        Remote = '/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/tymp_vascular_niche_v1/tables/boundary_distance_profiles_v1.tsv'
        Relative = 'results\part4_spatial\tymp_vascular_niche_v1\tables\boundary_distance_profiles_v1.tsv'
        Bytes = [int64]55737
        SHA256 = '99e131d874e388b9a1bc58191944b6bbfe5ae77cea7cdf83aef905fd9e0fa330'
    }
)

Write-Host 'project26210 supplemental minimal spatial tables only'
Write-Host ('Files: {0}' -f $Files.Count)
Write-Host ('Bytes: {0}' -f (($Files | Measure-Object -Property Bytes -Sum).Sum))

foreach ($File in $Files) {
    $LocalPath = Join-Path $DestinationRoot $File.Relative
    $LocalDir = Split-Path -Parent $LocalPath
    New-Item -ItemType Directory -Force -Path $LocalDir | Out-Null
    Write-Host ('Downloading: {0}' -f $File.Remote)
    & scp -- ('sgs10:' + $File.Remote) $LocalPath
    if ($LASTEXITCODE -ne 0) { throw ('scp failed: {0}' -f $File.Remote) }
    $ActualBytes = (Get-Item -LiteralPath $LocalPath).Length
    if ($ActualBytes -ne $File.Bytes) { throw ('size mismatch: {0}' -f $LocalPath) }
    $ActualHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $LocalPath).Hash.ToLowerInvariant()
    if ($ActualHash -ne $File.SHA256.ToLowerInvariant()) { throw ('SHA256 mismatch: {0}' -f $LocalPath) }
}
Write-Host 'DOWNLOAD_AND_SHA256_CHECK=PASS'
