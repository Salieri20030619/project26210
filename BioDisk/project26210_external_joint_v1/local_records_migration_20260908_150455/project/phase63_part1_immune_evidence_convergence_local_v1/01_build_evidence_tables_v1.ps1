$ErrorActionPreference = 'Stop'

$phaseDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$p60 = 'D:\project26210_external_joint_v1\phase60_part1_immunotherapy_f1_tymp_local_v1\results_v1'
$p61 = 'C:\Users\31937\_codex_phase61_part1_detailed_immune_local_v1\results_v2'
$p62 = 'C:\Users\31937\_codex_phase62_part1_immune_pseudobulk_local_v1\results_v5'

function Write-Tsv($rows, $path) {
    $rows | Export-Csv -LiteralPath $path -Delimiter "`t" -NoTypeInformation -Encoding utf8
}

$rules = @(
    [pscustomobject]@{grade='A'; label_cn='主结论/较强证据'; criteria='患者或配对为统计单位；多重校正通过；样本量基本可用；严格定义或独立队列/方法复现'; reporting='可进入主结果，但仍需写明队列和设计边界'},
    [pscustomobject]@{grade='B'; label_cn='有限支持'; criteria='多重校正通过，但样本量小、仅单队列、仅扩展定义或缺少复现；或多个分析方向一致'; reporting='可作为支持性结果，不称为已验证标志物'},
    [pscustomobject]@{grade='C1'; label_cn='FDR支持的探索性机制'; criteria='通路/机制层面多重校正通过，但模型样本小、组间不平衡或仅单队列'; reporting='报告方向、P、FDR与样本量；用于提出机制假设'},
    [pscustomobject]@{grade='C2'; label_cn='名义探索性信号'; criteria='P<0.05但FDR未通过；或列入每模型前15名义显著通路'; reporting='必须标为名义信号，不用于确认性结论'},
    [pscustomobject]@{grade='D'; label_cn='当前未支持'; criteria='分析可检验但多重校正未通过，且无足够的跨方法/跨队列支持'; reporting='写作“当前数据未支持”，不写“证明不存在”'},
    [pscustomobject]@{grade='U'; label_cn='不可判定'; criteria='样本量不足、缺少响应标签、无合格患者或无法执行推断'; reporting='仅描述缺口，不解释为阴性'}
)
Write-Tsv $rules (Join-Path $phaseDir 'IMMUNE_EVIDENCE_GRADING_RULES_v1.tsv')

$evidence = @(
    [pscustomobject]@{evidence_id='E01'; domain='F1_TYMP_vs_ICB_response'; dataset='GSE294109'; analysis='malignant F1–TYMP rho, benefit vs no-effect'; unit='patient'; n='expanded 5 vs 2; strict 4 vs 1'; effect='median difference benefit-minus-no-effect < 0 in both tiers'; p_value='0.5613 expanded; 0.7237 strict'; fdr='NA'; grade='D'; conclusion_cn='未观察到F1–TYMP与临床获益的支持性差异'; limitation='有效恶性层患者远少于完整9 vs 3响应队列'; source=(Join-Path $p60 'GSE294109_icb_response_comparison_v1.tsv')},
    [pscustomobject]@{evidence_id='E02'; domain='immune_features_vs_ICB_response'; dataset='GSE294109'; analysis='27 immune composition/state features, benefit vs no-effect'; unit='patient after region pooling'; n='9 vs 3'; effect='T/NK IFN response nominally lower in benefit'; p_value='min p=0.04196'; fdr='min q=0.46158'; grade='D'; conclusion_cn='没有患者层面免疫组成或状态特征通过家族内FDR'; limitation='no-effect仅3例，功效有限'; source=(Join-Path $p61 'GSE294109_icb_immune_response_comparison_v1.tsv')},
    [pscustomobject]@{evidence_id='E03'; domain='longitudinal_immune_features'; dataset='GSE264586'; analysis='54 paired immune composition/state tests'; unit='paired patient'; n='Cycle2 n=4 pairs; EOT n=3 pairs'; effect='immune fraction and broad composition showed nominal trends'; p_value='min p=0.10881'; fdr='min q=0.10881'; grade='D'; conclusion_cn='纵向免疫组成/状态未形成经FDR支持的变化'; limitation='仅3–4对患者'; source=(Join-Path $p61 'GSE264586_immune_longitudinal_summary_v1.tsv')},
    [pscustomobject]@{evidence_id='E04'; domain='F1_TYMP_vs_immune_composition'; dataset='GSE294109'; analysis='expanded malignant F1–TYMP rho vs Treg fraction'; unit='patient'; n='7'; effect='Spearman rho=0.92857'; p_value='0.002519'; fdr='q=0.027714'; grade='B'; conclusion_cn='扩展恶性定义下，F1–TYMP与Treg占比呈强正相关'; limitation='仅扩展层显著；严格层n=5、rho=0.9但q=0.411，未复现FDR'; source=(Join-Path $p61 'immune_feature_vs_f1_tymp_rho_associations_v1.tsv')},
    [pscustomobject]@{evidence_id='E05'; domain='F1_TYMP_vs_immune_state'; dataset='GSE294109'; analysis='strict malignant F1–TYMP rho vs T/NK dysfunction'; unit='patient'; n='5'; effect='Spearman rho=0.9'; p_value='0.03739'; fdr='q=0.4112'; grade='C2'; conclusion_cn='存在名义正相关，但未通过多重校正'; limitation='极小样本且无独立复现'; source=(Join-Path $p61 'immune_feature_vs_f1_tymp_rho_associations_v1.tsv')},
    [pscustomobject]@{evidence_id='E06'; domain='longitudinal_F1_TYMP'; dataset='GSE264586'; analysis='F1–TYMP longitudinal paired change'; unit='paired patient'; n='shown comparisons n_pairs=1'; effect='direction only'; p_value='NA'; fdr='NA'; grade='U'; conclusion_cn='无法进行纵向推断'; limitation='合格配对患者不足'; source=(Join-Path $p60 'GSE264586_longitudinal_delta_rho_v1.tsv')},
    [pscustomobject]@{evidence_id='E07'; domain='ICB_response_F1_TYMP'; dataset='GSE264586'; analysis='response association'; unit='patient'; n='NA'; effect='NA'; p_value='NA'; fdr='NA'; grade='U'; conclusion_cn='该队列无可用响应分析'; limitation='缺少适用的响应终点/映射'; source=(Join-Path $p60 'clinical_mapping_audit_v1.tsv')},
    [pscustomobject]@{evidence_id='E08'; domain='pseudobulk_gene_level'; dataset='GSE264586 + GSE294109'; analysis='12 edgeR QL pseudobulk models'; unit='paired patient or pooled patient'; n='3–4 pairs or 9–12 patients/model'; effect='8 gene-level FDR hits across 12 models'; p_value='model-specific'; fdr='q<0.05'; grade='C1'; conclusion_cn='基因层面差异稀疏，未形成稳定的单基因标志物集合'; limitation='小样本、组间不平衡、模型多'; source=(Join-Path $p62 'pseudobulk_model_summary_v1.tsv')},
    [pscustomobject]@{evidence_id='E09'; domain='pseudobulk_pathways'; dataset='GSE264586'; analysis='paired longitudinal Hallmark enrichment'; unit='paired patient pseudobulk'; n='3–4 pairs/model'; effect='TNFA/NFKB、IFN-gamma、mTORC1及allograft rejection等呈时间/细胞型特异方向'; p_value='see pathway register'; fdr='2–38 Hallmark q<0.05/model'; grade='C1'; conclusion_cn='治疗过程中存在广泛但高度情境依赖的免疫通路重塑'; limitation='仅3–4对，方向不能直接外推为疗效机制'; source=(Join-Path $p62 'pseudobulk_hallmark_fgsea_v1.tsv')},
    [pscustomobject]@{evidence_id='E10'; domain='pseudobulk_pathways'; dataset='GSE294109'; analysis='benefit vs no-effect Hallmark enrichment'; unit='patient pseudobulk after region pooling'; n='9–12 total, no-effect=3'; effect='T/NK氧化磷酸化降低、CD4/单核IFN-gamma降低、C1QC巨噬EMT升高等'; p_value='see pathway register'; fdr='7–23 Hallmark q<0.05/model'; grade='C1'; conclusion_cn='获益相关免疫状态更像多细胞型通路重排，而非稳定单基因差异'; limitation='no-effect仅3例；单队列；不能作预测标志物验证'; source=(Join-Path $p62 'pseudobulk_hallmark_fgsea_v1.tsv')}
)
Write-Tsv $evidence (Join-Path $phaseDir 'IMMUNE_EVIDENCE_REGISTER_v1.tsv')

$themeMap = @{
    'HALLMARK_INTERFERON_ALPHA_RESPONSE'='immune_interferon'
    'HALLMARK_INTERFERON_GAMMA_RESPONSE'='immune_interferon'
    'HALLMARK_TNFA_SIGNALING_VIA_NFKB'='immune_inflammation'
    'HALLMARK_INFLAMMATORY_RESPONSE'='immune_inflammation'
    'HALLMARK_IL2_STAT5_SIGNALING'='immune_cytokine'
    'HALLMARK_IL6_JAK_STAT3_SIGNALING'='immune_cytokine'
    'HALLMARK_COMPLEMENT'='immune_innate'
    'HALLMARK_ALLOGRAFT_REJECTION'='immune_activation'
    'HALLMARK_TGF_BETA_SIGNALING'='immune_suppression_tme'
    'HALLMARK_APOPTOSIS'='cell_stress_death'
    'HALLMARK_HYPOXIA'='tymp_tme'
    'HALLMARK_ANGIOGENESIS'='tymp_tme'
    'HALLMARK_COAGULATION'='tymp_tme'
    'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'='tumor_microenvironment'
    'HALLMARK_OXIDATIVE_PHOSPHORYLATION'='immune_metabolism'
    'HALLMARK_GLYCOLYSIS'='immune_metabolism'
    'HALLMARK_MTORC1_SIGNALING'='immune_metabolism'
    'HALLMARK_PI3K_AKT_MTOR_SIGNALING'='immune_metabolism'
    'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY'='immune_metabolism'
    'HALLMARK_FATTY_ACID_METABOLISM'='immune_metabolism'
    'HALLMARK_CHOLESTEROL_HOMEOSTASIS'='immune_metabolism'
    'HALLMARK_P53_PATHWAY'='tumor_stress'
}

$models = Import-Csv -Delimiter "`t" -LiteralPath (Join-Path $p62 'pseudobulk_model_summary_v1.tsv')
$modelByKey = @{}
foreach ($m in $models) {
    $key = "$($m.dataset)|$($m.contrast)|$($m.fine_annotation)"
    $modelByKey[$key] = $m
}

$fgsea = Import-Csv -Delimiter "`t" -LiteralPath (Join-Path $p62 'pseudobulk_hallmark_fgsea_v1.tsv')
$selected = New-Object System.Collections.Generic.List[object]
foreach ($g in ($fgsea | Group-Object dataset,contrast,fine_annotation)) {
    $valid = @($g.Group | Where-Object { $_.pval -ne 'NA' -and [double]$_.pval -lt 0.05 } | Sort-Object {[double]$_.pval})
    $rankByPath = @{}
    for ($i=0; $i -lt $valid.Count; $i++) { $rankByPath[$valid[$i].pathway] = $i + 1 }
    $topPaths = @($valid | Select-Object -First 15 | ForEach-Object {$_.pathway})
    foreach ($row in $valid) {
        $isTop = $topPaths -contains $row.pathway
        $theme = if ($themeMap.ContainsKey($row.pathway)) {$themeMap[$row.pathway]} else {'not_predefined_theme'}
        $isTheme = $theme -ne 'not_predefined_theme'
        if (-not ($isTop -or $isTheme)) { continue }
        $key = "$($row.dataset)|$($row.contrast)|$($row.fine_annotation)"
        $m = $modelByKey[$key]
        $padjNum = if ($row.padj -eq 'NA') {$null} else {[double]$row.padj}
        $reasons = @()
        if ($isTop) {$reasons += 'TOP15_NOMINAL_P'}
        if ($isTheme) {$reasons += 'THEME_RELATED_P_LT_0_05'}
        $selected.Add([pscustomobject]@{
            dataset=$row.dataset
            contrast=$row.contrast
            fine_annotation=$row.fine_annotation
            effect_definition=$row.effect_definition
            n_units=$m.n_units
            n_benefit=$m.n_benefit
            n_no_effect=$m.n_no_effect
            n_pairs=$m.n_pairs
            pathway=$row.pathway
            theme_category=$theme
            nominal_p_rank_within_model=$rankByPath[$row.pathway]
            p_value=$row.pval
            p_fdr_bh=$row.padj
            NES=$row.NES
            direction=if ([double]$row.NES -gt 0) {'positive_in_effect_direction'} else {'negative_in_effect_direction'}
            retention_reason=($reasons -join '|')
            significance_class=if ($null -ne $padjNum -and $padjNum -lt 0.05) {'FDR_LT_0_05'} else {'NOMINAL_P_LT_0_05_ONLY'}
            evidence_grade=if ($null -ne $padjNum -and $padjNum -lt 0.05) {'C1'} else {'C2'}
            source=(Join-Path $p62 'pseudobulk_hallmark_fgsea_v1.tsv')
        })
    }
}
$selected = @($selected | Sort-Object dataset,contrast,fine_annotation,{[int]$_.nominal_p_rank_within_model})
Write-Tsv $selected (Join-Path $phaseDir 'IMMUNE_PATHWAY_REPORTING_v1.tsv')

$index = foreach ($g in ($selected | Group-Object dataset,contrast,fine_annotation)) {
    $rows = @($g.Group)
    $first = $rows[0]
    [pscustomobject]@{
        dataset=$first.dataset
        contrast=$first.contrast
        fine_annotation=$first.fine_annotation
        n_units=$first.n_units
        n_benefit=$first.n_benefit
        n_no_effect=$first.n_no_effect
        n_pairs=$first.n_pairs
        retained_total=$rows.Count
        retained_top15=[int](@($rows | Where-Object {$_.retention_reason -match 'TOP15_NOMINAL_P'}).Count)
        retained_theme_p_lt_0_05=[int](@($rows | Where-Object {$_.retention_reason -match 'THEME_RELATED_P_LT_0_05'}).Count)
        retained_fdr_lt_0_05=[int](@($rows | Where-Object {$_.significance_class -eq 'FDR_LT_0_05'}).Count)
        retained_nominal_only=[int](@($rows | Where-Object {$_.significance_class -eq 'NOMINAL_P_LT_0_05_ONLY'}).Count)
    }
}
Write-Tsv $index (Join-Path $phaseDir 'IMMUNE_PATHWAY_MODEL_INDEX_v1.tsv')

$status = @(
    [pscustomobject]@{check='phase60_source_exists'; status=$(if (Test-Path $p60) {'PASS'} else {'FAIL'}); detail=$p60},
    [pscustomobject]@{check='phase61_source_exists'; status=$(if (Test-Path $p61) {'PASS'} else {'FAIL'}); detail=$p61},
    [pscustomobject]@{check='phase62_source_exists'; status=$(if (Test-Path $p62) {'PASS'} else {'FAIL'}); detail=$p62},
    [pscustomobject]@{check='evidence_rows'; status=$(if ($evidence.Count -ge 10) {'PASS'} else {'FAIL'}); detail=$evidence.Count},
    [pscustomobject]@{check='pathway_models'; status=$(if (@($index).Count -eq 12) {'PASS'} else {'FAIL'}); detail=@($index).Count},
    [pscustomobject]@{check='all_retained_pathways_nominal_p_lt_0_05'; status=$(if (@($selected | Where-Object {[double]$_.p_value -ge 0.05}).Count -eq 0) {'PASS'} else {'FAIL'}); detail='union(top15 nominal p<0.05, theme-related p<0.05)'},
    [pscustomobject]@{check='grading_preserves_nominal_only'; status=$(if (@($selected | Where-Object {$_.significance_class -eq 'NOMINAL_P_LT_0_05_ONLY'}).Count -gt 0) {'PASS'} else {'FAIL'}); detail=@($selected | Where-Object {$_.significance_class -eq 'NOMINAL_P_LT_0_05_ONLY'}).Count},
    [pscustomobject]@{check='workbook'; status='BLOCKED_RUNTIME'; detail='artifact-tool controlled runtime unavailable in current Windows environment; TSV/Markdown delivered'}
)
Write-Tsv $status (Join-Path $phaseDir 'PHASE63_VALIDATION_v1.tsv')

Write-Output "Evidence rows: $($evidence.Count)"
Write-Output "Pathway rows retained: $($selected.Count)"
Write-Output "Models covered: $(@($index).Count)"
