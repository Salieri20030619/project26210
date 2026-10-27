#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results_v1'
OUTPUT = RESULTS / 'PHASE66_INTERPRETATION_v1.md'
if OUTPUT.exists():
    raise SystemExit('Refusing to overwrite Phase66 interpretation')


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


rows = read_tsv(RESULTS / 'TRANSFER_GSEA_FULL_v1.tsv')
fdr = [row for row in rows if float(row['transfer_q_value']) < 0.05]
family_counts = Counter(row['family'] for row in fdr)


def selected(dataset: str, annotation: str, family: str) -> list[dict[str, str]]:
    values = [
        row for row in fdr
        if row['dataset'] == dataset and row['fine_annotation'] == annotation and row['family'] == family
    ]
    return sorted(values, key=lambda row: (float(row['transfer_q_value']), -abs(float(row['NES']))))


def compact(values: list[dict[str, str]]) -> str:
    if not values:
        return '无'
    return '；'.join(
        f"{row['label']}（NES={float(row['NES']):.3f}, q={float(row['transfer_q_value']):.3g}）"
        for row in values
    )


g264_mac_reg = selected('GSE264586', 'Macrophage_C1QC', 'regulon')
g264_mono_reg = selected('GSE264586', 'Monocyte_FCN1', 'regulon')
g264_mac_met = selected('GSE264586', 'Macrophage_C1QC', 'metabolic')
g294_mac_reg = selected('GSE294109', 'Macrophage_C1QC', 'regulon')
g294_mono_reg = selected('GSE294109', 'Monocyte_FCN1', 'regulon')
g294_mono_met = selected('GSE294109', 'Monocyte_FCN1', 'metabolic')
g294_treg_reg = selected('GSE294109', 'Treg', 'regulon')
g294_treg_met = selected('GSE294109', 'Treg', 'metabolic')

text = f"""# Phase66：免疫程序独立临床情境转移结论

## 一、总体结果

50 个冻结签名在 5 个独立患者级模型中形成 146 个预注册组合；123 个通过 10–500 基因覆盖门禁并进入 fgsea。共有 42 个家族内 FDR 转移结果，其中 regulon {family_counts['regulon']} 个、代谢通路 {family_counts['metabolic']} 个；另有 10 个仅名义显著结果。

结论不是“发现统一获益程序”，而是：Phase64 的 TYMP 相关髓系/Treg 程序可以在独立免疫治疗队列中被检测到，但其富集方向随治疗时间、细胞亚群和临床获益背景改变。

## 二、GSE264586 纵向治疗转移

- End-of-treatment Macrophage_C1QC：14 个可检验 regulon 全部正向 FDR 富集，核心包括 IRF1、NFKB1、RELB、STAT2、RELA、STAT6 和 CEBPB。完整结果：{compact(g264_mac_reg)}。
- End-of-treatment Monocyte_FCN1：9 个 regulon 正向 FDR 富集，包括 IRF1、NFKB1、RELB、RELA、STAT6、STAT2、ELF2、IRF7 和 CEBPB。完整结果：{compact(g264_mono_reg)}。
- Macrophage_C1QC 的两个呼吸链代谢签名在治疗末期负向富集：{compact(g264_mac_met)}。Monocyte_FCN1 无代谢签名通过 FDR。
- 该队列每个纵向模型只有 3 对患者，因此属于强方向一致、但样本很小的治疗阶段证据。

## 三、GSE294109 ICB 获益情境转移

- Macrophage_C1QC：IRF1、NFKB1 在获益组负向，而 CEBPB、ZBTB7A 正向，呈混合调控重排：{compact(g294_mac_reg)}。
- Monocyte_FCN1：SPI1、IRF1、IRF7、STAT2 均在获益组负向富集：{compact(g294_mono_reg)}。代谢层为 folate/pterine metabolism 正向、cholesterol biosynthesis 负向：{compact(g294_mono_met)}。
- Treg：呼吸链和 cellular response to hypoxia 在获益组负向富集：{compact(g294_treg_met)}。调控层 STAT1、IRF2 负向，REL、RELB 正向：{compact(g294_treg_reg)}。
- no-effect 组仅 3 名患者；这些信号是患者级程序富集，但不构成可直接应用的预测模型。

## 四、跨队列判断

GSE264586 治疗末期髓系 regulon 主要正向，而 GSE294109 获益相关 Monocyte regulon 主要负向。因此，Phase64 程序具有跨队列可检测性，但不存在统一、方向一致的“疗效获益程序”。这进一步支持 Phase65 的结论：F1/TYMP 更适合作为谱系和状态依赖的生物学轴，而不是单独的免疫治疗预测标志物。

## 五、证据等级与限制

- 可升级：Myeloid/DC‑TYMP regulon 从单队列方法收敛，增加了两个独立免疫治疗队列中的临床情境转移证据；仍应标为情境依赖的机制支持。
- 不升级：F1/TYMP 的直接 response/PFS/OS 预测等级保持不变。
- 本分析没有直接复算外部队列内 TYMP 连续表达与程序活性的相关性；不同临床对比的 NES 方向不能机械当作重复。
- regulon 使用每个冻结 TF 正向 regulon 的前 200 个加权靶基因；代谢签名已排除 TYMP，且所有显著性均按模型与签名家族进行 BH 校正。
"""
OUTPUT.write_text(text, encoding='utf-8')
print(f'WROTE {OUTPUT}')
