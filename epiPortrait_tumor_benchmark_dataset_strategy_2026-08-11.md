# epiPortrait 肿瘤优先 Benchmark 数据集与论文设计建议
## 生信专家 / 表观遗传学 / Reviewer 视角

**日期：2026-08-11**  
**适用对象：epiPortrait 方法学 / 软件论文 benchmark 设计**

---

# 1. 总体结论

如果 epiPortrait 最终定位为：

> **Replicate-Aware Epigenomic Domain Profiling**

并且论文希望突出：

```text
肿瘤表观遗传学
+
replicate-aware
+
Intensity × Breadth decomposition
+
condition remodeling
+
active / repressive mark generalization
```

那么 benchmark 数据应优先选择：

1. **明确的肿瘤背景；**
2. **至少两个 biological replicates，最好 ≥3；**
3. **有 H3K4me3 / H3K27ac / H3K27me3 等核心 histone marks；**
4. **存在明确生物学 perturbation / resistance / tumor-vs-reference 对照；**
5. **最好有 RNA-seq、TF ChIP、P300/MED1 或其他 orthogonal evidence；**
6. **能对应 epiPortrait 的核心 claim，而不是单纯“能跑”。**

不建议为了数量堆大量 unrelated datasets。

更推荐：

```text
1 个 flagship multi-mark tumor dataset
+
2 个机制/治疗响应数据集
+
2–3 个补充 validation dataset
```

---

# 2. 推荐数据集总表

| 优先级 | 数据集 | 肿瘤类型 | 关键设计 | 主要 marks / data | 最适合验证 |
|---|---|---|---|---|---|
| **主数据 1** | **GSE259262** | Glioblastoma | 10 GIC vs 10 matched/syngeneic iNSC | H3K4me3, H3K27ac, H3K27me3, H3K36me3 | **全框架 flagship benchmark** |
| **主数据 2** | **GSE136128** | Enzalutamide-resistant prostate cancer | C4-2 Control 3 vs ENZ-R 3 | H3K27ac, AR, FOXA1 | **therapy-resistance domain remodeling** |
| **主数据 3** | **GSE126004** | Breast cancer / metastatic signaling | T47D/MCF7, IL6 ± FOXA1 perturbation | H3K27ac + regulatory TFs | **enhancer-domain remodeling + mechanism** |
| Supplement | **GSE78206** | Multiple myeloma | MM1S DMSO n=5 vs KDM5-C70 n=3 | H3K4me3 | **Intensity vs Breadth perturbation** |
| Supplement / reference | **GSE151556** | MM / T-ALL | oncogene broad-H3K4me3 reference | H3K4me3, H3K27ac, multi-mark | **known broad-domain cancer loci** |
| Supplement / reference | **GSE145938** | Multiple myeloma | MM cell lines + primary MM + normal plasma cell | H3K27ac | **SE/ROSE cancer benchmark** |

---

# 3. Flagship 数据集：GSE259262
## Glioblastoma Initiating Cells vs matched iNSC

这是目前最推荐作为 epiPortrait **主旗舰数据集** 的肿瘤数据。

---

# 4. 为什么 GSE259262 特别适合 epiPortrait

该数据具有几个非常强的特点：

```text
10 个患者来源 GIC
+
10 个对应 iNSC
+
多种 histone marks
+
患者层面的 biological heterogeneity
```

包括：

```text
H3K4me3
H3K27ac
H3K27me3
H3K36me3
```

因此一个 dataset 就能同时支撑：

```text
H3K4me3 Breadth
H3K27ac Intensity/Breadth/Dual
H3K27me3 repressive-domain remodeling
```

---

# 5. 推荐分组方式

建议：

```text
Condition:
iNSC
GIC
```

如果 sample provenance 支持 paired origin，则再增加：

```text
PairID
```

用于连续指标 paired statistics。

---

# 6. 这个数据最适合证明什么

## H3K4me3

回答：

> GBM-associated promoter domains 是否出现 recurrent broadening / narrowing？

主要看：

```text
Breadth-Super
NativeMaxPeakWidth
Breadth expansion/contraction
```

---

## H3K27ac

回答：

> GBM regulatory domains 是否发生 Intensity、Breadth 或 Dual remodeling？

分类：

```text
Typical
Intensity-Super
Breadth-Super
Dual-Super
```

---

## H3K27me3

回答：

> GBM 是否存在 repressive-domain expansion/contraction 或 intensity extremes？

推荐使用：

```text
Extended-Domain
Intensity-Extreme
Dual-Extreme
```

等 mark-aware display terminology。

---

# 7. GSE259262 的一个重要优势：患者异质性

普通 cell-line replicate：

```text
3 replicates
```

更多是在衡量技术/实验重复性。

GSE259262 的：

```text
10 GIC
```

更接近：

> tumor biological heterogeneity

因此 replicate support 可以解释成：

```text
9/10 GIC
→ recurrent GBM-associated phenotype

5/10
→ heterogeneous/intermediate phenotype

2/10
→ patient-specific phenotype
```

这对肿瘤文章非常有价值。

---

# 8. GSE259262 推荐 Figure

## Figure A

```text
GIC vs iNSC
H3K4me3 / H3K27ac / H3K27me3
```

展示：

```text
domain class composition
replicate support
condition transitions
```

---

## Figure B

挑选典型癌症基因 / pathways：

```text
MYC-related
RTK signaling
stemness
mesenchymal programs
```

展示：

```text
browser track
Intensity
Breadth
SignalDispersion
```

但 browser snapshot 只作为说明，不应成为主要 benchmark。

---

# 9. 主数据 2：GSE136128
## Enzalutamide-resistant prostate cancer

这个数据非常适合做：

> **therapy-resistance epigenomic remodeling**

设计非常标准：

```text
C4-2 Control
n = 3

C4-2 ENZ-resistant
n = 3
```

并具有：

```text
H3K27ac
AR
FOXA1
```

---

# 10. 为什么这个数据非常适合 per-group benchmark

它基本就是 epiPortrait 最标准的使用场景：

```text
Control × 3
vs
Resistance × 3
```

可以直接验证：

```text
per-replicate calling
→ support
→ group phenotype
→ condition transition
```

---

# 11. 推荐重点 transition

重点分析：

```text
Typical
→ Intensity-Super

Typical
→ Breadth-Super

Intensity-Super
→ Dual-Super

Breadth-Super
→ Dual-Super

Dual-Super
→ Typical
```

。

---

# 12. 这个数据的 orthogonal validation 很强

同一研究还有：

```text
AR
FOXA1
```

因此可以问：

```text
新获得的 H3K27ac Intensity/Breadth/Dual domains
```

是否同时出现：

```text
AR occupancy gain
FOXA1 occupancy gain
```

。

这比单纯：

```text
H3K27ac 上升
```

更有机制解释力。

---

# 13. 论文中可以提出的 biological question

> Enzalutamide resistance 是否伴随不同类型的 active regulatory-domain remodeling，而不是单一 H3K27ac gain？

这是非常适合 epiPortrait 的故事。

---

# 14. 主数据 3：GSE126004
## Breast cancer IL6/STAT3/FOXA1 enhancer remodeling

这个数据适合做：

> **metastatic signaling-driven enhancer-domain remodeling**

比普通：

```text
MCF7 control vs treatment
```

更有癌症机制价值。

---

# 15. 推荐第一层分析

先只做：

```text
T47D
siControl + PBS
vs
siControl + IL6
```

分析：

```text
IL6-induced H3K27ac remodeling
```

---

# 16. 第二层机制验证

再加入：

```text
siFOXA1 + PBS
siFOXA1 + IL6
```

问：

> 哪些 IL6-induced domain transitions 依赖 FOXA1？

这样形成一个非常漂亮的：

```text
signal
×
chromatin domain phenotype
×
TF dependency
```

故事。

---

# 17. 推荐分析输出

例如：

```text
IL6 gained Intensity-Super
IL6 gained Breadth-Super
IL6 gained Dual-Super
```

然后比较：

```text
FOXA1-dependent
vs
FOXA1-independent
```

。

---

# 18. GSE126004 的优势

同 series 具有：

```text
STAT3 / pSTAT3
ER
FOXA1
```

等 regulatory factors。

因此可以做：

```text
epiPortrait domain phenotype
→ orthogonal TF occupancy enrichment
```

而不是只看 expression。

---

# 19. H3K4me3 专项：GSE78206
## MM1S + KDM5 inhibition

这个数据非常适合验证：

> **Intensity gain 不等于 Breadth expansion**

设计：

```text
MM1S DMSO
n = 5

MM1S KDM5-C70
n = 3
```

主要 mark：

```text
H3K4me3
```

。

---

# 20. 为什么 GSE78206 对 epiPortrait 很有意义

KDM5 inhibition 会增加 H3K4me3。

但是 epiPortrait 可以进一步区分：

```text
signal 增加
```

与：

```text
peak/domain 变宽
```

是否同步。

例如：

```text
Domain A:
Intensity ↑↑
PeakWidth ≈
→ Intensity-Super gain
```

而：

```text
Domain B:
Intensity ↑
PeakWidth ↑↑
→ Dual-Super / Breadth expansion
```

。

这正好证明：

> epiPortrait 的二维 decomposition 不是冗余。

---

# 21. GSE78206 推荐作为 Supplement

不一定需要放主正文。

但非常适合：

```text
Supplementary Figure:
H3K4me3 perturbation
```

证明：

```text
Intensity
≠
Breadth
```

。

---

# 22. Known-positive broad-domain reference：GSE151556

这个数据不一定是最佳 replicate benchmark，

但是它的最大价值在于：

> 已有明确 broad-H3K4me3 oncogene biology。

研究中涉及：

```text
CCND1
FGFR3
MAF
MYC
TAL1
LMO2
TLX3
```

等癌症相关 broad H3K4me3 / super-enhancer context。

---

# 23. 如何正确使用 GSE151556

不要用它证明：

```text
replicate support accuracy
```

。

更适合作为：

```text
known-positive loci validation
```

。

问：

```text
known oncogene-associated broad domains
```

是否被 epiPortrait 识别为：

```text
Breadth-Super
Breadth-expanded
Dual-Super
```

。

---

# 24. 为什么这种 validation 有价值

相比：

```text
MACS2 top 5% broadest peaks
```

known-positive cancer loci 提供的是：

```text
independent biological plausibility
```

。

因此很适合补强 Breadth 生物学故事。

---

# 25. H3K27ac / Super-enhancer cancer reference：GSE145938
## Multiple myeloma

这个数据适合：

```text
ROSE / super-enhancer benchmark
```

。

包含：

```text
MM cell lines
primary MM
normal plasma-cell reference
```

以及 H3K27ac profiles。

---

# 26. GSE145938 不适合做什么

不适合作为：

```text
replicate-aware primary benchmark
```

。

因为很多 primary samples 本质是：

```text
independent patients
```

而不是同一个 condition 的技术/生物重复。

---

# 27. 它最适合做什么

适合回答：

> ROSE-defined / H3K27ac-high cancer regulatory regions，在 epiPortrait 中分别属于什么 phenotype？

例如：

```text
ROSE SE
│
├─ Intensity-Super
├─ Breadth-Super
└─ Dual-Super
```

。

---

# 28. 推荐的论文数据组织

## Main Figure 1
Method framework

不对应具体 dataset。

---

## Main Figure 2
GSE259262

```text
GBM
10 GIC vs 10 iNSC
multi-mark flagship
```

核心：

```text
H3K4me3
H3K27ac
H3K27me3
```

。

---

## Main Figure 3
GSE136128

```text
Prostate cancer ENZ resistance
3 vs 3
H3K27ac
AR/FOXA1 orthogonal validation
```

。

---

## Main Figure 4
GSE126004

```text
Breast cancer
IL6 / FOXA1
H3K27ac regulatory remodeling
```

。

---

## Main Figure 5
Cross-dataset domain phenotype comparison

例如比较：

```text
GBM
Prostate resistance
Breast metastatic signaling
```

的：

```text
Intensity-only
Breadth-only
Dual
```

比例和 transition。

---

# 29. Supplementary datasets

## GSE78206

```text
H3K4me3
MM1S
KDM5 inhibition
```

验证 Intensity vs Breadth separation。

---

## GSE151556

```text
known cancer broad-H3K4me3 loci
```

positive-reference validation。

---

## GSE145938

```text
MM H3K27ac
ROSE / SE decomposition
```

。

---

# 30. 如果工作量要控制

我建议最小版本只做：

```text
GSE259262
+
GSE136128
+
GSE78206
```

这三个其实已经能分别证明：

```text
multi-mark tumor generality

therapy-resistance remodeling

H3K4me3 Intensity/Breadth separation
```

。

---

# 31. 如果追求更强的 paper

建议：

```text
GSE259262
GSE136128
GSE126004
GSE78206
```

主文/补充结合。

再用：

```text
GSE151556
GSE145938
```

做 reference validation。

---

# 32. 不建议的数据策略

不要：

```text
20 个癌症数据集
```

。

不要：

```text
每个 cancer type 只放一张 browser track
```

。

不要：

```text
所有数据都做同一种分析
```

。

正确方式：

```text
每个 dataset 对应一个不同 methodological claim
```

。

---

# 33. 每个数据集对应的 claim

## GSE259262

> epiPortrait 可以跨 active/repressive marks 描述 tumor-associated domain phenotypes。

---

## GSE136128

> epiPortrait 可以刻画 therapy-resistance-associated domain remodeling。

---

## GSE126004

> epiPortrait 可以分解 oncogenic signaling-induced enhancer-domain remodeling，并与 TF dependency 对应。

---

## GSE78206

> Intensity gain 与 Breadth expansion 是可分离的 epigenomic events。

---

## GSE151556

> epiPortrait Breadth phenotype 可以恢复已有 cancer broad-H3K4me3 biological reference。

---

## GSE145938

> epiPortrait 可以将传统 H3K27ac super-enhancer ranking 拆分为不同 domain phenotypes。

---

# 34. 推荐 benchmark 统一输出

所有 datasets 建议统一输出：

```text
replicate call concordance
support fraction
class composition
Intensity distribution
NativeMaxPeakWidth
SignalDispersion
Combined class
transition
```

。

对于 comparative datasets：

```text
bp Jaccard
Cohen's kappa
Spearman rho
odds ratio
95% CI
```

。

---

# 35. 肿瘤数据特别需要强调异质性

不要把所有：

```text
patient-derived samples
```

都简单叫 replicates。

建议区分：

```text
technical / biological replicate
```

与：

```text
independent tumor samples
```

。

例如 GSE259262：

```text
10 GIC
```

更适合描述为：

```text
independent patient-derived tumor samples
```

。

replicate-support 的统计解释可以扩展为：

```text
population recurrence / cross-patient consistency
```

而不是实验重复。

---

# 36. 论文中应避免过度 claim

例如：

```text
Dual-Super predicts causal oncogenes
```

不建议。

更合适：

```text
Dual-Super domains are enriched near / linked to
cancer-relevant regulatory programs.
```

。

---

# 37. 推荐优先级

当前建议排序：

```text
1. GSE259262
2. GSE136128
3. GSE126004
4. GSE78206
5. GSE145938
6. GSE151556
```

。

其中：

> **GSE259262 最适合作为旗舰 tumor dataset。**

---

# 38. 最终推荐论文框架

```text
epiPortrait
↓
Replicate-aware epigenomic domain profiling
↓
GBM multi-mark validation
↓
therapy-resistance remodeling
↓
oncogenic enhancer remodeling
↓
H3K4me3 perturbation
↓
known cancer broad-domain / SE reference validation
```

。

这比使用非肿瘤 differentiation datasets 更符合 epiPortrait 当前面向：

```text
cancer epigenomics
```

的应用场景。

---

# 39. 一句话结论

> **epiPortrait 的发表级 benchmark 应优先围绕肿瘤表观遗传重塑构建：以 GSE259262 作为多标记 GBM flagship，GSE136128 作为前列腺癌耐药 remodeling 主验证，GSE126004 作为乳腺癌 oncogenic signaling enhancer-remodeling 机制验证；再用 GSE78206、GSE145938、GSE151556 分别补强 H3K4me3 Intensity/Breadth 分离、H3K27ac SE decomposition 和 known broad-domain cancer biology。**
