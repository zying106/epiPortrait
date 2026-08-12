# epiPortrait 学术论文 Benchmark 与数据集设计蓝图
## 生物信息学 / 表观遗传学 / Reviewer 视角

**版本背景：epiPortrait 当前定位为 Replicate-Aware Epigenomic Domain Profiling**  
**日期：2026-08-11**

---

# 1. 总体策略

如果 epiPortrait 要发表软件/方法学论文，benchmark 不应该设计成：

```text
找很多公开数据
→ 跑 epiPortrait
→ 画很多漂亮图
→ 和几个工具做 overlap
```

而应该设计成：

```text
论文 Claim
↓
为每个 Claim 选择最合适的数据
↓
预先定义 benchmark metric
↓
同一 preprocessing
↓
独立 validation data
↓
effect size + robustness + biological interpretation
```

推荐只围绕四个主要 claim：

### Claim 1
epiPortrait 能将 epigenomic domain 的：

```text
signal magnitude
```

与：

```text
native breadth geometry
```

分开刻画。

### Claim 2
epiPortrait 的：

```text
replicate-level call
→ support aggregation
```

能够减少 single-replicate-driven calls。

### Claim 3
epiPortrait 可以检测：

```text
condition-specific domain remodeling
```

而不仅仅是静态 Super-domain。

### Claim 4
框架可以泛化到：

```text
active marks
repressive marks
ChIP-seq
CUT&Tag
```

但不同 mark 使用不同生物学解释。

---

# 2. 推荐 benchmark 总体架构

建议最终论文采用：

```text
Benchmark 0
Synthetic truth benchmark

Benchmark 1
H3K4me3 broad-domain benchmark

Benchmark 2
H3K27ac / ROSE super-enhancer benchmark

Benchmark 3
Condition-remodeling benchmark

Supplementary Benchmark 4
H3K27me3 / H3K9me3 generalization

Supplementary Benchmark 5
ChIP-seq vs CUT&Tag technology robustness

Supplementary Benchmark 6
Runtime / scaling / parameter robustness
```

这已经足够形成完整方法论文。

不建议为了“数据多”堆十几个 unrelated datasets。

---

# 3. Benchmark 0：Synthetic Ground-Truth Benchmark

这是论文里非常重要的一部分。

真实 epigenomic data 没有真正 ground truth：

```text
这个 domain 真的是 Breadth-Super 吗？
这个 domain 真的是 Dual-Super 吗？
```

没有绝对答案。

因此 synthetic benchmark 是唯一可以正式计算：

```text
Sensitivity
Specificity
Precision
Recall
F1
classification accuracy
```

的 benchmark。

---

# 4. Synthetic benchmark 推荐设计

建议生成：

```text
5,000–20,000 domains
2 conditions
3 biological replicates / condition
```

不要直接使用 package 内嵌 `example_se` 作为论文 benchmark。

package example data 用于：

```text
documentation
tests
vignette
```

论文 benchmark 使用独立 simulation pipeline。

---

# 5. Synthetic phenotype truth

至少四类：

```text
Typical

Intensity-only

Breadth-only

Dual
```

比例可以有两套 simulation：

## Balanced benchmark

```text
25%
25%
25%
25%
```

用于 classifier metrics。

## Realistic benchmark

例如：

```text
Typical       70%
Intensity     10%
Breadth       10%
Dual          10%
```

用于 realistic performance。

---

# 6. Intensity simulation

对每个 domain：

```text
baseline signal
+
replicate noise
+
condition effect
```

模拟：

```text
Intensity-Super:
signal × 2–5

Typical:
signal × ~1
```

增加：

```text
library scaling
background
replicate CV
```

例如：

```text
CV = 10%
20%
30%
```

。

---

# 7. Breadth simulation

生成 native peak boundaries：

```text
Typical:
1–3 kb

Broad:
5–20 kb
```

然后加入：

```text
boundary jitter
±5%
±10%
±20%

peak fragmentation
dropout
one-replicate outlier
```

。

这是非常关键的 robustness benchmark。

---

# 8. SignalDispersion simulation

构造同 width、同 AUC 但不同内部 signal architecture：

```text
single compact peak

two separated peaks

diffuse plateau

fragmented signal
```

证明：

```text
SignalDispersion
```

确实能提供：

> width 和 total intensity 无法表达的内部 architecture 信息。

但不要把它升级成第三个 Super axis。

---

# 9. Synthetic benchmark 应比较的 baseline

不要只报告 epiPortrait 自己。

建议至少比较：

### Baseline A
group mean → 一次 Intensity cutoff

### Baseline B
per-replicate call → 无 replicate support

### Baseline C
shared consensus interval width 直接作为 Breadth

### Baseline D
native peak width top fixed percentile

例如：

```text
top 5%
```

。

然后比较 epiPortrait：

```text
native peak width
+
replicate support
+
no-call semantics
```

。

---

# 10. Synthetic benchmark metrics

推荐：

```text
Intensity axis:
Sensitivity
Specificity
Precision
Recall
F1

Breadth axis:
Sensitivity
Specificity
Precision
Recall
F1

Combined classes:
Macro-F1
Cohen's kappa
overall accuracy
```

另外必须报告：

```text
replicate outlier robustness
boundary jitter robustness
dropout robustness
```

。

---

# 11. Benchmark 1：H3K4me3 Breadth Benchmark

这是 Breadth-Super 最重要的 biological benchmark。

经典研究显示 broad H3K4me3 domains 与 cell identity/function genes 和 transcriptional consistency 有关。

因此这条 benchmark 应用于证明：

> epiPortrait Breadth axis 捕捉的是有真实生物学意义的 broad chromatin architecture，而不仅仅是 peak-boundary artifact。

---

# 12. 推荐核心数据：GSE85158

GSE85158 是非常适合 epiPortrait 的公开 benchmark data。

它包含多个 breast epithelial / breast cancer cell lines，并针对多个 histone marks 做了 profiling。

其中包括：

```text
MCF10A
MCF7
MB361
UACC812
SKBR3
AU565
HCC1954
MDA-MB-231
MDA-MB-436
...
```

而且 H3K4me3 / H3K27ac / H3K27me3 等很多 mark 有：

```text
rep1
rep2
```

。

推荐不要一次分析全部。

主文：

```text
MCF10A
MCF7
MDA-MB-231
HCC1954
```

即可。

Supplement 可以扩展其他 cell lines。

---

# 13. 为什么 GSE85158 很适合

一个数据集同时可以测试：

```text
H3K4me3 Breadth
H3K27ac Intensity/Breadth
H3K27me3 generalization
replicate support
cross-cell-line phenotype specificity
```

而且 protocol 相对一致。

这比拼很多来源不同的 GEO datasets 更干净。

---

# 14. H3K4me3 comparator

推荐经典 literature baseline：

```text
MACS2 --broad
→ rank peak width
→ top 5% broadest H3K4me3 domains
```

Cao et al. 在 K562/MCF-7 的分析中使用 MACS2 `--broad` 后按长度排序并取 top 5% 作为 broad H3K4me3 domains。

这可以作为：

> literature-concordance baseline

但不能当绝对 truth。

---

# 15. 非常重要：避免 circular benchmark

如果 epiPortrait Breadth 也使用 MACS2 native peak width：

```text
epiPortrait Breadth
vs
top 5% MACS2 width
```

然后得到很高 overlap，

这本身不能证明 epiPortrait 更好。

因为两者共享同一个 width source。

因此主证据应该是：

```text
1. replicate reproducibility
2. parameter/caller robustness
3. biological enrichment
4. cross-cell-type specificity
```

而不是单纯：

```text
overlap with top 5%
```

。

---

# 16. H3K4me3 benchmark metrics

### A. Literature concordance

```text
region overlap
bp Jaccard
broad-reference recall
width-rank enrichment
```

### B. Replicate robustness

```text
3/3 或 2/2 support rate
replicate call concordance
Cohen's kappa
```

### C. Biological validation

Breadth-Super-associated genes 是否：

```text
更富集 cell identity genes
更具 cell-type specificity
与 lineage-relevant gene sets 相关
```

### D. Expression

如果有 matched RNA：

```text
expression level
expression consistency
```

都可以作为外部 biological validation。

---

# 17. 第二个 H3K4me3 benchmark：K562

推荐再放一个 K562。

原因不是数据量，而是：

```text
K562
```

是 super-enhancer 和 broad H3K4me3 文献中非常经典的模型。

ENCODE 有公开 H3K4me3 ChIP-seq，如：

```text
ENCSR000AKU
```

以及 H3K27ac：

```text
ENCSR000AKP
```

。

Cao et al. 2017 又同时分析 K562 和 MCF-7 的：

```text
super-enhancer
broad H3K4me3
chromatin interaction
```

所以 K562 是非常好的 historical reference。

---

# 18. Benchmark 2：H3K27ac / ROSE Benchmark

这是 Intensity-Super 和二维 phenotype 最重要的 benchmark。

但论文问题不能设计成：

> epiPortrait 是否比 ROSE call 更多 SE？

正确问题应该是：

> ROSE 的一维 ranking 是否混合了 signal magnitude 与 domain breadth，而 epiPortrait 能否把这些不同 phenotype 分解出来？

---

# 19. 推荐数据

## K562 ENCODE

```text
H3K27ac:
ENCSR000AKP
```

这是经典 SE benchmark 场景之一。

## MCF7

可以使用：

```text
GSE85158
```

中的 MCF7 H3K27ac：

```text
GSM2258722
GSM2258723
```

两个 replicate。

---

# 20. ROSE benchmark 公平性

必须尽量统一：

```text
same H3K27ac source
same genome
same promoter exclusion
same constituent peak set
same stitching distance
```

如果用标准 ROSE：

```text
12.5 kb stitching
promoter exclusion
H3K27ac ranking
tangent point
```

epiPortrait 的 H3K27ac preset 也用同一 candidate construction。

这样比较才公平。

---

# 21. H3K27ac benchmark 不应把 ROSE 当 ground truth

ROSE 是：

```text
reference method
```

不是：

```text
biological truth
```

因此不要报告：

```text
accuracy against ROSE
```

。

推荐报告：

```text
bp Jaccard
region overlap
rank correlation
ROSE SE recall
```

称：

```text
method concordance
```

。

---

# 22. 真正关键分析：ROSE SE 被 epiPortrait 如何分解

例如：

```text
ROSE SE
│
├─ Intensity-Super
├─ Breadth-Super
└─ Dual-Super
```

然后比较：

```text
signal
width
SignalDispersion
gene expression
regulatory contacts
```

。

如果：

```text
Dual-Super
```

具有更高：

```text
regulatory activity
contact support
gene-expression association
```

这会非常漂亮。

---

# 23. 推荐使用 2026 ENCODE-rE2G 作为 H3K27ac 的独立生物学辅助验证

2026 Nature ENCODE enhancer–gene encyclopedia 提供超过：

```text
92 million
```

enhancer–gene interactions，覆盖：

```text
1,458 biosamples
369 cell types/tissues
```

并建立了 CRISPR/eQTL/GWAS benchmark framework。

尤其 K562 有大规模 CRISPR enhancer–gene benchmark。

epiPortrait 可以只用它做：

```text
external validation
```

例如问：

> epiPortrait 的 Dual-Super / Intensity-Super H3K27ac domains 是否更容易包含有 CRISPR/rE2G 支持的 regulatory elements？

但是：

> 不要声称 epiPortrait 在预测 causal target genes。

因为这不是 epiPortrait 的核心任务。

---

# 24. H3K27ac benchmark metrics

建议：

```text
ROSE overlap / bp Jaccard
Spearman rank correlation
class-specific signal
class-specific width
class-specific SignalDispersion
```

外部 validation：

```text
rE2G / CRISPR regulatory-link enrichment
RNA expression
P300/MED1 occupancy
```

。

---

# 25. Benchmark 3：Dynamic Domain Remodeling

epiPortrait 不应该只证明：

```text
静态 Super-domain
```

。

否则 reviewer 会问：

> 相比已有 SE / broad-domain ranking，真正增加了什么？

condition remodeling 是非常好的回答。

---

# 26. 首选动态数据：GSE212265

这是我非常推荐的 dataset。

H3K27ac：

```text
TSC:
2 replicates

EVT D3:
3 replicates

EVT D8:
3 replicates
```

。

因此非常适合：

```text
replicate-aware
+
condition transition
```

。

同时数据来自同一个 trophoblast differentiation study。

---

# 27. GSE212265 可以回答什么

可以分析：

```text
TSC → EVT D3
EVT D3 → EVT D8
TSC → EVT D8
```

。

即使当前 epiPortrait 不提供正式 multi-time trajectory API，

pairwise transition 已经够论文使用。

---

# 28. Dynamic benchmark 重点

分析：

```text
Typical
→ Intensity-Super

Typical
→ Breadth-Super

Intensity-Super
→ Dual-Super

Breadth-Super
→ Dual-Super
```

等 transition。

同时报告：

```text
Intensity change
NativeMaxPeakWidth change
SignalDispersion change
```

。

---

# 29. 外部 biological validation

GSE212265 本身研究的就是：

```text
dynamic enhancer usage during EVT differentiation
```

，而且同系列还有：

```text
P300
MED1
H3K4me1
```

数据。

因此可以检验：

> epiPortrait H3K27ac phenotype gains 是否伴随 P300/MED1 enhancer evidence。

这个非常适合论文。

---

# 30. 第二个动态 validation：GSE78913

这是一个很好的较小独立 validation。

MCF7：

```text
Estradiol 45 min:
H3K27ac rep1
H3K27ac rep2

EtOH control:
H3K27ac rep1
H3K27ac rep2
```

。

适合测试：

```text
acute enhancer activation
```

。

它可以作为：

```text
independent perturbation validation
```

而不是主 benchmark。

---

# 31. Dynamic benchmark 加一个 orthogonal comparator

推荐：

```text
DiffBind
```

或：

```text
csaw
```

。

但注意定位：

> 它们不是 epiPortrait 的直接竞争算法。

DiffBind 的主要目标是：

```text
differential binding/enrichment
```

。

所以这里用它回答：

> epiPortrait 的 Intensity gain/loss 是否与传统 differential enrichment 相符？

而 epiPortrait 额外提供：

```text
Breadth remodeling
Combined phenotype
domain architecture
```

。

这是“互补性 benchmark”。

---

# 32. 推荐做法

在同一 domain universe：

```text
DiffBind/csaw:
differential signal

epiPortrait:
Intensity transition
Breadth transition
Combined transition
```

比较：

```text
Intensity gain
vs
significant differential enrichment
```

预期应该有较高 concordance。

但：

```text
Breadth-only remodeling
```

可能是 DiffBind 不能直接表达的。

这正是 epiPortrait 的增量价值。

---

# 33. Benchmark 4：Repressive Mark Generalization

这部分放：

```text
main final figure
```

或者 Supplement。

不要占过多正文。

---

# 34. 强烈推荐数据：GSE201262

这个数据非常适合 epiPortrait。

MCF7：

```text
WT
vs
MCM2-2A
```

。

每组有两个 CUT&Tag replicates，包括：

```text
H3K27me3
H3K9me3
H3K27ac
H3K4me3
H3K4me1
H3K36me3
```

并且还有 RNA-seq。

研究本身报道：

```text
dramatic epigenetic reprogramming
especially H3K27me3
```

。

---

# 35. 为什么 GSE201262 很强

一次 benchmark 同时证明：

```text
active mark
repressive mark
CUT&Tag input
condition comparison
RNA biological interpretation
```

。

因此它比找很多零散的 repressive datasets 更有价值。

---

# 36. repressive marks 的术语

不要在正文写：

```text
H3K27me3 Super-enhancer
```

。

应该用包当前 mark-aware display：

```text
Extended-Domain
Intensity-Extreme
Dual-Extreme
```

和：

```text
expansion
contraction
intensity gain
intensity loss
```

。

---

# 37. H3K27me3/H3K9me3 benchmark metrics

```text
replicate support
width expansion/contraction
Intensity change
combined remodeling state
```

。

外部 biological validation：

```text
RNA expression
Polycomb-associated genes
original-study reported reprogrammed regions
```

。

---

# 38. Supplementary Benchmark 5：ChIP-seq vs CUT&Tag

如果论文篇幅允许，可以做。

2025 Nature Communications 有一个很适合的数据设计：

```text
K562
H3K27ac
H3K27me3
CUT&Tag
vs
ENCODE ChIP-seq
```

。

该研究直接将 K562 CUT&Tag 与：

```text
ENCSR000AKP H3K27ac
ENCSR000EWB H3K27me3
```

进行 benchmark。

---

# 39. 这个 benchmark 的目的

不是要求：

```text
CUT&Tag epiPortrait calls
=
ChIP-seq epiPortrait calls 100%
```

。

因为 assay 本身存在系统差异。

正确问题：

> 在共享高置信 domains 中，epiPortrait 的 domain phenotype 是否具有合理跨 assay concordance？

报告：

```text
Intensity rank correlation
Breadth class concordance
Combined-class concordance
shared-domain Jaccard
```

。

---

# 40. 不要把 cross-assay benchmark 做成主 claim

因为 ChIP-seq 与 CUT&Tag 自身就有 peak recovery 差异。

如果结果不是特别漂亮，

放 Supplement 完全合理。

GSE201262 已经足以证明：

```text
CUT&Tag compatibility
```

。

---

# 41. Benchmark 6：Consensus / Parameter Robustness

这一部分是 reviewer 最容易追问 Breadth 的地方。

必须做。

---

# 42. Shared-domain universe sensitivity

至少两套：

```text
A. conservative/no-gap consensus

B. default/stitching consensus
```

报告：

```text
bp Jaccard
domain count
merge/split/lost domains
Intensity class concordance
Breadth class concordance
Dual class concordance
transition concordance
```

。

---

# 43. Native peak mapping threshold sensitivity

当前 default：

```text
min_peak_overlap_fraction = 0.5
```

建议 benchmark：

```text
0.25
0.50
0.75
```

报告：

```text
mapping retention
Broad call concordance
class concordance
```

。

---

# 44. Peak caller sensitivity

Breadth 最大 reviewer 风险之一：

> 你的 Breadth 是不是 peak caller boundary artifact？

所以至少选：

```text
primary caller
+
one alternative reasonable caller/setting
```

。

例如：

```text
MACS2 standard/broad settings
```

或适合 CUT&Tag 的 alternative caller。

不一定需要 5 个 caller。

两个合理方案就够。

---

# 45. Intensity robustness

建议做：

```text
leave-one-replicate-out
```

但只作为 validation script。

报告：

```text
rank Spearman
cutoff variability
Intensity-Super retention
Jaccard
```

。

不需要变成 exported API。

---

# 46. Replicate benchmark

Synthetic + real data 都应报告：

```text
per-replicate call table
support fraction
unanimous rate
majority-supported rate
no-call rate
```

。

这个非常符合 epiPortrait 的核心定位。

---

# 47. Runtime / Scalability Benchmark

软件论文最好有。

但不需要做得很复杂。

建议模拟：

```text
1,000 domains
10,000 domains
50,000 domains
100,000 domains
```

以及：

```text
4 samples
6 samples
12 samples
```

。

报告：

```text
wall-clock time
peak memory
```

。

---

# 48. BEDPE performance

如果论文强调 annotation：

可以 Supplement 测：

```text
10k loops
100k loops
1M loops
```

。

但如果 annotation 只是 secondary module，

不必放正文。

---

# 49. 不建议硬和所有软件比较 runtime

例如：

```text
epiPortrait vs ROSE vs DiffBind vs csaw
```

并不完全公平，

因为输入和任务不同。

更合理的是：

```text
epiPortrait internal scaling
```

。

ROSE 只比较：

```text
biological/method concordance
```

。

---

# 50. 一个非常重要的 reviewer 原则：
# Freeze parameters before final benchmark

不要：

```text
K562 调一套 threshold
MCF7 调一套
TSC 调一套
H3K27me3 再调一套
```

最后每个数据都很漂亮。

这会让 reviewer 怀疑 overfitting。

推荐：

## Development dataset

用：

```text
K562 / MCF7
```

完成方法开发和 preset freeze。

## External validation

然后完全冻结参数，用：

```text
GSE212265
GSE201262
GSE78913
```

做外部 validation。

这非常重要。

---

# 51. 哪些参数允许 mark-aware preset？

允许：

```text
H3K27ac promoter exclusion
H3K27ac stitching

H3K27me3 no re-stitching
H3K9me3 no re-stitching
```

因为这是事先定义的：

```text
mark-aware biological model
```

。

但是不要根据结果临时改变：

```text
elbow
support rule
overlap threshold
```

。

---

# 52. 推荐统一 preprocessing

所有自己重新处理的 raw datasets：

```text
FASTQ
↓
alignment
↓
MAPQ filtering
↓
duplicate policy
↓
blacklist
↓
per-replicate native peak calling
↓
normalized BigWig
```

。

尽量统一。

---

# 53. 不建议直接混用不同来源的 processed BigWig

如果：

```text
一个 dataset 用 fold-enrichment
一个用 CPM
一个用 RPKM
一个用 raw pileup
```

然后直接比较 absolute Intensity，

容易产生批次问题。

---

# 54. 推荐 normalized BigWig 策略

对于 condition 内比较：

优先使用：

```text
同一 study
同一 assay
同一 processing
```

。

如果没有 spike-in，

主要使用：

```text
relative calling
```

。

对：

```text
reference / pooled absolute gain-loss
```

应更谨慎，并明确 quantitative comparability assumption。

---

# 55. Native peak 的原则

**绝对不要为了方便使用 merged/consensus peak 代替 per-replicate native peaks。**

epiPortrait 的 Breadth 方法学价值恰恰依赖：

```text
per-replicate native peak geometry
```

。

---

# 56. Benchmark 的 statistical reporting

不要只：

```text
p < 2e-16
```

。

建议优先：

```text
effect size
95% CI
distribution
```

。

---

# 57. Genomic overlap metrics

推荐：

```text
bp Jaccard
Dice coefficient
region-level overlap
```

。

不要只报告：

```text
% overlap
```

，因为它不对称。

---

# 58. Class concordance

推荐：

```text
Cohen's kappa
confusion matrix
class-specific retention
```

。

---

# 59. Continuous feature

推荐：

```text
Spearman rho
```

比较：

```text
rank
Intensity
width
```

。

---

# 60. Biological enrichment

推荐报告：

```text
odds ratio
95% CI
Fisher exact
```

。

例如：

```text
Dual-Super 是否富集 rE2G-supported enhancers？
Breadth-Super 是否富集 cell-identity genes？
```

。

---

# 61. Expression

可以使用：

```text
median expression
log2 fold difference
effect size
```

。

不要把 expression 当作 domain class 的 circular validation。

它只是：

```text
external biological correlate
```

。

---

# 62. 推荐正文 Figure 结构

## Figure 1
epiPortrait conceptual framework

```text
shared domains
Intensity
native Breadth
SignalDispersion
replicate support
transition
```

---

## Figure 2
Synthetic truth + robustness

```text
classification accuracy
replicate outlier
boundary jitter
no-call
```

---

## Figure 3
H3K4me3 Breadth benchmark

```text
K562 / MCF7
literature broad-domain concordance
replicate support
cell identity biology
```

---

## Figure 4
H3K27ac ROSE benchmark

```text
ROSE vs epiPortrait
SE decomposition
Intensity/Breadth/Dual
rE2G / expression validation
```

---

## Figure 5
Dynamic remodeling

推荐：

```text
GSE212265
TSC → EVT D3 → EVT D8
```

。

---

## Figure 6
Generalization

可以是：

```text
GSE201262
H3K27me3 / H3K9me3
```

。

如果篇幅不足放 Supplement。

---

# 63. Supplementary Figures

建议：

```text
S1 preprocessing/QC
S2 parameter sensitivity
S3 consensus sensitivity
S4 mapping threshold sensitivity
S5 peak-caller sensitivity
S6 GSE78913 independent perturbation
S7 repressive mark benchmark
S8 CUT&Tag vs ChIP
S9 runtime/memory
S10 annotation/BEDPE validation
```

。

---

# 64. 推荐最终核心数据集清单

| Dataset | Mark | Replicates | Purpose | Priority |
|---|---|---:|---|---|
| Synthetic | multi-axis | 3/group | true accuracy / robustness | Essential |
| ENCODE K562 | H3K4me3/H3K27ac | biological replicates | historical benchmark | Essential |
| GSE85158 | H3K4me3/H3K27ac/H3K27me3 | 2/sample group | cross-cell-line / breadth | Essential |
| GSE212265 | H3K27ac | TSC 2, EVT D3 3, EVT D8 3 | dynamic remodeling | Essential |
| GSE78913 | H3K27ac | 2+2 | acute perturbation validation | Strong supplement |
| GSE201262 | H3K27me3/H3K9me3/H3K27ac/H3K4me3 | 2+2 | repressive / CUT&Tag generalization | Strong supplement |
| K562 CUT&Tag benchmark | H3K27ac/H3K27me3 | multiple | assay robustness | Optional supplement |

---

# 65. 如果时间有限，最小发表级 benchmark

如果不想把工程拉得过大：

## 必做

```text
1. Synthetic
2. H3K4me3: K562 + MCF7
3. H3K27ac: K562/ROSE
4. GSE212265 dynamic remodeling
```

## Supplement

```text
5. GSE201262 repressive marks
6. runtime / parameter robustness
```

这已经是一篇完整方法学软件论文的 benchmark 骨架。

---

# 66. 最不建议的 benchmark 设计

### 不要

```text
只展示很多 loci browser snapshots
```

。

### 不要

```text
只和 ROSE 做 Venn diagram
```

。

### 不要

```text
把 top5% broad H3K4me3 当 ground truth
```

。

### 不要

```text
每个 dataset 都重新调 parameter
```

。

### 不要

```text
只报告 p-value
```

。

### 不要

```text
为了显示泛化而塞十几个 histone marks
```

。

---

# 67. 最重要的论文故事

推荐最终 paper story：

```text
Existing epigenomic analyses typically quantify enrichment
or rank broad regulatory regions in one dimension.

epiPortrait separates domain signal magnitude
from native breadth geometry while retaining
within-domain signal architecture.

Replicate-level evidence yields robust domain phenotypes.

These phenotypes reveal biologically distinct
active-domain states and condition-specific remodeling,
and generalize to broad repressive chromatin domains.
```

然后 benchmark 每一部分都围绕这四句话服务。

---

# 68. 最终建议

如果我是 reviewer，我最想看到：

```text
Synthetic truth
+
H3K4me3 breadth
+
H3K27ac ROSE decomposition
+
dynamic condition remodeling
+
one independent generalization dataset
```

而不是：

```text
20 datasets
+
20 Venn diagrams
```

。

数据选择应该少而强，每个数据集解决一个清楚的问题。

---

# 69. 推荐数据优先级

## 第一批立刻开始

```text
GSE85158
ENCODE K562
GSE212265
```

。

## 第二批

```text
GSE201262
GSE78913
```

。

## 最后如果论文需要

```text
K562 CUT&Tag vs ENCODE ChIP
ENCODE-rE2G / CRISPR enhancer links
```

。

---

# 70. 一句话结论

> **epiPortrait 的发表级 benchmark 最好不是“和多少软件比较”，而是建立 Synthetic truth → H3K4me3 breadth validity → H3K27ac ROSE decomposition → dynamic remodeling → repressive/CUT&Tag generalization 的证据链。核心真实数据用 K562、GSE85158 和 GSE212265 即可搭出主论文，GSE201262 与 GSE78913 做独立外部验证会非常强。**
