# epiPortrait13 深度冻结审查
## 核心算法、文档重定位、Embedded Data、Validation 与发表准备

**检查对象：`epiPortrait13.zip`**  
**日期：2026-08-11**  
**视角：生物信息学 / 表观遗传学 / 软件方法学 Reviewer**

> 本轮为源码静态深审。当前执行环境没有 R，因此无法实际运行
> `R CMD build`、`R CMD check`、`BiocCheck` 或 `testthat`。
> 结论中的运行时部分必须最终在真实 R/Bioconductor 环境再确认。

---

# 1. 总体结论

`epiPortrait13` 比 `epiPortrait12` 又明显向最终冻结推进了一步。

## 当前判断

### 核心算法

**可以继续保持冻结。**

没有看到需要重新设计：

```text
Intensity
Breadth-Super
SignalDispersion
replicate support
combined taxonomy
condition transition
width remodeling
annotation
BEDPE
RNA prioritization
```

的 P0 级结构问题。

### 文档

新定位已经基本成功：

> **Replicate-Aware Epigenomic Domain Profiling**

DESCRIPTION、README、vignette、CITATION 已经基本统一，不再把整个包简单定义成
“super-domain caller”。

### Embedded data

上一轮主要问题已经基本修复：

- 500 domains 全部位于 hg38 chr21 范围内；
- `BioReplicate` 已加入；
- `bw_path / peak_path` 改为 `NA_character_`；
- `Simulation_Scenario` 替代容易误导的 static truth class；
- `simulation_provenance` 已加入；
- BED writer 已按标准 BED `0-based half-open` 输出；
- extdata README 已说明 fixture provenance。

### 当前状态

> **没有新的 P0 blocker。**

但是仍有 **5 个 P1 级问题** 值得在最终 submission / paper benchmark 前修掉，
其中最重要的是：

1. `compare_superdomains(reference/pooled)` 没有完整继承原始 replicate-support 设置；
2. multi-condition candidate-universe 的文档示例还不够严谨；
3. 新增 publication validation 脚本本身有几个会削弱 benchmark 的问题；
4. `negative_policy = "allow"` 的实际语义与文档不完全一致；
5. README/vignette Quick Start 仍有几处旧工作流残留。

修完这些后，我认为可以进入：

```text
R CMD check
→ BiocCheck
→ benchmark
→ manuscript
```

而不是继续扩展新功能。

---

# 2. v12 遗留问题：确认已经修复

## 2.1 Package identity 已重新定位

DESCRIPTION：

```text
Title: Replicate-Aware Epigenomic Domain Profiling
```

Description 现在先定义：

```text
quantitative phenotyping
shared epigenomic domains
Intensity
native breadth
SignalDispersion
replicate support
condition remodeling
```

最后才介绍：

```text
Intensity-Super
Breadth-Super
Dual-Super
```

这是正确的层级关系。

README headline：

```text
Replicate-Aware Epigenomic Domain Profiling and Remodeling
```

vignette：

```text
Replicate-Aware Profiling of Epigenomic Domain Architecture
```

CITATION：

```text
epiPortrait: Replicate-Aware Epigenomic Domain Profiling
```

整体已经统一。

---

# 3. v12 stitching sentinel 已正确修复

上一版最重要的 mark-preset 语义问题已经处理。

当前：

```r
H3K27me3:
stitch_distance = 0L

H3K9me3:
stitch_distance = 0L

unknown mark:
stitch_distance = 0L
```

文档也明确：

```text
0 / NULL is a sentinel
→ preset workflow should NOT call stitch_epi_peaks()
```

并特别说明：

```text
stitch_epi_peaks(x, 0)
```

本身仍然会：

```text
reduce overlaps / adjacency
```

因此不能把它理解为 no-op。

vignette 当前示例：

```r
if (preset$stitch_distance > 0) {
    domains <- stitch_epi_peaks(
        domains,
        stitch_distance = preset$stitch_distance
    )
}
```

这已经正确闭环。

---

# 4. Generic mark preset 也已经修复

未知 mark 不再偷偷继承：

```text
12.5 kb stitching
active taxonomy
```

当前：

```r
stitch_distance = 0L
mark_class = "generic"
taxonomy_style = "generic"
```

并 warning。

这是安全的默认行为。

---

# 5. Embedded example 数据已经大幅改善

当前：

```text
500 domains
chr21 hg38
3 Control
3 Treatment
```

坐标：

```r
start = seq(1e6, by = 80000, length.out = 500)
```

最大 interval 小于 chr21：

```text
46,709,983 bp
```

并有：

```r
stopifnot(
    all(end(gr) <= seqlengths(gr)[...])
)
```

这是上一轮 P0 数据错误的正确修法。

---

# 6. example_se replicate design 正确

当前：

```text
Control_1
Control_2
Control_3

Treatment_1
Treatment_2
Treatment_3
```

并新增：

```text
BioReplicate = 1,2,3
```

。

这是非常好的教学设计。

同时 synthetic 数据刻意加入：

```text
3/3
2/3
1/3
no-call
```

支持模式，

可以真实演示：

```text
majority
all
fraction
min_valid_replicates
get_replicate_calls()
```

。

---

# 7. BED 坐标问题已经修复

当前：

```r
start = start(gr) - 1L
end   = end(gr)
```

输出标准：

```text
0-based half-open BED
```

extdata README 也明确说明。

正确。

---

# 8. Z-score 文档已经同步

vignette 当前明确：

> Row-wise Z-score is intentionally not an analysis normalization method.

并说明：

```text
Z-score
→ destroys cross-domain magnitude ranking
```

。

DESCRIPTION / normalize_portrait docs 也已经一致。

这个问题关闭。

---

# 9. eval=FALSE 已经显著清理

当前 vignette 只发现：

```text
installation chunk:
eval = FALSE
```

其余核心分析 chunks 都是 executable。

相比之前版本已经非常好。

---

# 10. P1-1：compare_superdomains(reference/pooled)
# 没有完整继承 replicate-support 参数

这是当前最值得修的**代码逻辑问题**。

位置：

```text
R/05_utils.R
compare_superdomains()
```

在：

```r
cutoff_scope = "reference"
```

或：

```r
cutoff_scope = "pooled"
```

时，

当前已经正确继承：

```text
method
log_transform_used
min_quality
support_rule
```

但是没有完整继承：

```text
min_replicate_support
tie_policy
```

。

---

# 11. fraction support 被硬编码成 0.5

当前代码：

```r
r_type <- .replicate_support_call(
    .classify_vs_cutoff(...),
    feature,
    support_rule,
    0.5,
    m_valid
)$group_type
```

这里：

```text
0.5
```

是硬编码。

假设用户原始调用：

```r
call_super_domains(
    ...,
    support_rule = "fraction",
    min_replicate_support = 0.67
)
```

那么 relative mode 是：

```text
0.67
```

但 reference / pooled transition 会偷偷变成：

```text
0.50
```

。

这会造成：

```text
primary calling semantics
!=
common-cutoff transition semantics
```

。

---

# 12. 推荐最小修法

先取 provenance：

```r
support_fraction <- if (
    !is.null(prov$min_replicate_support)
) {
    prov$min_replicate_support
} else {
    0.5
}
```

然后：

```r
.replicate_support_call(
    ...,
    support_rule = support_rule,
    min_replicate_support = support_fraction,
    min_valid_replicates = m_valid
)
```

。

---

# 13. min_valid_replicates 的 fallback 也不完全一致

当前：

```r
m_valid <- if (!is.null(prov$min_valid_replicates)) {
    prov$min_valid_replicates
} else {
    1L
}
```

但是主 calling 中：

```r
min_valid_replicates = NULL
```

并不是：

```text
1
```

。

它表示：

```text
根据 support rule 和实际 replicate 数自动解析
```

例如 majority：

```text
n=2 → 2
n=3 → 2
n=4 → 3
```

。

因此这里 fallback 到：

```text
1L
```

会改变 semantics。

---

# 14. 更好的最小修法

不要转换成 1：

```r
m_valid <- if (!is.null(prov)) {
    prov$min_valid_replicates
} else {
    NULL
}
```

直接把：

```text
NULL
```

继续传给：

```r
.replicate_support_call()
```

让 helper 根据每组实际 replicate 数重新解析。

这对于：

```text
unbalanced design
3 vs 2
```

尤其正确。

---

# 15. tie_policy 也没有继承

`.classify_vs_cutoff()` 当前固定使用：

```r
v > cutoff
```

即：

```text
strict
```

。

但是主调用允许：

```text
tie_policy = "strict"
tie_policy = "inclusive"
```

。

因此：

```r
call_super_domains(..., tie_policy = "inclusive")
```

之后再：

```r
compare_superdomains(..., cutoff_scope = "reference")
```

会重新变回：

```text
strict
```

。

---

# 16. 推荐修法

修改：

```r
.classify_vs_cutoff <- function(
    mat,
    cutoff,
    feature,
    tie_policy = c("strict", "inclusive")
)
```

。

然后：

```r
hit <- if (tie_policy == "inclusive") {
    v >= cutoff
} else {
    v > cutoff
}
```

并从 provenance 继承：

```r
c_tie <- prov$tie_policy %||% "strict"
```

。

---

# 17. 建议新增 regression tests

至少加入：

```text
reference/pooled preserves fraction=0.67
```

和：

```text
reference/pooled preserves inclusive tie_policy
```

。

否则这一类 advanced semantics 很容易以后再次回归。

---

# 18. P1-2：fraction 参数本身建议增加合法范围验证

当前：

```text
min_replicate_support
```

没有明显的 [0,1] validation。

特别危险的是：

```r
support_rule = "fraction"
min_replicate_support = 0
```

。

此时：

```r
required = ceiling(0 * n_reps) = 0
```

理论上可能产生不合理的 group calls。

推荐在：

```r
call_super_domains()
```

进入后加入：

```r
if (
    support_rule == "fraction" &&
    (
        length(min_replicate_support) != 1L ||
        !is.numeric(min_replicate_support) ||
        !is.finite(min_replicate_support) ||
        min_replicate_support <= 0 ||
        min_replicate_support > 1
    )
) {
    stop(
        "min_replicate_support must be in (0, 1] ",
        "when support_rule = 'fraction'."
    )
}
```

。

这是一个小改动，但很值得。

---

# 19. P1-3：sample_sheet 缺少明确的 SampleID / group invariant

`build_portrait_matrix()` 当前检查：

```text
SampleID
Condition
bw_path
```

是否存在。

但是没有看到明确检查：

```text
SampleID non-NA
SampleID non-empty
SampleID unique
Condition non-NA
Condition non-empty
```

。

---

# 20. 为什么 SampleID uniqueness 很重要

后续大量代码依赖：

```r
names(native_peaks) <- sample_names
```

和：

```r
np[[s]]
```

以及：

```r
colnames(se)
replicate_call_matrix
```

。

如果：

```text
SampleID duplicated
```

就可能造成 ambiguous list indexing / output columns。

即使 SummarizedExperiment 最终某一步可能报错，

也应该在入口处给出明确错误。

---

# 21. 推荐入口 validation

加入：

```r
if (
    anyNA(sample_sheet$SampleID) ||
    any(!nzchar(as.character(sample_sheet$SampleID))) ||
    anyDuplicated(sample_sheet$SampleID)
) {
    stop("SampleID must be non-missing, non-empty, and unique.")
}
```

。

Condition：

```r
if (
    anyNA(sample_sheet$Condition) ||
    any(!nzchar(as.character(sample_sheet$Condition)))
) {
    stop("Condition must be non-missing and non-empty.")
}
```

。

---

# 22. P1-4：negative_policy = "allow" 的真实语义需要修正

当前文档：

```text
allow:
keep negatives
```

。

但实际代码中：

```r
.signal_dispersion()
```

内部仍然：

```r
xc[xc < 0] <- 0
```

。

因此：

```text
Intensity
```

可能保留负值，

但是：

```text
SignalDispersion
```

仍然把负值裁剪为 0。

---

# 23. 这在数学上其实有道理

SignalDispersion 使用：

```text
signal-weighted SD
```

：

```text
w_i = x_i / sum(x_i)
```

。

负 signal 不适合作为 probability-like weights。

因此：

```text
SignalDispersion 不接受 signed weight
```

在数学上合理。

真正的问题是：

> **当前文档没有明确这一点。**

---

# 24. 另外 default super calling 也不是真正 signed-signal aware

对于 Intensity：

```r
log_transform = NULL
```

默认会：

```r
log10(pmax(x, 0) + 1)
```

。

因此即使：

```r
negative_policy = "allow"
```

负 Intensity 在默认 calling 时还是被压成：

```text
0
```

。

`tangent` 方法也会把负值 clamp 到 0。

---

# 25. 推荐不要再扩大 signed-signal 支持

为了少改代码，推荐只修改文档：

```text
negative_policy = "allow"
```

明确写：

> Intended for advanced within-sample Intensity analyses with an explicitly
> compatible downstream transform. SignalDispersion always uses non-negative
> weights, and default log-transformed / tangent Super calling clips negative
> values to zero.

或者更简单：

> Signed tracks are not part of the recommended canonical workflow.

这样就够。

不建议为了 signed BigWig 再重构 SignalDispersion。

---

# 26. P1-5：multi-condition candidate universe 的推荐工作流还应该更明确

README 当前 Core Module 1 的示例仍接近：

```r
consensus_peaks <- get_consensus_peaks(
    peak_list,
    min_reps = 2
)
```

。

对单 condition replicate 很自然。

但是 epiPortrait 当前主要定位已经是：

```text
multi-condition remodeling
```

。

主推荐应该明确：

```text
Condition A replicate consensus
Condition B replicate consensus
↓
union of condition consensuses
↓
optional stitching
```

。

---

# 27. 为什么不能只强调“所有样本一起 min_reps=2”

例如：

```text
Control:
C1
C2
C3

Treatment:
T1
T2
T3
```

推荐：

```text
C consensus
T consensus
↓
union
```

这样：

```text
condition-specific domains
```

会保留。

这与：

```text
transition / gain / loss
```

的科学目标更一致。

---

# 28. 推荐 README/vignette 加一个明确范式

例如：

```r
ctrl_consensus <- get_consensus_peaks(
    ctrl_peak_list,
    min_reps = 2
)

treat_consensus <- get_consensus_peaks(
    treat_peak_list,
    min_reps = 2
)

candidate_domains <- GenomicRanges::reduce(
    c(ctrl_consensus, treat_consensus)
)
```

然后：

```r
if (preset$stitch_distance > 0) {
    candidate_domains <- stitch_epi_peaks(
        candidate_domains,
        preset$stitch_distance
    )
}
```

。

这与当前算法逻辑最匹配。

---

# 29. Publication validation 新问题 1：
# SignalDispersion 示例不是严格 controlled comparison

文件：

```text
validation/01_signal_dispersion_validation.R
```

当前：

```r
compact <- c(
    rep(0, 40),
    rep(1, 20),
    rep(0, 40)
)
```

长度：

```text
100
```

。

但是：

```r
bimodal <- c(
    rep(1, 10),
    rep(0, 60),
    rep(1, 10)
)
```

长度：

```text
80
```

。

---

# 30. 为什么这是 benchmark 设计问题

注释声称比较：

```text
same occupied width
different internal architecture
```

但 candidate interval length 同时从：

```text
100
```

变成：

```text
80
```

。

因此 reviewer 可以说：

> dispersion difference 部分来自 coordinate span change，而不只是
> internal signal redistribution。

---

# 31. 最小修复

改成：

```r
bimodal <- c(
    rep(1, 10),
    rep(0, 80),
    rep(1, 10)
)
```

长度同样：

```text
100
```

。

这样：

```text
compact:
20 occupied bp centered

bimodal:
20 occupied bp split at two edges

same interval width
same occupied bp
same signal mass
different architecture
```

。

这才是很干净的 validation。

---

# 32. Publication validation 新问题 2：
# consensus A/B 当前实际上很可能完全相同

文件：

```text
validation/02_breadth_super_validation.R
```

当前：

```r
dom_A <- rowRanges(example_se)

dom_B <- stitch_epi_peaks(
    dom_A,
    stitch_distance = 50000
)
```

。

但是 `example_se`：

```text
domains start every 80 kb
```

。

Typical domain width：

```text
1 kb
```

Broad scenario domain width：

```text
20 kb
```

因此间隔大约：

```text
Typical → next:
~79 kb

Broad → next:
~60 kb
```

。

`50 kb` stitching：

```text
50 kb < 60 kb
```

所以：

> **理论上几乎不会 merge 任何 domain。**

也就是说：

```text
dom_A
```

和：

```text
dom_B
```

可能完全相同。

---

# 33. 这会让当前 “consensus robustness” benchmark 失效

脚本本来想证明：

```text
different consensus construction
↓
same native peak broad calls
↓
domain-level robustness
```

。

但如果：

```text
A == B
```

得到：

```text
Jaccard = 1
PeakBroadCall concordance = 1
```

几乎没有验证价值。

---

# 34. 推荐最小修法

例如用：

```r
stitch_distance = 70000
```

这样：

```text
20 kb broad interval
→ 下一 domain gap ~60 kb
```

会产生部分 merge，

但：

```text
typical gap ~79 kb
```

通常仍不 merge。

比：

```text
100 kb
```

把大部分基因组全部 chain 起来更合理。

更重要的是加 invariant：

```r
stopifnot(
    length(dom_B) < length(dom_A)
)
```

。

只有真正得到：

```text
different candidate universes
```

才允许继续 robustness benchmark。

---

# 35. Publication validation 新问题 3：
# validation/02 可能缺少 SummarizedExperiment namespace

脚本目前加载：

```r
library(epiPortrait)
library(GenomicRanges)
```

但后面直接调用：

```r
rowData()
rowRanges()
```

。

这两个函数来自：

```text
SummarizedExperiment
```

。

`epiPortrait` importing SummarizedExperiment 并不意味着用户 global search path
自动 attach 这些函数。

因此独立执行这个脚本时很可能出现：

```text
could not find function "rowData"
```

或：

```text
could not find function "rowRanges"
```

。

---

# 36. 推荐修法

最简单：

```r
suppressMessages(
    library(SummarizedExperiment)
)
```

。

或者所有位置使用：

```r
SummarizedExperiment::rowData()
SummarizedExperiment::rowRanges()
```

。

对于 publication validation，我更推荐显式 namespace。

---

# 37. validation 目录整体方向是对的

这两个脚本本身非常值得保留。

尤其：

```text
01 SignalDispersion analytical validation
02 Breadth consensus independence validation
```

很好地对应：

```text
architecture descriptor
canonical breadth semantics
```

。

因此：

> 不需要删除 validation，只需要把上述 benchmark design 修严谨。

---

# 38. README/vignette Quick Start 仍有旧叙事残留

## 问题 A

vignette 写：

```text
500 domains
6 samples
and 3 dynamic assays
```

但 `example_se` 当前实际上有：

```text
5 assays
```

：

```text
Intensity
SignalDispersion
NativeMaxPeakWidth
NativeOccupiedWidth
NativePeakCount
```

。

建议直接改：

```text
5 quantitative assays
```

或者：

```text
two canonical measurement layers plus three native-geometry descriptors
```

。

---

# 39. 问题 B：Quick Start 仍以 tangent 为主

README Quick Start：

```r
call_super_domains(
    ...,
    method = "tangent",
    log_transform = FALSE
)
```

vignette Quick Start 也是类似。

但是当前文档已经明确：

```text
elbow = epiPortrait default/main analysis

tangent + raw
= optional H3K27ac / ROSE benchmark
```

。

因此 top-level quick start 继续使用 tangent 会给新用户错误信号。

---

# 40. 推荐改法

主 Quick Start：

```r
se <- call_super_domains(
    se,
    feature = "Intensity",
    mode = "per_group",
    group_var = "Condition",
    verbose = FALSE
)
```

即：

```text
default elbow
default log convention
per_group
```

。

然后专门在：

```text
ROSE benchmark section
```

展示：

```r
method = "tangent",
log_transform = FALSE
```

。

---

# 41. 问题 C：README Quick Start 仍然使用 default global_consensus

当前：

```r
se <- call_super_domains(
    se,
    feature = "Intensity",
    method = "tangent",
    log_transform = FALSE
)
```

没有：

```text
mode = "per_group"
```

。

Breadth 也是 default global mode。

但 README 前文已经说：

> `per_group` is recommended for the main analysis of multi-condition studies.

因此 Quick Start 应该与主 recommendation 一致。

---

# 42. 推荐 Quick Start

```r
se <- normalize_portrait(
    example_se,
    method = "None"
)

se <- call_super_domains(
    se,
    feature = "Intensity",
    mode = "per_group",
    group_var = "Condition",
    verbose = FALSE
)

se <- call_super_domains(
    se,
    feature = "Breadth",
    mode = "per_group",
    group_var = "Condition",
    verbose = FALSE
)

se <- combine_superdomain_calls(
    se,
    group_var = "Condition"
)

table(
    SummarizedExperiment::rowData(se)$
        Combined_Class__Control
)
```

。

这才是真正代表当前 package identity 的最小 workflow。

---

# 43. Vignette “complete core workflow” 的表述不准确

当前写：

> the rest of this vignette runs on a tiny synthetic BigWig dataset...
> so the complete core workflow is executed here

但是：

```r
build_portrait_matrix()
```

确实对 tiny BigWig 得到 `se` 后，

后续 normalization 又执行：

```r
se <- normalize_portrait(
    example_se,
    method = "None"
)
```

。

也就是说：

> tiny BigWig 主要实际验证的是 I/O + matrix construction，
> 后续 statistically meaningful calling 又切回 500-domain `example_se`。

---

# 44. 这不是算法问题，而是文档应该说实话

实际上这种设计是合理的：

```text
tiny 3-domain BigWig
→ I/O integration fixture

500-domain example_se
→ inflection / replicate calling demo
```

因为只有 3 个 domains 本来就不适合展示稳定 hockey-stick calling。

因此不要为了“complete pipeline”强行扩充 tiny BigWig。

推荐直接写：

> The tiny BigWigs exercise the real BigWig I/O and portrait-matrix
> construction path. Subsequent inflection-based calling is demonstrated on
> the 500-domain `example_se`, because a three-domain I/O fixture is not
> intended to provide a statistically meaningful ranked distribution.

这反而更专业。

---

# 45. Vignette log-transform 文案有一个小错误

当前：

```text
By default Intensity and width-like features are log10(x + 1)-transformed
```

实际上代码 auto：

```r
feature %in% c(
    "Intensity",
    "SignalDispersion"
)
```

才自动 log。

Breadth：

```text
peak-level native width
```

走独立 pipeline。

`IntervalWidth` 也不是 auto-log。

因此建议改为：

> By default `Intensity` and `SignalDispersion` use the package's automatic
> log10(x+1) ranking transform. Breadth uses the native peak-width calling
> pipeline described below.

---

# 46. P2：README vignette link

README 当前：

```text
articles/epiPortrait.html
```

。

如果没有 pkgdown site / `articles/` 目录，

GitHub README 上会成为 broken link。

建议当前 pre-submission 阶段直接：

```text
vignettes/epiPortrait.Rmd
```

。

以后有 pkgdown 后再换回：

```text
articles/epiPortrait.html
```

。

---

# 47. P2：PCA 文档仍有旧 terminology

当前：

```text
PCA of geometric features
```

以及：

```text
Performs Principal Component Analysis on a specific geometric feature
```

。

现在 Intensity 本身并不是“geometric feature”。

建议统一成：

```text
quantitative domain feature
```

。

例如：

```text
PCA of quantitative domain features
```

。

---

# 48. P2：LICENSE 可以进一步简化

DESCRIPTION：

```text
License: GPL (>= 3) + file LICENSE
```

而 `LICENSE` 本身又包含完整 GPLv3 license text。

对标准 GPL R/Bioconductor package 来说，这是冗余的。

最终真实：

```text
R CMD check
BiocCheck
```

如果对此没有 note，可以保持。

但从 package hygiene 看，更简洁的是：

```text
License: GPL (>= 3)
```

并按 R package convention 处理 LICENSE。

这不是当前 blocker。

---

# 49. H3K4me3 preset 的 5 kb stitching：
# 不建议现在立即重构，但发表时必须 benchmark / justify

当前：

```r
H3K4me3:
stitch_distance = 5000L
```

。

canonical Breadth-Super 本身使用：

```text
native PeakWidth
```

所以它不依赖 stitched-domain width。

这避免了最严重的 circularity。

但 H3K4me3 broad-domain 文献通常更强调：

```text
native broad peak/domain width
```

而不是类似 ROSE 的 stitching。

因此论文 benchmark 中建议至少测试：

```text
H3K4me3:
0 kb
vs
5 kb candidate stitching
```

看：

```text
Breadth class concordance
domain mapping
gene annotation
```

是否稳定。

如果稳定，可以保留 5 kb preset。

如果不稳定，再考虑：

```text
H3K4me3 default = 0
```

。

当前不建议为了这个问题立即再次改核心设计。

---

# 50. 当前最值得新增的 tests

## Test 1

```text
support_rule = fraction
min_replicate_support = 0.67
```

在：

```text
relative
reference
pooled
```

中保持一致。

---

## Test 2

```text
tie_policy = inclusive
```

在：

```text
reference / pooled
```

中保持一致。

---

## Test 3

```text
duplicate SampleID
```

必须 fail loudly。

---

## Test 4

```text
NA / empty Condition
```

必须 fail loudly。

---

## Test 5

```text
min_replicate_support = 0
1.1
NA
```

在 fraction mode 报错。

---

# 51. 当前不建议再做什么

不要再：

```text
新增 classifier
新增 statistical test
新增 report engine
新增 enrichment engine
新增 temporal module
新增 regulatory scoring
```

。

这些旧模块被删除是正确的。

当前 package 的聚焦度明显比早期版本更好。

---

# 52. 代码层面的最终冻结判断

## Intensity

保持：

```text
per replicate ranking
→ per replicate cutoff
→ support aggregation
```

正确。

---

## Breadth

保持：

```text
genome-wide native PeakWidth
→ per-replicate inflection
→ PeakBroadCall
→ unique peak-to-domain mapping
→ replicate support
```

正确。

---

## SignalDispersion

保持：

```text
secondary continuous architecture descriptor
```

正确。

不作为：

```text
third Super axis
```

正确。

---

## Combined class

保持：

```text
Typical
Intensity-Super
Breadth-Super
Dual-Super
Uncertain
```

正确。

---

## Transition

Relative mode 明确为：

```text
Relative_Prominence_Up / Down
```

而不是 absolute gain/loss。

正确。

Reference / pooled 也保留 replicate-aware classification。

仅需要把上面提到的：

```text
fraction threshold
tie policy
min-valid semantics
```

完全继承。

---

# 53. Annotation 层

本轮没有发现新的需要重构的 annotation 问题。

当前边界仍然合理：

```text
linear TxDb evidence
BEDPE promoter-contact evidence
RNA expression prioritization
```

并坚持：

```text
candidate gene
!=
causal target gene
```

。

这一部分保持冻结。

---

# 54. 文档定位判断

这次 repositioning 是成功的。

当前最好的 package identity 已经清楚：

> **replicate-aware epigenomic domain profiling**

而不是：

```text
4D peak shape
```

也不是：

```text
generalized ROSE
```

。

`Super-domain` 现在是：

```text
one analytical classification layer
```

而不是整个包的唯一身份。

这是正确方向。

---

# 55. Publication benchmark readiness

目前 package 已经具备进入肿瘤 benchmark 的条件。

但是正式跑：

```text
GBM
prostate resistance
breast cancer signaling
multiple myeloma
```

前，

建议先关闭：

```text
compare_superdomains common-cutoff semantics
```

。

否则论文中如果使用：

```text
reference cutoff
pooled cutoff
fraction support
```

可能出现主 calling 与 transition support semantics 不完全一致。

---

# 56. 推荐修复顺序

## P1-A — 必修，代码

```text
compare_superdomains:
inherit min_replicate_support
inherit tie_policy
preserve NULL min_valid_replicates semantics
```

---

## P1-B — 必修，入口安全

```text
SampleID unique/nonempty
Condition nonempty
fraction threshold validation
```

---

## P1-C — 必修，publication validation

```text
SignalDispersion controlled synthetic vectors
02 validation namespace
real A/B consensus perturbation
```

---

## P1-D — 文档

```text
Quick Start → per_group + elbow
3 dynamic assays → correct wording
tiny fixture ≠ full ranked-calling demonstration
log-transform wording
```

---

## P1-E — multi-condition workflow

```text
condition consensus
→ union
→ optional stitching
```

明确写入 README/vignette。

---

## P2

```text
README link
PCA wording
LICENSE hygiene
H3K4me3 stitching sensitivity
```

。

---

# 57. 最终评分

按当前源码静态状态：

| 模块 | 当前评价 |
|---|---:|
| Core conceptual design | **9.5 / 10** |
| Intensity calling | **9.4 / 10** |
| Breadth calling | **9.4 / 10** |
| Replicate semantics | **9.2 / 10** |
| Transition engine | **8.8 / 10** |
| SignalDispersion | **9.2 / 10** |
| Annotation | **9.2 / 10** |
| Embedded example data | **9.3 / 10** |
| README / vignette positioning | **9.0 / 10** |
| Publication validation scripts | **7.8 / 10** |
| Submission hygiene | **8.8 / 10** |

publication validation 分数较低不是因为算法本身差，

而是：

```text
validation 01 comparison 还不完全 controlled
validation 02 consensus A/B 实际可能没有改变 universe
validation 02 namespace execution risk
```

这些都属于非常容易修的小问题。

---

# 58. Freeze 判断

## 当前

```text
CORE ALGORITHM:
FREEZE

PUBLIC API:
FREEZE

ANNOTATION:
FREEZE

VISUALIZATION:
FREEZE unless bug

DOCUMENTATION:
small cleanup

VALIDATION:
fix before manuscript benchmark

SUBMISSION:
run actual R CMD check / BiocCheck after fixes
```

---

# 59. 一句话判断

> **epiPortrait13 已经不需要继续做算法重构；核心可以冻结。当前最值得修的是 `compare_superdomains()` 在 reference/pooled 模式下完整继承 replicate-support/tie semantics，以及把 publication validation 和 Quick Start 做到真正严谨。修完这些，再做真实 R CMD check/BiocCheck 和肿瘤 benchmark，比继续增加功能更重要。**
