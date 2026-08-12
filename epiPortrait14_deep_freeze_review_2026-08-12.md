# epiPortrait14 深度冻结审查
## 核心算法、Breadth 新增路径、Validation、文档与 Bioconductor 发布准备

**检查对象：** `epiPortrait14.zip`  
**日期：** 2026-08-12  
**审查视角：** 生物信息学 / 表观遗传学 / R-Bioconductor / 方法学 Reviewer

---

# 1. 总体结论

`epiPortrait14` 的**核心 Intensity × native Breadth 算法仍然可以冻结**。

相较 `epiPortrait13`，上一轮最重要的代码遗留问题已经基本修复：

- `compare_superdomains(reference/pooled)` 已继承：
  - `support_rule`
  - `min_replicate_support`
  - `tie_policy`
  - `min_valid_replicates`
- `support_rule = "fraction"` 已验证 `min_replicate_support ∈ (0,1]`
- `build_portrait_matrix()` 已验证：
  - SampleID 非 NA
  - SampleID 非空
  - SampleID 唯一
  - Condition 非 NA / 非空
- `negative_policy = "allow"` 文档已经明确：
  - Intensity 可保留 signed signal
  - SignalDispersion 内部仍使用非负权重
  - canonical Super workflow 并非 signed-signal aware
- README 已补上 multi-condition：
  - condition consensus → union → optional stitch
- Quick Start 已改成：
  - `per_group`
  - default elbow
  - default log convention

因此：

> **没有发现需要重新设计 Intensity、Breadth、Combined taxonomy、replicate support、annotation 或 transition engine 的新 P0 算法问题。**

但是，`epiPortrait14` 又新增了一个明显的“freeze 后扩功能”趋势：

```text
Breadth quantile_cutoff
+
Breadth n_bootstrap
+
大型 synthetic benchmark
```

其中部分实现有明确 edge-case / semantics 问题，而且你已经决定**不把 synthetic benchmark 作为论文核心证据**。

因此本轮最重要的建议不是继续扩展，而是：

> **收回不必要的新增面，修少量边界错误，重新冻结。**

---

# 2. 当前冻结判断

```text
Intensity core                  FREEZE
Breadth native-peak caller      FREEZE
replicate support               FREEZE
Combined taxonomy               FREEZE
condition transitions           FREEZE
annotation                      FREEZE
visualization                   FREEZE unless bug

Breadth quantile path           RECONSIDER / REMOVE preferred
Breadth bootstrap               KEEP only if useful; fix one bug
synthetic benchmark             REMOVE from formal manuscript workflow
validation/02 interpretation    REVISE
documentation                   SMALL CLEANUP
release/git state               MUST FIX
```

---

# 3. epiPortrait13 的关键遗留问题：哪些已经修掉？

## 3.1 compare_superdomains provenance inheritance —— 已修

当前 `R/05d_combine.R`：

```r
support_fraction <- if (!is.null(prov) &&
                        !is.null(prov$min_replicate_support)) {
    prov$min_replicate_support
} else {
    0.5
}

c_tie <- if (!is.null(prov) &&
             !is.null(prov$tie_policy)) {
    prov$tie_policy
} else {
    "strict"
}

m_valid <- if (!is.null(prov) &&
               !is.null(prov$min_valid_replicates)) {
    prov$min_valid_replicates
} else {
    NULL
}
```

然后继续传入：

```r
.replicate_support_call(...)
.classify_vs_cutoff(..., tie_policy = c_tie)
```

这是正确修复。

同时 tests 已加入：

```text
reference/pooled inherit fraction + tie_policy
.classify_vs_cutoff honors tie_policy
```

### 结论

**关闭。**

---

# 4. min_replicate_support validation —— 已修

当前：

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
    stop(...)
}
```

并有 test：

```text
0
1.1
NA
"half"
```

均要求报错。

### 结论

**关闭。**

---

# 5. SampleID / Condition invariant —— 已修

当前 `build_portrait_matrix()` 明确阻止：

```text
duplicate SampleID
empty SampleID
NA SampleID
NA/empty Condition
```

这是必要且正确的 object-entry contract。

### 结论

**关闭。**

---

# 6. negative_policy = allow —— 文档基本闭环

当前文档已经明确：

```text
allow
→ signed Intensity 可保留
```

但：

```text
SignalDispersion
→ negative clipped to 0 internally
```

且：

```text
default log transform
tangent
```

也并非真正 signed-aware。

这和当前代码语义基本一致。

### 建议

不要继续扩展 signed-signal calling。

### 结论

**关闭。**

---

# 7. multi-condition candidate universe —— README 已修，但 vignette 还缺一半

README 已加入：

```r
ctrl_consensus <- get_consensus_peaks(ctrl_peak_list, min_reps = 2)
treat_consensus <- get_consensus_peaks(treat_peak_list, min_reps = 2)

candidate_domains <- GenomicRanges::reduce(
    c(ctrl_consensus, treat_consensus)
)
```

这是正确范式：

```text
per-condition consensus
→ union
→ optional stitch
```

但是 vignette 主体目前仍没有同样清楚地加入这一范式。

### 建议

README 和 vignette 同步。

---

# 8. 本轮最重要的新设计问题：
# Breadth quantile_cutoff 是否真的有必要？

`epiPortrait14` 新增：

```r
call_super_domains(
    feature = "Breadth",
    quantile_cutoff = 0.95
)
```

当前解释：

```text
top 5% widest native peaks
```

。

从代码和 validation 内容看，这一新增路径与 balanced synthetic benchmark 直接绑定：

```text
balanced
→ quantile_cutoff = 0.50
→ Intensity + Breadth
```

。

但目前项目方向已经决定：

> synthetic ground-truth benchmark 不需要成为论文主线，甚至可以不做。

因此现在应重新问：

> **是否值得为了一个不准备主推的 synthetic scenario，扩大 canonical Breadth API？**

我的判断：

## 推荐

**优先删除 Breadth 的 `quantile_cutoff` 支持。**

保留：

```text
Breadth canonical:
native PeakWidth
→ data-driven elbow / tangent
→ no reliable inflection = no_call
```

。

理由：

1. 保持 Breadth 的生物学定义干净；
2. 不引入人为 top-X% broad threshold；
3. 避免 integer peak widths 的大量 ties；
4. 避免用户把 `top 5%` 当成“Broad”的标准定义；
5. 让 v1.0 API 更小、更稳定；
6. 之前已经明确过“不自动 top5% / fixed kb fallback”。

注意：

这不是说 quantile threshold 数学上错误，而是：

> **它对当前 epiPortrait 的核心定位没有必要。**

---

# 9. 如果决定保留 Breadth quantile：
# 当前至少有 3 个必须修的问题

---

## 9.1 `valid_chroms` 过滤后 0 eligible peaks 会报错

当前：

```r
eligible <- .eligible_native_peaks(peaks, valid_chroms)
w <- width(eligible)

if (!is.null(quantile_cutoff)) {
    cutoff_value <- quantile(w, ...)
}
```

如果：

```text
原 native peaks 存在
但 valid_chroms 把全部 peak 过滤掉
```

则：

```text
length(w) = 0
```

。

默认 elbow 路径可以通过 `find_hockey_inflection()` 返回 no_call。

但 quantile 路径会直接进入：

```r
stats::quantile(numeric(0), ...)
```

，产生运行错误。

### 必修

在 quantile/elbow 分支之前：

```r
if (length(eligible) == 0L) {
    rep_broad[[s]] <- logical(0)
    rep_call_valid[s] <- FALSE
    prov_replicates[[s]] <- list(
        replicate = s,
        n_eligible = 0L,
        call_status = "no_call",
        reason = "no eligible native peaks after filtering"
    )
    next
}
```

---

# 10. quantile 的 “top X%” 语义目前并不严格成立

generic Intensity quantile route：

```r
cutoff_value <- quantile(...)
inflection_idx <- which(Value >= cutoff_value)[1]
```

随后又：

```r
cutoff_value <- rank_df$Value[inflection_idx]
```

再根据默认 strict：

```r
Value > cutoff_value
```

。

因此，例如：

```text
n = 100
q = 0.95
```

最终不一定得到 5%。

连续值情况下可能少一个点；
大量 ties 时偏差可能更明显。

Breadth 的 peak width 又是整数，ties 更常见。

例如所有 peak 都同宽：

```text
strict
→ 0% Broad

inclusive
→ 100% Broad
```

而不是 5%。

### 结论

如果保留 quantile：

不要写：

> top `100*(1-q)%` peaks

更准确应写：

> **quantile-based width threshold**

并明确 tie policy 可使实际比例偏离名义比例。

或者代码中直接用原始 quantile value 作为 cutoff，不把 cutoff 改写成第一个 observed value。

---

# 11. Breadth bootstrap 有一个明确实现 bug

当前：

```r
br <- find_hockey_inflection(...)
cutoffs[b] <- br$cutoff_value
```

但 `find_hockey_inflection()` 在：

```text
quality < min_quality
```

时会返回：

```text
call_status = "no_call"
```

同时：

```text
cutoff_value 仍可能是 finite
```

。

于是 Breadth bootstrap：

```r
finite_idx <- is.finite(boot_cutoffs)
bootstrap_success_rate <- mean(finite_idx)
```

会把：

```text
no_call 但 finite cutoff
```

错误计为：

```text
bootstrap success
```

。

Intensity bootstrap 没这个问题，因为 `.call_super_domains_on_vector()` 在 no_call 时最终会把 cutoff 置为 NA。

### 正确修法

```r
cutoffs[b] <- if (
    identical(br$call_status, "called")
) {
    br$cutoff_value
} else {
    NA_real_
}
```

这是一个明确的 P1 correctness bug。

---

# 12. validation/02 的结果本身暴露了一个重要事实

当前保存的输出：

```text
Consensus A:
Breadth-Super = 54

Consensus B:
Breadth-Super = 51

PeakBroadCall concordance = 1.000
```

但是：

```text
bp-level genomic Jaccard = 0.243
```

。

这不是一个“robust domain geometry”的结果。

它真正说明的是：

> **peak-level Broad status 对 consensus construction 独立，但映射后的最终 shared-domain geometry 仍明显依赖 candidate universe。**

这其实并不违反 epiPortrait 的算法设计：

```text
Broad status
→ native peak-level
```

是 independent。

但：

```text
最终 Breadth-Super domain 坐标
```

仍然是：

```text
native Broad evidence
+
shared candidate universe
+
mapping
```

的结果。

所以：

## 不要再写

> consensus robustness

更适合改成：

> **candidate-universe sensitivity / peak-level call independence**

。

---

# 13. 0.243 Jaccard 应怎样解释？

你当前的 B universe 使用：

```text
70 kb stitching
```

，500 domains → 440 domains。

这会显著改变最终 genomic span。

因此低 genomic Jaccard 不一定代表 Breadth caller 不稳定，而是：

```text
共享坐标框架被主动改变
```

。

### 推荐 validation/02 改成两个层级

## Layer A：必须接近 1

```text
native PeakBroadCall concordance
```

这是算法设计真正声称 independent 的层。

## Layer B：敏感性结果

```text
mapped Breadth domain count
bp Jaccard
mapped/unmapped rate
gene-set concordance
```

明确承认：

> final domain representation can change with candidate-universe construction.

这样更科学。

---

# 14. validation/02 的脚本与保存 output 已经不同步

当前脚本后半段已经新增：

```text
Peak-set perturbation sensitivity
drop 20%
50 repetitions
```

。

但当前保存的：

```text
validation/02_breadth_super_validation_output.txt
```

并没有这部分输出。

因此：

```text
validation script
!=
saved validation result
```

。

### 必须

重新运行 validation/02 并覆盖 output，或者不要把 output 当正式证据。

---

# 15. synthetic benchmark：现在建议从正式 validation 叙事中撤掉

当前仓库中包含：

```text
03_synthetic_truth_benchmark.R
03_synthetic_truth_benchmark_revised.R
BENCHMARK0_*.csv
BENCHMARK0_*.pdf
BENCHMARK0_results.*
```

虽然 `.Rbuildignore` 已排除整个 `validation/`，所以不会进入 package tarball，

但从方法论文/repo维护角度：

> 项目现在已经决定不需要大型 synthetic benchmark。

### 推荐

两种方案均可：

## A. 最干净

删除：

```text
validation/03*
BENCHMARK0*
```

。

## B. 保留开发历史

移动到：

```text
validation/archive/synthetic_exploratory/
```

并明确：

```text
not part of manuscript validation
```

。

### 不建议

为了 synthetic balanced mode 继续扩 Breadth public API。

---

# 16. Quick Start 已修，但 vignette 中仍大量使用 tangent + raw

顶部 Quick Start 已正确变成：

```text
default elbow
per_group
default log convention
```

这是好事。

但后面的 main executable chunks 仍大量：

```r
method = "tangent"
log_transform = FALSE
```

例如：

```text
per-group
combined taxonomy
transitions
hockey plot
academic visualization
```

。

这与当前叙事：

```text
elbow = primary
tangent + raw = optional ROSE-style comparison
```

不一致。

### 推荐

只在：

```text
3.1 methods comparison
ROSE benchmark section
```

保留 tangent。

其余 main workflow 全部改成 default elbow。

---

# 17. README Core Module 2 首个代码块仍使用 global_consensus default

README：

```r
se <- call_super_domains(se, feature = "Intensity")
```

随后才解释：

```text
multi-condition studies 推荐 per_group
```

。

对于一个主要面向条件比较的工具，首个核心调用仍容易误导。

### 推荐

改为：

```r
se <- call_super_domains(
    se,
    feature = "Intensity",
    mode = "per_group",
    group_var = "Condition"
)
```

。

function default 不需要改。

---

# 18. README candidate-universe snippet 使用未定义 `preset`

当前 README：

```r
if (preset$stitch_distance > 0) {
    ...
}
```

但该代码块之前没有定义：

```r
preset
```

。

### 修法

加：

```r
preset <- get_mark_preset("H3K27ac")
```

或者直接写：

```r
stitch_distance = 12500
```

。

---

# 19. vignette 开头的一句分类仍不准确

当前写：

> two canonical measurement layers (Intensity, SignalDispersion) plus three native-geometry descriptors

但包当前 canonical Super taxonomy 是：

```text
Intensity
+
native Breadth
```

而：

```text
SignalDispersion
```

已经反复定义成：

```text
secondary architecture descriptor
```

。

### 推荐改为

> two continuous signal-derived assays (Intensity and SignalDispersion) plus three native-geometry descriptors

或者：

> one canonical signal-magnitude assay (Intensity), one secondary architecture assay (SignalDispersion), and three native-geometry descriptors.

不要再叫 SignalDispersion canonical measurement layer。

---

# 20. multi-condition universe 建议同步进入 vignette

README 已经正确说明：

```text
condition-specific consensus
→ union
→ stitch
```

。

vignette 主流程目前没有同等明确的 workflow。

这会导致：

```text
README 知道 condition-specific domains
vignette 用户却容易直接把所有 replicate 混起来 consensus
```

。

### 建议

vignette Module 1 增加一个简短 code block。

---

# 21. Test suite 很强，但还有两个“假覆盖”

当前约：

```text
163 test_that()
```

，整体覆盖已经明显足够。

但看到：

## A

```r
test_that("build_portrait_matrix enforces unique domain IDs", {
    ...
    expect_true(TRUE)
})
```

这个 test 实际并没有验证任何行为。

虽然正式函数里确实有 unique Domain_ID enforcement，

但这个 test 对未来 regression 没保护。

### 建议

用 tiny extdata BigWig 真正跑一次 duplicate domain names，并验证输出 rownames 唯一。

---

## B

```r
test_that("plot_peak_track returns ggplot", {
    skip("...")
})
```

是无条件 skip。

可以用 `inst/extdata` BigWig 增加真正 integration test。

这两个都不是 release blocker。

---

# 22. Release blocker：
# Git repository 当前并没有包含 epiPortrait14 的正式历史

这是当前最明确的发布工程风险。

当前 repo：

```text
HEAD = fcb1550
origin/main = fcb1550
```

但 working tree 有大量：

```text
modified
deleted
untracked
```

。

几乎整个 v1.0 redesign 都还没有成为 Git commit。

尤其：

```text
R/05a*
R/05b*
R/05c*
R/05d*
R/05e*
R/06_io_utils.R
R/07_qc_signal.R
R/08_visualization.R
R/09_object_contract.R
```

目前都是 Git untracked。

这意味着：

> **GitHub URL 当前仍可能展示旧的 4D package，而真正的 epiPortrait14 只存在于 working tree / zip 中。**

### 这是 Bioconductor / publication 前必须解决的

冻结后：

```bash
git status
git add -A
git commit
git push
```

然后确认：

```bash
git status --short
```

输出为空。

再从同一个 commit：

```bash
R CMD build .
R CMD check epiPortrait_0.99.x.tar.gz
BiocCheck(...)
```

。

---

# 23. 当前没有正式 R CMD check / BiocCheck log

包内没有看到：

```text
00check.log
BiocCheck output
```

。

虽然 validation benchmark 的 sessionInfo 证明：

```text
epiPortrait 0.99.0
R 4.5.1
```

曾经实际加载并运行主要 caller，

这不等于：

```text
R CMD check clean
BiocCheck clean
```

。

### 所以当前不能写

> submission ready

只能写：

> **algorithmically freeze-ready, pending clean package checks.**

---

# 24. H3K4me3 preset 5 kb

当前：

```r
H3K4me3:
stitch_distance = 5000
```

。

这仍然是之前已经知道的 sensitivity point。

结合当前正在看的 pan-cancer H3K4me3 broad-domain biology，

更建议真实 benchmark 中比较：

```text
0 kb
vs
5 kb
```

。

但：

> **不要现在因为这个再改核心 preset。**

先真实数据敏感性分析。

---

# 25. 当前最推荐的最小修复清单

为了避免 endless：

## 必修 1：发布工程

```text
把 v14 全部 commit/push
```

。

## 必修 2：Breadth optional path 二选一

### 推荐：

```text
删除 Breadth quantile_cutoff
```

。

如果保留，则必须修：

```text
0 eligible peaks
quantile tie / top-fraction wording
```

。

## 必修 3：Breadth bootstrap

```text
no_call resample
不能算 bootstrap success
```

。

## 必修 4：validation/02

```text
robustness → sensitivity
重新运行保存 output
不要把 Jaccard=0.243 写成 robust
```

。

## 必修 5：文档

```text
main workflow tangent/raw → elbow default
README undefined preset
vignette 增加 condition consensus → union
```

。

---

# 26. 可以不修 / 不建议再扩大的部分

当前不要继续动：

```text
Intensity definition
Breadth unique mapping
min_peak_overlap_fraction default
majority support rule
Uncertain semantics
Combined taxonomy
reference/pooled transition engine
annotation
RNA prioritization
BEDPE
mark-aware display
SignalDispersion formula
H3K27me3/H3K9me3 core pathway
```

。

这些已经足够稳定。

---

# 27. 本轮风险评级

## P0

### 没有新的核心算法 P0。

但是有一个**release P0**：

```text
current code not committed/pushed
```

。

---

## P1

1. Breadth bootstrap success-rate 错计 no_call。
2. Breadth quantile 0-eligible edge case。
3. Breadth quantile “top X%” 语义受 observed cutoff/ties 影响。
4. validation/02 Jaccard=0.243 不支持“domain robustness”措辞。
5. validation/02 script/output 不同步。
6. synthetic benchmark 与当前论文策略冲突，建议撤出正式 validation。
7. main vignette 仍有大量 tangent/raw old workflow。

---

## P2

1. README core first call 仍是 global_consensus default。
2. README candidate-universe snippet 中 `preset` 未定义。
3. vignette “canonical measurement layers” 对 SignalDispersion 措辞不准确。
4. vignette 未同步 multi-condition candidate-universe workflow。
5. 一个 `expect_true(TRUE)` test。
6. 一个无条件 skipped plot test。
7. H3K4me3 5 kb preset 需要真实数据 sensitivity，而不是代码重构。

---

# 28. 最终判断

与 `epiPortrait13` 相比：

> **+14 确实修掉了上一轮最关键的 provenance / input validation / documentation contract 问题。**

所以核心质量在提升。

但是 +14 同时又开始增加：

```text
Breadth top-fraction
Breadth bootstrap
synthetic benchmark
```

这已经有一点重新扩大 scope 的趋势。

结合当前项目方向，我建议：

> **不要继续往前加功能。把 Breadth quantile 这一新增面收回或修严谨，把 validation/02 解释正确，把 Git/release 状态收干净，然后冻结。**

---

# 29. 最终评分

| 模块 | 评价 |
|---|---:|
| Core Intensity logic | 9.3/10 |
| Canonical Breadth logic | 9.2/10 |
| Replicate support / Uncertain | 9.3/10 |
| Combined / transition logic | 9.1/10 |
| Annotation / object contract | 9.1/10 |
| Documentation positioning | 8.8/10 |
| Validation interpretation | 7.7/10 |
| Optional Breadth quantile/bootstrap | 7.4/10 |
| Release engineering | 6.5/10 until committed/check-clean |
| **Overall core scientific readiness** | **~9/10** |

---

# 30. 一句话结论

> **epiPortrait14 的核心算法已经可以冻结，没有必要再重新设计；当前剩余风险主要来自 freeze 后新增的 Breadth quantile/bootstrap 和 validation 解释，而不是 Intensity × Breadth 主框架本身。最值得做的是删减而不是继续扩展：撤出 synthetic benchmark、优先移除不必要的 Breadth top-fraction API，修 bootstrap/no-call 与 validation/02，同步文档，然后 commit/push 并跑完整 R CMD check + BiocCheck。**
