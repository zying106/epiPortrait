# epiPortrait12 冻结候选版深度复核
## 核心算法、Mark preset、Vignette、Examples 与 Bioconductor 提交准备

**日期：2026-08-11**  
**检查对象：`epiPortrait12.zip`**  
**评估视角：生物信息学 / 表观遗传学 / Bioconductor reviewer**  
**方式：源码静态深审。当前执行环境没有 R，因此本轮未实际运行 `testthat`、`R CMD build/check` 或 `BiocCheck`。**

---

# 1. 最终结论

`epiPortrait12` 已经进入真正的 **release-freeze / submission-cleanup** 阶段。

上一版 `epiPortrait11` 的主要冻结项已经确认修复：

- Z-score 已从 `normalize_portrait()` 移除，不再覆盖 canonical `Intensity`；
- `plot_domain_landscape()` 已支持并尊重 `group_var`；
- `filter_promoter_peaks()` 已加入 TxDb seqlevel compatibility hard guard；
- `filter_blacklist()` 已加入 zero-shared seqlevel hard guard；
- promoter window 参数已验证；
- `validate_epiportrait_object()` 已检查 replicate-call matrix contract；
- `stitch_epi_peaks()` 的 exact-distance off-by-one 已修；
- `get_replicate_calls(group=NULL)` 已恢复原始 sample order；
- `Rplots.pdf` 已删除，并写入 `.Rbuildignore`；
- Intensity/Breadth replicate accessor 已闭环；
- annotation / BEDPE / RNA-seq 模块没有看到新的结构性回归。

因此：

> **当前不建议再增加任何新的生物学模块。**

这一版真正值得收尾的是：

1. mark preset 的 stitching sentinel 语义；
2. vignette / man page 与当前 API 的文档同步；
3. Bioconductor 对可执行 vignette / runnable examples 的要求；
4. source tree 清理；
5. synthetic BED 坐标 convention；
6. 最终真实 `R CMD check / BiocCheck`。

---

# 2. 核心算法：可以冻结

当前没有看到需要重新设计：

```text
Intensity
Breadth
SignalDispersion
replicate support
Combined taxonomy
condition transitions
domain annotation
BEDPE evidence
RNA expression prioritization
```

的结构性问题。

尤其：

## Intensity

仍是：

```text
replicate-specific rank
→ replicate-specific cutoff
→ replicate-specific Super/Typical/no-call
→ replicate support aggregation
→ group call
```

不是 group-mean-only calling。

## Breadth

仍是：

```text
replicate native PeakWidth distribution
→ replicate-specific width inflection
→ native Broad/Typical/no-call
→ unique native-peak-to-shared-domain mapping
→ replicate support
→ Breadth-Super
```

并且现在可由：

```r
get_replicate_calls()
```

完整审计。

这两条 canonical axes 已经对称闭环。

---

# 3. +11 P1：Z-score canonical Intensity 风险已正确修复

当前：

```r
normalize_portrait()
```

只允许：

```text
None
TotalSignal
TMM
Quantile
```

不再允许：

```text
Z-score
```

并明确说明：

> row-wise Z-score 是 display/clustering transform，会破坏 Super calling 所需的 cross-domain magnitude ranking。

这是正确修法。

相关测试也已加入：

```r
expect_error(
    normalize_portrait(example_se, method = "Z-score")
)
```

因此上一版最重要的 normalization safety 问题已关闭。

---

# 4. +11 P1：plot_domain_landscape(group_var) 已修复

现在：

```r
plot_domain_landscape(
    se,
    group = "A",
    group_var = "Arm"
)
```

会把：

```text
group_var
```

传给：

```r
.epi_feature_vector()
```

并实际：

```r
meta[[group_var]]
```

取样本。

相应 test 已存在。

正确。

---

# 5. +11 P1：promoter / blacklist seqlevel guards 已修复

## promoter

现在：

```r
filter_promoter_peaks()
```

会调用：

```r
.check_seqlevel_compatibility(
    gr,
    txdb = txdb,
    enforce = TRUE
)
```

所以：

```text
chr1 vs 1
```

不会再被错误解释为：

```text
0 promoter overlaps
```

。

## blacklist

现在：

```text
zero shared seqlevels
→ stop
partial input seqlevels missing from blacklist
→ warning
```

方向正确。

---

# 6. +11 P1：object validator 已增加 replicate matrix invariant

现在如果：

```r
metadata(se)$superdomain_calls[[feature]]$groups[[g]]$
    replicate_call_matrix
```

存在，会检查：

```text
nrow(matrix) == nrow(se)
rownames(matrix) == rownames(se)
colnames(matrix) ⊆ colnames(se)
```

这是很好的对象 contract 防护。

---

# 7. +11 P2：stitch exact-distance off-by-one 已修复

现在：

```r
GenomicRanges::reduce(
    gr,
    min.gapwidth = stitch_distance + 1
)
```

因此文档：

> maximum distance = 12500 bp

真正表示：

```text
gap <= 12500
```

都会 stitch。

正确。

---

# 8. 当前最值得修的科学/API语义问题：
# `stitch_distance = 0` 仍不等于 “do not re-stitch”

这是 `epiPortrait12` 当前最重要的剩余设计一致性问题之一。

`get_mark_preset()` 文档对 broad repressive marks 写：

```text
stitch_distance = 0
→ do not re-stitch caller-defined domains
```

H3K27me3/H3K9me3 preset 也确实返回：

```r
stitch_distance = 0L
```

但是 vignette 仍然执行：

```r
domains <- stitch_epi_peaks(
    domains,
    stitch_distance = preset$stitch_distance
)
```

即：

```r
stitch_epi_peaks(domains, stitch_distance = 0)
```

---

# 9. 为什么 0 不是 no-op

当前：

```r
reduce(
    gr,
    min.gapwidth = 0 + 1
)
```

仍然会 `reduce()`。

因此：

```text
overlapping domains
adjacent domains
```

仍可能被 merge。

所以：

```text
stitch_distance = 0
```

的真实含义是：

> 不跨正距离 gap stitching，但仍 reduce overlapping/adjacent ranges。

它不等于：

> 原样保留 upstream caller 的 domains。

这与 preset 文档：

```text
0 / NULL = do not re-stitch
```

不完全一致。

---

# 10. 推荐最小修法：不要改 stitch_epi_peaks() 本身

不建议把：

```r
stitch_epi_peaks(gr, 0)
```

改成 return(gr)。

因为从一个通用 interval function 的语义看：

```text
gap <= 0
```

合并 overlap / adjacency 本身是合理的。

更好的做法是：

## workflow / preset 层将 0 当 sentinel

```r
if (!is.null(preset$stitch_distance) &&
    preset$stitch_distance > 0) {
    domains <- stitch_epi_peaks(
        domains,
        stitch_distance = preset$stitch_distance
    )
}
```

对于：

```text
H3K27me3
H3K9me3
```

：

```text
stitch_distance = 0
→ 根本不调用 stitch_epi_peaks()
```

这样才真正保留 upstream broad-domain caller geometry。

---

# 11. 同一问题：unknown/generic mark 仍然默认 12.5 kb stitching

当前 unknown mark：

```r
get_mark_preset("H4K20me3")
```

warning：

```text
using a generic preset
(no mark-specific stitching or terminology)
```

但实际返回：

```r
stitch_distance = 12500L
```

这是自相矛盾的。

更重要的是：

```text
12.5 kb
```

本质上是 active enhancer / ROSE-like 的强 mark-specific assumption。

对：

```text
H4K20me3
H3K36me3
其他未知 mark
```

未必合理。

---

# 12. 推荐 generic preset

建议：

```r
stitch_distance = 0L
```

同时文档：

```text
0 means no automatic re-stitching at the preset-workflow level.
```

用户如果知道自己的 mark 应该 stitch：

```r
stitch_epi_peaks(..., stitch_distance = user_value)
```

显式指定即可。

这是比 generic 12.5 kb 更安全的 fallback。

---

# 13. P1 文档残留：vignette 仍列出已删除的 Z-score

当前 vignette normalization 表仍有：

```text
Z-score | For clustering/heatmaps only
```

但是：

```r
normalize_portrait(method = "Z-score")
```

现在已经明确报错。

因此用户会看到：

```text
vignette 说有
代码说没有
```

。

必须删除这一行。

如果需要说明，可以在表后写：

> Row-wise Z-score is intentionally not an analysis normalization method in epiPortrait; plotting/PCA layers perform display scaling internally.

---

# 14. P1 文档残留：package-level documentation 仍说 Normalize SignalDispersion

`R/00_package.R`：

```text
normalize_portrait — Normalize Intensity and SignalDispersion assays
```

对应：

```text
man/epiPortrait.Rd
```

也仍然是相同表述。

但真实实现已经明确：

```text
Intensity adjusted
SignalDispersion left unchanged
```

。

建议改为：

```text
normalize_portrait — Normalize Intensity while preserving SignalDispersion
```

然后重新：

```r
roxygen2::roxygenise()
```

同步 `.Rd`。

---

# 15. 当前最值得重视的 Bioconductor 提交问题：
# vignette 还有多个非 installation 的 `eval = FALSE`

当前：

```text
install            eval = FALSE
custom             eval = FALSE
annotate           eval = FALSE
save-analysis      eval = FALSE
preset-workflow    eval = FALSE
```

其中 installation chunk 使用：

```text
eval = FALSE
```

是合理的。

但当前 Bioconductor development/review documentation 明确强调：

- vignette 应包含 executable code；
- shortcuts such as `eval=FALSE` generally undermine vignette reproducibility；
- 除 installation instructions 等有合理理由的情况外，应尽量避免。

因此这一项很值得在 submission 前处理。

官方文档：
https://contributions.bioconductor.org/docs.html

---

# 16. 推荐处理各个 vignette chunk

## installation

保留：

```r
eval = FALSE
```

这是官方明确允许/要求的场景。

---

## custom feature workflow

当前依赖：

```text
samples
macro_domains
BigWig I/O
```

而 Windows vignette fallback 不一定定义这些对象。

不建议硬执行。

更干净：

把：

````
```{r custom, eval = FALSE}
...
```
````

改成普通展示代码：

````
```r
...
```
````

并在文字中写：

> Example pattern for real BigWig inputs.

它就不再是“声称执行但关闭执行”的 knitr chunk。

---

## annotate

这是一个真实 exported module。

如果 build environment 的 TxDb Suggests 已安装，建议直接执行：

```r
se <- annotate_epi_domains(
    example_se,
    genome = "hg38"
)
```

不要 `eval=FALSE`。

如果担心 TxDb runtime，可构建一个 tiny custom TxDb 测试，但不建议为了 vignette过度复杂。

---

## save-analysis

可以改成真正可执行且不污染工作目录：

```r
f <- tempfile(fileext = ".rds")
saveRDS(se_viz, f)
se2 <- readRDS(f)

outdir <- tempfile("epiPortrait_results_")
dir.create(outdir)
export_epiportrait_results(se2, outdir = outdir)
```

必要时：

```r
unlink(...)
```

。

这样可以取消 `eval=FALSE`。

---

## preset-workflow

需要真实外部：

```text
H3K27me3_domains.bed
```

，因此建议改为普通 `r` fenced example，而不是 knitr eval-false chunk。

同时修正：

```r
if (preset$stitch_distance > 0) {
    domains <- stitch_epi_peaks(...)
}
```

。

---

# 17. Man page 中仍有 7 个 `\donttest{}`

当前：

```text
annotate_epi_domains
build_portrait_matrix
check_signal_compatibility
filter_promoter_peaks
get_domain_genes
plot_peak_track
save_epiportrait_figure
```

使用：

```text
\donttest{}
```

。

Bioconductor 当前开发文档明确说：

> `donttest` / `dontrun` are discouraged and generally not allowed unless justified.

因此，如果第二轮 review 之前已经遇到“未执行 chunk / example”类意见，这里建议顺手继续清理。

---

# 18. 哪些 `donttest` 很容易直接去掉

## save_epiportrait_figure

当前：

```r
p <- plot_portrait_pca(...)
save_epiportrait_figure(
    p,
    file = tempfile(fileext = ".png")
)
```

本身就是安全的 runnable example。

建议直接去掉 `donttest`。

---

## build_portrait_matrix / check_signal_compatibility / plot_peak_track

已经：

```r
if (.Platform$OS.type != "windows") {
    ...
}
```

并使用包内 tiny BigWigs。

可以考虑直接取消外层：

```text
donttest
```

。

R example 本身会执行；Windows 只是不进入 BigWig branch。

---

## annotation / promoter / candidate genes

可以：

```r
if (requireNamespace(
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    quietly = TRUE
)) {
    ...
}
```

然后取消 `donttest`。

---

# 19. P1 source-tree hygiene：vignettes 中仍有生成的 HTML

当前源码树实际包含：

```text
vignettes/epiPortrait.html
```

虽然：

```text
.Rbuildignore
```

已经排除：

```text
^vignettes/.*\.html$
```

，但 Bioconductor documentation 明确建议：

> `vignettes/` should contain only source vignette files and necessary static images; generated vignette products/intermediate files should not remain there.

所以建议直接从 repository 删除：

```text
vignettes/epiPortrait.html
```

。

---

# 20. 同样建议从 repository 删除开发/生成残留

当前源码树还包含：

```text
Meta/
README.html
.DS_Store
inst/.DS_Store
validation/.DS_Store
man/.DS_Store
```

虽然多数已被 `.Rbuildignore` 过滤，

但 release/source repository 更干净的做法仍是：

```text
直接删除
```

，而不是依赖 ignore。

`Meta/` 尤其像 installed-package 生成内容，不应作为 source package 内容维护。

---

# 21. zip 顶层还有 `__MACOSX`

上传的 zip 顶层存在：

```text
__MACOSX/
```

这不是 epiPortrait package directory 内的文件，

所以如果最终提交的是：

```r
R CMD build epiPortrait
```

产生的：

```text
epiPortrait_0.99.x.tar.gz
```

不会进入 source tarball。

因此：

> 不要把当前开发 zip 当提交 artifact。

最终只提交/推送标准 R package source，并以 `R CMD build` 生成 tarball。

---

# 22. P1/P2：synthetic BED generator 有 BED 0-based off-by-one

`data-raw/make_tiny_data.R` 当前：

```r
write_bed <- function(gr, path) {
    data.frame(
        chr   = ...,
        start = start(gr),
        end   = end(gr)
    )
}
```

但是 standard BED 是：

```text
0-based start
half-open end
```

而 GRanges 是：

```text
1-based closed
```

。

正确转换通常应该：

```r
start = start(gr) - 1L
end   = end(gr)
```

。

当前 tiny BED：

```text
peaks.bed
C1_peaks.bed
...
```

因此与原始 GRanges generator 设计相比会在 import 后发生 1-bp start shift。

---

# 23. 这个 off-by-one 会影响主算法吗？

不会影响真实用户数据或核心算法。

它只影响：

```text
package synthetic example/test BEDs
```

。

而且 1 bp 对当前几 kb domain benchmark 几乎不会改变结果。

但从：

```text
package correctness
coordinate convention
reviewer impression
```

来说，最好修干净。

推荐：

```r
start = GenomicRanges::start(gr) - 1L
```

然后重新生成：

```text
inst/extdata/peaks.bed
C1/C2/T1/T2_peaks.bed
```

并重新跑 tests/vignette。

---

# 24. P2：vignette 无 peak_path 的描述有一处逻辑不顺

当前写：

> Without peak_path, native breadth geometry is NA and Breadth-Super calling is unavailable ...; the imported native peaks are stored in metadata(se)$native_peaks ...

但：

```text
without peak_path
```

本来就没有 native peaks 被 import。

应该改成：

> Without `peak_path`, native breadth geometry and canonical Breadth-Super calling are unavailable. When `peak_path` is supplied, imported native peaks are retained in `metadata(se)$native_peaks` for auditable peak-level breadth calling.

另外不建议写：

```text
only Intensity / SignalDispersion super-domain calling
```

因为 SignalDispersion 是 secondary/exploratory，不是 canonical Super axis。

更准确：

> Intensity-Super remains available; SignalDispersion can still be explored as a secondary architecture descriptor.

---

# 25. P2：vignette 仍把 alias `consensus` 当推荐名称

前文已经正确使用：

```text
global_consensus
```

但 Parameter guidance 又写：

```r
mode = "consensus"
```

。

建议统一：

```r
mode = "global_consensus"
```

，并只在第一次介绍时注明：

```text
alias = "consensus"
```

。

---

# 26. P2：Troubleshooting 的 quantile_cutoff 建议应限定 Intensity

当前：

```text
No reliable inflection
→ Use quantile_cutoff
```

但 canonical Breadth：

```r
feature = "Breadth"
```

明确：

```text
quantile_cutoff is not supported
```

。

因此建议：

```text
Intensity: consider an explicit quantile_cutoff or inspect the rank curve.
Breadth: inspect eligible native peaks / width distribution / caller quality;
no-call should remain no-call rather than be forced by quantile.
```

这样不会让 Breadth 用户得到错误建议。

---

# 27. P2：package-level “occupied breadth” 表述可再精确一点

README / package docs 多处写：

```text
sample-specific occupied breadth
```

而 canonical Breadth-Super 真正调用依据是：

```text
single native peak PeakWidth
```

不是：

```text
NativeOccupiedWidth
```

。

因此从论文/审稿语言严谨度看，更建议：

> sample-specific native peak breadth geometry

而不是：

> occupied breadth

。

`NativeOccupiedWidth` 可以继续作为 domain-level descriptor。

这能避免 reviewer 误解 Breadth caller 是：

```text
sum of occupied widths
```

。

---

# 28. P2：README Quick Start 默认 global_consensus 与主场景略不一致

当前 Quick Start：

```r
call_super_domains(se, feature = "Intensity")
call_super_domains(se, feature = "Breadth")
combine_superdomain_calls(se)
```

使用默认：

```text
global_consensus
```

。

但 sample data 本身：

```text
Control
Treatment
```

而 README 又明确推荐：

```text
per_group
```

作为 multi-condition main analysis。

这不是算法 bug。

建议：

## 保持 function default 不变

避免 release 前改 API default。

但是 README Quick Start 可改成：

```r
call_super_domains(
    se,
    feature = "Intensity",
    mode = "per_group",
    group_var = "Condition"
)
```

Breadth 同样。

然后：

```r
combine_superdomain_calls(
    se,
    group_var = "Condition"
)
```

更贴合 epiPortrait 的主要 biological use case。

---

# 29. P2：NEWS 仍保留同一 0.99.0 开发期的 Z-score “Added” 记录

NEWS 顶部现在正确写：

```text
Z-score normalization REMOVED
```

但较老的同一：

```text
Changes in version 0.99.0
```

内容仍写：

```text
Added normalize_portrait() with ... Z-score ...
```

由于：

```text
0.99.0
```

尚是同一个 pre-release version，

release NEWS 最好反映最终用户可见 API，而不是完整开发过程。

建议删除或改成：

> Initial normalization API was added; row-wise Z-score was subsequently removed before release because it is incompatible with canonical Intensity ranking.

---

# 30. P2：TMM/Quantile hard dependency可以以后再优化，不建议当前再动

`limma` 当前在 Imports。

技术上 Quantile 是 optional method，可以将 limma 移到 Suggests + `requireNamespace()`。

但：

```text
现在不是必要修改
```

。

freeze 阶段不建议为了减少一个依赖再扩大改动面。

---

# 31. 当前 tests 覆盖已经明显增强

看到针对以下问题的 test：

```text
Z-score rejected
Breadth accessor
custom group_var accessor
plot_domain_landscape custom group_var
promoter chr/1 mismatch
blacklist chr/1 mismatch
mm9 blacklist
domain/TxDb mismatch
RNA wide matching
VST threshold
candidate rank
BEDPE evidence
domain_gene_links export
```

方向很好。

下一步比继续增加 tests 更重要的是：

```text
实际执行整套 tests
```

。

---

# 32. 当前 Bioconductor review 风险排序

## P1 — 建议 freeze 前修

### 1. Mark preset stitching sentinel

```text
H3K27me3/H3K9me3:
0 必须意味着 workflow 层不调用 stitch

generic unknown mark:
不要默认 12.5 kb
```

### 2. 删除 vignette stale Z-score

### 3. package-level normalization 文档同步

### 4. 减少非 installation 的 eval=FALSE chunks

### 5. 尽量移除可运行 example 的 donttest

### 6. 删除 vignettes/epiPortrait.html 等 source-tree generated artifacts

---

## P2 — 成本低，建议一起做

7. synthetic BED start 改为 `start - 1`;
8. 无 peak_path vignette 表述修正；
9. `consensus` 统一为 `global_consensus`;
10. Troubleshooting quantile 仅限定 Intensity；
11. “occupied breadth” 改为 “native peak breadth geometry”；
12. Quick Start 用 explicit `per_group`;
13. NEWS 清理 Z-score 开发历史矛盾；
14. 删除 `.DS_Store` / `Meta/` / `README.html` 等生成残留。

---

# 33. 不建议再增加的新功能

当前不建议加入：

```text
assess_call_stability()
Dispersion-Super
Height-Super
TargetScore
DESeq2 pipeline
GSEA engine
motif engine
ABC model
Hi-C normalization
loop caller
更多 annotation API
```

。

目前这些功能的边际收益都低于：

```text
validation
benchmark
submission hygiene
```

。

---

# 34. 冻结后真正值得投入的工作

## H3K4me3

验证：

```text
Breadth-Super
```

对 broad H3K4me3 的：

```text
replicate reproducibility
native geometry robustness
consensus robustness
biological enrichment
```

。

## H3K27ac

重点比较：

```text
ROSE-style 1D ranking
vs
epiPortrait Intensity / Breadth / Dual decomposition
```

。

不是只比较 overlap。

## Intensity robustness

benchmark script：

```text
leave-one-replicate-out
rank correlation
cutoff variability
Intensity-Super retention
```

即可，不必变成 public API。

## Consensus robustness

继续使用：

```text
true bp Jaccard
class concordance
transition concordance
```

。

---

# 35. 当前评分

| 模块 | epiPortrait12 静态评价 |
|---|---:|
| Intensity calling | 9.5/10 |
| Breadth calling | 9.4/10 |
| replicate support | 9.5/10 |
| replicate accessor | 9.3/10 |
| transitions | 9.1/10 |
| annotation | 9.2/10 |
| BEDPE/RNA optional evidence | 8.9/10 |
| genome/resource compatibility | 9.2/10 |
| normalization safety | 9.3/10 |
| visualization API | 9.0/10 |
| mark preset semantics | 8.2/10 |
| documentation consistency | 8.4/10 |
| Bioconductor packaging hygiene | 8.3/10 |
| core algorithm freeze readiness | **9.5/10** |
| overall submission freeze readiness | **约 9.1/10** |

---

# 36. 最终判断

`epiPortrait12` 已经不是“还要不要继续优化算法”的阶段。

从源码静态审查看：

> **核心算法可以冻结。**

真正还值得做的一轮修改是：

```text
mark preset 0/no-stitch 语义
→ generic preset 不再 12.5 kb
→ vignette/man 文档同步
→ eval=FALSE / donttest 清理
→ source-tree hygiene
→ synthetic BED coordinate fix
→ roxygenise
→ R CMD build
→ R CMD check
→ BiocCheck
→ benchmark
```

这轮完成后，如果：

```text
R CMD check = clean
BiocCheck = clean / only clearly justified notes
tests pass
vignette fully builds
```

就不建议继续修改 package architecture。

下一步应该转入论文 benchmark 和发布冻结，而不是继续添加功能。
