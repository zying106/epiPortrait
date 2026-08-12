# epiPortrait12 内嵌数据设计与优化总报告
## Example data / extdata / replicate / grouping / simulation truth / Bioconductor review

**版本：epiPortrait12**  
**日期：2026-08-11**  
**评估视角：生物信息学、表观遗传学、Bioconductor reviewer**

---

# 1. 总体结论

当前 epiPortrait 的包内数据设计总体上是合理的，而且已经足够支撑 v1.0。

目前实际上形成了一个很好的“两层 example-data 架构”：

```text
data/example_se.rda
        ↓
快速、预计算、适合算法和可视化示例

inst/extdata/*.bw + *.bed
        ↓
真实文件 I/O、BigWig/BED integration test
```

这是比单纯提供一个矩阵更完整、又比内嵌大型真实 ChIP-seq 数据更轻量的方案。

当前：

```text
data/ + inst/extdata
```

总大小约数百 KB，属于非常适合 Bioconductor source package 的规模。

因此核心建议是：

> **不要再增加大型真实数据，也不要增加多个冗余 example object。**

真正需要做的是：

```text
修正 synthetic 数据坐标
规范 BED coordinate convention
增强 replicate/group metadata
让 synthetic replicate disagreement 更有教学价值
规范 simulation truth / provenance
```

---

# 2. 推荐保留的最终数据层级

## 2.1 `data/`

只保留一个主要 processed example：

```text
example_se.rda
```

用途：

```text
Intensity calling
Breadth calling
replicate support
Combined class
condition transition
annotation
candidate gene extraction
PCA
plots
get_replicate_calls()
```

不建议增加：

```text
example_se2
example_breadth_se
example_annotation_se
example_2rep_se
example_3group_se
```

一个主 object 最清楚。

---

## 2.2 `inst/extdata/`

继续保留：

```text
tiny BigWigs
tiny candidate BED
per-sample native peak BEDs
blacklists
README/provenance
```

例如：

```text
C1.bw
C2.bw
T1.bw
T2.bw

peaks.bed

C1_peaks.bed
C2_peaks.bed
T1_peaks.bed
T2_peaks.bed

hg19-blacklist.v2.bed
hg38-blacklist.v2.bed
mm9-blacklist.v2.bed
mm10-blacklist.v2.bed
```

用途：

```text
BigWig import
BED import
build_portrait_matrix()
signal extraction
SignalDispersion
native peak geometry
blacklist filtering
I/O integration tests
```

---

# 3. P0：`example_se` 的 chr21 坐标必须修正

当前生成逻辑：

```r
n_dom <- 500

start = seq(
    10000,
    by = 200000,
    length.out = n_dom
)
```

最后一个 domain 起点约：

```text
99.8 Mb
```

但 hg38 chr21 长度约：

```text
46.7 Mb
```

因此大量 synthetic domains 已经超出 chr21。

这属于真实 coordinate consistency 问题。

---

# 4. 推荐修法

仍然使用 chr21 即可。

例如：

```r
gr <- GRanges(
    "chr21",
    IRanges(
        start = seq(
            1e6,
            by = 80000,
            length.out = 500
        ),
        width = widths
    ),
    seqinfo = Seqinfo(
        seqnames = "chr21",
        seqlengths = 46709983,
        genome = "hg38"
    )
)
```

这样最后一个 domain 仍位于 chr21 内。

---

# 5. 建议增加生成时 invariant

在保存 `example_se` 前：

```r
stopifnot(
    all(
        end(gr) <=
        GenomeInfoDb::seqlengths(gr)[
            as.character(seqnames(gr))
        ]
    )
)
```

避免以后修改生成参数时再次出现越界。

---

# 6. P1：tiny BED 必须使用标准 BED 坐标

当前 generator：

```r
start = start(gr)
end   = end(gr)
```

但：

```text
GRanges = 1-based closed
BED     = 0-based half-open
```

因此应改为：

```r
start = start(gr) - 1L
end   = end(gr)
```

---

# 7. 推荐 BED writer

```r
write_bed <- function(gr, path) {
    write.table(
        data.frame(
            chr   = as.character(seqnames(gr)),
            start = start(gr) - 1L,
            end   = end(gr)
        ),
        file = path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
}
```

然后重新生成：

```text
peaks.bed
C1_peaks.bed
C2_peaks.bed
T1_peaks.bed
T2_peaks.bed
```

---

# 8. 当前 replicate / 分组设计是否合理？

## 结论：非常合理，不建议改

当前：

```text
Control_1
Control_2
Control_3

Treatment_1
Treatment_2
Treatment_3
```

即：

```text
2 conditions × 3 biological replicates
```

正好适合 epiPortrait。

---

# 9. 为什么 3 replicates/group 很合适

因为 epiPortrait 的核心特点之一就是：

```text
replicate-aware calling
```

3 replicates 可以自然展示：

```text
3/3 → unanimous

2/3 → majority-supported

1/3 → unsupported

0/3 → unsupported
```

尤其：

```r
support_rule = "majority"
```

时：

```text
2/3
```

非常直观。

如果只有 2 replicates：

```text
2/2
1/2
0/2
```

majority 的教学价值反而弱很多。

---

# 10. 两个 condition 也已经足够

目前：

```text
Control
Treatment
```

就能完成完整示例：

```text
replicate call
↓
group support
↓
Control phenotype
Treatment phenotype
↓
transition
```

例如：

```text
Control:
Broad / Broad / Broad
→ Breadth-Super

Treatment:
Typical / Typical / Broad
→ Typical

Transition:
Breadth-Super → Typical
```

因此没有必要把主 example 改成：

```text
Control
Early
Late
Resistance
```

多组场景放到 unit tests 即可。

---

# 11. 建议增加 `BioReplicate` 字段

当前 `colData` 已有：

```text
SampleID
Condition
```

建议增加：

```text
BioReplicate
```

例如：

| SampleID | Condition | BioReplicate |
|---|---|---:|
| Control_1 | Control | 1 |
| Control_2 | Control | 2 |
| Control_3 | Control | 3 |
| Treatment_1 | Treatment | 1 |
| Treatment_2 | Treatment | 2 |
| Treatment_3 | Treatment | 3 |

推荐：

```r
BioReplicate = rep(1:3, 2)
```

---

# 12. 必须注明：BioReplicate 不表示 paired design

例如：

```text
Control_1
Treatment_1
```

不能因为都叫 replicate 1 就默认它们配对。

如果以后是真正 paired experiment，应增加：

```text
PairID
```

而不是靠 BioReplicate 编号推断。

---

# 13. 推荐让 synthetic data 故意包含不同 replicate-support 模式

这是当前 example data 最值得增强的地方。

不要让所有 intended Super domains 都完美：

```text
3/3
```

。

应刻意设计一部分：

```text
3/3
2/3
1/3
NA/no-call
```

这样可以真正展示：

```text
majority
all
fraction
min_valid_replicates
get_replicate_calls()
```

。

---

# 14. 推荐的 synthetic replicate evidence 结构

## 强 phenotype

多数 intended Super domains：

```text
Rep1 Super
Rep2 Super
Rep3 Super
→ 3/3
```

---

## 边界 phenotype

部分：

```text
Rep1 Super
Rep2 Super
Rep3 Typical
→ 2/3
```

---

## 弱/不稳定 phenotype

少量：

```text
Rep1 Super
Rep2 Typical
Rep3 Typical
→ 1/3
```

---

# 15. Breadth 也应有同样模式

例如：

```text
Broad
Broad
Broad
```

和：

```text
Broad
Broad
Typical
```

以及少量：

```text
Broad
Typical
Typical
```

。

这样：

```r
get_replicate_calls(
    example_se,
    feature = "Breadth",
    group = "Treatment"
)
```

才真正具有教学价值。

---

# 16. 建议保留少量 no-call / NA

可以设计约：

```text
5–10 domains
```

出现：

```text
Broad
Broad
NA
```

或：

```text
Broad
Typical
NA
```

。

用于演示：

```text
valid replicate count
min_valid_replicates
Uncertain
```

。

不建议让大量 domains 缺失，否则示例会变得嘈杂。

---

# 17. 推荐的整体比例

不需要精确固定，但可以大致：

```text
70–80% intended strong phenotype
→ 3/3

15–25%
→ 2/3

少量
→ 1/3 / no-call
```

这样 synthetic example 既稳定，又能展示 replicate-aware 边界。

---

# 18. n=2 / 多组 / unbalanced 不需要新 example object

这些情形应该放在 unit tests。

## 2 groups × 2 replicates

直接 subset：

```r
se2 <- example_se[, c(
    "Control_1",
    "Control_2",
    "Treatment_1",
    "Treatment_2"
)]
```

---

## 3 groups × 2 replicates

可以临时：

```r
colData(se)$Group3 <- c(
    "A", "A",
    "B", "B",
    "C", "C"
)
```

。

---

## unbalanced

例如：

```text
Control n=3
Treatment n=2
```

直接删掉：

```text
Treatment_3
```

测试即可。

---

# 19. 推荐测试矩阵

unit tests 至少覆盖：

```text
2 groups × 2 replicates
2 groups × 3 replicates
3 groups × 2 replicates
unbalanced 3 vs 2
custom group_var
missing replicate evidence
no-call
all/majority/fraction
```

。

主内嵌 example 仍保持：

```text
2 × 3 balanced
```

。

---

# 20. `True_Class` 当前语义建议修改

目前 synthetic object 有：

```text
True_Class
```

例如：

```text
Typical
Intensity-Super
Breadth-Super
Dual-Super
```

但 simulation 中：

```text
Breadth phenotype
```

可能是 condition-specific 的。

例如：

```text
Control:
broad

Treatment:
contracted
```

那么一个静态：

```text
True_Class = Breadth-Super
```

容易误导。

---

# 21. 推荐最小方案：改成 `Simulation_Scenario`

例如：

```text
Typical
Intensity-high
Breadth-contraction
Dual-with-breadth-contraction
```

。

这个字段表示：

> synthetic domain 的生成机制

而不是：

> classifier 在所有 condition 下都应该输出同一个 class。

---

# 22. 如果以后做 synthetic accuracy benchmark

再增加：

```text
True_Intensity_Class__Control
True_Intensity_Class__Treatment

True_Breadth_Class__Control
True_Breadth_Class__Treatment

True_Combined_Class__Control
True_Combined_Class__Treatment
```

这样才适合计算：

```text
sensitivity
specificity
precision
recall
class concordance
```

。

但 v1.0 不是必须。

---

# 23. 建议增加 simulation provenance

例如：

```r
metadata(example_se)$simulation_provenance <- list(
    dataset = "synthetic",
    seed = 7L,
    genome = "hg38",
    chromosome = "chr21",
    n_domains = 500L,
    n_control = 3L,
    n_treatment = 3L,
    generator = "data-raw/make_example_se.R"
)
```

同时：

```r
GenomeInfoDb::genome(
    rowRanges(example_se)
) <- "hg38"
```

。

---

# 24. `bw_path` / `peak_path` 建议使用 NA

对于 precomputed synthetic `example_se`：

```r
bw_path = ""
peak_path = ""
```

建议改：

```r
bw_path = NA_character_
peak_path = NA_character_
```

。

因为：

```text
空字符串
```

像是“漏填路径”。

而：

```text
NA
```

明确表示：

> 该 processed example object 没有真实 backing file。

---

# 25. tiny BigWig 的 domain 数量是否够？

## 够

tiny BigWig 的用途不是展示真实 hockey-stick curve，而是：

```text
I/O integration
signal extraction
build_portrait_matrix
BigWig compatibility
```

。

所以 3 个左右 tiny domains 就够。

---

# 26. 不要让 tiny BigWig 承担 rank/calling demonstration

真正展示：

```text
Intensity-Super
Breadth-Super
Combined class
transition
```

继续使用：

```text
500-domain example_se
```

。

这是正确分工。

---

# 27. 是否增加 BEDPE / RNA-seq 内嵌数据？

## 当前不是必须

已有 synthetic unit tests 已能覆盖：

```text
BEDPE symmetric contacts
multiple promoters
support count
invalid coordinates
RNA wide
RNA long
replicate aggregation
sample matching
gene ID matching
candidate prioritization
```

。

因此 correctness coverage 已经有。

---

# 28. 只有一个情况下建议增加 BEDPE/RNA tiny fixture

如果你决定让 vignette 真正执行：

```text
annotate_epi_domains()
+ BEDPE
+ expression
+ get_domain_genes()
```

那么可以加两个极小文件：

```text
example_contacts.bedpe
example_expression.tsv
```

。

每个几行即可。

否则：

> 不要为了“数据更丰富”而加。

---

# 29. 不建议增加 tiny TxDb 文件

annotation unit tests 现在可以动态：

```r
txdbmaker::makeTxDbFromGRanges()
```

生成 tiny TxDb。

这比在：

```text
inst/extdata
```

永久存一个 SQLite/TxDb artifact 更干净。

所以暂时不建议添加。

---

# 30. blacklist 数据当前设计合理

当前：

```text
hg19
hg38
mm9
mm10
```

blacklist 总大小很小。

而且：

```text
filter_blacklist()
```

直接需要。

因此建议继续内嵌。

---

# 31. extdata README 建议补齐

建议列出：

```text
C1_peaks.bed
C2_peaks.bed
T1_peaks.bed
T2_peaks.bed
```

并明确：

```text
BED coordinate convention:
0-based, half-open
```

。

同时对 synthetic files 说明：

```text
generated for package examples/testing
not biological reference data
```

。

---

# 32. 当前不建议增加真实大型数据

不建议把：

```text
ENCODE H3K27ac
real H3K4me3
real HiChIP
real RNA-seq
```

直接放 source package。

真实 benchmark 数据应通过：

```text
GEO
ENCODE
ExperimentHub
Zenodo
Figshare
benchmark repository
```

提供。

---

# 33. 为什么不建议内嵌大数据

因为会增加：

```text
package size
build time
test runtime
download burden
license complexity
maintenance cost
```

。

Bioconductor source package 的 example data 应：

```text
small
fast
deterministic
reproducible
```

。

---

# 34. 最推荐的最终 example data 设计

```text
example_se
│
├─ 500 domains
│
├─ hg38 / chr21
│
├─ coordinates 全部合法
│
├─ Control × 3
│
├─ Treatment × 3
│
├─ BioReplicate metadata
│
├─ Intensity
│
├─ SignalDispersion
│
├─ NativeMaxPeakWidth
│
├─ NativeOccupiedWidth
│
├─ NativePeakCount
│
├─ metadata(native_peaks)
│
├─ 3/3 replicate patterns
│
├─ 2/3 replicate patterns
│
├─ 少量 1/3
│
├─ 少量 no-call
│
└─ simulation provenance
```

。

---

# 35. 推荐最终 colData

```text
SampleID
Condition
BioReplicate
bw_path
peak_path
```

例如：

| SampleID | Condition | BioReplicate | bw_path | peak_path |
|---|---|---:|---|---|
| Control_1 | Control | 1 | NA | NA |
| Control_2 | Control | 2 | NA | NA |
| Control_3 | Control | 3 | NA | NA |
| Treatment_1 | Treatment | 1 | NA | NA |
| Treatment_2 | Treatment | 2 | NA | NA |
| Treatment_3 | Treatment | 3 | NA | NA |

---

# 36. 推荐优先级

## P0

1. 修复 `example_se` chr21 越界。

---

## P1

2. BED writer 使用标准 0-based half-open；
3. 重新生成所有 tiny BED；
4. 增加 `BioReplicate`；
5. synthetic data 刻意包含 3/3、2/3、1/3；
6. 增加少量 no-call；
7. `True_Class` → `Simulation_Scenario`；
8. 增加 simulation/genome provenance；
9. `bw_path/peak_path` 改 NA；
10. extdata README 补 peak files + coordinate convention。

---

## P2 / optional

11. 如果 vignette 需要完整 BEDPE/RNA runnable workflow，再增加 tiny BEDPE/RNA fixture；
12. 不增加第二个 example SE；
13. 不增加大型真实 dataset；
14. 不增加 permanent tiny TxDb artifact。

---

# 37. 从 Bioconductor reviewer 角度的最终评价

如果完成上述 P0/P1，我会认为 epiPortrait 的 example-data 体系非常合理：

```text
processed synthetic object
+
real file-format integration fixtures
+
small package footprint
+
explicit replicate-aware design
+
clear provenance
```

。

这已经足够支持：

```text
examples
vignette
unit tests
API demonstration
```

。

不需要更多数据。

---

# 38. 一句话结论

> **epiPortrait 当前内嵌数据的结构、规模、replicate 数和分组设计都已经合适；主 example 保持 2 groups × 3 biological replicates，不需要扩充。真正应该做的是修正 chr21 越界和 BED 坐标、增加 BioReplicate / provenance，并让 synthetic 数据包含少量 3/3、2/3、1/3 和 no-call 模式，以更充分展示 replicate-aware calling。**
