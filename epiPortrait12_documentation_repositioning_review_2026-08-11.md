# epiPortrait 文档重构建议
## 生信专家 / 表观遗传学 / Bioconductor reviewer 视角

**评估版本：epiPortrait12**  
**日期：2026-08-11**

---

# 1. 总结结论

## 需要重新写，但不是“所有文档推倒重来”

当前代码和方法学已经明显从早期：

```text
4D peak shape analysis
↓
Super-domain calling
```

演化为：

```text
replicate-aware epigenomic domain profiling
↓
Intensity
native peak breadth
SignalDispersion
↓
replicate-supported phenotype
↓
condition remodeling / transitions
↓
domain annotation
```

因此当前真正需要的是：

> **重写 package-level narrative / positioning，而不是重写所有函数 Rd。**

推荐强度：

| 文档 | 修改程度 |
|---|---:|
| DESCRIPTION Title / Description | **重写** |
| README.Rmd / README.md | **重写前半部分与工作流结构** |
| 主 vignette | **较大幅重构** |
| package-level `?epiPortrait` | **重写定位部分** |
| Citation / BibTeX title | **同步修改** |
| NEWS | 清理历史定位矛盾 |
| individual function Rd | **不需要重写，只做术语一致性修订** |
| tests / API | 不因文档重写而改变 |

---

# 2. 当前最大问题：包身份仍然过度等同于 “Super-Domain”

当前 DESCRIPTION：

```text
Title:
Quantitative Profiling of Epigenomic Signal Super-Domains
```

README：

```text
Quantitative Profiling of Epigenomic Super-Domains from Continuous Signal Tracks
```

vignette：

```text
epiPortrait: Profiling Epigenomic Super-Domains from Continuous Signal Tracks
```

这几个标题都会让 reviewer / 用户第一眼认为：

> epiPortrait 是一个新的 super-enhancer / super-domain caller。

但现在的实际方法已经明显更广：

```text
domain-level quantitative phenotype
replicate-aware calling
native peak breadth geometry
SignalDispersion
condition remodeling
heterochromatin domain expansion/contraction
domain annotation
```

因此当前 headline 已经**低估了包的真正方法学范围，同时错误强化了与 ROSE 的直接竞争关系**。

---

# 3. 推荐新的包身份

## 首选定位

> **Replicate-aware epigenomic domain profiling and remodeling**

更完整的一句话：

> **epiPortrait is a replicate-aware framework for quantitative phenotyping of epigenomic domains by signal magnitude, native peak breadth, and within-domain signal architecture, with support for condition-specific classification and remodeling analysis.**

这句话比：

```text
Super-domain caller
```

更准确。

---

# 4. 推荐 DESCRIPTION Title

## 首选

```text
Replicate-Aware Profiling of Epigenomic Domain Architecture
```

或者：

```text
Quantitative Profiling of Epigenomic Domain Architecture
```

如果希望与 logo：

```text
Epigenomic Domain Profiling
```

高度一致，推荐：

```text
Replicate-Aware Epigenomic Domain Profiling
```

---

# 5. 不建议 DESCRIPTION Title 继续出现 “Super-Domains”

原因：

1. `Super-domain` 只是 package 的一个 downstream phenotype / classification outcome；
2. H3K27me3 / H3K9me3 使用的是：
   - Intensity-Extreme
   - Extended-Domain
   - Dual-Extreme
   等 display aliases；
3. SignalDispersion 是 architecture descriptor，而非 Super axis；
4. condition transition / width remodeling 已经超出单纯 Super calling；
5. annotation / candidate genes 也是 domain-level downstream interpretation。

所以：

> **Super-domain 应从“package identity”降为“one analytical module”。**

---

# 6. 推荐 DESCRIPTION Description 结构

当前 DESCRIPTION 前几句已经比较接近正确设计，但顺序仍以 super-domain 为中心。

推荐改成四段逻辑：

```text
1. What problem?
2. What domain features?
3. How are replicates handled?
4. What downstream analyses are supported?
```

建议文本：

```text
epiPortrait is a replicate-aware computational framework for quantitative
profiling and comparison of epigenomic domains from normalized continuous
signal tracks and per-sample native peak/domain calls.

It characterizes each shared genomic domain by integrated signal magnitude
(Intensity), sample-specific native peak breadth geometry, and continuous
within-domain signal architecture (SignalDispersion).

Replicate-level evidence is retained explicitly and aggregated using
transparent support rules, enabling condition-specific classification of
Intensity-Super, Breadth-Super, Dual-Super, Typical, or Uncertain states
without using consensus-domain width as the Breadth caller.

The package further supports condition transitions, continuous domain-width
remodeling, signal quality control, genome-aware domain annotation, optional
BEDPE contact evidence, expression-aware candidate-gene prioritization, and
publication-oriented visualization.
```

重点：

> 第一段不要先说 Super-domain。

---

# 7. README 需要重写的程度：约前 40–50%

README 后半部分 API/example 基本可以保留。

真正应该重写：

```text
Title
What is epiPortrait?
feature table
package positioning
workflow diagram
Quick Start narrative
comparison / scope
```

---

# 8. README 推荐新标题

当前：

```text
Quantitative Profiling of Epigenomic Super-Domains from Continuous Signal Tracks
```

建议：

```text
Replicate-Aware Epigenomic Domain Profiling and Remodeling
```

或者更简洁：

```text
Replicate-Aware Epigenomic Domain Profiling
```

logo 已经写：

```text
Epigenomic Domain Profiling
```

两者统一会非常好。

---

# 9. README “What is epiPortrait?” 应重新组织

当前：

```text
super-domains
signal magnitude
sample-specific occupied breadth
```

建议改成：

```text
epiPortrait profiles the quantitative state of predefined/shared
epigenomic domains across biological replicates and conditions.
```

然后立即定义三个 measurement layers：

```text
1. Intensity
2. Native breadth geometry
3. SignalDispersion
```

再说明：

```text
Intensity + Breadth
→ canonical domain phenotype taxonomy

SignalDispersion
→ secondary architecture descriptor
```

最后才介绍：

```text
Intensity-Super
Breadth-Super
Dual-Super
Typical
Uncertain
```

这样读者不会把：

```text
Super
```

误认为输入或包唯一目标。

---

# 10. “occupied breadth” 应减少使用

README 当前写：

```text
sample-specific occupied breadth
```

但 canonical Breadth-Super 实际依据：

```text
single native PeakWidth distribution
→ replicate-specific width inflection
```

不是：

```text
NativeOccupiedWidth
```

所以建议统一为：

```text
native peak breadth geometry
```

或：

```text
sample-specific native peak breadth
```

`NativeOccupiedWidth` 仍保留为 descriptive domain-level metric。

---

# 11. 推荐 README Feature Table

建议把表格改成：

| Dimension | Definition | Role |
|---|---|---|
| **Intensity** | Integrated normalized signal across the shared domain | Canonical phenotype axis |
| **PeakWidth** | Width of each replicate's native peak/domain | Canonical Breadth-Super evidence |
| **NativeMaxPeakWidth** | Maximum mapped native peak width per shared domain/sample | Domain descriptor / transitions |
| **NativeOccupiedWidth** | Sum of reduced native occupied segments within the shared domain | Domain descriptor |
| **SignalDispersion** | Signal-weighted genomic SD | Secondary architecture descriptor |
| **IntervalWidth** | Width of the shared coordinate frame | Coordinate descriptor only |

这样 reviewer 一眼能看清：

> 哪个才是 Breadth caller。

---

# 12. README 应增加一个非常清楚的 workflow figure / text diagram

建议：

```text
Normalized BigWig + native peak/domain calls
                    │
                    ▼
        shared candidate-domain universe
                    │
        ┌───────────┼───────────┐
        ▼           ▼           ▼
    Intensity   Peak breadth   SignalDispersion
        │           │           │
        └──────┬────┘           │
               ▼                │
       replicate-supported      │
         domain phenotype       │
               │                │
               ▼                ▼
     condition transitions   architecture
               │
               ▼
      domain annotation /
      candidate genes
```

这会比 “Core Module 1/2/3/4” 更能说明方法学逻辑。

---

# 13. vignette 需要比 README 更大幅度重构

Bioconductor 明确要求 vignette Introduction：

- 有 package motivation；
- 说明目标；
- 描述关键独特功能；
- 与相近工具比较；
- 解释为什么该 package 适合 Bioconductor。

因此 epiPortrait 当前 vignette 不应该只是：

```text
What the assays are
→ run functions
```

而应该把“方法学边界”讲清楚。

---

# 14. 推荐 vignette 新结构

## 1. Introduction

回答：

```text
为什么需要 epiPortrait？
```

核心问题：

> Conventional peak-based analyses emphasize differential enrichment or
> peak presence, whereas broad epigenomic domains also vary in signal
> magnitude, native breadth geometry, and internal signal organization.

然后引出：

```text
epiPortrait = domain phenotyping
```

---

## 2. Conceptual model

专门解释：

```text
shared coordinate frame
vs
sample-specific native peak geometry
```

这是 epiPortrait 最重要、也是最容易被 reviewer 误解的地方。

---

## 3. The three quantitative layers

### Intensity

### Native breadth

### SignalDispersion

明确：

```text
Intensity and Breadth are canonical classification axes.

SignalDispersion is a secondary architecture descriptor.
```

---

## 4. Why replicate-level calling?

解释：

```text
each replicate called independently
→ support aggregation
→ no-call remains no-call
```

这是包的核心差异化价值。

---

## 5. Building the domain universe

解释：

```text
union/reduce
consensus support
optional stitching
mark presets
```

以及：

> consensus width is not the Breadth-Super metric.

这句必须显眼。

---

## 6. Domain phenotype calling

### Intensity

### Breadth

### Combined taxonomy

---

## 7. Condition remodeling

包括：

```text
relative prominence
reference/pooled cutoff
width expansion/contraction
```

---

## 8. Architecture descriptors

这里再讲：

```text
SignalDispersion
NativeOccupiedWidth
NativePeakCount
```

不要和 canonical classifier 混在一起。

---

## 9. Annotation

```text
TxDb linear annotation
BEDPE evidence
RNA expression prioritization
```

明确：

```text
candidate / associated genes
not causal target genes
```

---

## 10. Mark-specific interpretation

### H3K4me3

broad active promoters

### H3K27ac

enhancer/super-enhancer architecture

### H3K27me3 / H3K9me3

domain expansion/contraction / intensity extremes

这样 package 的泛化价值更明显。

---

## 11. Relationship to existing methods

必须有一小节。

不需要攻击其他工具。

推荐写：

```text
DiffBind / csaw:
differential enrichment

ChIPseeker:
genomic annotation

ROSE:
1D enhancer ranking / super-enhancer identification

deepTools:
signal visualization / profiling

epiPortrait:
replicate-aware quantitative domain phenotyping
```

这会非常帮助 reviewer 理解生态位。

---

# 15. README/vignette 里 ROSE 的位置应调整

当前 README 对：

```text
ROSE-style comparison
```

写得很突出。

建议保留，但把它降为：

```text
H3K27ac benchmark / optional comparison
```

而不是让用户感觉：

```text
epiPortrait = generalized ROSE
```

epiPortrait 的真正价值：

```text
Intensity × native Breadth decomposition
replicate support
condition remodeling
```

比 ROSE 模仿本身重要。

---

# 16. package-level documentation 需要同步

`R/00_package.R` / `?epiPortrait` 应重新写 opening。

推荐：

```text
epiPortrait provides replicate-aware quantitative phenotyping of epigenomic
domains across biological samples and conditions.
```

然后：

```text
Core analyses:
- domain matrix construction
- Intensity and native Breadth classification
- SignalDispersion characterization
- replicate support
- condition transitions
- domain annotation
```

不要开头就定义为：

```text
Super-domain package
```

。

---

# 17. individual Rd 不需要重写

当前绝大多数 function docs 已经经历多轮 correctness 修复。

所以：

```text
call_super_domains
build_portrait_matrix
get_replicate_calls
annotate_epi_domains
get_domain_genes
compare_superdomains
```

不需要为了新 positioning 再推倒。

只需要做 terminology sweep：

```text
occupied breadth
super-domain package
peak-shape
4D
physical conformation
```

这些旧叙事清掉即可。

---

# 18. Citation title 必须同步

当前 README citation：

```text
epiPortrait: Quantitative Profiling of Epigenomic Super-Domains
```

建议改成：

```text
epiPortrait: Replicate-Aware Profiling of Epigenomic Domain Architecture
```

或者：

```text
epiPortrait: Replicate-Aware Epigenomic Domain Profiling
```

这和新 logo：

```text
Epigenomic Domain Profiling
```

形成统一品牌。

---

# 19. Logo tagline 与 package title 建议形成三层结构

推荐：

## Logo

```text
Epigenomic Domain Profiling
```

## DESCRIPTION title

```text
Replicate-Aware Epigenomic Domain Profiling
```

## manuscript title

可以更强调方法学：

```text
epiPortrait: Replicate-Aware Quantitative Phenotyping of Epigenomic Domains
```

三者互相一致，但不需要完全重复。

---

# 20. Heterochromatin 叙事应该正式进入 vignette

现在代码已经支持：

```text
H3K27me3
H3K9me3
```

并有：

```text
Extended-Domain
Intensity-Extreme
Dual-Extreme
```

display terminology。

如果 vignette 仍然主要写：

```text
Super-domain
```

会显得这个扩展很突兀。

建议用：

> epiPortrait's internal taxonomy is mark-independent, while biological
> interpretation is mark-aware.

这是一个很好的方法学表述。

---

# 21. ATAC-seq 应降低宣传强度

当前 DESCRIPTION：

```text
ChIP-seq, ATAC-seq, CUT&Tag
```

技术上支持。

但 canonical Breadth / domain phenotype 的生物学验证主要还是 histone marks。

所以推荐写：

```text
histone-mark ChIP-seq / CUT&Tag as the primary use case,
with continuous-track workflows also applicable to ATAC-seq-derived domains.
```

不要让 reviewer 误以为已经对 ATAC 全面 benchmark。

---

# 22. 建议 README 增加 “What epiPortrait is not”

保留且加强当前 scope box：

```text
epiPortrait does not:
- call peaks from BAM/FASTQ
- replace differential peak callers
- infer causal enhancer–gene regulation
- call chromatin loops
- perform RNA-seq differential expression
```

这样 BEDPE/RNA 模块不会让包显得 scope creep。

---

# 23. README Quick Start 应以 `per_group` 为主要多条件示例

如果 example data 是：

```text
Control × 3
Treatment × 3
```

README quick start 应明确：

```r
call_super_domains(
    se,
    feature = "Intensity",
    mode = "per_group",
    group_var = "Condition"
)

call_super_domains(
    se,
    feature = "Breadth",
    mode = "per_group",
    group_var = "Condition"
)
```

而不是让初学者默认先跑：

```text
global_consensus
```

。

function default 可以不改。

---

# 24. README 与 vignette 不应完全重复

建议：

## README

目标：

```text
2–5 分钟看懂 package
```

包括：

```text
定位
核心概念
quick start
主要图
installation
citation
```

## Vignette

目标：

```text
完整方法学与 workflow 教程
```

包括：

```text
概念模型
输入要求
replicate semantics
domain construction
calling
transition
annotation
mark-specific interpretation
method comparison
```

这样文档层次更专业。

---

# 25. Bioconductor reviewer 视角：为什么现在值得重写

当前 Bioconductor review expectations 强调：

```text
ease of use
documentation
well-written code
interoperability
```

而 vignette Introduction 要明确：

```text
objective
unique functions
key distinguishing points
comparison to related packages
motivation for inclusion in Bioconductor
```

现在代码已经成熟，但旧的 package-level narrative 会：

```text
低估方法学价值
混淆 Super-domain 与 domain profiling
让 reviewer 错误拿 ROSE 当唯一 benchmark
使 heterochromatin generalization 看起来像后加功能
```

因此：

> **现在重写文档的收益已经高于继续加新功能。**

---

# 26. 推荐重写优先级

## P0

### DESCRIPTION

- Title
- Description

### README

- headline
- What is epiPortrait?
- feature terminology
- workflow
- quick start
- citation

### vignette

- title
- Introduction
- conceptual architecture
- relationship to other tools
- mark-specific interpretation

---

## P1

### package-level documentation

- `R/00_package.R`
- `man/epiPortrait.Rd`

### terminology sweep

搜索：

```text
4D
peak shape
occupied breadth
physical conformation
Super-domain package
consensus peak calling
```

逐项确认。

### citation/BibTeX

同步。

---

## P2

### NEWS

保留历史 change log，但顶部写清当前 final scope。

### function Rd

只做术语一致性，不改 API。

---

# 27. 不建议的做法

不要为了新定位：

```text
重命名所有 functions
改 class 名
改 Intensity/Breadth taxonomy
改 API default
增加新功能
重新设计 annotation
```

。

文档重构应该：

> **反映已经稳定的代码，而不是再次驱动代码变化。**

---

# 28. 推荐最终 narrative

一句话：

> **epiPortrait is a replicate-aware framework for quantitative phenotyping and remodeling analysis of epigenomic domains.**

第二层：

> It separates integrated signal magnitude from native peak breadth while
> retaining continuous within-domain signal architecture as an additional
> descriptor.

第三层：

> Replicate-supported domain states can be compared across biological
> conditions and interpreted through genome-aware annotation, optional
> chromatin-contact evidence, and expression-aware candidate-gene
> prioritization.

这三句话基本就构成整个 package 的新叙事主干。

---

# 29. 最终判断

## 是否需要重新写？

**需要。**

但准确说：

> **需要重新写 package-level narrative，大约 README/vignette 顶层内容的 40–50%；不需要重新写算法文档和所有函数手册。**

当前最需要修的不是代码，而是：

```text
代码已经是“domain phenotyping framework”
文档仍部分写成“Super-domain caller”
```

把两者同步后，epiPortrait 的方法学定位会明显更清楚，也更符合目前新 logo：

```text
Epigenomic Domain Profiling
```

的品牌方向。
