# epiPortrait README / Vignette 叙事重构方案
## 从“API manual”升级为“科学问题 → 方法创新 → 应用场景 → 工作流”

**目标版本：** epiPortrait14+  
**日期：** 2026-08-12

---

# 1. 总体判断

当前 README / vignette 已经在技术定义上明显改善：

- package identity 已从旧的 shape / 4D / physical-conformation 叙事转向
  `Replicate-Aware Epigenomic Domain Profiling and Remodeling`
- Intensity / native Breadth / SignalDispersion 已基本区分
- Breadth-Super 已说明来自 replicate-specific native PeakWidth
- `IntervalWidth` 已明确只是 shared coordinate frame
- replicate support / Uncertain / transition 已进入主文档
- 与 DiffBind / csaw / ChIPseeker / ROSE / deepTools 的边界已开始说明

但当前文档仍存在一个明显问题：

> **内容是正确的，但价值没有被“讲出来”。**

现在 README 的阅读逻辑仍接近：

```text
What is epiPortrait?
→ feature definition
→ workflow
→ install
→ module 1
→ module 2
→ module 3
```

这更像：

```text
software manual
```

而不是：

```text
methods package landing page
```

。

对于一个准备投稿 / Bioconductor 的包，前 1–2 个屏幕应该首先让用户和 reviewer 理解：

```text
What biological/computational problem exists?
Why are existing one-dimensional views insufficient?
What exactly does epiPortrait add?
Which biological questions can I answer with it?
```

然后才进入 API。

---

# 2. README 最核心的问题：
# “What is epiPortrait?” 太快进入 feature table

当前 opening 很快进入：

```text
Intensity
Native peak breadth
SignalDispersion
```

技术上没有错。

但用户还没有被说服：

> 为什么我要同时看这三个东西？

更好的 opening 应该先给出“问题”。

---

# 3. README 顶部建议的新主命题

推荐 hero statement：

> **Epigenomic domains can remodel by becoming stronger, broader, or both. epiPortrait separates these modes while retaining replicate-level evidence.**

中文概念：

> **表观遗传域的改变不只是信号增强或减弱；它还可以发生宽度扩展、收缩以及二者耦合。epiPortrait 将这些变化拆解，并保留 biological replicate 层面的证据。**

这是比：

```text
profiles quantitative state
```

更能让人立即理解价值的一句话。

---

# 4. 第二句应该回答“为什么现有分析不够”

推荐：

> Conventional peak-based and differential-enrichment workflows are highly effective for identifying enriched regions and changes in signal abundance, but they do not explicitly distinguish whether a regulatory domain is remodeled through signal magnitude, native genomic breadth, or a combination of both.

注意措辞：

```text
complements
```

而不是：

```text
replaces
```

。

不要攻击 DiffBind / csaw / ROSE。

---

# 5. README 应增加一个明确的 “Why epiPortrait?” section

放在 feature table 之前。

建议：

## Why epiPortrait?

Many epigenomic analyses collapse a region into one principal quantity:

```text
How much signal is present?
```

But biologically distinct domains can have similar total signal:

```text
Domain A
high, focal signal

Domain B
moderate signal spread across a broad enriched territory
```

Likewise, two broad domains can differ strongly in total activity.

epiPortrait therefore separates:

```text
signal magnitude
from
native breadth geometry
```

and asks whether each state is reproducibly supported across biological replicates.

---

# 6. 最重要的创新表达：
# 不要写“我们发明了 breadth”

Breadth biology 已经有成熟 H3K4me3 / SE 文献基础。

所以不能把创新写成：

> epiPortrait introduces domain breadth.

真正创新应写成：

> **epiPortrait operationalizes signal magnitude and conventional native-domain breadth as separable, replicate-aware phenotype axes within one shared-domain framework.**

再补一句：

> **The package then converts these axes into interpretable domain states and tracks how those states remodel across conditions.**

---

# 7. README 建议增加 “What does epiPortrait add?” 四点

建议不要写十几个 feature。

只写四件真正具有识别度的东西：

### 1. Two separable canonical phenotype axes

```text
Intensity
×
native Breadth
```

而不是把 width 和 signal 混成一个 rank。

### 2. Replicate-aware evidence

```text
per-replicate call
→ explicit support rule
→ group phenotype
→ Uncertain when evidence is insufficient
```

### 3. Condition remodeling

```text
Typical
↔ Intensity-Super
↔ Breadth-Super
↔ Dual-Super
```

以及 continuous width expansion/contraction。

### 4. Mark-aware interpretation

同一 quantitative framework 可用于：

```text
H3K27ac active enhancer domains
H3K4me3 promoter domains
H3K27me3 / H3K9me3 repressive domains
```

但 biological terminology 按 mark 改变。

---

# 8. README 应增加一个最直观的科学图示

当前 workflow diagram 是工程流程图。

还需要一个“科学概念图”。

推荐：

```text
                   NATIVE BREADTH
                        ↑

      Breadth-Super     |      Dual-Super
                        |
                        |
------------------------+------------------------→ INTENSITY
                        |
          Typical       |    Intensity-Super
                        |
```

旁边解释：

```text
Intensity-Super
= high signal without necessarily broad native geometry

Breadth-Super
= unusually broad native enrichment without necessarily extreme signal

Dual-Super
= both high magnitude and broad geometry
```

这个图比 pipeline diagram 更能体现创新。

---

# 9. README 必须增加 “What biological questions can I ask?”

这是当前文档最缺的一部分。

建议直接放 5 个应用场景。

---

# 10. Application 1 — H3K4me3 promoter-domain remodeling

推荐文案：

### Broad promoter remodeling

H3K4me3 breadth has been associated with stable transcriptional programs,
cell identity, tumor-suppressor regulation and cancer-associated
reprogramming.

epiPortrait can distinguish:

```text
Intensity gain
Breadth expansion
Dual gain
Breadth contraction
```

因此适用于：

```text
normal → cancer
differentiation
drug perturbation
therapy resistance
oncogene / tumor-suppressor-associated remodeling
```

关键卖点：

> **Signal gain and promoter broadening do not need to be treated as the same event.**

---

# 11. Application 2 — H3K27ac enhancer / super-enhancer architecture

推荐：

### Active enhancer-domain remodeling

传统 SE workflow 强调：

```text
stitched enhancer signal ranking
```

epiPortrait 可以进一步问：

```text
Is the domain intensity-dominant?
breadth-dominant?
or dual-extreme?
```

适用于：

```text
oncogenic enhancer acquisition
drug-response enhancer remodeling
lineage-state switching
super-enhancer decomposition
```

可以明确写：

> ROSE comparison is supported as an optional reference analysis; epiPortrait is not intended as a bit-exact ROSE replacement.

---

# 12. Application 3 — Repressive chromatin domains

H3K27me3 / H3K9me3 不要继续使用很强的 “Super” biology 语言。

推荐：

### Repressive-domain remodeling

关注：

```text
intensity increase/decrease
domain expansion/contraction
coupled remodeling
```

surface terminology：

```text
Intensity-Extreme
Extended-Domain
Dual-Extreme
```

核心句：

> **The quantitative engine is mark-independent, while biological interpretation is mark-aware.**

这句话应成为文档中的 signature sentence。

---

# 13. Application 4 — Perturbation and resistance

这是 condition-transition module 的真正应用场景。

例如：

```text
Control → drug
Sensitive → resistant
Normal → tumor
WT → KO
before → after differentiation
```

epiPortrait 不只输出：

```text
signal logFC
```

而是描述：

```text
Typical → Intensity-Super
Typical → Breadth-Super
Intensity-Super → Dual-Super
Breadth-Super → Typical
```

。

推荐文案：

> **This state-based representation is intended to complement conventional differential-enrichment statistics by describing how the quantitative organization of a domain changes.**

---

# 14. Application 5 — Cohort recurrence and heterogeneity

这一点非常适合癌症。

但文档必须谨慎区分：

```text
biological replicates
```

和：

```text
independent tumor models / patients
```

。

可以写：

> In cohort-style datasets, the same support machinery can also summarize cross-sample recurrence, provided that independent patient samples are interpreted as population-level recurrence rather than experimental replicates.

应用：

```text
recurrent tumor-associated domain phenotype
patient-specific domain state
pan-cancer recurrence
```

。

这可以让包从普通 3 vs 3 workflow 扩展到 tumor atlas 使用场景。

---

# 15. README 推荐新增 “Choose epiPortrait when…”

这是非常有用的一块。

推荐：

### Use epiPortrait when you want to know:

| Question | epiPortrait output |
|---|---|
| Which domains carry extreme integrated signal? | Intensity-Super |
| Which domains are unusually broad in native peak geometry? | Breadth-Super |
| Which are both strong and broad? | Dual-Super |
| Is the call reproducible across replicates? | replicate support / Uncertain |
| How does a domain change between conditions? | class transition |
| Is a domain expanding or contracting? | width transition |
| How is signal spatially distributed inside the domain? | SignalDispersion |
| Which genes are spatially associated with the domain? | annotation / candidate links |

这张表非常建议进入 README。

---

# 16. README 同时应该加入 “Do not use epiPortrait when…”

建议很短：

epiPortrait is not designed to:

```text
align FASTQ/BAM
call peaks from reads
replace DiffBind/csaw differential testing
infer causal enhancer–gene regulation
call chromatin loops
perform RNA-seq differential expression
```

这一段可以减少 reviewer 对 scope creep 的担忧。

---

# 17. Relationship to existing methods 不应只有 bullets

目前：

```text
DiffBind
ChIPseeker
ROSE
deepTools
```

已经有了，但太简短。

建议改成表格：

| Tool / family | Primary question | Where epiPortrait differs |
|---|---|---|
| MACS2 / SICER / epic2 | Where are enriched peaks/domains? | epiPortrait starts after domain calling |
| DiffBind / csaw | Where does enrichment differ statistically? | epiPortrait describes domain phenotype and remodeling mode |
| ROSE | Which stitched enhancers are signal-extreme? | epiPortrait separates Intensity from native Breadth and adds replicate support |
| ChIPseeker | Where is a peak relative to genes/features? | epiPortrait uses annotation after quantitative phenotyping |
| deepTools | How do signal tracks look across regions/samples? | epiPortrait converts track measurements into domain-level phenotypes |

重点：

> 不写 better / superior。

写：

```text
different analytical target
```

。

---

# 18. README 推荐最终目录

建议从现在：

```text
What is
Installation
Core Module 1
Core Module 2
Core Module 3
...
```

改成：

```text
1. Hero statement
2. Why epiPortrait?
3. What does epiPortrait add?
4. Conceptual domain phenotype
5. What biological questions can it address?
6. Relationship to existing methods
7. Input and data model
8. Quick Start
9. Core workflow
10. Interpreting classes and transitions
11. Mark-specific application patterns
12. Scope and limitations
13. Installation
14. Citation
```

其中：

```text
Core Module 1/2/3
```

可以保留，但降到 Quick Start 后面。

---

# 19. README 长度不用大幅增加

当前问题不是“字少”。

而是前半部分：

```text
feature → parameter → function
```

占比过高。

建议：

```text
前 35–40%
= scientific identity + applications

后 60–65%
= API + executable examples
```

。

---

# 20. Vignette 与 README 不能只是长短区别

Vignette 应从：

```text
long README
```

升级成：

> **scientific tutorial / analysis reasoning guide**

用户看完 vignette 应知道：

```text
为什么这样构建 universe
为什么 Breadth 不能用 consensus width
为什么 replicate 要独立 call
什么时候用 relative vs pooled/reference transition
不同 marks 如何解释
哪些 conclusions 不能下
```

。

---

# 21. Vignette 推荐新标题

首选：

> **From Signal Tracks to Replicate-Supported Epigenomic Domain Phenotypes**

副标题：

> **A practical guide to Intensity, native Breadth, domain architecture and condition remodeling with epiPortrait**

相比：

```text
Replicate-Aware Profiling...
```

更像 tutorial。

---

# 22. Vignette Introduction 建议用“矛盾案例”开场

不要先定义三个 assays。

可以这样开：

```text
Two epigenomic domains may carry similar integrated signal while differing
substantially in their native genomic extent. Conversely, two equally broad
domains may differ by an order of magnitude in signal abundance.
```

然后：

> These scenarios represent different regulatory states but can be obscured when a region is summarized by a single ranking statistic.

再引出：

```text
Intensity × native Breadth
```

。

这样 reader 会自然理解为什么需要 package。

---

# 23. Vignette 应明确 5 个 design principles

建议在 Introduction 后加入：

## Design principles

### Principle 1
Shared coordinates are for comparison, not Breadth truth.

### Principle 2
Breadth comes from replicate-specific native geometry.

### Principle 3
Intensity and Breadth are separable canonical axes.

### Principle 4
Replicate evidence is retained before group aggregation.

### Principle 5
Biological interpretation depends on the histone mark.

这五条几乎就是 epiPortrait 的方法学 identity。

---

# 24. Vignette 增加 “A domain can change in different ways”

强烈建议增加概念表：

| Ref | Target | Interpretation |
|---|---|---|
| Typical | Intensity-Super | magnitude-dominant activation |
| Typical | Breadth-Super | breadth expansion |
| Typical | Dual-Super | coupled magnitude + breadth acquisition |
| Intensity-Super | Dual-Super | acquisition of extreme breadth |
| Breadth-Super | Dual-Super | acquisition of extreme magnitude |
| Breadth-Super | Typical | breadth contraction / loss of broad state |

并明确：

```text
relative transition
≠ absolute biochemical gain/loss
```

。

这个 section 会非常体现 transition module 的价值。

---

# 25. Vignette 应增加 “Mark-specific interpretation” 主章节

不能只放在 Best Practices 后半段。

建议独立成一个主 Module。

---

# 26. H3K4me3 章节建议这样讲

## H3K4me3: promoter breadth and transcriptional state

强调：

```text
broad promoter domains
cell-identity programs
cancer-associated redistribution
```

然后解释 epiPortrait 的问题：

```text
Is a broad promoter also intensity-extreme?
Does cancer-associated broadening occur without signal gain?
Does shortening precede or accompany intensity loss?
```

不要 claim：

```text
epiPortrait discovers broad H3K4me3
```

。

而是：

> **epiPortrait decomposes known breadth biology into separable quantitative remodeling modes.**

---

# 27. H3K27ac 章节建议这样讲

## H3K27ac: enhancer-domain magnitude versus regulatory territory

问题：

```text
Is a high-ranking enhancer strong because of concentrated signal,
because it spans a broad regulatory territory,
or because both occur together?
```

然后解释：

```text
Intensity-Super
Breadth-Super
Dual-Super
```

作为 enhancer domain phenotypes。

ROSE 放在这里做 optional comparison，而不是前面频繁出现。

---

# 28. H3K27me3 / H3K9me3 章节

## Broad repressive domains

解释：

```text
“Super” is an internal class label
```

但 biological interpretation：

```text
Intensity-Extreme
Extended-Domain
Dual-Extreme
```

。

并明确：

```text
upstream broad caller should be mark-appropriate
```

。

这会显著加强包不是 H3K27ac-only 的印象。

---

# 29. Vignette 应增加 “Application patterns” section

用问题而不是 mark 来组织：

### Pattern A — Disease state remodeling

```text
Normal → Cancer
```

### Pattern B — Drug / perturbation response

```text
Vehicle → Drug
WT → KO
```

### Pattern C — Therapy resistance

```text
Sensitive → Resistant
```

### Pattern D — Differentiation / lineage transition

```text
State A → State B
```

### Pattern E — Cohort recurrence

```text
Independent patient / cell-line models
```

。

每个只用 3–5 句，不需要真实 benchmark 数据塞进 vignette。

---

# 30. Vignette 应增加一个“interpretation guardrails”章节

这是 Bioconductor reviewer 很喜欢的。

例如：

## Interpretation guardrails

### Intensity
requires quantitatively comparable BigWigs for between-condition magnitude claims.

### Breadth
depends on upstream native peak/domain caller and boundary quality.

### Dual
is a conjunction of the two axes, not a new independent measurement.

### SignalDispersion
describes spatial signal organization, not physical chromatin conformation.

### Annotation
nearest / overlapping genes are candidates, not causal targets.

### Cohort support
independent patients are recurrence evidence, not conventional experimental replicates.

这些会明显提升方法成熟度。

---

# 31. README 与 vignette 最值得加入的 signature phrases

推荐固定使用以下 5 句。

### Phrase 1

> **epiPortrait asks not only whether an epigenomic domain changes, but how it changes.**

### Phrase 2

> **Signal magnitude and native genomic breadth are treated as separable quantitative properties.**

### Phrase 3

> **Replicate-level evidence is retained before group-level phenotypes are assigned.**

### Phrase 4

> **Shared domain coordinates provide a comparison frame; they do not define sample-specific breadth.**

### Phrase 5

> **The quantitative framework is mark-independent, while biological interpretation is mark-aware.**

这五句话应该贯穿：

```text
README
vignette
DESCRIPTION
manuscript Introduction
```

形成统一品牌语言。

---

# 32. 推荐 README opening 完整示例

可以直接考虑这样的 opening：

> **epiPortrait asks not only whether an epigenomic domain changes, but how it changes.**
>
> Conventional epigenomic workflows are highly effective at identifying peaks and differential enrichment, but a regulatory domain can remodel through distinct quantitative modes: its signal can become stronger, its native enriched territory can become broader, or both can occur together.
>
> epiPortrait provides a replicate-aware framework that separates **integrated signal magnitude (Intensity)** from **replicate-specific native peak breadth**, while retaining **SignalDispersion** as a complementary descriptor of within-domain signal organization. Replicate-level evidence is aggregated using explicit support rules to assign interpretable domain phenotypes — `Intensity-Super`, `Breadth-Super`, `Dual-Super`, `Typical`, or `Uncertain` — and to follow their remodeling across biological conditions.
>
> The framework is designed primarily for histone-mark ChIP-seq and CUT&Tag domains, including active enhancer domains (H3K27ac), broad promoter domains (H3K4me3), and broad repressive chromatin (H3K27me3/H3K9me3). The quantitative engine is mark-independent, while biological interpretation is mark-aware.

这个 opening 比直接 feature table 更有吸引力。

---

# 33. 推荐 vignette Introduction opening 示例

> Epigenomic domains are often summarized by signal enrichment, yet signal magnitude is only one aspect of domain organization. A focal, high-amplitude H3K27ac region and a broad, moderately enriched regulatory domain may carry similar integrated signal while representing different chromatin states. Likewise, broad H3K4me3 or H3K27me3 domains can expand or contract without proportional changes in signal abundance.
>
> epiPortrait was designed to represent these changes explicitly. It separates integrated signal magnitude from replicate-specific native domain breadth, retains within-domain signal dispersion as a complementary structural descriptor, and aggregates evidence only after each biological replicate has been evaluated independently.
>
> The resulting domain portrait is not intended to replace peak calling or differential-enrichment statistics. Instead, it provides a downstream phenotype layer for asking whether a domain is predominantly intensity-extreme, breadth-extreme, jointly extreme, or undergoing a transition between these states.

---

# 34. README 应减少哪些内容的 prominence？

当前需要降级：

```text
ROSE/tangent
bootstrap details
custom_features
HDF5/cache implementation details
```

不是删除。

但这些不应该占据 landing page 的高位置。

可以移到：

```text
Advanced usage
Performance
Reference
```

。

---

# 35. 哪些功能反而应该提前？

需要提前：

```text
Intensity × Breadth concept
replicate support
Uncertain
condition transitions
application scenarios
mark-aware interpretation
```

。

这些才是用户决定“要不要用这个包”的依据。

---

# 36. README / vignette 中 “Super-domain” 的使用建议

不要删除 taxonomy。

但 package identity 继续使用：

```text
domain phenotyping
domain profiling
domain remodeling
```

。

“Super” 只出现在 classification 层：

```text
Intensity-Super
Breadth-Super
Dual-Super
```

。

对于 repressive mark：

surface labels 用：

```text
Intensity-Extreme
Extended-Domain
Dual-Extreme
```

。

---

# 37. Innovation section 推荐最终写法

建议 README 加一个很短的：

## Methodological focus

> epiPortrait does not propose that domain breadth itself is a new biological concept. Its methodological contribution is to place conventional native-domain breadth and continuous signal magnitude into a unified replicate-aware framework, preserve their distinct evidence sources, and convert them into interpretable condition-specific domain phenotypes and remodeling states.

这句话很重要。

它能提前化解 reviewer：

> broad H3K4me3 / super-enhancer width 已经有人研究过了。

你的创新不是：

```text
breadth exists
```

而是：

```text
joint + replicate-aware + transition framework
```

。

---

# 38. 最终文档价值主线

README / vignette 都应该围绕：

```text
A domain can be:

strong but not broad
broad but not strong
both strong and broad
neither
uncertain because evidence is insufficient
```

然后：

```text
these states can remodel across conditions
```

。

这是 epiPortrait 最容易被记住的概念。

---

# 39. 最终建议

不建议继续在现有 README 结构上小修小补。

建议做一次明确的：

```text
narrative refactor
```

而不是：

```text
API rewrite
```

。

代码示例大部分可以继续保留。

真正重写的主要是：

```text
README 前 35–40%
Vignette Introduction
Conceptual model
Relationship to methods
新增 Application scenarios
新增 Mark-specific interpretation
新增 Interpretation guardrails
```

。

---

# 40. 一句话结论

> **现在 epiPortrait 的代码身份已经比文档身份成熟。下一步文档不应继续从“有哪些 assay / 参数 / module”出发，而应该从“epigenomic domains 可以通过强度、宽度或二者共同发生重塑”这一科学问题出发，把 Intensity × native Breadth、replicate support、condition transition 和 mark-aware interpretation 作为四个核心价值，再用 H3K4me3、H3K27ac、repressive domains、药物/耐药和肿瘤异质性场景告诉用户为什么以及什么时候使用 epiPortrait。**
