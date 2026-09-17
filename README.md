# epiPortrait <img src="man/figures/logo.jpg" align="right" width="160" alt="epiPortrait Logo" />

**Replicate-Aware Epigenomic Domain Profiling and Remodeling**

[![R-CMD-check](https://github.com/zying106/epiPortrait/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zying106/epiPortrait/actions/workflows/R-CMD-check.yaml)
[![BiocCheck](https://github.com/zying106/epiPortrait/actions/workflows/bioccheck.yaml/badge.svg)](https://github.com/zying106/epiPortrait/actions/workflows/bioccheck.yaml)
[![codecov](https://codecov.io/gh/zying106/epiPortrait/branch/main/graph/badge.svg)](https://app.codecov.io/gh/zying106/epiPortrait)
[![License: GPL (\>=
3)](https://img.shields.io/badge/License-GPL%20(%3E%3D%203)-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)

------------------------------------------------------------------------

**epiPortrait describes signal and native peak breadth on a shared set of
genomic domains.**

Conventional epigenomic workflows are highly effective at identifying peaks and
differential enrichment, but a regulatory domain can remodel through distinct
quantitative modes: its signal can become stronger, its native enriched
territory can become broader, or both can occur together.

epiPortrait measures **integrated signal (Intensity)** and
**replicate-specific native peak breadth** separately, then summarizes their
support and transitions across biological conditions.

Here, *replicate-aware* means that samples are called separately before their
calls are aggregated with a stated support rule. It does not mean that the
rank-based calls model between-replicate variance. Continuous-signal inference
is available separately through `analyze_differential_domains()`.

------------------------------------------------------------------------

## Why epiPortrait?

Many epigenomic analyses collapse a region into one principal quantity:
*how much signal is present?* But biologically distinct domains can carry
similar integrated signal:

```text
Domain A: high, focal signal
Domain B: moderate signal spread across a broad enriched territory
```

Likewise, two equally broad domains can differ strongly in total signal.
When a region is summarized by a single ranking statistic, these distinct
signal geometries are conflated.

epiPortrait therefore separates **signal magnitude** from **native breadth
geometry**, and asks whether each state is reproducibly supported across
biological replicates.

**Methodological focus.** epiPortrait does not propose that domain breadth is
itself a new biological concept—broad H3K4me3 domains
([Benayoun *et al.*, 2014](https://doi.org/10.1016/j.cell.2014.06.027)) and
stitched H3K27ac super-enhancers
([Whyte *et al.*, 2013](https://doi.org/10.1016/j.cell.2013.03.035)) have an
established literature. Its contribution is to
place conventional native-domain breadth and continuous signal magnitude into
a **unified replicate-aware framework**, preserve their distinct evidence
sources, and convert them into interpretable condition-specific domain
phenotypes and remodeling states.

## What does epiPortrait add?

### 1. Two separable canonical phenotype axes

`Intensity` and native `Breadth` are computed and thresholded separately rather
than merged into one rank. They are not assumed to be statistically
independent: integrated signal can increase with both signal amplitude and the
length of the analysed interval.

| Axis | Measures | Canonical class |
|:-----|:---------|:----------------|
| **Intensity** | Integrated normalized signal across the shared domain | `Intensity-Super` |
| **Native breadth** | Width of each replicate's native peak/domain call | `Breadth-Super` |
| **SignalDispersion** | Signal-weighted genomic SD within the domain (secondary architecture descriptor) | — |

### 2. Replicate-aware evidence

Each biological replicate is called independently; replicate calls are
aggregated with an explicit support rule before a group phenotype is assigned.
`Uncertain` is assigned when evidence is insufficient — it is an abstention,
never a silent relabel.

### 3. Condition remodeling

Domains transition between states across conditions
(`Typical ↔ Intensity-Super ↔ Breadth-Super ↔ Dual-Super`), and continuous
width expansion / contraction is tracked separately.

### 4. Mark-aware interpretation

The same calculations can be applied to active enhancer domains (H3K27ac), broad
promoter domains (H3K4me3), and broad repressive chromatin
(H3K27me3 / H3K9me3), with biological terminology adapted per mark.

The calculations are mark-independent, but their biological interpretation is
not.

## How do domains combine the two axes?

Each domain is classified separately on the two axes; their combination gives
four architecture states:

```text
                       BREADTH-SUPER (broad native geometry)
                              │
      Breadth-Super           │          Dual-Super
      (broad, not            │          (broad and
       signal-extreme)       │           signal-extreme)
──────────────────────────────┼──────────────────────────────  INTENSITY-SUPER
      Typical                 │          Intensity-Super
      (neither)               │          (signal-extreme,
                              │           not necessarily broad)
                              │
```

- `Intensity-Super` = extreme integrated signal, without necessarily broad native geometry.
- `Breadth-Super` = unusually broad native enrichment, without necessarily extreme signal.
- `Dual-Super` = both high magnitude and broad geometry.

`Uncertain` is reported when the evidence in a condition is not reliable
(e.g. no native peaks, or an unstable width inflection).

## What biological questions can epiPortrait address?

| Question | epiPortrait output |
|:---------|:-------------------|
| Which domains carry extreme integrated signal? | `Intensity-Super` |
| Which domains are unusually broad in native peak geometry? | `Breadth-Super` |
| Which are both strong and broad? | `Dual-Super` |
| Is the call reproducible across replicates? | replicate support / `Uncertain` |
| How does a domain change between conditions? | class transition |
| Does measured native breadth increase or decrease? | width transition |
| How is signal spatially organized inside the domain? | `SignalDispersion` |
| Which genes are spatially associated with the domain? | annotation / candidate links |

**Representative applications.** Broad-promoter remodeling (H3K4me3:
intensity gain vs breadth expansion vs coupled change; normal → cancer,
differentiation, drug perturbation). Active enhancer / super-enhancer
architecture (H3K27ac: intensity-extreme vs breadth-extreme vs dual-extreme;
oncogenic enhancer acquisition, drug response, lineage switching). Repressive
domain remodeling (H3K27me3 / H3K9me3: expansion/contraction, with
mark-appropriate surface labels). Perturbation and therapy resistance
(Control → drug, Sensitive → Resistant, WT → KO), where epiPortrait describes
the form of the state change in addition to a signal log-fold-change.

**Cohort-style studies.** When a condition comprises independent patients, the
support fraction describes cross-patient recurrence rather than technical or
within-subject replication. Choose the support rule to match the recurrence
question and report the number of contributing patients. With one sample per
condition, use `mode = "per_sample"`; a group-level support claim is not
available.

## Relationship to existing methods

| Tool / family | Primary question | Where epiPortrait differs |
|:--------------|:-----------------|:--------------------------|
| MACS2 / SICER / epic2 | Where are enriched peaks/domains? | epiPortrait starts after domain calling |
| DiffBind / csaw | Where does enrichment differ statistically? | epiPortrait describes domain phenotype and remodeling mode |
| ROSE | Which stitched enhancers are signal-extreme? | epiPortrait separates Intensity from native Breadth and adds replicate support |
| ChIPseeker | Where is a peak relative to genes/features? | epiPortrait uses annotation after quantitative phenotyping |
| deepTools | How do signal tracks look across regions/samples? | epiPortrait converts track measurements into domain-level phenotypes |

These tools address different analytical targets; epiPortrait is a downstream
phenotype layer, not a replacement for peak calling or differential-enrichment
testing.

## Scope and limitations

epiPortrait is **not** designed to:

- align FASTQ/BAM or call peaks from reads;
- replace DiffBind / csaw differential testing;
- infer causal enhancer–gene regulation;
- call chromatin loops;
- perform RNA-seq differential expression;
- treat sharp/focal peaks (e.g. TF ChIP-seq) as broad domains: the `Breadth`
  axis is protected by the `min_broad_width_bp` sharp-peak guard, and
  sharp-peak data should use the `Intensity` axis only.

For quantitative claims, the input BigWigs must be comparable (CPM / RPGC /
spike-in normalization). Native breadth calls depend on the upstream peak /
domain caller and boundary quality. `SignalDispersion` describes spatial
signal organization, not physical chromatin conformation. Annotated
nearest/overlapping genes are candidates, not causal targets.

All `Super` labels are relative to the analysed candidate-domain universe and
the selected cutoff. They do not by themselves establish a functional
super-enhancer. Integrated Intensity is an area under the signal track and is
therefore affected by interval length. CPM or RPGC corrects library scale but
does not establish absolute comparability when a perturbation causes a global
occupancy shift; spike-in normalization or appropriately qualified relative
interpretation is then required. Use the same peak caller and parameters across
samples because sequencing depth and boundary calling affect native breadth.
The candidate-domain construction rule is also part of the analysis: report
`min_reps`, promoter filtering, and stitching, and assess condition-specific
domains when they are central to the biological question.

Direct local BigWig import uses the `rtracklayer` UCSC backend, which is not
available on Windows. Windows users can build the `SummarizedExperiment` on
Linux or macOS, save it with `saveRDS()`, and run downstream calling,
annotation, and visualization on Windows.

------------------------------------------------------------------------

## Installation

``` r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "GenomicRanges", "SummarizedExperiment", "GenomicFeatures",
  "org.Hs.eg.db", "org.Mm.eg.db"
))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("zying106/epiPortrait", build_vignettes = TRUE)
```

------------------------------------------------------------------------

## Core workflow

### Module 1: Build the Portrait Matrix

**Required inputs.**

| Input | Needed for | How to obtain |
|:------|:-----------|:--------------|
| `bw_path` (BigWig) | Intensity, SignalDispersion, all signal assays | Normalized, comparable tracks (CPM/RPGC/spike-in) |
| `consensus_peaks` (GRanges) | **the shared candidate-domain universe** — Intensity-Super is computed *on these intervals* | From peak-caller bed/narrowPeak via `rtracklayer::import()` + `get_consensus_peaks()` |

**Two different "peak files" are involved — do not confuse them.**

1. **Candidate-domain universe** (`consensus_peaks`, one shared set for all
   samples): defines *which regions you analyze*. It is derived from your peak
   caller's output (MACS2, SICER, epic2, ... bed/narrowPeak), the same upstream
   input ROSE consumes. **Every analysis needs this — including single-sample
   Intensity-only** — because without it there are no intervals to integrate
   signal over.
   - *Single sample*: `consensus_peaks` is simply that sample's own peak file.
   - *Multiple samples*: `consensus_peaks` is the consensus/union of all
     samples' peak calls (`get_consensus_peaks()`).
   - Exception: you may import an externally defined domain set (a fixed
     tiling, a database interval table, or published domains) instead.

2. **Per-sample native peaks** (`peak_path`, one file per sample): each
   sample's *own* peak calls, used to measure how wide the native enrichment is
   per sample. **Only the `Breadth` axis needs this** — Intensity-only
   analyses use the candidate universe alone and omit `peak_path`.

**Bottom line: peak files are always required** (to build the candidate
universe) — the only difference is how many and what for. Intensity-only uses
them to define the domains; Breadth additionally uses each sample's own peaks
for width measurement. There is no analysis mode that runs from BigWig alone.

**Single-sample, Intensity-only.** One `bw_path` + one candidate-domain
`GRanges` (read from that sample's peak file, or any external domain set).
Note this still requires the peak file — it just serves as the single-sample
candidate universe:

``` r
my_domains <- rtracklayer::import("sample1_narrowPeak.bed")  # candidate universe (required)

samples1 <- data.frame(
  SampleID  = "S1",
  Condition = "Control",
  bw_path   = "sample1.bw"
)
se1 <- build_portrait_matrix(samples1, consensus_peaks = my_domains)
se1 <- call_super_domains(se1, feature = "Intensity", mode = "per_sample")
```

With one sample there is no replicate support, so `Intensity-Super` reflects a
single-sample rank state (use `mode = "per_sample"`).

``` r
library(epiPortrait)

# 1. Prepare sample sheet
samples <- data.frame(
  SampleID  = c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Treat_1", "Treat_2", "Treat_3"),
  Condition = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"),
  bw_path   = c("ctrl1.bw", "ctrl2.bw", "ctrl3.bw", "treat1.bw", "treat2.bw", "treat3.bw"),
  peak_path = c("ctrl1.bed", "ctrl2.bed", "ctrl3.bed",
                "treat1.bed", "treat2.bed", "treat3.bed")
)

# 2. Get consensus peaks from your peak files
peak_list <- lapply(samples$peak_path, rtracklayer::import)
consensus_peaks <- get_consensus_peaks(peak_list, min_reps = 2)

# 3. Optionally stitch proximal peaks into macro-domains
# macro_domains <- stitch_epi_peaks(consensus_peaks, stitch_distance = 12500)

# 4. Build the portrait matrix (peak_path is optional, for Breadth-Super)
se <- build_portrait_matrix(
  sample_sheet    = samples,
  consensus_peaks = consensus_peaks,
  workers         = 4
)

# 5. Normalization is OFF by default: input BigWigs should already be
# quantitatively comparable (CPM/RPGC/spike-in). See vignette for options.
se <- normalize_portrait(se, method = "None")
```

`normalize_portrait()` intentionally does not offer count-based TMM for
continuous BigWig signal. `TotalSignal` and Quantile normalization are
explicitly warned alternatives for sensitivity analysis; applied methods,
assumptions, factors and sample totals are recorded in
`metadata(se)$normalization`, and repeated post-hoc normalization is rejected.

**Key Parameters:**

- **`workers`**: Parallel threads via `BiocParallel`.
- **`on_disk`**: Store assays as HDF5 for mammalian-scale processing.
- **`cache_dir`**: Cache per-sample coverage views as RDS to skip BigWig I/O.
- **`custom_features`**: Named list of functions for user-defined dimensions.
  Each function receives one domain's sanitized per-base numeric coverage and
  returns one finite numeric value or `NA_real_`; names must be unique and must
  not collide with built-in features.
- **`negative_policy`**: `"error"` (default), `"clip_zero"`, or `"allow"` for
  negative signal handling.

**Output assays (dynamic, per-sample):**

| Assay | Description |
|:---|:---|
| `Intensity` | Total signal abundance (area under BigWig curve) |
| `SignalDispersion` | Signal-weighted genomic SD within the domain (architecture descriptor) |
| `NativeMaxPeakWidth` | Max native peak width mapped to the domain (needs `peak_path`) |
| `NativeOccupiedWidth` | Sum of reduced native peak widths inside the domain (needs `peak_path`) |
| `NativePeakCount` | Number of native peaks overlapping the domain (needs `peak_path`) |

The static genomic interval length is stored in `rowData(se)$IntervalWidth`.
Per-sample native peak files (via the optional `peak_path` column) enable
Breadth-Super calling; they are stored in `metadata(se)$native_peaks`.

### Candidate domain universe for multi-condition studies

For multi-condition analyses (e.g. Control vs Treatment), build the candidate
universe as the **union of per-condition replicate consensuses** so
condition-specific domains are retained for transition / gain-loss analysis,
then optionally stitch:

``` r
preset <- get_mark_preset("H3K27ac")   # mark-aware defaults, e.g. stitch_distance

ctrl_consensus <- get_consensus_peaks(ctrl_peak_list, min_reps = 2)
treat_consensus <- get_consensus_peaks(treat_peak_list, min_reps = 2)
candidate_domains <- GenomicRanges::reduce(c(ctrl_consensus, treat_consensus))

if (preset$stitch_distance > 0) {
  candidate_domains <- stitch_epi_peaks(candidate_domains,
                                        stitch_distance = preset$stitch_distance)
}
```

A single global `min_reps = 2` across all samples would let rare
condition-specific loci fall below the replicate threshold; the
per-condition-then-union pattern keeps them.

### Module 2: Super-Domain Calling

The canonical workflow for condition comparison is per-group, replicate-aware
calling (each replicate called independently, then aggregated per condition):

``` r
se <- call_super_domains(
    se,
    feature = "Intensity",
    mode = "per_group",
    group_var = "Condition"
)
table(rowData(se)$Intensity_Domain_Type)
```

For a single-condition study (no `Condition` grouping), the default
`mode = "global_consensus"` aggregates all replicates into one call.

**Key Parameters:**

- **`feature`**: canonical Super axes — `"Intensity"` → Super-Element,
  `"Breadth"` → peak-level Breadth-Super. Secondary/exploratory (not canonical
  Super axes): `"SignalDispersion"` (within-domain architecture descriptor),
  `"IntervalWidth"` (one static ranking of the shared coordinate frame; no
  replicate support and no `per_group` / `per_sample` mode).
- **`method`**: `"elbow"` (default, max perpendicular distance) or `"tangent"`
  (ROSE-inspired tangent-optimization inflection; Whyte et al., 2013 *Cell*).
  Note that the tangent implementation is a geometric variant, not a bit-exact
  reproduction of ROSE's `calculate_cutoff()`.
- **`log_transform`**: `NULL` (auto per-feature, default), `TRUE`, or `FALSE`
  (rank on the raw feature scale).
- **`mode`**: `"global_consensus"` (default; alias `"consensus"`),
  `"per_group"`, or `"per_sample"`. In `"per_group"` (recommended for
  multi-condition studies) each replicate is called independently on its own
  ranked feature distribution, then replicate calls are aggregated within each
  condition using the selected support rule; group-mean ranks are stored for
  visualization only and never determine the group call.
- **`n_bootstrap`**: Cutoff stability interval + success rate (Intensity:
  candidate domains; Breadth: per-replicate native PeakWidth distribution).
  A resampling-stability diagnostic, not a classical CI.
- **Breadth only**: `min_peak_overlap_fraction` (default `0.5`, unique
  peak-to-domain mapping), `min_broad_width_bp` (default `500`; sharp-peak
  guard, see below) and `valid_chroms` (allowed chromosomes).
- **Sharp-peak guard (`min_broad_width_bp`)**: if a replicate's widest
  eligible native peak is narrower than this floor, the width distribution is
  in the sharp-peak regime and the replicate provides no Broad evidence
  (`Uncertain`, with a warning) rather than confidently labelling ~200 bp
  peaks as "broad". Applied to both the inflection and `quantile_cutoff`
  paths; set `min_broad_width_bp = NULL` to disable. Sharp-peak data (TF
  ChIP-seq, narrow marks) should use `feature = "Intensity"` only and skip
  `Breadth` and stitching.

**Optional H3K27ac benchmark against ROSE.** For a head-to-head comparison
with the ROSE super-enhancer pipeline (Whyte et al., 2013), use
`method = "tangent"` with `log_transform = FALSE` (raw signal, ROSE-style
scale). This is an optional benchmark setting, not the package's primary
workflow. Use the reference ROSE implementation for claims about numerical
agreement with ROSE.

**Breadth-Super is a peak-level call.** Each replicate's genome-wide eligible
native PeakWidth distribution is cut by an elbow/inflection; broad peaks are
mapped to the shared domains by unique assignment and aggregated across
replicates. It is decoupled from consensus construction.

Breadth calling also records peak presence independently of the canonical
Super/Typical taxonomy:

``` r
evidence <- get_breadth_evidence(se, group = "Control", long = TRUE)
head(evidence)
```

The replicate states are `Broad`, `Typical`, `PeakAbsent`, and `NoCall`.
`PeakAbsent` requires a valid replicate-level width call and no eligible native
peak overlapping the domain. A domain touched only by ambiguous or
threshold-failing peaks remains `NoCall`. Therefore `PeakAbsent` means absence
from the supplied peak calls under the stated pipeline—not biological proof
that the chromatin domain disappeared. Group summaries are stored as
`Breadth_PresenceStatus__<group>`, `Breadth_PresenceFraction__<group>`,
`Breadth_AbsenceFraction__<group>`, and
`Breadth_N_Assessable__<group>`; fractions use all replicates as the
denominator, so technical no-calls cannot inflate evidence.

### Combined taxonomy

For multi-condition studies (the main biological use case), use the
`per_group` mode so each replicate is called independently and replicate calls
are aggregated within each condition:

``` r
se <- call_super_domains(se, feature = "Intensity",
                         mode = "per_group", group_var = "Condition")
se <- call_super_domains(se, feature = "Breadth",
                         mode = "per_group", group_var = "Condition")
se <- combine_superdomain_calls(se, group_var = "Condition")
table(rowData(se)$Combined_Class__Control)  # Intensity-Super / Breadth-Super / Dual-Super / Typical / Uncertain
```

### Module 3: Condition Transitions

``` r
se <- call_super_domains(se, feature = "Intensity",
                         mode = "per_group", group_var = "Condition")
se <- compare_superdomains(se, group_var = "Condition",
                           ref_group = "Control", target_group = "Treatment")
table(rowData(se)[[
  "Intensity_Transition__relative__Control_vs_Treatment"]])
```

The default `cutoff_scope = "relative"` calls each group with its own cutoff;
the output classes `Relative_Prominence_Up` / `Relative_Prominence_Down`
express a *relative rank-state* change and must **not** be interpreted as
absolute signal gain/loss. Use `cutoff_scope = "reference"` or `"pooled"` for
a common-scale `Gain` / `Loss` label when the BigWigs are quantitatively
comparable. These labels still do not establish an absolute biochemical
change. `Uncertain` is assigned when a group call was not reliable. Pair and
cutoff scope are encoded in every canonical transition column; no mutable
latest-result alias is created.

### Module 4: Annotation & Visualization

``` r
# Domain-aware annotation (nearest TSS / promoter / gene-body / contained)
se <- annotate_epi_domains(se, genome = "hg38")
candidates <- get_domain_genes(se, group = "Control")

# One row per domain, one row per domain-gene pair, and raw relationship detail
head(metadata(se)$annotation_summary)
head(metadata(se)$domain_gene_links_dedup)
head(metadata(se)$domain_gene_links)

# External annotation tools are supported through a standard long table.
# Convert software-specific labels before import; epiPortrait does not guess
# ABC or loop-workflow output schemas.
external_links <- data.frame(
  domain_id = rownames(se)[1],
  gene_id = "1017",
  gene_symbol = "CDK2",
  relation_type = "bedpe_promoter_contact",
  bedpe_record_id = "loop_0001",
  external_relation = "ABC enhancer-gene link",
  abc_score = 0.037
)
se_external <- import_domain_annotations(
  se, external_links, source = "ABC_loop_links", mode = "append")

# Hockey-stick ranking plot
se <- call_super_domains(se, feature = "Intensity")
plot_hockey_stick(se, feature = "Intensity")

# Single-domain feature profile (raw values across conditions)
plot_domain_feature_profile(se, peak_id = rownames(se)[1], group_var = "Condition")

# QC
plot_portrait_pca(se, feature = "Intensity", group_var = "Condition")
plot_portrait_correlation(se, feature = "Intensity")
```

### Module 4.5: Differential Domain Analysis

`compare_superdomains()` reports changes between rank-based states.
`analyze_differential_domains()` instead fits a limma model to log-transformed
continuous signal:

``` r
se <- analyze_differential_domains(
  se, feature = "Intensity",
  ref_group = "Control", target_group = "Treatment")
table(rowData(se)[[
  "Intensity_Diff__Control_vs_Treatment__DiffStatus"]])

plot_differential_volcano(
  se, feature = "Intensity", result_name = "Control_vs_Treatment", label_n = 10)
```

For a paired or blocked design, include the blocking variable in the model and
specify a valid design coefficient:

``` r
se <- analyze_differential_domains(
  se, feature = "Intensity",
  design = ~ Patient + Condition,
  contrast = "ConditionTreatment")
```

This is a limma analysis of normalized continuous track values, not a
fragment-count model. Use per-domain counts with DiffBind, csaw, edgeR, or
DESeq2 when count-based inference is required. Simple two-group models need
replication to estimate residual variance.

### Module 5: Functional Interpretation (GO ORA / GSEA)

Links domain phenotypes and remodeling modes to biological programs. The
background is always the genes linked to the object's own domains
(`get_domain_gene_universe()`), never all organism genes; every tested term
gets an explicit log2 odds ratio (ORA) or NES (GSEA) so groups of different
sizes remain comparable. Requires `clusterProfiler`, `GO.db` and an OrgDb
(all in `Suggests`).

``` r
# Discrete phenotype: which programs characterize super domains?
se <- call_super_domains(se, feature = "Intensity", verbose = FALSE)
sel <- rowData(se)$Intensity_Domain_Type == "Intensity_Super_Element"
se <- enrich_epi_domains(se, domains = sel, mark = "H3K27ac")
get_epi_enrichment(se)

# Continuous remodeling: rank genes by the statistic, do not cut a phenotype.
se <- analyze_differential_domains(se, ref_group = "Control",
                                   target_group = "Treatment")
se <- enrich_epi_domains(
  se, method = "GSEA",
  score_col = "Intensity_Diff__Control_vs_Treatment__t")
get_epi_enrichment(se)

# Compare phenotypes on ONE shared universe (log2 OR +/- 95% CI)
se <- compare_epi_enrichment(se, sets = list(
  Super   = sel,
  Typical = rowData(se)$Intensity_Domain_Type == "Intensity_Typical"))
plot_epi_enrichment(se)
```

`as_enrich_result(res)` returns the stored clusterProfiler object, so the
whole enrichplot ecosystem (dotplot, cnetplot, emapplot) works without
reimplementation; `plot_epi_enrichment()` adds only the comparative
effect-size heatmap.

Only GO is wrapped on purpose. For KEGG / Reactome, pass the same domain-derived
foreground and background to the native functions:

``` r
genes <- get_domain_genes(se, domains = sel, unique_genes = TRUE)
universe <- get_domain_gene_universe(se, mark = "H3K27ac")
clusterProfiler::enrichKEGG(gene = genes$gene_id, universe = universe,
                            organism = "hsa", pvalueCutoff = 1)
```

------------------------------------------------------------------------

## Quick Start with Example Data

``` r
library(epiPortrait)
data(example_se)

se <- normalize_portrait(example_se, method = "None")
se <- call_super_domains(se, feature = "Intensity",
                         mode = "per_group", group_var = "Condition",
                         verbose = FALSE)
se <- call_super_domains(se, feature = "Breadth",
                         mode = "per_group", group_var = "Condition",
                         verbose = FALSE)
se <- combine_superdomain_calls(se, group_var = "Condition")
table(rowData(se)$Combined_Class__Control)
```

For an optional H3K27ac / ROSE-style benchmark, use
`method = "tangent", log_transform = FALSE` (see the parameter notes above);
the default `elbow` method is the main analysis.

See the [vignette](vignettes/epiPortrait.Rmd) for the full workflow.

------------------------------------------------------------------------

## Contact

**Ying ZHANG**\
Zhejiang University

Issues and feature requests:
<https://github.com/zying106/epiPortrait/issues>

------------------------------------------------------------------------

## Citation

> ZHANG Y. (2026). *epiPortrait: Replicate-Aware Epigenomic Domain Profiling*.
> R package version 0.99.4.
> <https://github.com/zying106/epiPortrait>

The package citation is generated from `DESCRIPTION` and can be retrieved
programmatically with `citation("epiPortrait")`; the BibTeX entry below mirrors
the current development version.

``` bibtex
@Manual{epiPortrait,
  title  = {epiPortrait: Replicate-Aware Epigenomic Domain Profiling},
  author = {Ying ZHANG},
  year   = {2026},
  note   = {R package version 0.99.4},
  url    = {https://github.com/zying106/epiPortrait}
}
```

------------------------------------------------------------------------
