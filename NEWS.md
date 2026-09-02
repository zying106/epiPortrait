# epiPortrait 0.99.1

* Made elbow calling orientation-aware so inverse/left-loaded ranked curves
  produce `no_call` instead of an implausibly large Super set.
* Propagated both `NA` and the explicit `Uncertain` class through combined and
  per-sample transitions.
* Fixed two-group limma analysis for grouping factors with unused levels.
* Added genome/seqlevel compatibility checks for signal extraction and
  cross-mark integration; epigenomic cross-mark overlaps are now explicitly
  strand-agnostic.
* Modernized the package citation and strengthened caller input validation.

Changes in version 0.99.0 (2026-08-10)
---------------------------------------

Cross-mark integration and per-sample transitions (2026-08-26)

* New `integrate_cross_mark(se_a, se_b, ...)`: overlays the per-sample
  super-domain calls of object B onto the (possibly DIFFERENT) domain
  universe of object A. Each object A domain is classified from the
  object B calls of the object B domains that overlap it (`any_super` or
  `majority` aggregation, optional `min_overlap_bp`). Enables the two-axis
  design (e.g. enhancer-mark Intensity on the promoter mark's domains, or
  boundary-rescue of promoter-proximal enhancers). Adds
  `<mark_b>_Call__<tp>` + `<mark_b>_NOverlaps__<tp>` rowData columns, a
  long-format table in `metadata(se_a)$cross_mark_integration`, and
  provenance in `metadata(se_a)$cross_mark_provenance`.
* New `transition_matrix_per_sample(se, feature/call_fmt, timepoints,
  ref)`: pairwise transitions between per-SAMPLE call columns (mode
  `per_sample`, one sample per time point), which
  `compare_superdomains()`/`compare_superdomain_classes()` do not support
  (they need per-GROUP columns). Adds
  `<prefix>_SampleTransition__<ref>_vs_<target>` rowData columns with
  `compare_superdomain_classes()`-consistent labels (Persistent_<value> /
  <from>_to_<to> / Uncertain) and per-pair counts tables in
  `metadata(se)$transitions`. With `ref = NULL` all ordered time-point
  pairs are computed.

BEDPE contact-score aggregation in annotate_epi_domains() (2026-08-26)

* New `bedpe_score_col` argument: an optional BEDPE column (integer index
  >= 7 or column name) holding a per-record contact STRENGTH, e.g. the
  FitHiChIP TMM-normalized contact frequency. Motivated by real datasets
  (GSE251898) where the published TMM BEDPE loop SET is identical across
  conditions and only the score varies — count-based contact integration
  alone is invariant there by construction.
* When supplied: `metadata(se)$domain_gene_links$contact_score` carries the
  per-evidence-row score; the dedup pair table and the annotation summary
  gain `bedpe_contact_score` (summed over UNIQUE supporting records; domain
  level = sum over all contacted records); `rowData(se)$bedpe_contact_score`
  mirrors the domain-level value. `bedpe_provenance$score_column` records
  what was used.
* Semantics: no BEDPE evidence -> 0; evidence without usable scores (or runs
  without `bedpe_score_col`) -> NA, so "no 3D evidence" and "unscored 3D
  evidence" remain distinguishable.
* `get_domain_genes()` gains a `bedpe_contact_score` column and a new
  `rank_by = c("tier", "bedpe_score")` argument: "bedpe_score" prioritizes
  candidate genes primarily by descending 3D contact strength
  (dominant-loop-style target assignment), with the evidence tier and
  expression as tie-breakers.
* Default behaviour (`bedpe_score_col = NULL`) is unchanged apart from the
  two new NA-valued columns.

Statistical hardening of analyze_differential_domains() (2026-08-23)

* `analyze_differential_domains()` now fits the limma mean-variance trend by
  default (`trend = TRUE`): integrated intensity exhibits a strong
  mean-variance dependence across domains, and without the trend low-signal
  domains were over-called as significant. The trend is disabled automatically
  with a warning when fewer than 20 domains pass `min_signal` filtering.
* New `robust = FALSE` argument exposes limma's robust empirical-Bayes
  moderation for outlier-resistant inference.
* Both settings are recorded in `metadata(se)$differential_domains`
  provenance; invalid values are rejected.
* Documented the scale-dependence of the log(x + 1) pseudo-count offset
  (CPM vs RPGC vs spike-in units) and its interaction with `min_signal`.
* `eBayes(trend = TRUE)` changes P-values slightly relative to previous
  releases; re-run differential analyses after upgrading.
* Bug fix: a two-element character `contrast` (coefficient-name pair,
  documented as `coef1 - coef2`) was silently ignored and the uncontrasted
  fit was returned. It is now converted to numeric contrast weights, with a
  clear error when a name does not match a design coefficient. Supplying
  `design` without any `contrast` is now also rejected instead of silently
  returning an uncontrasted fit.
* Bug fix: a user-supplied formula `design` (e.g. `~ Condition` or
  `~ 0 + Condition`) crashed inside limma with an obscure internal error,
  because the formula object was passed to `lmFit()` unconverted. Formulas
  are now materialized with `model.matrix()` on colData before fitting, so
  custom paired/batched designs are actually usable.

Expert review fixes (2026-08-14)

* normalize_portrait(): the completion message no longer uses a sprintf()
  format string without a placeholder (this produced a runtime warning
  "argument out of range" and truncated the message to
  "Normalization complete! Intensity adjusted; SignalDispersion "). The
  "(bp-scale spatial descriptor) is left unchanged" note is now emitted
  correctly.
* plot_hockey_stick(): when the object was called with mode = "per_group"
  and `group` is omitted, the error now states that the calls are
  per-group and lists the available condition groups, instead of the
  misleading "Run call_super_domains first" message.
* annotate_epi_domains(): the strand-aware signed nearest-TSS distance is
  now fully vectorized (one match() over the gene model instead of a
  per-domain O(genes) subsetting loop); output is unchanged.
* build_portrait_matrix(): the batch import path accumulates Intensity into
  a double (`numeric()`) vector instead of `integer()` (values were never
  truncated because R promotes on assignment, but the type was misleading
  for continuous CPM/RPGC signal).
* plot_peak_track(): replaced scale_fill_brewer("Set1") (capped at 9
  groups) with the package's own condition palette, consistent with the
  rest of the visualization suite.
* call_super_domains(): documented the intentional asymmetry of
  min_valid_replicates under support_rule = "fraction" (a low threshold,
  e.g. 0.5 with n=2, resolves to 1 and can call a "Super + no_call" pair
  Super because the user explicitly chose how much support suffices).
* Added regression tests: fractional Intensity preservation through
  build_portrait_matrix() (double precision, cache and non-cache paths
  identical) and the per-group hockey-stick error message.

Documentation repositioning (2026-08-11)

* Package identity repositioned from "super-domain caller" to
  "replicate-aware epigenomic domain profiling": DESCRIPTION title,
  README/vignette narrative, and CITATION now describe the framework as
  quantitative phenotyping of shared epigenomic domains (Intensity x native
  peak breadth + within-domain signal architecture), with super-domain calling
  as one analytical module. ROSE is presented as an optional H3K27ac benchmark;
  the tangent inflection is documented as a geometric variant, not a bit-exact
  ROSE reproduction.

Breadth top-fraction & stability calling (2026-08-11)

* `call_super_domains(feature = "Breadth")` now supports `quantile_cutoff`
  as an explicit top-fraction opt-in, symmetric with Intensity: the cutoff is
  the `quantile_cutoff` quantile of each replicate's genome-wide eligible
  native PeakWidth distribution (top `100*(1-quantile_cutoff)`% widest peaks
  are Broad; e.g. `0.95` = top 5%). It never auto-falls-back to, or from, the
  data-driven elbow/tangent inflection, so no-call semantics of the default
  path are unchanged. `quantile_cutoff` and `min_replicate_support` are now
  validated.
* `n_bootstrap` is supported for Breadth (symmetric with Intensity): a
  per-replicate within-set cutoff stability interval and success rate are
  stored in `get_call_provenance(se, "Breadth")$replicates[[s]]`. As with
  Intensity, this is a resampling-stability diagnostic whose information
  content shrinks with the size of the ranked set; the dominant uncertainty
  for Breadth (upstream native-peak set composition) is assessed separately
  via the peak-set perturbation check in `validation/02`.

Freeze pre-check fixes (deep review 2026-08-10)

* Z-score normalization REMOVED from normalize_portrait(): it is a row-wise
  display/clustering transform that destroys the cross-domain magnitude ranking
  required by Super calling and must never overwrite the canonical Intensity
  assay. plot_portrait_pca() and plotting layers already apply their own
  scaling for display. normalize_portrait() description now states that
  SignalDispersion is preserved in native bp units.
* plot_domain_landscape() now accepts and respects group_var (was hard-coded to
  "Condition"), matching call_super_domains() / get_replicate_calls() /
  combine_superdomain_calls() for custom group variables.
* filter_promoter_peaks() now enforces domain/TxDb seqlevel compatibility
  (zero shared -> error, partial -> warning) and validates upstream/downstream
  parameters, consistent with annotate_epi_domains().
* filter_blacklist() now fails loudly on zero shared seqlevels (chr1-vs-1
  mismatch) instead of silently removing 0 peaks, and warns on partial overlap.
* validate_epiportrait_object() now verifies the replicate_call_matrix contract
  (rows = domains, rownames = domain IDs, columns a subset of samples).
* stitch_epi_peaks(): min.gapwidth = stitch_distance + 1 so a gap exactly equal
  to the documented "maximum distance" is merged (off-by-one fix).
* get_replicate_calls(group = NULL) preserves the original sample order instead
  of the alphabetical group order used during column binding.
* Removed package-root Rplots.pdf; added ^Rplots[.]pdf$ to .Rbuildignore.
* Updated stale comments (annotation resolver, NEWS font wording).

Nature-style publication theming (nature-figure skill)

* All plots now use a Nature-style publication theme: no background grid lines,
  thin black axis lines, platform-safe sans-serif, compact base sizes, and
  frameless
  legends. Affected: .epi_theme_publication() (hockey-stick, PCA, landscape,
  transition, replicate support, class composition, feature profile),
  plot_peak_track() and plot_portrait_correlation() (both switched from
  theme_minimal to the grid-free classic theme).
* Export via save_epiportrait_figure() produces editable-text SVG / PDF
  (vector fonts, no outlined paths).


### Intensity replicate-calling transparency (design 2026-08-10)


* New get_replicate_calls(): exposes the per-domain x per-replicate super-domain
  call matrix underlying the group-level replicate-support call (matrix or long
  format). For mode = "per_group" it reads the stored provenance
  (replicate_call_matrix is now saved per group); for mode = "per_sample" it
  reconstructs from the <feature>_Call__<sample> rowData columns.
* call_super_domains(mode = "per_group") now stores the full domain x replicate
  call matrix in metadata(se)$superdomain_calls[[feature]]$groups[[g]]$
  replicate_call_matrix, making 2/3-style support fully auditable.
* Documentation clarified (mode descriptions): per_group calls EACH replicate
  independently on its own ranked feature distribution, then aggregates
  replicate calls by the support rule; group-mean ranks are stored for
  visualization only and never determine the group call.
* Fixed min_valid_replicates documentation: NULL (default) resolves to the
  support-rule-required number (majority: floor(n/2)+1; all: n; fraction:
  ceiling(fraction*n)), not "default 1".


### Species/genome & annotation audit fixes (deep audit 2026-08-10)


* P0-1: filter_blacklist() no longer depends on the TxDb/OrgDb resolver.
  Blacklist resolution is independent: built-in hg38/hg19/mm10/mm9 blacklists
  load directly from inst/extdata, so mm9 works and blacklist filtering never
  requires a TxDb Suggests package. Custom BEDs via blacklist_path.
* P0-2: domain-TxDb seqlevel compatibility is now ENFORCED in
  annotate_epi_domains(): zero shared seqlevels -> error (with domain-only
  seqlevels reported); partial overlap -> warning. A genome-build / chr-vs-1
  naming mismatch is never silently reinterpreted as intergenic domains.
  Compatibility details (shared, domain_without_txdb, txdb_only) are stored in
  annotation_provenance.
* P1-1: min_expression is now genuinely NULL by default. expressed_first on
  TPM/CPM uses a documented convenience threshold of 1; on VST/rlog/
  normalized_counts it requires an explicit user threshold. high_expression_first
  needs no threshold and leaves expression_status NA.
* P1-2: wide RNA expression matrices with ANY unmatched sample column now
  error (previously silently dropped); a subset of ChIP samples is allowed.
* P1-3: wide RNA validates unique non-empty rownames and numeric values.
* P1-4: candidate_priority is re-ranked to a contiguous 1..N ordinal after
  unique_genes collapse (no gaps from dropped rows).
* P1-5: added a real evidence-vs-expression test (BEDPE-supported low-expression
  gene outranks nearest-only high-expression gene).
* P1-6: BEDPE provenance/comment wording unified: records are not modified or
  re-scored; bedpe_record_id links derived evidence back to the source.
* P1-7: gene_symbol is NA when no symbol mapping is available (never a gene_id
  masquerading as a symbol).
* P1-8: a user-provided gene_id_keytype whose mapping fails now warns (not
  silent); gene_id_keytype is recorded in annotation_provenance.
* P1-9: BEDPE seqlevel_style inferred from the anchors (UCSC or custom/unknown),
  not hard-coded.
* P1-10: promoter_upstream / promoter_downstream validated (finite, >= 0).
* P1-11: get_domain_genes() validates domains length/NA, group existence, and
  requires group + expression data for expression_priority != "none".
* P1-12/13: documentation clarifies the gene-level representative TSS and the
  difference between filter_promoter_peaks() (transcript-level promoters) and
  annotate_epi_domains() (gene-level promoter window).
* P1-14: export_epiportrait_results() writes annotation/domain_gene_links.tsv
  and annotation/provenance.txt when domain annotation is present.
* Added 11 audit tests (mm9 blacklist, blacklist/TxDb independence, zero/partial
  seqlevels, wide unmatched samples, rownames/numeric validation, VST threshold,
  invalid group, unique rerank, evidence-vs-expression, annotation export).


### Annotation correctness fixes (deep review 2026-08-10)


* P0-1: .gene_model_from_txdb() no longer reduce()s TxDb genes. One TxDb gene
  is now always one gene-level feature; overlapping genes are never merged into
  a pseudo-gene ("geneA,geneB").
* P0-2: nearest_tss_distance_bp is now the strand-aware SHORTEST interval-to-TSS
  distance: 0 when the TSS lies inside the domain (previously the midpoint
  offset could give several kb for an overlapping TSS).
* P0-3: nearest-TSS hits are backfilled by queryHits() so domains on contigs
  absent from the TxDb stay NA instead of misaligning.
* P0-4: primary_genomic_context simplified to Promoter-associated /
  Gene-body-associated / Intergenic, consistent with the count columns (no more
  over-produced "Mixed" or unbounded "Downstream").
* P0-5: long-format RNA aggregation now uses stats::aggregate (correct gene x
  replicate collapse; tapply length mismatch fixed).
* P0-6: wide RNA matrix columns are matched to epiPortrait samples BY NAME, not
  by position; independent RNA sample IDs require the long-format input.
* P0-7: candidate_priority is now a true integer rank (order() permutation was
  being assigned directly).
* P0-8: BEDPE evidence score uses pmax() so the strongest 3D evidence is never
  overwritten by a weaker relation on the same domain-gene pair.
* P1-1: gene symbols are resolved once over ALL linked gene ids and backfilled
  into every relation row (previously only nearest_tss got symbols).
* P1-2: gene_id_keytype parameter; custom TxDbs require an explicit keytype
  instead of hard-coded ENTREZID.
* P1-3: unique_genes = TRUE now retains the BEST-SUPPORTED domain-gene
  association per gene (sorted by priority first).
* P1-4: bedpe_support_count counts UNIQUE BEDPE records (a record matching via
  both directions is no longer double counted).
* P1-5: BEDPE seqlevels are actually checked against the domains (zero shared
  -> error; partial -> warning).
* P1-6: BEDPE input validation (non-empty chrom, numeric coordinates, start>=0,
  end>start).
* P1-7: BEDPE provenance notes that records are not modified/re-scored and that
  bedpe_record_id links derived evidence back to the source (accurate promise).
* P1-8: BEDPE promoter-contact linking is range-vectorized (one findOverlaps per
  direction + join by record id), scaling to 10^4-10^6 interactions.
* P1-9: expression provenance records the RNA sample IDs (not the ChIP ones).
* P1-10: RNA gene-ID match QC (n_matched / match_fraction in provenance; stop on
  zero matches, warning on <50%).
* P1-11: min_expression = NULL default semantics; expressed_first on non-TPM/CPM
  types requires an explicit user threshold.
* P1-13: ChIPseeker removed from README/vignette install commands and vignette
  section title (DESCRIPTION no longer imports it).
* Genome resolver: custom TxDb seqstyle inferred from seqlevels (no longer
  hard-coded "UCSC"); filter_blacklist() now uses the shared resolver.
* Added 17 edge-case tests (overlapping genes, TSS-inside-domain distance=0,
  contig mismatch, RNA long/wide alignment, candidate rank, BEDPE validation).


### Domain-aware annotation & candidate-gene prioritization (integrated design

2026-08-10)

* New annotate_epi_domains(): genome-flexible, domain-aware annotation of the
  final shared domains (not individual peaks). Adds 7 compact rowData columns:
  primary_genomic_context, nearest_tss_gene_id, nearest_tss_gene_symbol,
  nearest_tss_distance_bp (strand-aware), promoter_overlap_gene_count,
  gene_body_overlap_gene_count, fully_contained_gene_count. Stores a long-format
  domain-gene link table in metadata(se)$domain_gene_links.
* New get_domain_genes(): transparent candidate-gene extraction / ordinal
  prioritization from the stored evidence (evidence first, expression as a
  secondary tie-breaker; no black-box composite score).
* New internal .resolve_genome_resources() / .check_seqlevel_compatibility():
  unified genome resource contract (hg38/hg19/mm10 built-in TxDb/OrgDb, mm9
  blacklist-only, custom TxDb/OrgDb/blacklist) shared by annotate_epi_domains(),
  filter_promoter_peaks() and filter_blacklist().
* BEDPE 3D evidence: optional bedpe_promoter_contact (symmetric / bidirectional
  domain-anchor to gene-promoter contact). BEDPE is read-only: no calling,
  filtering, merging, or FDR re-scoring; bedpe_record_id links derived evidence back to the source; per-record
  bedpe_record_id; provenance in metadata(se)$bedpe_provenance.
* RNA-seq expression integration: user-declared expression_type (TPM/CPM/VST/
  rlog), per-condition median aggregation over replicates, and expression-aware
  prioritization (none / expressed_first / high_expression_first). Expression
  alone never creates a domain-gene link (evidence required first).
* Removed annotate_epi_peaks(): it was a thin ChIPseeker::annotatePeak()
  wrapper and epiPortrait's analysis object is the shared domain, not the
  individual peak. Peak-level annotation remains available via ChIPseeker.
* Provenance: metadata(se)$annotation_provenance, $bedpe_provenance,
  $expression_provenance.
* Added 23 tests (genome resolution, linear annotation, BEDPE contacts, RNA-seq
  integration, candidate prioritization).

Changes in version 0.99.0 (2026-08-09)
---------------------------------------

Pre-submission release for Bioconductor. Metric redesign (v1.0 design docs,
2026-08-09): the taxonomy is simplified to the two-axis architecture
Intensity x Breadth, and Breadth-Super is moved to a peak-level, consensus-free
calling paradigm.

Heterochromatin / broad-repressive mark compatibility (plan 2026-08-10)

* get_mark_preset() now returns mark_class ("active" / "broad_repressive"),
  taxonomy_style ("super_domain" / "repressive_remodeling") and class_labels
  (display aliases). H3K27me3 / H3K9me3 default to stitch_distance = 0: their
  domains are already spatially clustered by the upstream broad-domain caller
  (SICER / epic2 / RECOGNICER / MACS2_broad), so epiPortrait no longer applies
  an extra 5 kb stitching by default (users can pass stitch_distance = 5000
  explicitly).
* Display-only terminology: with mark = "H3K27me3" (or H3K9me3), plot_domain_
  landscape(), plot_class_transition() and plot_domain_class_composition()
  render Intensity-Extreme / Extended-Domain / Dual-Extreme. The internal
  canonical classes (Intensity-Super / Breadth-Super / Dual-Super) are never
  changed, so the API and transition engine stay mark-independent.
* New compute_width_transition(): continuous domain-width expansion/contraction
  descriptor (WidthDelta_bp, log2WidthRatio, WidthDirection) between two
  condition groups, computed from per-domain median native width. This is a
  continuous effect descriptor, not a significance test (no new p-value
  framework at n=2/n=3).
* build_portrait_matrix() records domain-source / caller provenance in
  metadata(se)$domain_provenance when the sample sheet carries optional
  domain_source / domain_caller / caller_version / caller_mode columns.
* Added tests: mark presets, display alias mapping, mark-aware plots,
  compute_width_transition semantics, domain provenance.

Correctness fixes from the deep review (2026-08-09)

* Performance: .native_peak_geometry() is now fully vectorized. The previous
  per-domain R loop (one intersect/reduce per domain) took ~100 s for 3,000
  domains and ~17 min for 26,000 domains, making build_portrait_matrix()
  unusable at mammalian scale. The rewrite clips all peaks to their domains in
  one pintersect() call and reduces only multi-peak domains (single-peak
  domains need no reduce), giving an ~10,000x speedup with identical output.
  Verified on real GBM10 H3K4me3 (26,483 domains): 98.7 s -> 0.01 s for the
  geometry step; full build_portrait_matrix() incl. a 450 MB BigWig runs in
  ~6 s.
* Breadth-Super: a replicate whose native PeakWidth elbow is NOT reliable
  (no_call) now contributes NO evidence (its domains become Uncertain) instead
  of being coerced to "all Typical". Previously no_call was encoded as an
  all-FALSE logical vector, silently downgrading super support.
* Breadth-Super mapping: Unmapped / Ambiguous peaks now never produce Typical
  evidence. Only uniquely-mapped peaks contribute evidence (broad wins over
  typical inside a domain). Previously any-overlap domains could be labelled
  Typical even when their only peak was below the overlap threshold or
  ambiguous.
* build_portrait_matrix(): has_peaks is now always a per-sample logical vector.
  Previously, when no peak_path column existed, has_peaks was a scalar FALSE
  and has_peaks[i] for i>1 evaluated to NA, breaking BigWig-only multi-sample
  mode. has_peaks is also correctly re-subset after fail_action = "drop".
* normalize_portrait(): SignalDispersion (a spatial descriptor in bp) is no
  longer rescaled by library-size factors, quantile or z-score transforms.
  Only Intensity is normalized; SignalDispersion stays in bp and is
  scale-invariant by construction.
* compare_superdomains(): fixed an inverted condition that forced
  support_rule = "majority" for reference/pooled transitions even when the
  stored call used "all" or "fraction".
* call_super_domains(feature = "Breadth"): quantile_cutoff is supported as an
  explicit top-fraction opt-in on the genome-wide native PeakWidth
  distribution (see "Breadth top-fraction calling" above); n_bootstrap is
  supported per replicate (within-set cutoff stability, symmetric with
  Intensity); and min_peak_overlap_fraction is validated to lie in [0, 1].
* Breadth peak-to-domain mapping provenance now records OverlapBp,
  PeakOverlapFraction, DomainOverlapFraction, the row index AND the actual
  Domain_ID (epiDomain_xxxxxx).
* plot_hockey_stick(feature = "Breadth"): new peak-level branch that plots the
  replicate-specific native PeakWidth rank curve from
  metadata(se)$breadth_peak_calls; the replicate is passed via `group`.
* validate_epiportrait_object() now requires Intensity and SignalDispersion
  (native geometry assays remain optional).
* Added unit tests for Breadth no-call propagation, below-threshold Unmapped,
  ambiguous ties, unique broad/typical mapping, broad-over-typical precedence,
  and overlap-fraction provenance; consensus A/B robustness benchmark in
  validation/02_breadth_super_validation.R (Jaccard / PeakBroadCall
  concordance).

Breaking changes (v1.0 metric simplification)

* Assays: RobustHeight and EffectiveWidth (alias Spreading / Central80Width)
  are REMOVED from build_portrait_matrix(). Core assays are now Intensity,
  SignalDispersion, and (when peak files are provided) NativeMaxPeakWidth,
  NativeOccupiedWidth, NativePeakCount. The static interval length stays in
  rowData(se)$IntervalWidth.
* signal_contract metadata no longer carries width_baseline_quantile (the
  central-80% span was removed); SignalSpan80 / .central80_width() are gone.
* call_super_domains(): feature = "Breadth" routes to a NEW peak-level
  Breadth-Super pipeline (see below). RobustHeight / EffectiveWidth /
  Spreading are no longer valid features; use Intensity, SignalDispersion or
  Breadth.
* combine_superdomain_calls(): combines Intensity x Breadth and produces
  Intensity-Super / Breadth-Super / Dual-Super / Typical / Uncertain
  (Width-Super is renamed Breadth-Super).
* get_mark_preset(): broad marks (H3K4me3 / H3K27me3 / H3K9me3) now default to
  primary_feature = "Breadth".

New: peak-level Breadth-Super calling

* Breadth-Super is called at the NATIVE PEAK level, before and independently
  of consensus construction: each replicate's genome-wide eligible native
  PeakWidth distribution is ranked and cut by an elbow/inflection
  (method / min_quality / tie_policy), producing PeakBroadCall per peak.
* Broad peaks are mapped to the shared domain universe by UNIQUE assignment
  (maximum overlap + min_peak_overlap_fraction, default 0.5); a peak cannot
  double-count into neighbouring domains (MappingStatus = Unique / Ambiguous /
  Unmapped).
* Replicate-aware aggregation (support_rule / min_replicate_support /
  min_valid_replicates) yields group-level / global Breadth-Super.
* Sample-specific native peak files are OPTIONAL: provide peak_path in the
  sample_sheet to build_portrait_matrix(). Imported peaks are stored in
  metadata(se)$native_peaks; peak-level calls and mapping provenance are stored
  in metadata(se)$breadth_peak_calls / breadth_peak_mapping.
* New params on call_super_domains(): min_peak_overlap_fraction, valid_chroms.

New: SignalDispersion

* Signal-weighted genomic SD within each domain. Secondary within-domain
  signal architecture descriptor; scale-invariant, spatially sensitive, and
  requires no background baseline. Replaces the removed EffectiveWidth as the
  continuous-profile descriptor, but does NOT define Breadth-Super.

Visualization

* plot_domain_landscape(): y-axis defaults to NativeMaxPeakWidth (falls back to
  NativeOccupiedWidth / SignalDispersion).
* plot_peak_track(): show_effective_width replaced by show_native_occupancy
  (NativeOccupiedWidth band).
* plot_peak_portrait(): features default to Intensity, SignalDispersion,
  NativeMaxPeakWidth.
* Example data (example_se) regenerated with the 5-assay structure and native
  peaks; validation scripts updated (01_signal_dispersion_validation.R,
  02_breadth_super_validation.R).

Changes in version 0.99.0 (2026-08-08)
---------------------------------------

Pre-submission release for Bioconductor. Complete refocus of the package from
the legacy "4D Shape Shift / MANOVA" framework to a BigWig-first
super-domain profiling framework.

Algorithm & science

* P0-1: negative_policy is now applied once via .sanitize_signal() BEFORE any
  feature is computed; clip_zero actually clips negatives (Intensity,
  RobustHeight, EffectiveWidth, custom features all share the sanitized
  vector). NA positions are treated as 0 to preserve the genomic coordinate
  axis. The negative fraction is now computed over signal positions (not
  domain count) and stored in colData(se)$NegativeFraction.
* P0-2: bootstrap resampling preserves the user-selected log scale. The
  effective transform is resolved once; log_transform=FALSE stays FALSE inside
  the bootstrap (previously it silently re-resolved to TRUE for Intensity).
  metadata records log_transform_requested and log_transform_used.
* P0-3: replicate support uses a configurable rule (support_rule =
  "majority" | "all" | "fraction") with a min_valid_replicates requirement.
  no-call replicates are never counted as Super and never inflate support
  (previously Super + no_call gave 100% support). n=2 majority now requires
  2/2.
* P0-4: compare_superdomains() gained cutoff_scope = "relative" (default) |
  "reference" | "pooled". Relative transitions are reported as
  Relative_Prominence_Up/Down and must not be interpreted as absolute
  Gain/Loss; only reference/pooled common-cutoff modes report Gain/Loss.
* P0-5: Spreading (EffectiveWidth) is sensitive to low-level distributed
  background; added width_baseline_quantile to build_portrait_matrix() to
  subtract a per-profile baseline before computing the span. Simulation shows
  a 100 bp peak in a 1000 bp domain spans ~81 bp for background 0-2 with
  q10 correction vs 81-700 bp without (validation/01_signal_width_validation.R).
* Added min_quality argument to call_super_domains() (previously implicit).
* Removed legacy MANOVA and Shape Shift modules: test_global_shape_shift(),
  test_shape_shift(), test_differential_feature(), classify_shape_shift(),
  detect_conformational_outliers(), simulate_power(), plot_power_curve(),
  and their dependent visualization/export/report functions.
* Removed the Skewness assay from build_portrait_matrix() (was signal-value
  skewness, not genomic spatial directionality). Portrait now produces
  Intensity, RobustHeight, and Spreading.
* Added Spreading (per-sample effective signal span = central 80% signal-mass
  span) as a core assay, aliased as EffectiveWidth / Central80Width.
* Moved the static genomic interval length out of the assays into
  rowData(se)$IntervalWidth (alias "Width").
* Height is now RobustHeight (95th percentile of the signal profile), robust
  to single-bin spikes.
* Added find_hockey_inflection() with two inflection heuristics:
  - "elbow" (default): maximum perpendicular distance on [0,1]-normalized
    coordinates; dimensionless and scale-invariant.
  - "tangent": tangent-optimization cutoff ported from the ROSE algorithm
    (Whyte et al., 2013, Cell 153:307-319).
* Added log_transform argument to call_super_domains(): NULL (auto),
  TRUE, or FALSE (rank on raw scale for ROSE-style comparison).
* Added combine_superdomain_calls() for the Intensity-Super / Width-Super /
  Dual-Super / Typical / Uncertain taxonomy. Default width feature is now
  EffectiveWidth (dynamic breadth), not static interval length.
* call_super_domains() now stores full call provenance in
  metadata(se)$superdomain_calls (cutoff, cutoff_CI, quality, status, ...).
* Consensus mode is condition-aware: replicate support is stratified by
  group. per_group mode is replicate-aware (per-replicate calls + group
  support) rather than a single group-mean ranking.
* normalize_portrait() now defaults to method = "None"; TotalSignal and
  Quantile emit warnings that they may remove genuine global biological
  shifts.
* build_portrait_matrix() gained negative_policy ("error" | "clip_zero" |
  "allow") and fail_action ("stop" | "drop") arguments.
* Added check_signal_compatibility() for pre-analysis BigWig QC.
* build_portrait_matrix() now imports BigWig coverage in region batches
  (memory-bounded for mammalian-scale domains).
* Cache keys include full path + size + mtime to avoid collisions and stale
  reads.
* Rewrote vignette and README around the super-domain workflow.
* Updated example_se dataset (Spreading assay, RobustHeight, IntervalWidth).

Bioconductor packaging

* Version bumped to 0.99.0 for first submission.
* Removed LazyData; removed unused Suggests; added TxDb and rmarkdown
  Suggests; unified maintainer metadata via Authors@R.
* R CMD check passes with 0 errors / 0 warnings / 0 notes (including
  --run-donttest).

Changes in version 0.1.0 (2026-05-12)
-----------------------------------------

Initial Release

* Added build_portrait_matrix() core engine for 4D geometric feature extraction
  (Width, Intensity, Height, Skewness) from BigWig coverage with BiocParallel
  support, HDF5 on-disk storage option, and per-BigWig caching.
* Added get_consensus_peaks() and stitch_epi_peaks() for consensus peak
  calling and macro-domain stitching.
* Added call_super_domains() for ROSE-style single-dimension super-element /
  broad-domain / steep-peak identification.
* Added normalize_portrait() with TotalSignal, Quantile, and TMM methods.
  (A row-wise Z-score method was prototyped but removed before release: it is a
  display/clustering transform that would destroy the cross-domain magnitude
  ranking used by Super calling.)
* Added visualization suite: plot_peak_portrait(), plot_peak_track(),
  plot_hockey_stick(), plot_portrait_pca(), plot_portrait_correlation().
* Added annotate_epi_peaks() for ChIPseeker-based gene mapping with lazy
  TxDb loading for hg38, hg19, and mm10.
* Added filter_promoter_peaks() for promoter-proximal peak exclusion.
