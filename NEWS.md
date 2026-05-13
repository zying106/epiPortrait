Changes in version 0.1.0 (2026-05-12)
-----------------------------------------

Initial Release

* Added build_portrait_matrix() core engine for 4D geometric feature extraction
  (Width, Intensity, Height, Skewness) from BigWig coverage with BiocParallel
  support, HDF5 on-disk storage option, and per-BigWig caching.
* Added test_global_shape_shift() for MANOVA-based multivariate shape shift
  testing with support for simple two-group, paired/blocked, and complex
  multi-factor formula designs. Returns Pillai trace, partial eta-squared,
  and BH-adjusted P-values.
* Added test_differential_feature() for per-dimension differential testing
  with limma (empirical Bayes), t.test, and Wilcoxon backends.
* Added test_temporal_shape_shift() for polynomial-contrast MANOVA across
  ordered conditions (timepoints, doses).
* Added classify_shape_shift() translating abstract Pillai scores into six
  biologically interpretable categories (Concentration, Flattening, Polarity
  Shift, Global Gain, Global Loss, Complex).
* Added detect_conformational_outliers() for 4D Mahalanobis distance-based
  outlier detection with robust MCD and classic covariance estimators.
  Supports single-sample query-vs-reference analysis.
* Added plot_outlier_manhattan() for genomic visualization of conformational
  outlier significance across chromosomes.
* Added get_consensus_peaks() and stitch_epi_peaks() for consensus peak
  calling and macro-domain stitching.
* Added call_super_domains() for ROSE-style single-dimension super-element /
  broad-domain / steep-peak identification.
* Added normalize_portrait() with TotalSignal, Quantile, Z-score, and TMM
  methods.
* Added visualization suite: plot_portrait_volcano() (4D volcano),
  plot_peak_portrait() (radar chart), plot_shift_comparison() (group
  overlay tracks), plot_peak_track() (faceted coverage), plot_hockey_stick()
  (ROSE ranking), plot_portrait_pca(), plot_portrait_correlation(),
  plot_power_curve().
* Added annotate_epi_peaks() for ChIPseeker-based gene mapping with lazy
  TxDb loading for hg38, hg19, and mm10.
* Added enrich_shape_shifted() for per-shape-class GO/KEGG enrichment via
  clusterProfiler.
* Added summarize_findings() for automated structured biological narrative
  synthesis from classified shape-shift results.
* Added simulate_power() and plot_power_curve() for MANOVA power analysis
  and prospective experimental design.
* Added export_shifted_bed() and export_igv_session() for result export.
* Added render_report() for one-click parameterised HTML report generation.
* Added filter_promoter_peaks() for promoter-proximal peak exclusion.
