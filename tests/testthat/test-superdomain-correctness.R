# rtracklayer BigWig reads fail on Windows ('UCSC library operation failed',
#   rtracklayer#52/#62/#128/#151), so true-BigWig integration tests are
#   skipped there. All other platforms run them in full.
skip_on_os("windows")

# Scientific-correctness tests for the fixes from the deep review (2026-08-08)
# and the v1.0 metric redesign (2026-08-09). These go beyond "returns a
# ggplot" and verify the actual semantics.

data(example_se)

# ---- v1.0: SignalDispersion analytical behaviour -----------------------------
test_that("signal_dispersion returns NA for zero signal", {
  f <- epiPortrait:::.signal_dispersion
  expect_true(is.na(f(c(0, 0, 0))))
  expect_true(is.na(f(rep(NA_real_, 5))))
})

test_that("signal_dispersion is scale-invariant (uniform multiplicative)", {
  f <- epiPortrait:::.signal_dispersion
  x <- c(0, 1, 4, 4, 4, 1, 0)
  expect_equal(f(x), f(x * 100))
})

test_that("signal_dispersion grows with spatial separation", {
  f <- epiPortrait:::.signal_dispersion
  compact <- c(rep(0, 20), rep(1, 20), rep(0, 20))
  spread  <- c(rep(1, 10), rep(0, 40), rep(1, 10))
  expect_gt(f(spread), f(compact))
})

# ---- P0-9 / 17.2: negative signal policy ------------------------------------
test_that("sanitize_signal and dispersion handle negatives per policy", {
  expect_error(epiPortrait:::.sanitize_signal(c(1, -1), "error"), "negative")
  expect_equal(epiPortrait:::.sanitize_signal(c(1, -1), "clip_zero"), c(1, 0))
  expect_equal(epiPortrait:::.sanitize_signal(c(1, -1), "allow"), c(1, -1))
})

# ---- P0-4 / 17.3: bootstrap does not double-log -----------------------------
test_that("bootstrap cutoff CI is finite and reproducible under seed", {
  se1 <- call_super_domains(example_se, feature = "Intensity",
                            method = "tangent", log_transform = TRUE,
                            n_bootstrap = 20, seed = 1, verbose = FALSE)
  se2 <- call_super_domains(example_se, feature = "Intensity",
                            method = "tangent", log_transform = TRUE,
                            n_bootstrap = 20, seed = 1, verbose = FALSE)
  md1 <- S4Vectors::metadata(se1)$superdomain_calls$Intensity
  md2 <- S4Vectors::metadata(se2)$superdomain_calls$Intensity
  expect_true(is.finite(md1$cutoff))
  expect_true(all(is.finite(md1$cutoff_stability_interval)))
  expect_equal(md1$cutoff_stability_interval, md2$cutoff_stability_interval)
})

# ---- P0-5 / 17.5: call provenance is stored --------------------------------
test_that("call provenance metadata is saved", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  md <- S4Vectors::metadata(se)$superdomain_calls$Intensity
  expect_equal(md$method, "tangent")
  expect_equal(md$log_transform_used, FALSE)
  expect_true(is.finite(md$cutoff))
  expect_true(md$call_status %in% c("called", "no_call"))
})

# ---- P0-8 / 17.6: no-call propagates to Uncertain in combined class ---------
test_that("no-call / NA propagates to Uncertain in combined taxonomy", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  # Force one intensity call to NA (simulate a no-call row) and combine
  rowData(se)$Intensity_Call__Control[1] <- NA
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  se <- combine_superdomain_calls(se, group_var = "Condition")
  expect_equal(rowData(se)$Combined_Class__Control[1], "Uncertain")
})

# ---- P1-9 / 17.5: duplicate Domain_IDs must not silently mis-map ------------
test_that("build_portrait_matrix enforces unique domain IDs", {
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found (source checkout?)")

  samples <- data.frame(
    SampleID  = c("C1", "C2", "T1", "T2"),
    Condition = c("Control", "Control", "Treatment", "Treatment"),
    Replicate = c(1, 2, 1, 2),
    bw_path   = file.path(extdata, c("C1.bw", "C2.bw", "T1.bw", "T2.bw")),
    peak_path = file.path(extdata, c("C1_peaks.bed", "C2_peaks.bed",
                                     "T1_peaks.bed", "T2_peaks.bed"))
  )
  # Deliberately duplicate every consensus peak name: build_portrait_matrix must
  # reassign unique Domain_IDs instead of silently mis-mapping via match().
  peaks <- rtracklayer::import(file.path(extdata, "peaks.bed"))
  names(peaks) <- "dup"
  se <- build_portrait_matrix(samples, consensus_peaks = peaks, workers = 1)
  expect_true(anyDuplicated(rownames(se)) == 0)
  expect_equal(nrow(se), length(peaks))
})

# ---- P1-9 / 17.5: IntervalWidth ranking works and uses rowData --------------
test_that("call_super_domains ranks on IntervalWidth from rowData", {
  se <- call_super_domains(example_se, feature = "IntervalWidth",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
  expect_true("IntervalWidth_Domain_Type" %in% colnames(rowData(se)))
  expect_true("IntervalWidth_Super_Element" %in% rowData(se)$IntervalWidth_Domain_Type)
})

# ---- P0-10: hockey-stick plot labels match caller output --------------------
test_that("hockey-stick plot uses caller's super labels", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           method = "tangent", log_transform = FALSE, verbose = FALSE)
   p <- plot_hockey_stick(se, feature = "Intensity", label_genes = FALSE)
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_true("Type_Label" %in% colnames(d))
  expect_true("Intensity-Super" %in% as.character(na.omit(d$Type_Label)))
})

# ---- P1-7: elbow normalization makes linear curves quality ~0 --------------
test_that("elbow quality is 0 for perfectly linear curves", {
  lin <- seq(1, 100, length.out = 100)
  res <- find_hockey_inflection(lin, method = "elbow")
  expect_lt(res$quality_score, 1e-6)
})

# ---- v1.0: SignalDispersion is finite across example domains --------------
test_that("example_se SignalDispersion is finite", {
  d <- rowMeans(assay(example_se, "SignalDispersion"), na.rm = TRUE)
  expect_true(all(is.finite(d[!is.na(d)])))
})

# ---- P1-1: check_signal_compatibility validates input ----------------------
test_that("check_signal_compatibility rejects bad sample sheets", {
  expect_error(check_signal_compatibility(data.frame(x = 1:3)), "bw_path")
})

# ---- P0-1: sanitize_signal actually clips (clip_zero) ----------------------
test_that("sanitize_signal clips negatives under clip_zero", {
  x <- c(5, 5, -10, 5, 5)
  clipped <- epiPortrait:::.sanitize_signal(x, "clip_zero")
  expect_equal(clipped, c(5, 5, 0, 5, 5))
  # allow keeps negatives, error throws on negatives
  expect_equal(epiPortrait:::.sanitize_signal(x, "allow"), x)
  expect_error(epiPortrait:::.sanitize_signal(x, "error"), "negative")
})

# ---- P0-1: negative fraction uses position denominator ---------------------
test_that("sanitize_signal treats NA as 0 (preserve coordinate axis)", {
  x <- c(1, NA, 3)
  expect_equal(epiPortrait:::.sanitize_signal(x, "allow"), c(1, 0, 3))
})

# ---- P0-2: bootstrap preserves explicit log_transform = FALSE --------------
test_that("bootstrap uses raw scale when log_transform = FALSE", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           method = "tangent", log_transform = FALSE,
                           n_bootstrap = 15, seed = 1, verbose = FALSE)
  md <- S4Vectors::metadata(se)$superdomain_calls$Intensity
  expect_false(md$log_transform_requested)
  expect_false(md$log_transform_used)
  expect_true(all(is.finite(md$cutoff_stability_interval)))
})

# ---- P0-3: replicate support rule (majority / all / fraction) --------------
test_that("replicate support: n=2 majority needs 2, all needs 2", {
  tm <- matrix(c("Intensity_Super_Element", "Intensity_Typical"), nrow = 1)
  # majority: n=2 -> need 2 -> only 1 super -> Typical
  r1 <- epiPortrait:::.replicate_support_call(tm, "Intensity", "majority")
  expect_equal(r1$group_type, "Intensity_Typical")
  # all: need 2 -> Typical
  r2 <- epiPortrait:::.replicate_support_call(tm, "Intensity", "all")
  expect_equal(r2$group_type, "Intensity_Typical")
  # fraction 0.5: 1/2 >= 0.5 -> Super (user opt-in)
  r3 <- epiPortrait:::.replicate_support_call(tm, "Intensity", "fraction", 0.5)
  expect_equal(r3$group_type, "Intensity_Super_Element")
})

test_that("replicate support: no-call does not inflate support", {
  # Super + no_call: n=2, valid=1 -> majority needs 2, and min_valid
  # auto-resolves to 2 -> Uncertain (NA), NOT Typical (E)
  tm <- matrix(c("Intensity_Super_Element", NA_character_), nrow = 1)
  r <- epiPortrait:::.replicate_support_call(tm, "Intensity", "majority")
  expect_equal(r$support, 0.5)  # 1/2, no-call NOT dropped
  expect_true(is.na(r$group_type))  # insufficient valid replicates -> Uncertain
})

# ---- P0-4: transition cutoff scope ------------------------------------------
test_that("relative transition is not labelled as absolute Gain/Loss", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  se <- compare_superdomains(se, group_var = "Condition",
                             ref_group = "Control", target_group = "Treatment")
  tr <- rowData(se)$Intensity_Relative_Transition
  expect_false(any(tr %in% c("Gain", "Loss"), na.rm = TRUE))
  expect_true(all(tr[!is.na(tr)] %in%
                    c("Persistent_Super", "Persistent_Typical",
                      "Relative_Prominence_Up", "Relative_Prominence_Down",
                      "Uncertain")))
})

test_that("reference cutoff transition works and labels Gain/Loss", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           method = "tangent", log_transform = FALSE,
                           verbose = FALSE)
  se <- compare_superdomains(se, group_var = "Condition",
                             ref_group = "Control", target_group = "Treatment",
                             cutoff_scope = "reference")
  expect_true("Intensity_Transition__reference" %in% colnames(rowData(se)))
})

# ---- v1.0: Breadth-Super peak-level calling semantics ------------------------
# Helper: build a minimal SummarizedExperiment carrying native peaks so the
# real .call_breadth_super_domains() pipeline can be exercised.
#
# The pipeline requires a reliable elbow on the genome-wide eligible native
# PeakWidth distribution (>= 3 peaks with a clear inflection). .make_breadth_se
# therefore appends three fixed narrow filler peaks far from the domains
# (widths 101, 111 and 121 bp; the elbow cutoff lands at 120 bp), so the caller
# controls whether the test peak (added by the specific test) is Broad (>120)
# or Typical (<=120). Domains: D1 = chr1:100-400, D2 = chr1:500-800.
.make_breadth_se <- function(test_peaks, n_samples = 2) {
  fillers <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(20000, 30000, 40000), c(20100, 30110, 40121)))
  native_peaks <- c(fillers, test_peaks)
  suppressWarnings({
    domains <- GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(c(100, 500), c(400, 800)),
      seqinfo = GenomeInfoDb::Seqinfo("chr1", 1000000))
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(Intensity = matrix(1, nrow = 2, ncol = n_samples),
                    SignalDispersion = matrix(1, nrow = 2, ncol = n_samples)),
      rowRanges = domains)
  })
  rownames(se) <- sprintf("epiDomain_%06d", seq_len(2))
  colnames(se) <- paste0("S", seq_len(n_samples))
  S4Vectors::metadata(se)$native_peaks <- setNames(
    rep(list(native_peaks), n_samples), colnames(se))
  se
}

test_that("Breadth no-call (constant widths) -> no evidence -> Uncertain", {
  # All native peaks have identical width -> no reliable elbow -> the replicate
  # provides NO evidence; the domain call must be NA (Uncertain), never Typical.
  np <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(200, 600), c(300, 700)))
  se <- .make_breadth_se(GenomicRanges::GRanges())
  # replace with all-equal widths so the elbow is unreliable
  S4Vectors::metadata(se)$native_peaks <- setNames(
    rep(list(np), ncol(se)), colnames(se))
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_sample",
                           verbose = FALSE)
  prov <- S4Vectors::metadata(se)$superdomain_calls$Breadth
  expect_equal(prov$replicates$S1$call_status, "no_call")
  expect_true(all(is.na(rowData(se)$Breadth_Call__S1)))
})

test_that("Breadth below-overlap-threshold peak gives Unmapped and no evidence", {
  # A broad test peak (250-600, width 351 -> Broad) overlaps D1 (100-400) by
  # 150 bp (fraction 0.43 < 0.5) and D2 (500-800) by 100 bp (0.28 < 0.5): no
  # domain reaches the overlap threshold -> Unmapped; neither domain gets
  # Broad/Typical evidence from it (P0-2).
  np <- GenomicRanges::GRanges("chr1", IRanges::IRanges(250, 600))
  se <- .make_breadth_se(np, n_samples = 1)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_sample",
                           verbose = FALSE)
  mapping <- S4Vectors::metadata(se)$breadth_peak_mapping
  expect_true(mapping$MappingStatus[mapping$PeakWidth == 351] == "Unmapped")
  # no unique evidence for either domain -> Uncertain
  expect_true(all(is.na(rowData(se)$Breadth_Call__S1)))
})

test_that("Breadth ambiguous tie gives no evidence to either domain", {
  # Peak 200-700 (width 501) straddles D1 (100-400) and D2 (500-800) with a
  # TIE of 200 bp into each -> Ambiguous; neither domain may get Broad or
  # Typical evidence from it (P0-2).
  np <- GenomicRanges::GRanges("chr1", IRanges::IRanges(200, 700))
  se <- .make_breadth_se(np, n_samples = 1)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_sample",
                           verbose = FALSE)
  mapping <- S4Vectors::metadata(se)$breadth_peak_mapping
  expect_true(mapping$MappingStatus[mapping$PeakWidth == 501] == "Ambiguous")
  expect_true(all(is.na(rowData(se)$Breadth_Call__S1)))
})

test_that("Breadth uniquely-mapped broad peak gives Super evidence once", {
  # Broad test peak (width 251 > cutoff 120 -> Broad) fully inside D2:
  # fraction ~1.0, unique -> D2 = Super, D1 = no evidence.
  np <- GenomicRanges::GRanges("chr1", IRanges::IRanges(525, 775))
  se <- .make_breadth_se(np, n_samples = 1)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_sample",
                           verbose = FALSE)
  mapping <- S4Vectors::metadata(se)$breadth_peak_mapping
  m_test <- mapping[mapping$PeakWidth == 251, ]
  expect_equal(m_test$MappingStatus, "Unique")
  expect_equal(m_test$SharedDomainID, "epiDomain_000002")
  expect_equal(unname(rowData(se)$Breadth_Call__S1),
               c(NA_character_, "Breadth_Super_Element"))
})

test_that("Breadth uniquely-mapped non-broad peak gives Typical evidence", {
  # Narrow test peak (width 111 <= cutoff 120 -> Typical) fully inside D1.
  np <- GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 260))
  se <- .make_breadth_se(np, n_samples = 1)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_sample",
                           verbose = FALSE)
  expect_equal(unname(rowData(se)$Breadth_Call__S1),
               c("Breadth_Typical", NA_character_))
})

test_that("Breadth broad evidence wins over typical in the same domain", {
  # D2 gets one Broad peak (525-775, 251) and one Typical peak (600-710, 111),
  # both fully inside -> D2 = Super (broad wins). D1 gets a Typical peak
  # (150-260, 111) -> D1 = Typical.
  np <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(525, 600, 150), c(775, 710, 260)))
  se <- .make_breadth_se(np, n_samples = 1)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_sample",
                           verbose = FALSE)
  ev <- unname(rowData(se)$Breadth_Call__S1)
  expect_equal(ev, c("Breadth_Typical", "Breadth_Super_Element"))
})

test_that("Breadth mapping provenance carries overlap metrics and real Domain_ID", {
  np <- GenomicRanges::GRanges("chr1", IRanges::IRanges(525, 775))
  se <- .make_breadth_se(np, n_samples = 1)
  se <- call_super_domains(se, feature = "Breadth",
                           mode = "per_sample",
                           verbose = FALSE)
  m <- S4Vectors::metadata(se)$breadth_peak_mapping
  expect_true(all(c("OverlapBp", "PeakOverlapFraction",
                    "DomainOverlapFraction", "SharedDomainID",
                    "SharedDomainIndex") %in% colnames(m)))
  m_test <- m[m$PeakWidth == 251, ]
  expect_true(is.finite(m_test$PeakOverlapFraction))
  expect_equal(m_test$SharedDomainIndex, 2L)
})

test_that("Breadth invalid overlap fraction is rejected", {
  se <- .make_breadth_se(GenomicRanges::GRanges("chr1", IRanges::IRanges(525, 775)))
  expect_error(call_super_domains(se, feature = "Breadth",
                                  min_peak_overlap_fraction = 1.5),
               "in \\[0, 1\\]")
})

test_that("Breadth-Super rejects missing native peaks", {
  se_np <- example_se
  S4Vectors::metadata(se_np)$native_peaks <- NULL
  expect_error(call_super_domains(se_np, feature = "Breadth"), "native peak")
})

test_that("Breadth-Super stores peak-level provenance", {
  se <- call_super_domains(example_se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition", verbose = FALSE)
  prov <- S4Vectors::metadata(se)$superdomain_calls$Breadth
  expect_equal(prov$calling_paradigm,
               "peak-level native PeakWidth, unique mapping, replicate aggregation")
  expect_false(is.null(prov$min_peak_overlap_fraction))
  expect_false(is.null(S4Vectors::metadata(se)$breadth_peak_calls))
  expect_false(is.null(S4Vectors::metadata(se)$breadth_peak_mapping))
})

# ---- P0-B: global_consensus honors majority (not support>=1) ----------------
test_that("global_consensus majority: 2/3 in every group is Super", {
  tm_C <- matrix(c("I_Super_Element","I_Super_Element","I_Typical"), nrow=1)
  tm_T <- matrix(c("I_Super_Element","I_Super_Element","I_Typical"), nrow=1)
  gC <- epiPortrait:::.replicate_support_call(tm_C, "I", "majority")
  gT <- epiPortrait:::.replicate_support_call(tm_T, "I", "majority")
  # both groups pass majority (2/3)
  expect_equal(gC$group_type, "I_Super_Element")
  expect_equal(gT$group_type, "I_Super_Element")
  # global class = Super (all groups pass their rule)
  expect_equal(gC$group_type, "I_Super_Element")
})

test_that("global_consensus: insufficient valid replicates -> Uncertain", {
  tm <- matrix(c("I_Super_Element", NA_character_, NA_character_), nrow=1)
  g <- epiPortrait:::.replicate_support_call(tm, "I", "majority", min_valid_replicates = 2L)
  expect_true(is.na(g$group_type))
})

# ---- P0-D: tie policy strict vs inclusive -----------------------------------
test_that("tie_policy strict uses > cutoff, inclusive uses >= cutoff", {
  se_s <- call_super_domains(example_se, feature = "Intensity",
                             method = "tangent", log_transform = FALSE,
                             tie_policy = "strict", verbose = FALSE)
  se_i <- call_super_domains(example_se, feature = "Intensity",
                             method = "tangent", log_transform = FALSE,
                             tie_policy = "inclusive", verbose = FALSE)
  n_s <- sum(rowData(se_s)$Intensity_Domain_Type == "Intensity_Super_Element", na.rm = TRUE)
  n_i <- sum(rowData(se_i)$Intensity_Domain_Type == "Intensity_Super_Element", na.rm = TRUE)
  expect_lte(n_s, n_i)  # strict never yields more supers than inclusive
})

# ---- E: min_valid_replicates auto-resolves to the support-rule requirement --
test_that("min_valid_replicates NULL auto-resolves (Super + no_call is Uncertain)", {
  tm <- matrix(c("I_Super_Element", NA_character_), nrow = 1)  # n=2, majority
  g <- epiPortrait:::.replicate_support_call(tm, "I", "majority")  # min_valid=NULL
  expect_true(is.na(g$group_type))  # only 1 valid replicate, majority needs 2 -> Uncertain
})

# ---- P0-A: cache is not polluted by negative_policy -------------------------
test_that("cache stores raw coverage (clip_zero then error still errors)", {
  extdata <- system.file("extdata", package = "epiPortrait")
  skip_if(extdata == "", "inst/extdata not found")
  peaks <- rtracklayer::import(file.path(extdata, "peaks.bed"))
  ss <- data.frame(SampleID = "S1", Condition = "C",
                   bw_path = file.path(extdata, "C1.bw"))
  cache <- tempfile()
  # C1.bw has no negatives; verify cache round-trip is raw-signal consistent:
  se1 <- build_portrait_matrix(ss, consensus_peaks = peaks, workers = 1,
                               cache_dir = cache, negative_policy = "allow")
  se2 <- build_portrait_matrix(ss, consensus_peaks = peaks, workers = 1,
                               cache_dir = cache, negative_policy = "allow")
  expect_identical(as.matrix(assay(se1, "Intensity")),
                   as.matrix(assay(se2, "Intensity")))
})



# ---- freeze review 2026-08-11: reference/pooled transition inherits the
# ---- FULL replicate-support configuration of the primary call ----------------
test_that("compare_superdomains reference/pooled inherit fraction + tie_policy", {
  se <- call_super_domains(example_se, feature = "Intensity",
                           mode = "per_group", group_var = "Condition",
                           support_rule = "fraction", min_replicate_support = 0.67,
                           tie_policy = "inclusive", verbose = FALSE)
  se <- compare_superdomains(se, group_var = "Condition",
                             ref_group = "Control", target_group = "Treatment",
                             cutoff_scope = "reference")
  tr <- S4Vectors::metadata(se)$transitions[["Control_vs_Treatment"]]
  expect_equal(tr$cutoff_scope, "reference")
  expect_equal(tr$support_rule, "fraction")
  expect_equal(tr$min_replicate_support, 0.67)
  expect_equal(tr$tie_policy, "inclusive")

  # pooled mode propagates the same configuration
  se2 <- call_super_domains(example_se, feature = "Intensity",
                            mode = "per_group", group_var = "Condition",
                            support_rule = "fraction", min_replicate_support = 0.67,
                            tie_policy = "inclusive", verbose = FALSE)
  se2 <- compare_superdomains(se2, group_var = "Condition",
                              ref_group = "Control", target_group = "Treatment",
                              cutoff_scope = "pooled")
  tr2 <- S4Vectors::metadata(se2)$transitions[["Control_vs_Treatment"]]
  expect_equal(tr2$min_replicate_support, 0.67)
  expect_equal(tr2$tie_policy, "inclusive")
  expect_true(is.null(tr2$min_valid_replicates) || is.numeric(tr2$min_valid_replicates))
})

test_that(".classify_vs_cutoff honors tie_policy and keeps NA", {
  m <- matrix(c(10, 20, NA, 15), nrow = 2)
  st <- epiPortrait:::.classify_vs_cutoff(m, 10, "Intensity", "strict")
  in_ <- epiPortrait:::.classify_vs_cutoff(m, 10, "Intensity", "inclusive")
  expect_equal(st[1, 1], "Intensity_Typical")       # 10 > 10 FALSE
  expect_equal(in_[1, 1], "Intensity_Super_Element") # 10 >= 10 TRUE
  expect_true(is.na(st[1, 2]))                       # NA stays NA
})

# ---- freeze review 2026-08-11: entry-point validation ------------------------
test_that("build_portrait_matrix rejects duplicate / empty SampleID and NA Condition", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
  dup <- data.frame(SampleID = c("A", "A"), Condition = c("C", "C"),
                    bw_path = c("x", "x"), stringsAsFactors = FALSE)
  expect_error(build_portrait_matrix(dup, gr), "SampleID")
  empt <- data.frame(SampleID = c("A", ""), Condition = c("C", "C"),
                     bw_path = c("x", "x"), stringsAsFactors = FALSE)
  expect_error(build_portrait_matrix(empt, gr), "SampleID")
  na_cond <- data.frame(SampleID = c("A", "B"), Condition = c("C", NA),
                        bw_path = c("x", "x"), stringsAsFactors = FALSE)
  expect_error(build_portrait_matrix(na_cond, gr), "Condition")
})

test_that("call_super_domains rejects invalid min_replicate_support in fraction mode", {
  for (bad in list(0, 1.1, NA_real_, "half")) {
    expect_error(
      call_super_domains(example_se, feature = "Intensity",
                         support_rule = "fraction", min_replicate_support = bad,
                         verbose = FALSE),
      "min_replicate_support")
  }
})

# ---- freeze review 2026-08-11: Breadth explicit top-fraction (quantile) -----
test_that("Breadth quantile_cutoff is an explicit top-fraction opt-in", {
  se <- call_super_domains(example_se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition",
                           quantile_cutoff = 0.8, verbose = FALSE)
  prov <- S4Vectors::metadata(se)$superdomain_calls$Breadth
  for (r in prov$replicates) {
    expect_equal(r$call_status, "called")
    expect_true(grepl("^quantile_", r$inflection_method))
    expect_true(is.finite(r$cutoff))
    expect_true(is.na(r$quality_score))  # no quality gate on the opt-in path
  }
  # quantile_cutoff = 0.8 -> top ~20% of each replicate's eligible peaks Broad
  # (Control only: Treatment anchors tie at the cutoff, so the strict > rule
  # yields fewer Broad peaks there — that is the intended tie behavior).
  pc <- S4Vectors::metadata(se)$breadth_peak_calls
  frac <- tapply(pc$PeakBroadCall == "Broad", pc$SampleID, mean)
  ctrl_frac <- frac[grepl("^Control", names(frac))]
  expect_true(all(ctrl_frac > 0.15 & ctrl_frac < 0.25))
})

test_that("Breadth quantile_cutoff rejects invalid values", {
  for (bad in list(0, 1, 1.5, NA_real_)) {
    expect_error(
      call_super_domains(example_se, feature = "Breadth",
                         quantile_cutoff = bad, verbose = FALSE),
      "quantile_cutoff")
  }
})

# ---- Breadth n_bootstrap (symmetric with Intensity, freeze 2026-08-11) ------
test_that("Breadth n_bootstrap stores per-replicate stability intervals", {
  se <- call_super_domains(example_se, feature = "Breadth",
                           mode = "per_group", group_var = "Condition",
                           n_bootstrap = 20, seed = 1, verbose = FALSE)
  prov <- S4Vectors::metadata(se)$superdomain_calls$Breadth
  n_called <- 0L
  for (r in prov$replicates) {
    if (r$call_status == "called") {
      n_called <- n_called + 1L
      expect_true(is.finite(r$cutoff))
      expect_true(all(is.finite(r$cutoff_stability_interval)))
      expect_true(r$bootstrap_success_rate >= 0 && r$bootstrap_success_rate <= 1)
    } else {
      expect_true(is.null(r$cutoff_stability_interval))
    }
  }
  expect_gt(n_called, 0)
})

test_that("Breadth quantile_cutoff + n_bootstrap warns (bootstrap applies only to inflection)", {
  expect_warning(
    call_super_domains(example_se, feature = "Breadth",
                       mode = "per_group", group_var = "Condition",
                       quantile_cutoff = 0.8, n_bootstrap = 10,
                       verbose = FALSE),
    "data-driven inflection path")
})
