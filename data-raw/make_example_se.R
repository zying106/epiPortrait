# Recreate data/example_se.rda with the v1.0 assay structure.
# 500 domains x 6 samples (3 Control, 3 Treatment).
# Assays: Intensity, SignalDispersion, NativeMaxPeakWidth, NativeOccupiedWidth,
# NativePeakCount. Native peaks are stored in metadata(se)$native_peaks so that
# Breadth-Super calling works out of the box.
#
# Design (2026-08-11, embedded-data review):
#   * All domains lie inside chr21 (hg38, 46,709,983 bp) — no out-of-bounds
#     ranges; the genome is declared on the seqinfo.
#   * colData carries SampleID / Condition / BioReplicate / bw_path / peak_path
#     (bw_path / peak_path are NA: this is a precomputed object with no backing
#     files).
#   * Synthetic data deliberately mixes replicate-support patterns so
#     replicate-aware calling is instructive:
#       - Intensity-high domains: ~70% strong (8/8/8 -> 3/3 support),
#         ~20% boundary (8/8/0.5 -> 2/3), ~10% weak (8/0.5/0.5 -> 1/3).
#       - Breadth domains: per-replicate broad peaks in Control (same
#         strong/boundary/weak mix), contracted in Treatment; 6 no-call
#         domains have no native peak in 2 Control replicates (NA evidence ->
#         Uncertain via min_valid_replicates).
#   * rowData$Simulation_Scenario records the generative mechanism, not a
#     static "true class": Intensity-high / Breadth-contraction /
#     Dual-with-breadth-contraction / Typical.
#   * metadata$simulation_provenance records seed / genome / generator.
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(GenomicRanges))

set.seed(7)
n_dom <- 500
n_ctrl <- 3
n_treat <- 3

sl <- Seqinfo(seqnames = "chr21", seqlengths = 46709983, isCircular = FALSE,
              genome = "hg38")

# Ground-truth generative mechanisms: ~8% Intensity-high, ~8% Breadth-only
# (broad in Control only), ~4% Dual (intensity-high + breadth-contraction),
# rest Typical.
int_super <- sample(n_dom, round(0.08 * n_dom))
brd_super <- sample(setdiff(seq_len(n_dom), int_super), round(0.08 * n_dom))
dual_super <- sample(setdiff(seq_len(n_dom), c(int_super, brd_super)),
                     round(0.04 * n_dom))
scenario <- rep("Typical", n_dom)
scenario[int_super] <- "Intensity-high"
scenario[brd_super] <- "Breadth-contraction"
scenario[dual_super] <- "Dual-with-breadth-contraction"

# Shared candidate domain intervals: Breadth/Dual-Super domains are WIDE
# macro-domain intervals (typical of stitched super-enhancer regions) so that
# native broad peaks fit inside them (peak-overlap fraction >= 0.5). Typical /
# Intensity-Super domains are ~1 kb. 500 domains spaced every 80 kb from 1 Mb
# keep every interval inside chr21 (max end ~40.9 Mb < 46.7 Mb).
widths <- rep(1000L, n_dom)
widths[c(brd_super, dual_super)] <- 20000L
gr <- GRanges("chr21",
              IRanges(start = seq(1e6, by = 80000, length.out = n_dom),
                      width = widths),
              seqinfo = sl)
stopifnot(all(end(gr) <= seqlengths(gr)[as.character(seqnames(gr))]))

# ---- per-replicate Intensity support patterns --------------------------------
int_high <- unique(c(int_super, dual_super))  # disjoint: int_super + dual_super
set.seed(77)
int_lvl <- sample(c("strong", "boundary", "weak"), length(int_high),
                  replace = TRUE, prob = c(0.7, 0.2, 0.1))
int_rep <- lapply(int_lvl, function(p) switch(p,
  strong   = c(50, 50, 50),
  boundary = c(50, 50, 0.3),
  weak     = c(50, 0.3, 0.3)))
names(int_rep) <- int_high

int_mat <- matrix(rlnorm(n_dom * (n_ctrl + n_treat), meanlog = 3, sdlog = 0.15),
                  nrow = n_dom)
for (j in seq_len(n_ctrl + n_treat)) {
  rep_j <- if (j <= n_ctrl) j else j - n_ctrl
  mult <- rep(1, n_dom)
  for (d in int_high) mult[d] <- int_rep[[as.character(d)]][rep_j]
  int_mat[, j] <- int_mat[, j] * mult
}

# ---- per-replicate Breadth support patterns (Control) ------------------------
brd_set <- unique(c(brd_super, dual_super))
set.seed(123)
brd_lvl <- sample(c("strong", "boundary", "weak"), length(brd_set),
                  replace = TRUE, prob = c(0.7, 0.2, 0.1))
brd_rep <- lapply(brd_lvl, function(p) switch(p,
  strong   = c(TRUE, TRUE, TRUE),
  boundary = c(TRUE, TRUE, FALSE),
  weak     = c(TRUE, FALSE, FALSE)))
names(brd_rep) <- brd_set

# ---- no-call domains: no native peak in 2 Control replicates (NA evidence) ---
set.seed(321)
no_call_doms <- sample(setdiff(seq_len(n_dom), brd_set), 6)
no_call_reps <- lapply(no_call_doms, function(d) sort(sample(1:3, 2)))
names(no_call_reps) <- no_call_doms

# Native peaks per sample. Breadth domains are broad (10-30 kb) in Control and
# contract (2-4 kb) in Treatment (breadth contraction); the per-replicate broad
# pattern decides which Control replicates carry the wide peak.
native_peaks <- vector("list", n_ctrl + n_treat)
nwmax_mat <- matrix(rep(NA_real_, n_dom * (n_ctrl + n_treat)), nrow = n_dom)
nwocc_mat <- matrix(rep(NA_real_, n_dom * (n_ctrl + n_treat)), nrow = n_dom)
nwcnt_mat <- matrix(rep(0, n_dom * (n_ctrl + n_treat)), nrow = n_dom)

make_native <- function(dom_idx, w_scale, amp) {
  # Exactly ONE peak per domain, centered, sized so that the peak-overlap
  # fraction against the domain is >= 0.6 (>= min_peak_overlap_fraction). This
  # keeps the unique peak->domain mapping clean (no Ambiguous ties, no
  # Unmapped): a typical 1 kb domain gets a ~600 bp peak, a broad 20 kb domain
  # gets a ~12 kb peak in Control. Width still varies per sample so the
  # replicate-aware support patterns are meaningful.
  n <- length(dom_idx)
  d <- gr[dom_idx]
  w0 <- pmax(200, rnorm(n, 2000, 400))
  w <- pmax(200, w0 * w_scale)
  w <- pmin(w, round(0.6 * width(d)))
  start1 <- start(d) + (width(d) - w) / 2
  start1 <- as.integer(pmax(start1, start(d)))
  GRanges(seqnames(d), IRanges(start = start1, width = as.integer(w)),
          score = amp)
}

# Fixed-width single peak per domain (used for the contracted Treatment
# replicates so the width distribution has no noisy tail).
make_native_fixed <- function(dom_idx, w, amp) {
  d <- gr[dom_idx]
  start1 <- start(d) + (width(d) - w) / 2
  start1 <- as.integer(pmax(start1, start(d)))
  GRanges(seqnames(d), IRanges(start = start1, width = as.integer(w)),
          score = amp)
}

for (j in seq_len(n_ctrl + n_treat)) {
  is_treat <- j > n_ctrl
  base_doms <- setdiff(seq_len(n_dom), c(brd_set, no_call_doms))
  if (is_treat) {
    # Treatment = breadth contraction: every domain gets a fixed ~400 bp peak
    # (brd domains contracted to typical width). A few residual broad peaks are
    # placed in the inter-domain GAPS (no domain overlap -> no evidence) purely
    # to anchor the per-sample width distribution, so the width inflection is
    # valid and the contracted brd domains are called Typical (not no-call).
    pk <- make_native_fixed(base_doms, 400L, 5)
    pk <- c(pk, make_native_fixed(no_call_doms, 400L, 5))
    pk <- c(pk, make_native_fixed(brd_set, 400L, 5))
    gap_centers <- as.integer((end(gr)[-length(gr)] + start(gr)[-1]) / 2)
    anchors <- GRanges("chr21",
                       IRanges(start = gap_centers - 6000L, width = 12000L))
    set.seed(999 + j)
    pk <- c(pk, anchors[sample(length(anchors), 30)])
  } else {
    # Control baseline (typical) peaks: ~400 bp, clearly below the width
    # inflection. Breadth domains are EXCLUDED here (their evidence is generated
    # explicitly below: a wide ~12 kb peak when the replicate pattern is broad,
    # otherwise a ~400 bp typical peak).
    pk <- make_native(base_doms, w_scale = 0.2, amp = 5)
    # no-call domains: present only in the Control replicates not marked missing.
    for (d in no_call_doms) {
      if (!j %in% no_call_reps[[as.character(d)]]) {
        pk <- c(pk, make_native(d, w_scale = 0.2, amp = 5))
      }
    }
    # Breadth domains: wide (x10 -> ~12 kb, Broad) or typical (~400 bp) per the
    # replicate pattern (strong 3/3 / boundary 2/3 / weak 1/3).
    for (d in brd_set) {
      wide <- brd_rep[[as.character(d)]][j]
      pk <- c(pk, make_native(d, w_scale = if (wide) 10 else 0.2, amp = 5))
    }
  }
  pk <- sort(pk)
  native_peaks[[j]] <- pk
  hits <- findOverlaps(gr, pk)
  qh <- queryHits(hits)
  nwcnt_mat[, j] <- tabulate(qh, nbins = n_dom)
  for (d in unique(qh)) {
    ww <- width(pk[subjectHits(hits)[qh == d]])
    nwmax_mat[d, j] <- max(ww)
    red <- reduce(intersect(pk[subjectHits(hits)[qh == d]], gr[d],
                            ignore.strand = TRUE))
    nwocc_mat[d, j] <- sum(width(red))
  }
}

names(native_peaks) <- c(paste0("Control_", 1:n_ctrl), paste0("Treatment_", 1:n_treat))

# SignalDispersion: proportional to native breadth (wider -> larger dispersion)
# plus noise; typical narrow domains ~150-250 bp dispersion.
disp_mat <- matrix(NA_real_, n_dom, n_ctrl + n_treat)
for (j in seq_len(n_ctrl + n_treat)) {
  disp_mat[, j] <- nwocc_mat[, j] / 20 + rnorm(n_dom, 100, 30)
  disp_mat[, j] <- pmax(disp_mat[, j], 20)
}

colnames(int_mat) <- names(native_peaks)
colnames(disp_mat) <- names(native_peaks)
colnames(nwmax_mat) <- names(native_peaks)
colnames(nwocc_mat) <- names(native_peaks)
colnames(nwcnt_mat) <- names(native_peaks)

col_data <- data.frame(
  SampleID = names(native_peaks),
  Condition = rep(c("Control", "Treatment"), c(n_ctrl, n_treat)),
  BioReplicate = rep(1:3, 2),
  # NA: this is a precomputed synthetic object with no backing BigWig/peak
  # files (an empty string would look like a forgotten path).
  bw_path = NA_character_,
  peak_path = NA_character_,
  row.names = names(native_peaks)
)

example_se <- SummarizedExperiment(
  assays = list(Intensity = int_mat,
                SignalDispersion = disp_mat,
                NativeMaxPeakWidth = nwmax_mat,
                NativeOccupiedWidth = nwocc_mat,
                NativePeakCount = nwcnt_mat),
  rowRanges = gr,
  colData = col_data
)
rownames(example_se) <- sprintf("epiDomain_%06d", seq_len(n_dom))
rowData(example_se)$Simulation_Scenario <- scenario
rowData(example_se)$IntervalWidth <- width(gr)
S4Vectors::metadata(example_se)$native_peaks <- native_peaks
S4Vectors::metadata(example_se)$native_peak_available <- rep(TRUE, ncol(example_se))
S4Vectors::metadata(example_se)$feature_definitions <- list(
  Intensity = list(definition = "Integrated BigWig signal within candidate domain"),
  SignalDispersion = list(definition = "Signal-weighted genomic SD within domain"),
  NativeMaxPeakWidth = list(definition = "Max width of native peaks mapped to the domain"),
  NativeOccupiedWidth = list(definition = "Sum of reduced native peak widths inside the domain"),
  NativePeakCount = list(definition = "Number of native peaks overlapping the domain"))
S4Vectors::metadata(example_se)$simulation_provenance <- list(
  dataset = "synthetic",
  seed = 7L,
  genome = "hg38",
  chromosome = "chr21",
  n_domains = n_dom,
  n_control = n_ctrl,
  n_treatment = n_treat,
  generator = "data-raw/make_example_se.R",
  intensity_support_patterns =
    "Intensity-high domains: ~70% strong (3/3), ~20% boundary (2/3), ~10% weak (1/3)",
  breadth_support_patterns =
    "Breadth domains: per-replicate broad peaks in Control (strong/boundary/weak), contracted in Treatment; 6 no-call domains with no native peak in 2 Control replicates")

save(example_se, file = "data/example_se.rda", compress = "xz")
message("example_se written: ", n_dom, " x ", ncol(example_se),
        ", assays: ", paste(assayNames(example_se), collapse = ", "))
