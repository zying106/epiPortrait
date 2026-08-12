# Recreate inst/extdata tiny synthetic dataset (C1/C2/T1/T2 .bw + consensus
# peaks.bed + per-sample native peak .bed files).
# Run once; outputs are committed to inst/extdata.
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))

set.seed(42)
sl <- Seqinfo(seqnames = "chr21", seqlengths = 46709983, isCircular = FALSE)

# Signal is spread across the FULL peak width (60 sample points, per-base bw).
make_bw <- function(path, peaks, amp = 10) {
  sig <- GRanges()
  for (i in seq_len(length(peaks))) {
    w <- width(peaks[i])
    x <- seq(0, w, length.out = 60)
    bump <- amp * exp(-((x - w / 2) / (w / 8))^2) + runif(60, 0, 0.5)
    pos <- round(seq(start(peaks[i]), end(peaks[i]), length.out = 60))
    sig <- c(sig, GRanges(seqnames(peaks)[i],
                          IRanges(pos, width = 1),
                          score = bump))
  }
  seqlengths(sig) <- seqlengths(sl)["chr21"]
  export.bw(sig, path)
}

write_bed <- function(gr, path) {
  # BED is 0-based half-open; GRanges is 1-based closed.
  write.table(
    data.frame(chr = as.character(seqnames(gr)),
               start = start(gr) - 1L,
               end = end(gr)),
    file = path, row.names = FALSE, col.names = FALSE,
    sep = "\t", quote = FALSE
  )
}

d <- "inst/extdata"
# Control (C1/C2): 2 narrow (~1-2 kb) + 1 broad (~10 kb) native peaks.
peaksA <- GRanges("chr21",
                  IRanges(c(44000000, 44020000, 44030000),
                          c(44001000, 44022000, 44040000)))
# Treatment (T1/T2): same loci but breadth contraction (broad -> ~4 kb);
# higher amplitude so Intensity stays comparable.
peaksB <- GRanges("chr21",
                  IRanges(c(44000500, 44020500, 44033000),
                          c(44001200, 44022200, 44037000)))

make_bw(file.path(d, "C1.bw"), peaksA, amp = 10)
make_bw(file.path(d, "C2.bw"), peaksA, amp = 10)
make_bw(file.path(d, "T1.bw"), peaksB, amp = 25)
make_bw(file.path(d, "T2.bw"), peaksB, amp = 25)

# Consensus (shared) domain universe: union of all native peaks (fixed width
# bed used by the build example).
write_bed(GRanges("chr21",
                  IRanges(c(43999000, 44019000, 44029000),
                          c(44002000, 44023000, 44044000))),
          file.path(d, "peaks.bed"))

# Per-sample native peak files (used by Breadth-Super calling). Slight
# per-sample variation so replicate-aware support is meaningful.
write_bed(peaksA, file.path(d, "C1_peaks.bed"))
write_bed(peaksA, file.path(d, "C2_peaks.bed"))
write_bed(peaksB, file.path(d, "T1_peaks.bed"))
write_bed(peaksB, file.path(d, "T2_peaks.bed"))

message("done")
