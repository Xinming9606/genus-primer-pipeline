# =============================================================================
# design_primers.R
#
# Input  (via snakemake@ object):
#   snakemake@input[[1]]          : aligned FASTA (.aln)
#
# Output (via snakemake@ object):
#   snakemake@output[["tsv"]]     : primer pairs TSV
#   snakemake@output[["plot"]]    : diversity PNG
#
# Params (via snakemake@params):
#   primer_len        : primer length in bp          (default 20)
#   amplicon_min_len  : minimum amplicon length       (default 300)
#   amplicon_max_len  : maximum amplicon length       (default 1000)
#   div_cut           : Shannon entropy cutoff        (default 0.5)
#   GC_tol            : max GC% difference in a pair (default 0.1)
#
# TSV columns:
#   primer_id, fwd, rev, fwd_pos, rev_pos, amplicon_len,
#   fwd_GC, rev_GC, pair_diversity, delta_GC, combined_score
# =============================================================================

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

library(seqinr)
library(zoo)

# ---------------------------------------------------------------------------
# 0. Parameters
# ---------------------------------------------------------------------------
aln_file         <- snakemake@input[[1]]
out_tsv          <- snakemake@output[["tsv"]]
out_plot         <- snakemake@output[["plot"]]

primer_len       <- as.integer(snakemake@params[["primer_len"]])
amplicon_min_len <- as.integer(snakemake@params[["amplicon_min_len"]])
amplicon_max_len <- as.integer(snakemake@params[["amplicon_max_len"]])
div_cut          <- as.numeric(snakemake@params[["div_cut"]])
GC_tol           <- as.numeric(snakemake@params[["GC_tol"]])

message("Parameters:")
message("  aln_file         = ", aln_file)
message("  primer_len       = ", primer_len)
message("  amplicon_min_len = ", amplicon_min_len)
message("  amplicon_max_len = ", amplicon_max_len)
message("  div_cut          = ", div_cut)
message("  GC_tol           = ", GC_tol)

# ---------------------------------------------------------------------------
# 1. Load alignment
# ---------------------------------------------------------------------------
fasta   <- read.fasta(aln_file, seqtype = "DNA", forceDNAtolower = TRUE)
DNA_mat <- do.call(rbind, fasta)

n_seqs  <- nrow(DNA_mat)
aln_len <- ncol(DNA_mat)

message("Loaded ", n_seqs, " sequences, alignment length ", aln_len, " bp")

if (n_seqs < 2) {
  stop("Need at least 2 sequences to design primers; only ", n_seqs, " found.")
}

# ---------------------------------------------------------------------------
# 2. Per-position Shannon entropy
# ---------------------------------------------------------------------------
CONSENSUS <- apply(DNA_mat, 2, function(x) names(sort(table(x), decreasing = TRUE))[1])

divs <- numeric(aln_len)
for (i in seq_len(aln_len)) {
  per_pos <- factor(DNA_mat[, i], levels = c("a", "g", "c", "t", "-"))
  probs   <- table(per_pos) / n_seqs
  divs[i] <- -sum(probs * ifelse(probs == 0, 0, log(probs)))
}

message("Mean per-position entropy: ", round(mean(divs), 4))

# ---------------------------------------------------------------------------
# 3. Diversity plot (points + rolling mean, primer sites added later)
# ---------------------------------------------------------------------------
roll_k       <- min(10, aln_len)
roll_means   <- rollmean(divs, k = roll_k, fill = NA)

png(out_plot, width = 1800, height = 900, res = 150)
par(mar = c(4, 4, 2, 1), family = "sans")
plot(
  divs, pch = 16, cex = 0.3,
  ylim = c(0, max(divs) * 1.1),
  type = "p",
  main = paste0("Sequence diversity — ", basename(aln_file)),
  xlab = "Alignment position (bp)",
  ylab = "Shannon entropy"
)
lines(roll_means, col = "#2ca25f", lwd = 1.5)
legend("topright",
  legend = c("Per-position entropy", paste0("Rolling mean (k=", roll_k, ")")),
  col    = c("black", "#2ca25f"),
  pch    = c(16, NA), lty = c(NA, 1), pt.cex = 0.6, lwd = c(NA, 1.5),
  bty    = "n"
)
# primer rectangles are added below, after we know the positions
plot_env <- environment()   # capture so we can add rects after dev is open
dev.off()

# ---------------------------------------------------------------------------
# 4. Helper: reverse complement of a sequence string
# ---------------------------------------------------------------------------
rev_comp <- function(seq_str) {
  bases <- rev(strsplit(toupper(seq_str), "")[[1]])
  comp  <- chartr("ACGTRYMKSWHBVDN", "TGCAYRKMSWDVBHN", paste(bases, collapse = ""))
  comp
}

# ---------------------------------------------------------------------------
# 5. Build kmer table (consensus sequence + per-kmer diversity + GC)
# ---------------------------------------------------------------------------
n_kmers <- aln_len - primer_len
if (n_kmers < 1) {
  stop("Alignment (", aln_len, " bp) is shorter than primer_len (", primer_len, " bp).")
}

kmers <- data.frame(
  pos  = integer(n_kmers),
  kmer = character(n_kmers),
  divs = numeric(n_kmers),
  GC   = numeric(n_kmers),
  stringsAsFactors = FALSE
)

for (j in seq_len(n_kmers)) {
  idx          <- j:(j + primer_len - 1)
  kmers$pos[j] <- j
  kmers$kmer[j] <- paste(CONSENSUS[idx], collapse = "")
  kmers$divs[j] <- sum(divs[idx])
  gc_count      <- sum(CONSENSUS[idx] %in% c("g", "c"))
  kmers$GC[j]   <- gc_count / primer_len
}

# ---------------------------------------------------------------------------
# 6. Filter candidates by diversity cutoff — no auto-increment
# ---------------------------------------------------------------------------
candidates <- kmers[kmers$divs <= div_cut, ]

if (nrow(candidates) == 0) {
  message("No primer candidates pass div_cut = ", div_cut,
          ". Writing empty TSV.")
  empty <- data.frame(
    primer_id = character(), fwd = character(), rev = character(),
    fwd_pos = integer(), rev_pos = integer(), amplicon_len = integer(),
    fwd_GC = numeric(), rev_GC = numeric(),
    pair_diversity = numeric(), delta_GC = numeric(),
    combined_score = numeric()
  )
  write.table(empty, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  sink(type = "output"); sink(type = "message"); close(log_con)
  quit(save = "no", status = 0)
}

message(nrow(candidates), " candidate kmers pass div_cut filter")

# ---------------------------------------------------------------------------
# 7. Evaluate all candidate pairs
# ---------------------------------------------------------------------------
nc      <- nrow(candidates)
results <- list()

for (i in seq_len(nc - 1)) {
  for (j in (i + 1):nc) {
    amp_len <- candidates$pos[j] - candidates$pos[i]

    if (amp_len < amplicon_min_len || amp_len > amplicon_max_len) next

    delta_gc <- abs(candidates$GC[i] - candidates$GC[j])
    if (delta_gc >= GC_tol) next

    pair_div <- candidates$divs[i] + candidates$divs[j]
    score    <- 1 / (abs(pair_div) + 10 * delta_gc^2 + 0.01)

    fwd_seq  <- toupper(candidates$kmer[i])
    rev_seq  <- rev_comp(candidates$kmer[j])

    results[[length(results) + 1]] <- data.frame(
      primer_id      = NA_character_,
      fwd            = fwd_seq,
      rev            = rev_seq,
      fwd_pos        = candidates$pos[i],
      rev_pos        = candidates$pos[j],
      amplicon_len   = amp_len,
      fwd_GC         = round(candidates$GC[i], 4),
      rev_GC         = round(candidates$GC[j], 4),
      pair_diversity = round(pair_div, 4),
      delta_GC       = round(delta_gc, 4),
      combined_score = round(score, 6),
      stringsAsFactors = FALSE
    )
  }
}

if (length(results) == 0) {
  message("No valid primer pairs found within amplicon length and GC constraints.",
          " Writing empty TSV.")
  empty <- data.frame(
    primer_id = character(), fwd = character(), rev = character(),
    fwd_pos = integer(), rev_pos = integer(), amplicon_len = integer(),
    fwd_GC = numeric(), rev_GC = numeric(),
    pair_diversity = numeric(), delta_GC = numeric(),
    combined_score = numeric()
  )
  write.table(empty, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  sink(type = "output"); sink(type = "message"); close(log_con)
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# 8. Sort, assign IDs, write TSV
# ---------------------------------------------------------------------------
out_df <- do.call(rbind, results)
out_df <- out_df[order(-out_df$combined_score), ]
out_df$primer_id <- paste0("primer_pair_", seq_len(nrow(out_df)))

write.table(out_df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", nrow(out_df), " primer pairs to ", out_tsv)

# ---------------------------------------------------------------------------
# 9. Re-draw diversity plot with primer site rectangles overlaid
# ---------------------------------------------------------------------------
# Take top 5 pairs and highlight their positions
top_n   <- min(5, nrow(out_df))
top_df  <- out_df[seq_len(top_n), ]

png(out_plot, width = 1800, height = 900, res = 150)
par(mar = c(4, 4, 2, 1), family = "sans")
plot(
  divs, pch = 16, cex = 0.3,
  ylim = c(0, max(divs) * 1.1),
  type = "p",
  main = paste0("Sequence diversity — ", basename(aln_file)),
  xlab = "Alignment position (bp)",
  ylab = "Shannon entropy"
)
lines(roll_means, col = "#2ca25f", lwd = 1.5)

# Draw primer site rectangles (semi-transparent)
rect_col <- adjustcolor("#e34a33", alpha.f = 0.25)
for (k in seq_len(top_n)) {
  # forward primer
  rect(top_df$fwd_pos[k], -0.05,
       top_df$fwd_pos[k] + primer_len, max(divs) * 1.05,
       col = rect_col, border = NA)
  # reverse primer
  rect(top_df$rev_pos[k], -0.05,
       top_df$rev_pos[k] + primer_len, max(divs) * 1.05,
       col = rect_col, border = NA)
}

legend("topright",
  legend = c("Per-position entropy",
             paste0("Rolling mean (k=", roll_k, ")"),
             paste0("Top ", top_n, " primer sites")),
  col    = c("black", "#2ca25f", rect_col),
  pch    = c(16, NA, 15), lty = c(NA, 1, NA),
  pt.cex = c(0.6, NA, 1.5), lwd = c(NA, 1.5, NA),
  bty    = "n"
)
dev.off()

message("Wrote diversity plot to ", out_plot)

# ---------------------------------------------------------------------------
# 10. Done
# ---------------------------------------------------------------------------
sink(type = "output")
sink(type = "message")
close(log_con)
