#!/usr/bin/env Rscript
# Design primer candidates from a multiple sequence alignment.
#
# Inputs (via Snakemake `snakemake` object OR commandArgs fallback):
#   - alignment FASTA (input)
#   - primer params (length range, Tm range, GC range, amplicon range, etc.)
#
# Outputs:
#   - primers.tsv: one row per primer pair, ranked by combined score
#   - diversity.png: Shannon diversity plot with primer sites highlighted

suppressPackageStartupMessages({
  for (pkg in c("seqinr", "zoo")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  library(seqinr)
  library(zoo)
})

# ---------- argument handling: snakemake or CLI ----------
if (exists("snakemake")) {
  ALN_FILE   <- snakemake@input[["alignment"]]
  PRIMER_TSV <- snakemake@output[["primers"]]
  DIV_PLOT   <- snakemake@output[["diversity_plot"]]
  cfg        <- snakemake@params[["cfg"]]
} else {
  args <- commandArgs(trailingOnly = TRUE)
  ALN_FILE   <- args[1]
  PRIMER_TSV <- args[2]
  DIV_PLOT   <- args[3]
  cfg <- list(
    length_min     = as.numeric(args[4]),
    length_max     = as.numeric(args[5]),
    tm_min         = as.numeric(args[6]),
    tm_max         = as.numeric(args[7]),
    gc_min         = as.numeric(args[8]),
    gc_max         = as.numeric(args[9]),
    amplicon_min   = as.numeric(args[10]),
    amplicon_max   = as.numeric(args[11]),
    max_degeneracy = as.numeric(args[12]),
    top_n          = as.numeric(args[13])
  )
}

LEN_MIN  <- as.integer(cfg$length_min)
LEN_MAX  <- as.integer(cfg$length_max)
TM_MIN   <- cfg$tm_min
TM_MAX   <- cfg$tm_max
GC_MIN   <- cfg$gc_min / 100
GC_MAX   <- cfg$gc_max / 100
AMP_MIN  <- cfg$amplicon_min
AMP_MAX  <- cfg$amplicon_max
MAX_DEG  <- as.integer(cfg$max_degeneracy)
TOP_N    <- as.integer(cfg$top_n)

# ---------- helpers ----------

IUPAC <- list(
  "a"="a","c"="c","g"="g","t"="t",
  "ag"="r","ct"="y","cg"="s","at"="w","gt"="k","ac"="m",
  "cgt"="b","agt"="d","act"="h","acg"="v",
  "acgt"="n"
)
iupac_code <- function(bases) {
  bases <- sort(unique(tolower(bases)))
  bases <- bases[bases %in% c("a","c","g","t")]
  if (length(bases) == 0) return("n")
  key <- paste(bases, collapse = "")
  code <- IUPAC[[key]]
  if (is.null(code)) "n" else code
}

column_degeneracy <- function(col) {
  bases <- col[col %in% c("a","c","g","t")]
  length(unique(bases))
}

primer_degeneracy_fold <- function(mat_slice) {
  folds <- apply(mat_slice, 2, function(col) {
    d <- column_degeneracy(col)
    max(d, 1)
  })
  prod(folds)
}

build_degenerate_primer <- function(mat_slice) {
  apply(mat_slice, 2, function(col) {
    bases <- col[col %in% c("a","c","g","t")]
    if (length(bases) == 0) return("n")
    iupac_code(bases)
  })
}

gc_fraction <- function(seq) {
  s <- tolower(seq)
  gc_strong <- sum(s %in% c("g","c","s"))
  half      <- sum(s %in% c("r","y","k","m","b","d","h","v","n"))
  total     <- length(s)
  if (total == 0) return(NA)
  (gc_strong + 0.5 * half) / total
}

# Tm via Wallace (short oligos) or salt-adjusted formula (longer).
# Quick estimate sufficient for ranking; for downstream design verify with
# a dedicated tool (e.g. Primer3 or Bio.SeqUtils.MeltingTemp).
tm_estimate <- function(seq) {
  s <- tolower(seq)
  n <- length(s)
  gc_count <- sum(s %in% c("g","c","s")) +
              0.5 * sum(s %in% c("r","y","k","m","b","d","h","v","n"))
  at_count <- n - gc_count
  if (n < 14) return(2 * at_count + 4 * gc_count)
  64.9 + 41 * (gc_count - 16.4) / n
}

revcomp <- function(seq) {
  comp <- c(a="t", t="a", g="c", c="g",
            r="y", y="r", s="s", w="w", k="m", m="k",
            b="v", v="b", d="h", h="d", n="n", "-"="-")
  chars <- strsplit(tolower(seq), "")[[1]]
  rc <- comp[chars]
  rc[is.na(rc)] <- "n"
  paste(rev(rc), collapse = "")
}

# ---------- read alignment ----------
fasta   <- read.fasta(ALN_FILE, forceDNAtolower = TRUE)
DNA_mat <- do.call(rbind, fasta)
max_len <- ncol(DNA_mat)
n_seqs  <- nrow(DNA_mat)

cat(sprintf("[design_primers] %d sequences, alignment length %d\n",
            n_seqs, max_len))

# ---------- per-position Shannon diversity ----------
divs <- numeric(max_len)
for (i in seq_len(max_len)) {
  per_pos <- factor(DNA_mat[, i], levels = c("a","g","c","t","-"))
  probs   <- table(per_pos) / n_seqs
  probs_nz <- probs[probs > 0]
  divs[i] <- -sum(probs_nz * log(probs_nz))
}
roll_div <- rollmean(divs, k = 10, fill = NA)

# ---------- scan primer-length k-mers across all allowed lengths ----------
cat("[design_primers] scanning primer candidates...\n")
candidates <- list()
for (L in LEN_MIN:LEN_MAX) {
  for (start in seq_len(max_len - L + 1)) {
    end <- start + L - 1
    slice <- DNA_mat[, start:end, drop = FALSE]

    if (any(slice == "-")) next

    fold <- primer_degeneracy_fold(slice)
    if (fold > MAX_DEG) next

    primer_seq <- build_degenerate_primer(slice)
    primer_str <- paste(primer_seq, collapse = "")

    gc <- gc_fraction(primer_seq)
    if (is.na(gc) || gc < GC_MIN || gc > GC_MAX) next

    tm <- tm_estimate(primer_seq)
    if (tm < TM_MIN || tm > TM_MAX) next

    div_sum <- sum(divs[start:end])

    candidates[[length(candidates) + 1]] <- data.frame(
      pos = start, end = end, length = L, seq = primer_str,
      gc = gc, tm = tm, degeneracy = fold, div_sum = div_sum,
      stringsAsFactors = FALSE
    )
  }
}

empty_tsv <- function() {
  write.table(
    data.frame(
      primer_id = character(0), fwd = character(0), rev = character(0),
      fwd_pos = integer(0), rev_pos = integer(0),
      amplicon_len = integer(0),
      fwd_tm = numeric(0), rev_tm = numeric(0), tm_diff = numeric(0),
      fwd_gc = numeric(0), rev_gc = numeric(0),
      fwd_deg = integer(0), rev_deg = integer(0),
      conservation_score = numeric(0), internal_diversity = numeric(0),
      score = numeric(0)
    ),
    PRIMER_TSV, sep = "\t", quote = FALSE, row.names = FALSE
  )
}

empty_plot <- function() {
  png(DIV_PLOT, width = 2400, height = 1200, res = 200)
  par(mar = c(4, 4, 2, 1))
  plot(divs, type = "p", pch = 16, cex = 0.3,
       xlab = "Alignment position", ylab = "Shannon diversity",
       main = "Per-position diversity")
  lines(roll_div, col = "red")
  legend("topright",
         legend = c("per-position", "rolling mean (k=10)"),
         col    = c("black", "red"),
         pch    = c(16, NA),
         lty    = c(NA, 1))
  dev.off()
}

if (length(candidates) == 0) {
  warning("No primer candidates passed single-primer filters.")
  empty_tsv()
  empty_plot()
  quit(save = "no", status = 0)
}

cands <- do.call(rbind, candidates)
cat(sprintf("[design_primers] %d candidates passed single-primer filters\n",
            nrow(cands)))

# ---------- pair candidates -> primer pairs ----------
cat("[design_primers] pairing candidates...\n")
pairs_list <- list()
for (i in seq_len(nrow(cands))) {
  fwd <- cands[i, ]
  rev_candidates <- cands[
    cands$pos > fwd$end &
    (cands$pos + cands$length - 1 - fwd$pos + 1) >= AMP_MIN &
    (cands$pos + cands$length - 1 - fwd$pos + 1) <= AMP_MAX, ,
    drop = FALSE
  ]
  if (nrow(rev_candidates) == 0) next

  for (j in seq_len(nrow(rev_candidates))) {
    rv <- rev_candidates[j, ]
    tm_diff <- abs(fwd$tm - rv$tm)
    if (tm_diff > 5) next

    amplicon_len <- rv$pos + rv$length - 1 - fwd$pos + 1
    inner_start  <- fwd$end + 1
    inner_end    <- rv$pos - 1
    internal_div <- if (inner_end >= inner_start)
                    mean(divs[inner_start:inner_end]) else 0

    conservation <- 1 / (1 + fwd$div_sum + rv$div_sum)
    score <- conservation * (1 + internal_div) / (1 + tm_diff)

    pairs_list[[length(pairs_list) + 1]] <- data.frame(
      fwd_pos = fwd$pos, rev_pos = rv$pos,
      fwd = fwd$seq, rev = rv$seq,
      amplicon_len = amplicon_len,
      fwd_tm = fwd$tm, rev_tm = rv$tm, tm_diff = tm_diff,
      fwd_gc = fwd$gc, rev_gc = rv$gc,
      fwd_deg = fwd$degeneracy, rev_deg = rv$degeneracy,
      conservation_score = conservation,
      internal_diversity = internal_div,
      score = score,
      stringsAsFactors = FALSE
    )
  }
}

if (length(pairs_list) == 0) {
  warning("No valid primer pairs found.")
  empty_tsv()
  empty_plot()
  quit(save = "no", status = 0)
}

pairs_df <- do.call(rbind, pairs_list)

# reverse-complement the reverse primer for actual PCR ordering
pairs_df$rev <- vapply(pairs_df$rev, revcomp, character(1))

pairs_df <- pairs_df[order(-pairs_df$score), ]
pairs_df$primer_id <- sprintf("pair_%03d", seq_len(nrow(pairs_df)))
pairs_df <- pairs_df[seq_len(min(TOP_N, nrow(pairs_df))), ]

out_df <- pairs_df[, c(
  "primer_id", "fwd", "rev",
  "fwd_pos", "rev_pos", "amplicon_len",
  "fwd_tm", "rev_tm", "tm_diff",
  "fwd_gc", "rev_gc",
  "fwd_deg", "rev_deg",
  "conservation_score", "internal_diversity", "score"
)]
for (col in c("fwd_tm","rev_tm","tm_diff","fwd_gc","rev_gc",
              "conservation_score","internal_diversity","score")) {
  out_df[[col]] <- round(out_df[[col]], 3)
}

write.table(out_df, PRIMER_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[design_primers] wrote %d primer pairs to %s\n",
            nrow(out_df), PRIMER_TSV))

# ---------- diversity plot with top primer sites highlighted ----------
png(DIV_PLOT, width = 2400, height = 1200, res = 200)
par(mar = c(4, 4, 2, 1))
plot(divs, type = "p", pch = 16, cex = 0.3, col = "grey40",
     xlab = "Alignment position", ylab = "Shannon diversity",
     main = sprintf("Diversity profile and top %d primer sites",
                    nrow(out_df)))
lines(roll_div, col = "red", lwd = 1.2)

palette_n <- max(nrow(out_df), 1)
cols <- rainbow(palette_n, alpha = 0.5)
for (i in seq_len(nrow(out_df))) {
  fp <- pairs_df$fwd_pos[i]; rp <- pairs_df$rev_pos[i]
  flen <- nchar(pairs_df$fwd[i]); rlen <- nchar(pairs_df$rev[i])
  rect(fp, -0.05, fp + flen, max(divs) + 0.05, col = cols[i], border = NA)
  rect(rp, -0.05, rp + rlen, max(divs) + 0.05, col = cols[i], border = NA)
}
legend("topright",
       legend = c("per-position diversity", "rolling mean (k=10)",
                  "primer sites"),
       col    = c("grey40", "red", cols[1]),
       pch    = c(16, NA, 15),
       lty    = c(NA, 1, NA),
       bty    = "n")
dev.off()
cat(sprintf("[design_primers] diversity plot saved to %s\n", DIV_PLOT))
