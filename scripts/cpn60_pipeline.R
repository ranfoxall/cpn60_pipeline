# ==========================================
# cpn60_pipeline: End-to-end cpn60 universal target microbiome pipeline
#
# Developer:  Randi Foxall
# PI:         Dr. Cheryl Whistler
# Lab:        Whistler Lab, University of New Hampshire
# Funding:    New Hampshire Agricultural Experiment Station
#             NHAES CREATE (Collaborative Research Enhancement
#             and Team Exploration) program
# ==========================================

suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
  library(dplyr)
  library(magrittr)
  library(tidyr)
  library(ggplot2)
  library(optparse)
  library(phyloseq)
})

# ----------------------
# 1. Parse command line arguments
# ----------------------
option_list <- list(
  make_option(c("--reads_path"), type="character", help="Path to primer-trimmed FASTQ reads"),
  make_option(c("--output_prefix"), type="character", help="Prefix for output files"),
  make_option(c("--error_model"), type="character", default="loess",
              help="Error model: 'loess' (default), 'default', or 'compare'")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$reads_path) | is.null(opt$output_prefix)) {
  stop("Please provide --reads_path and --output_prefix")
}

if (!opt$error_model %in% c("loess", "default", "compare")) {
  stop("--error_model must be one of: loess, default, compare")
}

reads_path    <- opt$reads_path
output_prefix <- opt$output_prefix
error_model   <- opt$error_model

cat("Reading FASTQ files from:", reads_path, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("Error model:", error_model, "\n")

# ----------------------
# 2. List forward and reverse reads
# ----------------------
fnFs <- sort(list.files(reads_path, pattern="_R1_001\\.fastq\\.gz$", full.names=TRUE))
fnRs <- sub("_R1_001\\.fastq\\.gz$", "_R2_001.fastq.gz", fnFs)

paired <- file.exists(fnRs)
fnFs <- fnFs[paired]
fnRs <- fnRs[paired]

stopifnot(length(fnFs) == length(fnRs))

# Robust sample name extraction: strip _S##_L###_R1_001.fastq.gz
sample_names <- sub("_S[0-9]+_L[0-9]+_R1_001\\.fastq\\.gz$", "", basename(fnFs))
names(fnFs) <- sample_names
names(fnRs) <- sample_names

cat("Found", length(fnFs), "samples:\n")
print(sample_names)

# ----------------------
# 3. Calculate run-level truncation lengths with safeguards
# ----------------------
# NOTE: cpn60 UT amplicon is ~560 bp. With 300 bp paired reads the overlap is
# already tight (~40 bp), so hard truncation is NOT used — truncating reads to
# even 200 bp each (400 bp total) would make merging impossible.
# Quality filtering relies entirely on maxEE instead.
# truncLen=0 tells filterAndTrim not to truncate.
cat("Truncation disabled for cpn60 UT amplicon (relying on maxEE filtering).\n")

repFs <- fnFs[1:min(3, length(fnFs))]
repRs <- fnRs[1:min(3, length(fnRs))]

pdf(file = paste0(output_prefix, "_quality_profiles.pdf"))
plotQualityProfile(repFs)
plotQualityProfile(repRs)
dev.off()

# ----------------------
# 4. Setup filtered FASTQ file paths
# ----------------------
filt_path <- paste0(output_prefix, "_filtered")
if (!dir.exists(filt_path)) {
  dir.create(filt_path)
  cat("Created filtered directory at:", filt_path, "\n")
} else {
  cat("Filtered directory already exists at:", filt_path, "\n")
}

filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

cat("Example filtered forward files:\n")
print(head(filtFs, 3))
cat("Example filtered reverse files:\n")
print(head(filtRs, 3))

# ----------------------
# 4b. Filter and trim
# ----------------------
cat("Filtering and trimming reads...\n")

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen    = c(0, 0),
  maxN        = 0,
  maxEE       = c(2, 5),
  truncQ      = 2,
  rm.phix     = TRUE,
  compress    = TRUE,
  multithread = TRUE
)
cat("Read counts after filtering:\n")
print(out)

# Drop samples with no reads after filtering
exists <- file.exists(filtFs) & file.size(filtFs) > 0
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]
sample_names <- sample_names[exists]

if (length(filtFs) == 0) stop("No reads survived filtering. Check truncation lengths and input quality.")
cat(length(filtFs), "samples remain after filtering.\n")

# ----------------------
# 5. Error model learning
# ----------------------

# LOESS error function (defined here for use by any model option)
loessErrfun_mod <- function(trans) {
  qq  <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=16, ncol=length(qq))
  rownames(est) <- paste0(rep(c("A","C","G","T"), each=4), "2",
                          rep(c("A","C","G","T"), 4))
  colnames(est) <- colnames(trans)
  for (nti in c("A","C","G","T")) {
    for (ntj in c("A","C","G","T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot  <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        df   <- data.frame(q=qq, rlogp=log10((errs+1)/(tot+1)))
        mod.lo <- loess(rlogp ~ q, df, span=0.95, degree=1)
        pred   <- predict(mod.lo, qq)
        pred[is.na(pred)] <- 0
        est[paste0(nti,"2",ntj),] <- 10^pred
      }
    }
  }
  est[est > 0.25] <- 0.25
  est[est < 1e-7] <- 1e-7
  return(est)
}

# Model fit metric (used only for 'compare' mode)
error_fit_metric <- function(errObj) {
  trans      <- errObj$trans
  trans_norm <- sweep(trans, 2, colSums(trans), "/")
  err_mat    <- errObj$err_out
  err_mat    <- err_mat[rownames(trans_norm), colnames(trans_norm)]
  return(sum((err_mat - trans_norm)^2))
}

err_cache_file <- paste0(output_prefix, "_error_model.rds")

if (file.exists(err_cache_file)) {
  cat("Loading saved error model from:", err_cache_file, "\n")
  err_cache <- readRDS(err_cache_file)
  errF <- err_cache$errF
  errR <- err_cache$errR
  cat("Error model loaded.\n\n")

} else {

  if (error_model == "loess") {
    cat("Learning error rates using LOESS model...\n")
    errF <- learnErrors(filtFs, multithread=TRUE, nbases=1e9,
                        errorEstimationFunction=loessErrfun_mod)
    errR <- learnErrors(filtRs, multithread=TRUE, nbases=1e9,
                        errorEstimationFunction=loessErrfun_mod)
    cat("LOESS error model complete.\n\n")

  } else if (error_model == "default") {
    cat("Learning error rates using DADA2 default model...\n")
    errF <- learnErrors(filtFs, multithread=TRUE, nbases=1e9)
    errR <- learnErrors(filtRs, multithread=TRUE, nbases=1e9)
    cat("Default error model complete.\n\n")

  } else if (error_model == "compare") {
    cat("Learning error rates using both default and LOESS models (this will take longer)...\n")

    cat("  Running default model...\n")
    errF_default <- learnErrors(filtFs, multithread=TRUE, nbases=1e9)
    errR_default <- learnErrors(filtRs, multithread=TRUE, nbases=1e9)

    cat("  Running LOESS model...\n")
    errF_loess <- learnErrors(filtFs, multithread=TRUE, nbases=1e9,
                              errorEstimationFunction=loessErrfun_mod)
    errR_loess <- learnErrors(filtRs, multithread=TRUE, nbases=1e9,
                              errorEstimationFunction=loessErrfun_mod)

    fit_default_F <- error_fit_metric(errF_default)
    fit_loess_F   <- error_fit_metric(errF_loess)
    fit_default_R <- error_fit_metric(errR_default)
    fit_loess_R   <- error_fit_metric(errR_loess)

    cat("Forward fit scores:  Default =", fit_default_F, " | LOESS =", fit_loess_F, "\n")
    cat("Reverse fit scores:  Default =", fit_default_R, " | LOESS =", fit_loess_R, "\n")

    errF <- if (fit_loess_F < fit_default_F) {
      cat("Selected LOESS for Forward.\n"); errF_loess
    } else {
      cat("Selected DEFAULT for Forward.\n"); errF_default
    }
    errR <- if (fit_loess_R < fit_default_R) {
      cat("Selected LOESS for Reverse.\n"); errR_loess
    } else {
      cat("Selected DEFAULT for Reverse.\n"); errR_default
    }
    cat("Error model selection complete.\n\n")
  }

  cat("Saving error model to:", err_cache_file, "\n")
  saveRDS(list(errF=errF, errR=errR), err_cache_file)
}

# ----------------------
# 6. Dereplicate and DADA2 denoise
# ----------------------
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")

# justConcatenate=TRUE is required for cpn60 UT amplicon sequenced with 2x300 bp reads.
# The ~549-567 bp UT region (after primer removal) exceeds the ~550 bp span of two
# 275 bp trimmed reads, so R1 and R2 do not overlap. Reads are concatenated with
# a string of Ns in the gap. The cpn60 classifier handles concatenated reads correctly.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                      justConcatenate=TRUE, verbose=TRUE)
seqtab  <- makeSequenceTable(mergers)

# ----------------------
# 7. Filter by expected cpn60 UT amplicon length
# ----------------------
# With justConcatenate=TRUE, sequences are R1 + 10 Ns + R2.
# Expected length: ~275 + 10 + 275 = ~560 bp, but can vary by sample.
# Using a generous range to capture all legitimate concatenated amplicons.

# Diagnostic: print actual length distribution before filtering
cat("Sequence length distribution (before length filter):\n")
print(table(nchar(colnames(seqtab))))

target_range <- 500:620
seqtab <- seqtab[, nchar(colnames(seqtab)) %in% target_range]
cat("Sequences in target length range (500-620 bp):", ncol(seqtab), "\n")

# ----------------------
# 8. Remove chimeras
# ----------------------
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
cat("Sequences after chimera removal:", ncol(seqtab_nochim), "\n")

# ----------------------
# 9. Create ASV tables
# ----------------------
asv_seqs <- colnames(seqtab_nochim)
rownames(seqtab_nochim) <- sample_names

# Sequence-based ASVs
asv_tab_seq <- t(seqtab_nochim)
rownames(asv_tab_seq) <- asv_seqs
write.table(asv_tab_seq,
            paste0(output_prefix, "_Counts_seqASV_b.tsv"),
            sep="\t", quote=FALSE, col.names=NA)

# Numbered ASVs
asv_ids     <- paste0("ASV_", seq_along(asv_seqs))
asv_tab_num <- asv_tab_seq
rownames(asv_tab_num) <- asv_ids
write.table(asv_tab_num,
            paste0(output_prefix, "_Counts_numASV.tsv"),
            sep="\t", quote=FALSE, col.names=NA)

# ----------------------
# 10. Write FASTA of numbered ASVs
# ----------------------
asv_fasta_path <- paste0(output_prefix, "_ASVs.fa")
asv_fasta      <- c(rbind(paste0(">", asv_ids), asv_seqs))
write(asv_fasta, asv_fasta_path)
cat("FASTA written to:", asv_fasta_path, "\n")

# ----------------------
# 11. Create phyloseq object (counts only, no taxonomy yet)
# ----------------------
sample_metadata <- data.frame(
  SampleID   = colnames(asv_tab_num),
  TotalReads = colSums(asv_tab_num),
  row.names  = colnames(asv_tab_num),
  stringsAsFactors = FALSE
)

ps <- phyloseq(
  otu_table(as.matrix(asv_tab_num), taxa_are_rows=TRUE),
  sample_data(sample_metadata)
)

ps_file <- paste0(output_prefix, "_phyloseq.rds")
saveRDS(ps, ps_file)
cat("Phyloseq object written to:", ps_file, "\n")

# ----------------------
# 12. Read tracking table
# ----------------------
getN <- function(x) sum(getUniques(x))

track <- data.frame(
  input     = out[, "reads.in"],
  filtered  = out[, "reads.out"],
  denoised  = sapply(dadaFs, getN),
  merged    = sapply(mergers, getN),
  tabled    = rowSums(seqtab),
  nonchim   = rowSums(seqtab_nochim),
  row.names = sample_names
)
track$pct_merged  <- round(track$merged  / track$input * 100, 1)
track$pct_nonchim <- round(track$nonchim / track$input * 100, 1)

track_file <- paste0(output_prefix, "_read_tracking.csv")
write.csv(track, track_file, quote=FALSE)
cat("Read tracking table written to:", track_file, "\n")

# ----------------------
# 13. Summary
# ----------------------
cat("\nPipeline finished. Outputs:\n")
cat("1. Sequence-based ASV table:", paste0(output_prefix, "_Counts_seqASV_b.tsv"), "\n")
cat("2. Numbered ASV table:      ", paste0(output_prefix, "_Counts_numASV.tsv"), "\n")
cat("3. FASTA for classifier:    ", asv_fasta_path, "\n")
cat("4. Phyloseq object:         ", ps_file, "\n")
cat("5. Quality profiles PDF:    ", paste0(output_prefix, "_quality_profiles.pdf"), "\n")
cat("6. Read tracking table:     ", track_file, "\n")
cat("\nNext step -- run taxonomy classification:\n")
cat("  conda activate <your_qiime2_env>\n")
cat("  Rscript cpn60_classify.R \\\n")
cat("    --output_prefix", output_prefix, "\\\n")
cat("    --classifier /path/to/cpn60_classifier_v11_sklearn142.qza\n")
