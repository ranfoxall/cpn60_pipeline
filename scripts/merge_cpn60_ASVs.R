# ==========================================
# merge_cpn60_ASVs.R
# Merge ASV tables and FASTA files from multiple cpn60 pipeline runs
#
# Usage:
#   Rscript merge_cpn60_ASVs.R <output_prefix> \
#       <run1_Counts_seqASV_b.tsv> <run1_ASVs.fa> \
#       <run2_Counts_seqASV_b.tsv> <run2_ASVs.fa> \
#       [run3 ...]
#
# Example:
#   Rscript scripts/merge_cpn60_ASVs.R merged_cpn60 \
#       results/plate1_Counts_seqASV_b.tsv results/plate1_ASVs.fa \
#       results/plate2_Counts_seqASV_b.tsv results/plate2_ASVs.fa
#
# NOTE: Only sequence-based ASV tables (_b) are mergeable.
#       Do not merge numbered ASV tables directly.
#
# Developer:  Randi Foxall
# PI:         Dr. Cheryl Whistler
# Lab:        Whistler Lab, University of New Hampshire
# Funding:    NHAES CREATE program
# License:    CC BY-NC 4.0
# ==========================================

suppressPackageStartupMessages({
  library(dplyr)
  library(phyloseq)
})

# ----------------------
# 1. Parse arguments
# ----------------------
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3 || (length(args) - 1) %% 2 != 0) {
  stop(
    "Usage: Rscript merge_cpn60_ASVs.R <output_prefix> ",
    "<counts1.tsv> <asvs1.fa> [<counts2.tsv> <asvs2.fa> ...]\n",
    "At least one counts/FASTA pair is required."
  )
}

output_prefix <- args[1]
run_args      <- args[-1]

# Pair up counts and FASTA files
counts_files <- run_args[seq(1, length(run_args), by=2)]
fasta_files  <- run_args[seq(2, length(run_args), by=2)]

cat("Output prefix:", output_prefix, "\n")
cat("Merging", length(counts_files), "run(s):\n")
for (i in seq_along(counts_files)) {
  cat("  Run", i, "-- counts:", counts_files[i], "| fasta:", fasta_files[i], "\n")
}

# Check all files exist
for (f in c(counts_files, fasta_files)) {
  if (!file.exists(f)) stop("File not found: ", f)
}

# ----------------------
# 2. Load and merge counts tables
# ----------------------
# Sequence-based ASV tables use full sequences as row names,
# so identical ASVs across runs automatically collapse on merge.

cat("\nLoading counts tables...\n")
counts_list <- lapply(counts_files, function(f) {
  tab <- read.table(f, sep="\t", header=TRUE, row.names=1,
                    stringsAsFactors=FALSE, check.names=FALSE)
  cat("  Loaded:", f, "--", nrow(tab), "ASVs,", ncol(tab), "samples\n")
  tab
})

# Check for duplicate sample names across runs
all_samples <- unlist(lapply(counts_list, colnames))
dupes <- all_samples[duplicated(all_samples)]
if (length(dupes) > 0) {
  stop(
    "Duplicate sample names found across runs: ", paste(dupes, collapse=", "), "\n",
    "Sample names must be unique across all runs being merged."
  )
}

# Merge: union of all ASV sequences, fill missing with 0
cat("Merging counts tables...\n")
all_asvs <- unique(unlist(lapply(counts_list, rownames)))
merged_counts <- matrix(0L,
                        nrow = length(all_asvs),
                        ncol = length(all_samples),
                        dimnames = list(all_asvs, all_samples))

col_start <- 1
for (tab in counts_list) {
  col_end <- col_start + ncol(tab) - 1
  merged_counts[rownames(tab), colnames(tab)] <- as.matrix(tab)
  col_start <- col_end + 1
}

cat("Merged counts table:", nrow(merged_counts), "ASVs,",
    ncol(merged_counts), "samples\n")

# ----------------------
# 3. Write sequence-based merged counts table
# ----------------------
counts_seq_file <- paste0(output_prefix, "_merged_ASVs_b.tsv")
write.table(merged_counts, counts_seq_file,
            sep="\t", quote=FALSE, col.names=NA)
cat("Sequence-based counts written to:", counts_seq_file, "\n")

# ----------------------
# 4. Write numbered ASV counts table
# ----------------------
asv_seqs <- rownames(merged_counts)
asv_ids  <- paste0("ASV_", seq_along(asv_seqs))

merged_counts_num <- merged_counts
rownames(merged_counts_num) <- asv_ids

counts_num_file <- paste0(output_prefix, "_merged_ASVs_num.tsv")
write.table(merged_counts_num, counts_num_file,
            sep="\t", quote=FALSE, col.names=NA)
cat("Numbered counts written to:", counts_num_file, "\n")

# ----------------------
# 5. Write merged FASTA
# ----------------------
fasta_out  <- paste0(output_prefix, "_merged_ASVs_b.fa")
fasta_lines <- c(rbind(paste0(">", asv_ids), asv_seqs))
write(fasta_lines, fasta_out)
cat("FASTA written to:", fasta_out, "\n")

# ----------------------
# 6. Create phyloseq object
# ----------------------
sample_metadata <- data.frame(
  SampleID   = colnames(merged_counts_num),
  TotalReads = colSums(merged_counts_num),
  row.names  = colnames(merged_counts_num),
  stringsAsFactors = FALSE
)

ps <- phyloseq(
  otu_table(merged_counts_num, taxa_are_rows=TRUE),
  sample_data(sample_metadata)
)

ps_file <- paste0(output_prefix, "_merged_phyloseq.rds")
saveRDS(ps, ps_file)
cat("Phyloseq object written to:", ps_file, "\n")

# ----------------------
# 7. Summary
# ----------------------
cat("\nMerge complete. Outputs:\n")
cat("1. Sequence-based counts:", counts_seq_file, "\n")
cat("2. Numbered counts:      ", counts_num_file, "\n")
cat("3. Merged FASTA:         ", fasta_out, "\n")
cat("4. Phyloseq object:      ", ps_file, "\n")
cat("\nNext step -- run taxonomy classification:\n")
cat("  conda activate qiime2-2024.10.1\n")
cat("  Rscript scripts/cpn60_classify.R \\\n")
cat("    --output_prefix", output_prefix, "\\\n")
cat("    --classifier /path/to/cpn60_classifier_v11.qza\n")
