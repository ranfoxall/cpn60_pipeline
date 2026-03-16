# cpn60_classify.R
# QIIME2 taxonomy classification for cpn60 universal target ASVs
#
# Run AFTER cpn60_pipeline.R with your QIIME2 environment active:
#   conda activate qiime2-2020.2
#   Rscript cpn60_classify.R \
#     --output_prefix cpn60_run1 \
#     --classifier    /path/to/cpn60_classifier_v11.qza
#
# Note: this script uses base R only (no optparse) so it runs in
# the QIIME2 conda environment without extra R package installation.
#
# UNH Premise classifier path:
#   /mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
#
# Developer:  Randi Foxall
# PI:         Dr. Cheryl Whistler
# Lab:        Whistler Lab, University of New Hampshire
# Funding:    NHAES CREATE program, NH EPSCoR
# License:    CC BY-NC 4.0

# phyloseq is loaded conditionally below -- it is not available in all
# QIIME2 conda environments. The core classification steps do not require it.

# [1] arguments -------------------------------------------------------
# Uses base R commandArgs() -- no optparse dependency required.
# This allows the script to run in the QIIME2 conda environment (R 3.5+)
# without needing additional R package installations.
#
# Usage:
#   Rscript cpn60_classify.R \
#     --output_prefix cpn60_run1 \
#     --classifier /path/to/cpn60_classifier_v11.qza \
#     [--phyloseq_rds /path/to/file.rds]

args <- commandArgs(trailingOnly=TRUE)

parse_args_simple <- function(args) {
  result <- list(output_prefix=NULL, classifier=NULL, phyloseq_rds=NULL)
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--output_prefix" && i < length(args)) {
      result$output_prefix <- args[i+1]; i <- i + 2
    } else if (args[i] == "--classifier" && i < length(args)) {
      result$classifier <- args[i+1]; i <- i + 2
    } else if (args[i] == "--phyloseq_rds" && i < length(args)) {
      result$phyloseq_rds <- args[i+1]; i <- i + 2
    } else {
      i <- i + 1
    }
  }
  result
}

opt <- parse_args_simple(args)

for (arg in c("output_prefix", "classifier")) {
  if (is.null(opt[[arg]])) stop("Please provide --", arg)
}

output_prefix    <- opt$output_prefix
cpn60_classifier <- opt$classifier
phyloseq_rds     <- if (!is.null(opt$phyloseq_rds)) opt$phyloseq_rds else paste0(output_prefix, "_phyloseq.rds")

asv_fasta_path <- paste0(output_prefix, "_ASVs.fa")
counts_table   <- paste0(output_prefix, "_Counts_numASV.tsv")

for (f in c(asv_fasta_path, counts_table, cpn60_classifier)) {
  if (!file.exists(f)) stop("File not found: ", f)
}

cat("output_prefix: ", output_prefix, "\n")
cat("classifier:    ", cpn60_classifier, "\n")

# [2] check QIIME2 ----------------------------------------------------
qiime_check <- system2("which", args="qiime", stdout=TRUE, stderr=TRUE)
if (length(qiime_check) == 0 || !grepl("qiime", qiime_check[1])) {
  stop(
    "QIIME2 not found in PATH.\n",
    "Activate your QIIME2 environment before running:\n",
    "  conda activate qiime2-2024.10.1"
  )
}
cat("QIIME2 found at:", qiime_check[1], "\n")

# [3] classify --------------------------------------------------------
queries_qza         <- paste0(output_prefix, "_ASVs.qza")
taxonomy_qza        <- paste0(output_prefix, "_taxonomy.qza")
taxonomy_tsv        <- paste0(output_prefix, "_taxonomy.tsv")
taxonomy_export_dir <- paste0(output_prefix, "_taxonomy_export")

cat("Importing ASVs to QIIME2...\n")
ret <- system2("qiime", args=c(
  "tools", "import",
  "--type", "FeatureData[Sequence]",
  "--input-path",  asv_fasta_path,
  "--output-path", queries_qza
))
if (ret != 0) stop("QIIME2 import failed. Check that your QIIME2 environment is activated.")

cat("Classifying ASVs...\n")
ret <- system2("qiime", args=c(
  "feature-classifier", "classify-sklearn",
  "--i-classifier",     cpn60_classifier,
  "--i-reads",          queries_qza,
  "--o-classification", taxonomy_qza
))
if (ret != 0) stop(
  "Classification failed.\n",
  "Common cause: scikit-learn version mismatch.\n",
  "cpn60_classifier_v11.qza was trained with scikit-learn 0.24.1 (QIIME2 v2022.11).\n",
  "If your QIIME2 is newer, retrain the classifier -- see README."
)

cat("Exporting taxonomy...\n")
if (!dir.exists(taxonomy_export_dir)) dir.create(taxonomy_export_dir)
system2("qiime", args=c(
  "tools", "export",
  "--input-path",  taxonomy_qza,
  "--output-path", taxonomy_export_dir
))
file.rename(file.path(taxonomy_export_dir, "taxonomy.tsv"), taxonomy_tsv)
cat("Taxonomy written to:", taxonomy_tsv, "\n")

# [4] build taxonomy table --------------------------------------------
cat("Building taxonomy table...\n")

taxonomy_raw <- read.table(taxonomy_tsv, sep="\t", header=TRUE,
                            stringsAsFactors=FALSE, row.names=1)

# Split semicolon-delimited taxonomy string into individual rank columns.
# phyloseq expects a matrix with one column per rank.
# Ranks present in cpn60 UT classifier: Kingdom, Phylum, Class, Order, Family, Genus, Species
# rank prefixes used to fill unclassified ranks (standard phyloseq convention)
tax_ranks   <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_prefix <- c("k__",     "p__",    "c__",   "o__",   "f__",    "g__",   "s__")

tax_split <- strsplit(taxonomy_raw$Taxon, ";\\s*")
tax_mat <- t(sapply(tax_split, function(x) {
  length(x) <- length(tax_ranks)
  x
}))
rownames(tax_mat) <- rownames(taxonomy_raw)
colnames(tax_mat) <- tax_ranks

# replace missing/ambiguous values with the rank prefix (e.g. c__, s__)
# this is the standard QIIME2/phyloseq convention for unclassified ranks
for (i in seq_along(tax_ranks)) {
  empty <- is.na(tax_mat[, i]) | tax_mat[, i] == "" | tax_mat[, i] == "-" |
           tax_mat[, i] == rank_prefix[i]
  tax_mat[empty, i] <- rank_prefix[i]
}

# taxonomy table (phyloseq-compatible: ASV x rank, prefixes retained, no confidence)
taxonomy_table      <- as.data.frame(tax_mat, stringsAsFactors=FALSE)
taxonomy_table_file <- paste0(output_prefix, "_taxonomy_table.csv")
write.csv(taxonomy_table, taxonomy_table_file, quote=FALSE)
cat("Taxonomy table written to:", taxonomy_table_file, "\n")

# confidence scores (separate CSV)
confidence_file <- paste0(output_prefix, "_taxonomy_confidence.csv")
write.csv(data.frame(Confidence=taxonomy_raw$Confidence, row.names=rownames(taxonomy_raw)),
          confidence_file, quote=FALSE)
cat("Confidence scores written to:", confidence_file, "\n")

# [5] update phyloseq -------------------------------------------------
# phyloseq is loaded here conditionally -- skipped if not available
# (e.g. when running in a QIIME2 conda environment without phyloseq installed)
if (!requireNamespace("phyloseq", quietly=TRUE)) {
  cat("phyloseq not available -- skipping phyloseq update.\n")
  cat("To add taxonomy to your phyloseq object, run in your R environment:\n")
  cat("  library(phyloseq)\n")
  cat("  ps <- readRDS(\"", paste0(output_prefix, "_phyloseq.rds"), "\")\n", sep="")
  cat("  tax_mat <- as.matrix(read.csv(\"", taxonomy_table_file, "\", row.names=1))\n", sep="")
  cat("  ps <- merge_phyloseq(ps, tax_table(tax_mat))\n")
  cat("  saveRDS(ps, \"", paste0(output_prefix, "_phyloseq_taxonomy.rds"), "\")\n", sep="")
} else {
  library(phyloseq)
  if (!file.exists(phyloseq_rds)) {
    warning("Phyloseq RDS not found: ", phyloseq_rds, " -- skipping.")
  } else {
    cat("Adding taxonomy to phyloseq object...\n")
    ps          <- readRDS(phyloseq_rds)
    ps_tax_mat  <- as.matrix(taxonomy_table)
    ps          <- merge_phyloseq(ps, tax_table(ps_tax_mat))
    ps_tax_file <- paste0(output_prefix, "_phyloseq_taxonomy.rds")
    saveRDS(ps, ps_tax_file)
    cat("Phyloseq with taxonomy written to:", ps_tax_file, "\n")
  }
}

# [6] summary ---------------------------------------------------------
cat("\nClassification complete. Outputs:\n")
cat("  ", taxonomy_table_file, "\n")
cat("  ", confidence_file, "\n")
if (file.exists(paste0(output_prefix, "_phyloseq_taxonomy.rds")))
  cat("  ", paste0(output_prefix, "_phyloseq_taxonomy.rds"), "\n")
