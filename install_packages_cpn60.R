# install_packages.R
# cpn60_pipeline – R/Bioconductor Package Installation
# R 4.5 / Bioconductor 3.22
#
# Run this script ONCE after creating and activating the
# cpn60_pipeline conda environment:
#
#   conda activate cpn60_pipeline
#   Rscript install_packages.R
#
# Installation order matters -- do not rearrange sections.

# Helper: install only if missing
install_if_missing <- function(pkg, bioc=FALSE, version=NULL) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    if (!is.null(version)) {
      message("Installing ", pkg, " version ", version, " ...")
      remotes::install_version(pkg, version=version, repos="https://cran.r-project.org")
    } else if (bioc) {
      message("Installing ", pkg, " from Bioconductor ...")
      BiocManager::install(pkg, ask=FALSE, update=FALSE)
    } else {
      message("Installing ", pkg, " from CRAN ...")
      install.packages(pkg, repos="https://cran.r-project.org")
    }
  } else {
    message(pkg, " already installed -- skipping.")
  }
}

# 1. BiocManager -- pin to Bioconductor 3.22 for R 4.5

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cran.r-project.org")

BiocManager::install(version="3.22", ask=FALSE)

# 2. remotes -- needed for pinned CRAN versions

install_if_missing("remotes")
library(remotes)


# 3. Core CRAN dependencies (order matters)


install_if_missing("Matrix")
library(Matrix)
message("Matrix: ", packageVersion("Matrix"))

install_if_missing("MASS")
library(MASS)
message("MASS: ", packageVersion("MASS"))

install_if_missing("mgcv")
library(mgcv)
message("mgcv: ", packageVersion("mgcv"))

install_if_missing("png")
library(png)
message("png: ", packageVersion("png"))

install_if_missing("data.table")
library(data.table)
message("data.table: ", packageVersion("data.table"))

install_if_missing("igraph")
library(igraph)
message("igraph: ", packageVersion("igraph"))

install_if_missing("latticeExtra")
library(latticeExtra)
message("latticeExtra: ", packageVersion("latticeExtra"))

install_if_missing("dplyr")
library(dplyr)
message("dplyr: ", packageVersion("dplyr"))

install_if_missing("magrittr")
library(magrittr)
message("magrittr: ", packageVersion("magrittr"))

install_if_missing("tidyr")
library(tidyr)
message("tidyr: ", packageVersion("tidyr"))

install_if_missing("ggplot2")
library(ggplot2)
message("ggplot2: ", packageVersion("ggplot2"))

install_if_missing("optparse")
library(optparse)
message("optparse: ", packageVersion("optparse"))

# 4. Bioconductor genomic packages (order matters)

# Rhtslib requires xz/zlib/bzip2/curl/openssl system libs
# If this fails with "lzma.h not found", run:
#   conda install -c conda-forge xz zlib bzip2 curl openssl
install_if_missing("Rhtslib", bioc=TRUE)
library(Rhtslib)
message("Rhtslib: ", packageVersion("Rhtslib"))

# rhdf5filters/rhdf5/biomformat -- install via conda if compilation fails:
#   mamba install -c bioconda bioconductor-rhdf5filters bioconductor-rhdf5 bioconductor-biomformat
install_if_missing("rhdf5filters", bioc=TRUE)
install_if_missing("rhdf5", bioc=TRUE)
install_if_missing("biomformat", bioc=TRUE)

install_if_missing("DelayedArray", bioc=TRUE)
library(DelayedArray)
message("DelayedArray: ", packageVersion("DelayedArray"))

install_if_missing("Biostrings", bioc=TRUE)
install_if_missing("XVector", bioc=TRUE)
library(Biostrings)
library(XVector)
message("Biostrings: ", packageVersion("Biostrings"))
message("XVector: ", packageVersion("XVector"))

install_if_missing("GenomicRanges", bioc=TRUE)
install_if_missing("GenomicAlignments", bioc=TRUE)
install_if_missing("Rsamtools", bioc=TRUE)
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)
message("GenomicRanges: ", packageVersion("GenomicRanges"))
message("GenomicAlignments: ", packageVersion("GenomicAlignments"))
message("Rsamtools: ", packageVersion("Rsamtools"))

install_if_missing("ShortRead", bioc=TRUE)
library(ShortRead)
message("ShortRead: ", packageVersion("ShortRead"))


# 5. Pipeline-specific packages

install_if_missing("multtest", bioc=TRUE)
library(multtest)
message("multtest: ", packageVersion("multtest"))

install_if_missing("dada2", bioc=TRUE)
install_if_missing("phyloseq", bioc=TRUE)
install_if_missing("vegan")
install_if_missing("ade4")
library(dada2)
library(phyloseq)
library(biomformat)
library(vegan)
library(ade4)
message("dada2: ", packageVersion("dada2"))
message("phyloseq: ", packageVersion("phyloseq"))
message("biomformat: ", packageVersion("biomformat"))
message("vegan: ", packageVersion("vegan"))
message("ade4: ", packageVersion("ade4"))

# 6. Final check -- print all installed versions

pkgs <- c("Matrix", "MASS", "mgcv", "png", "data.table", "igraph",
          "latticeExtra", "dplyr", "magrittr", "tidyr", "ggplot2",
          "optparse", "DelayedArray", "Rhtslib", "rhdf5filters",
          "rhdf5", "Biostrings", "XVector", "GenomicRanges",
          "GenomicAlignments", "Rsamtools", "ShortRead", "multtest",
          "dada2", "phyloseq", "biomformat", "vegan", "ade4")

message("\n============================================================")
message("Installed package versions:")
message("============================================================")
for (p in pkgs) {
  if (requireNamespace(p, quietly=TRUE)) {
    message(sprintf("  %-25s %s", p, packageVersion(p)))
  } else {
    message(sprintf("  %-25s *** NOT INSTALLED ***", p))
  }
}
message("============================================================")
message("Installation complete.")
