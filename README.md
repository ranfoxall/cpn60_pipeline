# cpn60_pipeline README

## Overview
This pipeline processes **cpn60 universal target (UT) amplicon sequencing data** using DADA2, producing ASVs, counts tables, FASTA files, taxonomy files, and a phyloseq object ready for downstream analysis. It is designed for the cpn60 UT primers described in Links et al. 2012 and Schellenberg et al. 2009, which target a ~560 bp hypervariable region of the cpn60 gene and capture broad bacterial diversity across environmental, food, and host-associated communities. Primer trimming uses QIIME2/cutadapt with four degenerate primers containing inosine bases (H279, H280, H1612, H1613). Taxonomy classification using the cpn60 QIIME2 classifier is run as a separate second script, allowing DADA2 and QIIME2 to run in their own conda environments. The pipeline can also merge multiple runs or plates into a single dataset. There is a special section (9) for running on the Premise HPC at UNH.

---

## Step 1: Clone the repository
```bash
git clone https://github.com/yourusername/cpn60_pipeline.git
cd cpn60_pipeline
```
This will create a directory called `cpn60_pipeline` containing the `scripts/` folder and any test data.

---

## Step 2: Install Miniconda (if not already installed)
Download Miniconda from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html) and follow instructions for your OS.

---

## Step 3: Create the conda environment

The `environment.yml` file sets up R 4.5 and all required base dependencies, including system libraries and HDF5 Bioconductor packages that must be installed via conda to avoid compilation errors:

```bash
conda env create -f environment.yml
conda activate cpn60_pipeline
```

> **Note:** QIIME2 is **not included** in this environment. See Step 4.

---

## Step 3b: Install R/Bioconductor packages

Most Bioconductor packages cannot be reliably installed via conda due to version conflicts. Run the provided installation script **once** after activating the environment:

```bash
conda activate cpn60_pipeline
Rscript install_packages.R
```

This script installs all packages in the correct order, skips anything already installed, and prints a version summary at the end. It uses Bioconductor 3.22 for compatibility with R 4.5.

> **If you hit errors:** The most common issues are:
> - `Rhtslib` failing with `fatal error: lzma.h: No such file or directory` — the `environment.yml` includes `xz` to prevent this, but if it occurs run:
>   ```
>   conda install -c conda-forge xz
>   Rscript install_packages.R
>   ```
> - `rhdf5filters` failing with a C compiler conflict — the `environment.yml` installs this via conda to prevent it, but if it occurs run:
>   ```
>   mamba install -c bioconda bioconductor-rhdf5filters bioconductor-rhdf5 bioconductor-biomformat bioconductor-phyloseq
>   Rscript install_packages.R
>   ```
> - `png` failing on Linux — you may need `sudo apt install libpng-dev` (or `libpng-devel` on RHEL/Rocky)

---

## Step 4: Install QIIME2 separately (optional)
QIIME2 is only required if you want to run primer trimming (Step 7) and taxonomy classification (Step 8). It must be installed in its own separate conda environment — do **not** install it into the `cpn60_pipeline` environment. Follow the official instructions:
[https://docs.qiime2.org/2024.10/install/](https://docs.qiime2.org/2024.10/install/)

Once installed, note the name of your QIIME2 environment (e.g. `qiime2-2020.2`) — you will need to set it as `QIIME2_ENV` in the slurm scripts and activate it manually when running `cpn60_classify.R` in Step 8.

> **On UNH Premise:** QIIME2 is already available as a conda environment via the `anaconda/colsa` module. See Step 9 for Premise-specific instructions — you do not need to install it yourself.

---

## Step 5: Get the cpn60 classifier (optional)
The cpn60 classifier is required only for the classification step (Step 8, `cpn60_classify.R`). There are two options:

**Option A — Download from GitHub (v11.1):**
```bash
wget https://github.com/HillLabSask/cpn60-Classifier/releases/download/v11.1/cpn60-q2-feature-classifier-v11.tar.gz
tar -xzf cpn60-q2-feature-classifier-v11.tar.gz
# Classifier will be at: cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
```

**Option B — Use shared cluster path (UNH Premise HPC users):**
```
/mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
```

> **Important:** The pre-trained classifier was built with QIIME2 v2022.11 (scikit-learn v0.24.1). Verify your QIIME2 version matches before running:
> ```bash
> conda activate <your_qiime2_env>
> qiime info | grep scikit
> ```
> If there is a version mismatch, you will need to retrain the classifier. See Step 9 (UNH Premise) for retraining instructions, or the [cpn60-Classifier GitHub page](https://github.com/HillLabSask/cpn60-Classifier/releases/tag/v11.1) for general instructions.

---

## Step 6: Prepare your reads folder
Organize your demultiplexed raw FASTQ files in a single folder. Files must follow standard Illumina naming:
```
raw_reads/
  Sample1_S1_L001_R1_001.fastq.gz
  Sample1_S1_L001_R2_001.fastq.gz
  Sample2_S2_L001_R1_001.fastq.gz
  Sample2_S2_L001_R2_001.fastq.gz
  ...
```
The pipeline extracts sample names by stripping `_S##_L###_R1_001.fastq.gz` from the filename, so `2B_2_new_S1_L001_R1_001.fastq.gz` becomes `2B_2_new`.

> Files **must** be paired — every `_R1_` file must have a matching `_R2_` file.

---

## Step 7: Run the cpn60 pipeline (primer trimming + DADA2)

The pipeline runs in two parts handled by the slurm scripts: first QIIME2/cutadapt strips the cpn60 UT primers, then DADA2 performs quality filtering, error modeling, denoising, merging, and ASV table generation. Taxonomy classification is handled separately in Step 8.

### About the cpn60 UT primers and inosines

This pipeline uses four degenerate cpn60 UT primers: two forward (H279, H1612) and two reverse (H280, H1613). All four contain inosine bases (I), which pair with any nucleotide and are represented as `N` in the cutadapt command. The primers are anchored to the 5' end of each read using `--p-front-f`/`--p-front-r`, and `--p-discard-untrimmed` removes any read pair where a primer is not found — ensuring only genuine cpn60 amplicons reach DADA2.

Primer sequences (inosines shown as I):

| Primer | Direction | Sequence (5'→3') |
|--------|-----------|-----------------|
| H279   | Forward   | `GAIIIIGCIGGIGAYGGIACIACIAC` |
| H1612  | Forward   | `GAIIIIGCIGGYGACGGYACSACSAC` |
| H280   | Reverse   | `YKIYKITCICCRAAICCIGGIGCYTT` |
| H1613  | Reverse   | `CGRCGRTCRCCGAAGCCSGGIGCCTT` |

> **Before running:** make sure the `cpn60_pipeline` conda environment is active for DADA2, and that QIIME2 is available for primer trimming. The slurm scripts handle environment switching automatically.

```bash
Rscript scripts/cpn60_pipeline.R \
  --reads_path /full/path/to/primer_trimmed_reads \
  --output_prefix cpn60_run1
```

**With custom error model:**
```bash
Rscript scripts/cpn60_pipeline.R \
  --reads_path /full/path/to/primer_trimmed_reads \
  --output_prefix cpn60_run1 \
  --error_model compare
```

> **Error model guidance:** `loess` is the default and recommended for cpn60 amplicon data. Use `compare` to empirically verify which model fits your data better — note this roughly doubles runtime.

**Outputs:**

| File | Description |
|------|-------------|
| `{output_prefix}_Counts_seqASV_b.tsv` | Counts table with ASV sequences as row names |
| `{output_prefix}_Counts_numASV.tsv` | Counts table with numbered ASVs (ASV_1, ASV_2, …) |
| `{output_prefix}_ASVs.fa` | FASTA of ASVs (input for classification script) |
| `{output_prefix}_phyloseq.rds` | Phyloseq object with numbered ASVs and sample metadata |
| `{output_prefix}_quality_profiles.pdf` | Quality profile plots for representative samples |
| `{output_prefix}_read_tracking.csv` | Per-sample read counts at each step: input → filtered → denoised → merged → nonchim, with % merged and % retained |
| `{output_prefix}_error_model.rds` | Cached DADA2 error model — reused automatically on restart to avoid relearning |

---

## Step 8: Run taxonomy classification (QIIME2)

This script takes the FASTA output from Step 7 and runs the cpn60 classifier. It must be run in a **separate environment** from Step 7 — you need to deactivate the DADA2 environment and activate your QIIME2 environment first.

> **Before running:** switch conda environments:
> ```bash
> conda deactivate                        # deactivate cpn60_pipeline
> conda activate <your_qiime2_env>        # name set by QIIME2_ENV in slurm scripts, e.g. qiime2-2020.2
> ```

```bash
Rscript scripts/cpn60_classify.R \
  --output_prefix cpn60_run1 \
  --classifier    /path/to/cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
```

> **UNH Premise users** — use this path:
> ```
> --classifier /mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
> ```

The script automatically finds `{output_prefix}_ASVs.fa`, `{output_prefix}_Counts_numASV.tsv`, and `{output_prefix}_phyloseq.rds` from the prefix. If you need to point to a phyloseq file in a different location you can override it with `--phyloseq_rds /path/to/file.rds`.

**Outputs:**

| File | Description |
|------|-------------|
| `{output_prefix}_taxonomy.tsv` | Raw taxonomy assignments from cpn60 classifier (QIIME2 format) |
| `{output_prefix}_taxonomy_table.csv` | Taxonomy table with one column per rank (Kingdom–Species), phyloseq-compatible |
| `{output_prefix}_taxonomy_confidence.csv` | Per-ASV classifier confidence scores |
| `{output_prefix}_phyloseq_taxonomy.rds` | Complete phyloseq object with otu_table, tax_table, refseq, and sample_data |

> **If classification fails** with a scikit-learn version error, you need to retrain the classifier for your QIIME2 version. See the UNH Premise section for retraining instructions using the reference files in the shared directory.

---

## Step 9: Running on the UNH Premise Cluster

Premise is UNH's HPC cluster. This section covers everything you need to run the pipeline there.

### Loading software on Premise

Premise uses environment modules to manage software. Two module collections are relevant:

- `linuxbrew/colsa` — provides R and bioinformatics tools
- `anaconda/colsa` — provides conda environments including QIIME2

> **Note:** You cannot have both `linuxbrew/colsa` and `anaconda/colsa` loaded at the same time. Unload one before loading the other.

Load R via linuxbrew:
```bash
module purge
module load linuxbrew/colsa
Rscript --version  # confirm R is available
```

Load QIIME2 via anaconda:
```bash
module purge
module load anaconda/colsa
conda activate <your_qiime2_env>   # e.g. qiime2-2020.2 — set QIIME2_ENV in slurm scripts
qiime info
```

> **⚠️ Classifier compatibility warning:** The pre-trained cpn60 classifier (`cpn60_classifier_v11.qza`) was built with QIIME2 v2022.11 (scikit-learn v0.24.1). Premise currently has QIIME2 2024.10.1, which uses a newer version of scikit-learn. This **will cause a version mismatch error** when running classification. To fix this, retrain the classifier in the newer QIIME2 environment:
> ```bash
> module load anaconda/colsa
> conda activate <your_qiime2_env>   # e.g. qiime2-2020.2
>
> # Import reference sequences and taxonomy
> qiime tools import \
>   --type 'FeatureData[Sequence]' \
>   --input-path /mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/cpn60_v11_seqs.fasta \
>   --output-path cpn60_v11_seqs.qza
>
> qiime tools import \
>   --type 'FeatureData[Taxonomy]' \
>   --input-format HeaderlessTSVTaxonomyFormat \
>   --input-path /mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/cpn60_v11_taxonomy_table.txt \
>   --output-path cpn60_v11_taxonomy.qza
>
> # Train the classifier
> qiime feature-classifier fit-classifier-naive-bayes \
>   --i-reference-reads cpn60_v11_seqs.qza \
>   --i-reference-taxonomy cpn60_v11_taxonomy.qza \
>   --o-classifier cpn60_classifier_v11_sklearn142.qza
> ```
> This only needs to be done once. The reference files are already available in the shared directory.

### cpn60 classifier location on Premise

```
/mnt/home/whistler/shared/cpn60-Classifier/
  cpn60-q2-feature-classifier-v11/
    cpn60_classifier_v11.qza          # original classifier (QIIME2 v2022.11 / sklearn 0.24.1)
    cpn60_v11_seqs.fasta              # reference sequences (for retraining)
    cpn60_v11_taxonomy_table.txt      # taxonomy table (for retraining)
  cpn60_classifier_v11_sklearn142.qza # retrained classifier (QIIME2 2024.10.26 / sklearn 1.4.2)
```

### Submitting as a SLURM job

There are two approaches — pick whichever suits you.

Before running either option, edit the variables block at the top of the script:

| Variable | Description |
|----------|-------------|
| `READS_PATH` | Path to the folder containing your raw paired FASTQ files |
| `OUTPUT_PREFIX` | A short name for your run — all output files will be named with this prefix |
| `ERROR_MODEL` | DADA2 error model — use `loess` unless you have a reason to change it |
| `CLASSIFIER` | Full path to the cpn60 QIIME2 classifier `.qza` file |
| `WORKDIR` | Working directory where outputs will be written — does not need to contain the scripts |
| `SCRIPTS_DIR` | Path to the `scripts/` folder containing the R scripts — can be different from `WORKDIR` |
| `CONDA_ENV` | Name of your R/DADA2 conda environment (default: `cpn60_pipeline`) |
| `QIIME2_ENV` | Name of your QIIME2 conda environment (default: `qiime2-2020.2`) |

---

#### Option A: Two separate jobs (recommended)

Submit DADA2 and classification as separate jobs, with the second depending on the first finishing successfully.

```bash
JOB1=$(sbatch --parsable run_cpn60_dada2.slurm)
echo "Submitted DADA2 job: $JOB1"

JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 run_cpn60_classify.slurm)
echo "Submitted classify job: $JOB2 (depends on $JOB1)"
```

---

#### Option B: Single combined script

Submit everything at once with `run_cpn60_full.slurm`:

```bash
sbatch run_cpn60_full.slurm
```

> **Note:** Each step checks for success before continuing. If primer trimming or DADA2 fails, the job exits cleanly rather than running subsequent steps on missing files.

---

### Restarting a failed job

Both `run_cpn60_dada2.slurm` and `run_cpn60_full.slurm` have restart flags at the top so you can skip steps that already completed successfully.

| Flag | Default | What it skips |
|------|---------|---------------|
| `SKIP_IMPORT` | `0` | QIIME2 import of raw FASTQs |
| `SKIP_CUTADAPT` | `0` | Primer trimming |
| `SKIP_DADA2` | `0` | DADA2 ASV inference |
| `SKIP_CLASSIFY` | `0` | Taxonomy classification (`run_cpn60_full.slurm` only) |

Set any flag to `1` to skip that step.

**`TRIMMED_QZA` — required when skipping cutadapt:**

When `SKIP_CUTADAPT=1`, set `TRIMMED_QZA` to the trimmed reads artifact from the previous run. The script re-extracts FASTQs from it automatically:
```bash
TRIMMED_QZA=cpn60_run1_trimmed-reads.qza
```

**Common restart scenarios:**

*Trimming succeeded but DADA2 failed:*
```bash
SKIP_IMPORT=1
SKIP_CUTADAPT=1
TRIMMED_QZA=cpn60_run1_trimmed-reads.qza
```
Then resubmit with `sbatch run_cpn60_dada2.slurm` or `sbatch run_cpn60_full.slurm`.

> **Note:** If the error model cache (`{output_prefix}_error_model.rds`) exists from a previous run, DADA2 will load it automatically and skip the error learning step — saving significant time.

*DADA2 succeeded but classification failed (`run_cpn60_full.slurm` only):*
```bash
SKIP_IMPORT=1
SKIP_CUTADAPT=1
SKIP_DADA2=1
TRIMMED_QZA=cpn60_run1_trimmed-reads.qza
```
Then resubmit with `sbatch run_cpn60_full.slurm`, or run classification directly with `sbatch run_cpn60_classify.slurm`.

*Classification succeeded but taxonomy merge into phyloseq failed:*

Run `cpn60_add_taxonomy.R` directly from the `cpn60_pipeline` conda environment:
```bash
conda activate cpn60_pipeline
Rscript scripts/cpn60_add_taxonomy.R --output_prefix cpn60_run1
```
This reads `{output_prefix}_phyloseq.rds`, `{output_prefix}_taxonomy_table.csv`, and `{output_prefix}_ASVs.fa` and produces `{output_prefix}_phyloseq_taxonomy.rds`.

---

Monitor your jobs:
```bash
squeue -u $USER                          # check job status
sacct -u $USER                           # view completed/failed jobs
slurm-monitor <job_id>                   # monitor CPU and memory usage live
tail -f cpn60_full_<job_id>.log          # stream combined log
scancel <job_id>                         # cancel a job if needed
```

> For questions about Premise, contact [Toni Westbrook](mailto:anthony.westbrook@unh.edu) or see the [Premise software page](https://premise.sr.unh.edu/software.html).

---

## Step 10: (Optional) Merge Multiple Plates / Runs
If you have multiple sequencing runs or plates, merge them using the merge script:

```bash
Rscript scripts/merge_cpn60_ASVs.R merged_run \
    results/plate1_Counts_seqASV_b.tsv results/plate1_ASVs.fa \
    results/plate2_Counts_seqASV_b.tsv results/plate2_ASVs.fa
```

**Outputs:**

| File | Description |
|------|-------------|
| `merged_run_merged_ASVs_b.tsv` | Counts table with ASV sequences as row names |
| `merged_run_merged_ASVs_b.fa` | Merged ASV FASTA |
| `merged_run_merged_ASVs_num.tsv` | Counts table with numbered ASVs |
| `merged_run_merged_phyloseq.rds` | Phyloseq object with sample metadata |

> **Note:** Only sequence-based ASV tables (`_b`) are mergeable. Do not merge numbered ASV tables directly.

---

## Step 11: Use the phyloseq object in R
```R
library(phyloseq)
# complete object (after classification):
ps <- readRDS("{output_prefix}_phyloseq_taxonomy.rds")
otu_table(ps)      # ASV counts (numbered ASVs)
sample_data(ps)    # SampleID and TotalReads
tax_table(ps)      # Kingdom through Species, rank prefixes (k__, p__, etc.)
refseq(ps)         # ASV sequences (DNAStringSet)

# counts-only object (after DADA2, before classification):
ps_counts <- readRDS("{output_prefix}_phyloseq.rds")
```

---

## Directory structure
```
cpn60_pipeline/
  scripts/
    cpn60_pipeline.R          # DADA2 pipeline (Step 7)
    cpn60_classify.R          # QIIME2 classification (Step 8)
    cpn60_add_taxonomy.R      # merges taxonomy + refseq into phyloseq RDS (Step 8)
    merge_cpn60_ASVs.R        # Multi-run merge (Step 10)
  slurm/
    run_cpn60_dada2.slurm      # SLURM job for primer trimming + DADA2
    run_cpn60_classify.slurm   # SLURM job for classification + phyloseq update
    run_cpn60_full.slurm       # SLURM job for full pipeline
  tests/
    example FASTQ files (optional)
  environment.yml           # conda environment (R base + dependencies)
  install_packages.R        # Bioconductor package installer (run once)
  README.md
  LICENSE
```

---

## Arguments reference

### cpn60_pipeline.R

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--reads_path` | Yes | — | Path to folder containing primer-trimmed paired FASTQ files |
| `--output_prefix` | Yes | — | Prefix for all output files |
| `--error_model` | No | `loess` | Error model: `loess` (recommended), `default`, or `compare` (run both, pick best — slowest) |

### cpn60_classify.R

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--output_prefix` | Yes | — | Output prefix used in cpn60_pipeline.R — all input files are derived from this |
| `--classifier` | Yes | — | Path to cpn60 QIIME2 classifier `.qza` |
| `--phyloseq_rds` | No | `{prefix}_phyloseq.rds` | Override the phyloseq file path if it is not in the same directory |

### cpn60_add_taxonomy.R

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--output_prefix` | Yes | — | Output prefix used in cpn60_pipeline.R — used to locate phyloseq RDS, taxonomy CSV, and FASTA |

> This script is called automatically by `run_cpn60_full.slurm` and `run_cpn60_classify.slurm`. It merges the taxonomy table and ASV sequences (refseq) into the phyloseq object, producing `{output_prefix}_phyloseq_taxonomy.rds`.

---

## Citation

If you use this pipeline in your research, please cite it as:

> Foxall, R., Whistler, C. and Jones, S. (2026). *cpn60_pipeline: A cpn60 universal target amplicon sequencing pipeline*. Whistler Lab, University of New Hampshire. https://github.com/yourusername/cpn60_pipeline

A `CITATION.cff` file is included in this repository. GitHub will automatically display a **"Cite this repository"** button on the repository page.

Please also cite the cpn60 UT primer publications on which this pipeline is based:

> Schellenberg J, Links MG, Hill JE, Dumonceaux TJ, Peters GA, Tyler S, Kimpe SR, Severini A and Plummer FA (2009). Pyrosequencing of the chaperonin-60 universal target as a tool for determining microbial community composition. *Applied and Environmental Microbiology* 75(9):2889–2898. https://doi.org/10.1128/AEM.02645-08

> Links MG, Dumonceaux TJ, Hemmingsen SM and Hill JE (2012). The chaperonin-60 universal target is a barcode for bacteria and mitochondria. *PLoS ONE* 7(9):e49755. https://doi.org/10.1371/journal.pone.0049755

### Tool Citations

Please also cite the following tools used by this pipeline:

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA and Holmes SP (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods* 13:581–583. https://doi.org/10.1038/nmeth.3869

> Bolyen E, Rideout JR, Dillon MR et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology* 37:852–857. https://doi.org/10.1038/s41587-019-0209-9

> Martin M (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal* 17(1):10–12. https://doi.org/10.14806/ej.17.1.200

> McMurdie PJ and Holmes S (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE* 8(4):e61217. https://doi.org/10.1371/journal.pone.0061217

> R Core Team (2023). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

---

## Funding and Acknowledgements

This work was developed by Randi Foxall in the Whistler Lab at the University of New Hampshire under the supervision of Dr. Cheryl Whistler and Dr. Stephen Jones.

This work was supported by the New Hampshire Agricultural Experiment Station through the NHAES CREATE (Collaborative Research Enhancement and Team Exploration) program, which is part of the New Hampshire Agricultural Experiment Station at the University of New Hampshire.

Thanks to Dr. Cheryl Whistler (UNH, Molecular Cellular and Biomedical Sciences), Dr. Stephen Jones (UNH, Natural Resources and the Environment), and Dr. Ashley Busco (UNH, Biological Sciences) for their support and funding contributions to this work.

---

## Authors

| Name | Role | Affiliation |
|------|------|-------------|
| Randi Foxall | Developer | Whistler Lab, University of New Hampshire |
| Cheryl Whistler | Principal Investigator & Funder | Molecular Cellular and Biomedical Sciences, University of New Hampshire |
| Stephen Jones | Co-Investigator & Funder | Department of Natural Resources and the Environment, University of New Hampshire |
| Ashley Busco | Co-Investigator & Funder | Biological Sciences, University of New Hampshire |

---

## License

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).

You are free to share and adapt this material for non-commercial purposes, provided appropriate credit is given. Commercial use requires written permission from the authors.
