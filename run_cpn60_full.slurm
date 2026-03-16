#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=512000
#SBATCH -p shared
#SBATCH --job-name="cpn60_full"
#SBATCH --output=cpn60_full_%j.log
#SBATCH --error=cpn60_full_%j.err
#SBATCH --exclude=node117,node118

# [variables] ----------------------------------------------------------
# edit these before submitting
READS_PATH=/mnt/home/whistler/foxall/Scripts/CREATE/cpn60/reads
OUTPUT_PREFIX=cpn60_run1
ERROR_MODEL=loess
CLASSIFIER=/mnt/home/whistler/shared/cpn60-Classifier/cpn60_classifier_v11_sklearn142.qza
WORKDIR=/mnt/home/whistler/foxall/Scripts/CREATE/cpn60
SCRIPTS_DIR=/mnt/home/whistler/foxall/Scripts/CREATE/cpn60_pipeline/scripts
CONDA_ENV=cpn60_pipeline
QIIME2_ENV=qiime2-2020.2

# [restart] ------------------------------------------------------------
# set to 1 to skip steps already completed in a previous run
SKIP_IMPORT=0
SKIP_CUTADAPT=0
SKIP_DADA2=0
SKIP_CLASSIFY=0

# if skipping cutadapt, point to the trimmed .qza from the previous run
TRIMMED_QZA=

cd $WORKDIR

# [1] import -----------------------------------------------------------
if [ $SKIP_IMPORT -eq 0 ] && [ $SKIP_CUTADAPT -eq 0 ]; then
  module purge
  module load anaconda/colsa
  source activate $QIIME2_ENV

  srun qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path $READS_PATH \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path ${OUTPUT_PREFIX}_demux-paired-reads.qza

  if [ $? -ne 0 ]; then echo "ERROR: QIIME2 import failed."; exit 1; fi
fi

# [2] primer trimming --------------------------------------------------
# Trims cpn60 UT primers (H279, H280, H1612, H1613) from 5' end of reads.
# Inosines (I) in primer sequences are represented as N.
# --p-discard-untrimmed removes read pairs where no primer is found,
# ensuring only genuine cpn60 amplicons are passed to DADA2.
#
# Primers (5'->3', I converted to N):
#   H279  (fwd): GANNNNGCNGGNGAYGGNACNACNAC
#   H1612 (fwd): GANNNNGCNGGYGACGGYACSACSAC
#   H280  (rev): YKNYKNTCNCCRAANCCNGGNGCYTT
#   H1613 (rev): CGRCGRTCRCCGAAGCCSGGNGCCTT
if [ $SKIP_CUTADAPT -eq 0 ]; then
  module purge
  module load anaconda/colsa
  source activate $QIIME2_ENV

  srun qiime cutadapt trim-paired \
    --i-demultiplexed-sequences ${OUTPUT_PREFIX}_demux-paired-reads.qza \
    --p-cores 24 \
    --p-front-f GANNNNGCNGGNGAYGGNACNACNAC \
    --p-front-f GANNNNGCNGGYGACGGYACSACSAC \
    --p-front-r YKNYKNTCNCCRAANCCNGGNGCYTT \
    --p-front-r CGRCGRTCRCCGAAGCCSGGNGCCTT \
    --p-discard-untrimmed \
    --p-match-read-wildcards \
    --p-error-rate 0.15 \
    --o-trimmed-sequences ${OUTPUT_PREFIX}_trimmed-reads.qza \
    --verbose

  if [ $? -ne 0 ]; then echo "ERROR: Primer trimming failed."; exit 1; fi

  srun qiime demux summarize \
    --i-data ${OUTPUT_PREFIX}_trimmed-reads.qza \
    --o-visualization ${OUTPUT_PREFIX}_trimmed-reads.qzv

  rm -rf ${OUTPUT_PREFIX}_Exported_Reads

  srun qiime tools extract \
    --input-path ${OUTPUT_PREFIX}_trimmed-reads.qza \
    --output-path ${OUTPUT_PREFIX}_Exported_Reads

  TRIMMED_PATH=$(find ${OUTPUT_PREFIX}_Exported_Reads -name "*.fastq.gz" -printf '%h\n' | sort -u | head -1)
  echo "Trimmed reads exported to: $TRIMMED_PATH"
  if [ -z "$TRIMMED_PATH" ]; then echo "ERROR: Could not find exported FASTQ files."; exit 1; fi

else
  if [ -z "$TRIMMED_QZA" ]; then
    echo "ERROR: SKIP_CUTADAPT=1 but TRIMMED_QZA is not set."
    echo "Set TRIMMED_QZA to your trimmed reads .qza (e.g. cpn60_run1_trimmed-reads.qza)"
    exit 1
  fi

  echo "Skipping cutadapt. Re-extracting FASTQs from: $TRIMMED_QZA"
  module purge
  module load anaconda/colsa
  source activate $QIIME2_ENV

  rm -rf ${OUTPUT_PREFIX}_Exported_Reads

  srun qiime tools extract \
    --input-path $TRIMMED_QZA \
    --output-path ${OUTPUT_PREFIX}_Exported_Reads

  TRIMMED_PATH=$(find ${OUTPUT_PREFIX}_Exported_Reads -name "*.fastq.gz" -printf '%h\n' | sort -u | head -1)
  echo "Trimmed reads exported to: $TRIMMED_PATH"
  if [ -z "$TRIMMED_PATH" ]; then echo "ERROR: Could not find exported FASTQ files."; exit 1; fi
fi

# [3] DADA2 ------------------------------------------------------------
if [ $SKIP_DADA2 -eq 0 ]; then
  module purge
  module load anaconda/colsa
  source activate $CONDA_ENV

  srun Rscript $SCRIPTS_DIR/cpn60_pipeline.R \
    --reads_path $TRIMMED_PATH \
    --output_prefix $OUTPUT_PREFIX \
    --error_model $ERROR_MODEL

  if [ $? -ne 0 ]; then echo "ERROR: DADA2 pipeline failed."; exit 1; fi
else
  echo "Skipping DADA2."
fi

# [4] classify ---------------------------------------------------------
if [ $SKIP_CLASSIFY -eq 0 ]; then
  module purge
  module load anaconda/colsa
  source activate $QIIME2_ENV

  srun Rscript $SCRIPTS_DIR/cpn60_classify.R \
    --output_prefix $OUTPUT_PREFIX \
    --classifier $CLASSIFIER

  if [ $? -ne 0 ]; then echo "ERROR: Classification failed."; exit 1; fi
else
  echo "Skipping classification."
fi

# [5] add taxonomy to phyloseq RDS ------------------------------------
# switches back to cpn60_pipeline R environment where phyloseq is installed
module purge
module load anaconda/colsa
source activate $CONDA_ENV

export OUTPUT_PREFIX
srun Rscript $SCRIPTS_DIR/cpn60_add_taxonomy.R \
  --output_prefix $OUTPUT_PREFIX

if [ $? -ne 0 ]; then echo "ERROR: Failed to add taxonomy to phyloseq object."; exit 1; fi

echo "Done."
