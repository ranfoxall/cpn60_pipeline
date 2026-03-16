#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=512000
#SBATCH -p shared
#SBATCH --job-name="cpn60_classify"
#SBATCH --output=cpn60_classify_%j.log
#SBATCH --error=cpn60_classify_%j.err
#SBATCH --exclude=node117,node118

# [variables] ----------------------------------------------------------
# edit these before submitting
OUTPUT_PREFIX=cpn60_run1
CLASSIFIER=/mnt/home/whistler/shared/cpn60-Classifier/cpn60-q2-feature-classifier-v11/cpn60_classifier_v11.qza
WORKDIR=/mnt/home/whistler/foxall/Scripts/CREATE/cpn60
SCRIPTS_DIR=/mnt/home/whistler/foxall/Scripts/CREATE/cpn60_pipeline/scripts
CONDA_ENV=cpn60_pipeline
QIIME2_ENV=qiime2-2020.2

cd $WORKDIR

# [1] classify ---------------------------------------------------------
module purge
module load anaconda/colsa
source activate $QIIME2_ENV

srun Rscript $SCRIPTS_DIR/cpn60_classify.R \
  --output_prefix $OUTPUT_PREFIX \
  --classifier $CLASSIFIER

if [ $? -ne 0 ]; then
  echo "ERROR: Classification failed."
  exit 1
fi

# [2] add taxonomy to phyloseq RDS ------------------------------------
# switches back to cpn60_pipeline R environment where phyloseq is installed
module purge
module load anaconda/colsa
source activate $CONDA_ENV

srun Rscript $SCRIPTS_DIR/cpn60_add_taxonomy.R \
  --output_prefix $OUTPUT_PREFIX

if [ $? -ne 0 ]; then echo "ERROR: Failed to add taxonomy to phyloseq object."; exit 1; fi

echo "Done."
