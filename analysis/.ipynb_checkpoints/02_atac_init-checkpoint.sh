#!/bin/bash

#SBATCH -n 30
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --job-name=archr_init
#SBATCH --output=02_atac_init_%j.out
#SBATCH --error=02_atac_init_%j.err

source ~/.bashrc

conda activate multiome_winner

cd /data/peer/lpaddock/repos/multiome_QC_preprocessing/analysis/

# execute program
Rscript 02_atac_init.r /data/peer/lpaddock/repos/multiome_QC_preprocessing/data /data/peer/lpaddock/repos/multiome_QC_preprocessing/results

conda deactivate