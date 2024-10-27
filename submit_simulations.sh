#!/bin/bash -l

#SBATCH -p memory
#SBATCH -A naiss2023-5-517
#SBATCH --time=48:00:00
#SBATCH --mem=300GB

module load PDC/23.12
module load R/4.4.1-cpeGNU-23.12

Rscript test_run.R
