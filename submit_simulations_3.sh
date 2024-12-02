#!/bin/bash -l

#SBATCH -p shared
#SBATCH -A naiss2023-5-517
#SBATCH --time=1:00:00
#SBATCH --mem=30GB
#SBATCH -c 1

module load PDC/23.12
module load R/4.4.1-cpeGNU-23.12

Rscript test_run_3.R 0.5 1 all 2
