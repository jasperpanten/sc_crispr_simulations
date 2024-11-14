#!/bin/python -l

# run as 
# snakemake -s snakefile --use-envmodules --executor slurm --default-resources slurm_account=naiss2023-5-517 slurm_partition=shared --jobs 1

import pandas as pd
import pdb

effect_sizes = [0.5]
reps = [1, 2]

cres_to_run = pd.read_csv("../processed_data/cre_names_gasperini.txt")["x"]
cres_to_run = cres_to_run.to_numpy()
cres_to_run = cres_to_run[[1]]

rule all:
    input:
        expand("../results/gasperini_results/sim_res_{cre}_{effect_size}_{rep}.rds", cre = cres_to_run, effect_size = effect_sizes, rep = reps)

rule run_simulation:
    output: "../results/gasperini_results/sim_res_{cre}_{effect_size}_{rep}.rds",
    envmodules:
        "PDC/23.12",
        "R/4.4.1-cpeGNU-23.12",
    resources:
        #mem_mb=100000,
        #disk_mb=100000,
        slurm_partition="shared",
    shell: "Rscript test_run_3.R {wildcards.effect_size} 20 {wildcards.cre} {wildcards.rep}"
