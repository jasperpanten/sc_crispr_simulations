#!/bin/python -l

# run as 
# snakemake -s snakefile --use-envmodules --executor slurm --default-resources slurm_account=naiss2024-5-581 slurm_partition=shared --jobs 1

import pandas as pd
import pdb

effect_sizes = [0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5]
reps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

#cres_to_run = pd.read_csv("../processed_data/gasperini_target_sites_Jasper.txt")["x"]
#cres_to_run = cres_to_run.to_numpy()
#cres_to_run = cres_to_run[[1]]
cres_to_run = "all"

rule all:
    input:
        expand("../results/gasperini_results/sim_sceptre_res_{cre}_{effect_size}_{rep}.rds", cre = cres_to_run, effect_size = effect_sizes, rep = reps)

rule run_simulation:
    output: "../results/gasperini_results/sim_sceptre_res_{cre}_{effect_size}_{rep}.rds",
    envmodules:
        "PDC/23.12",
        "R/4.4.1-cpeGNU-23.12",
    resources:
        mem_mb=80000,
        runtime=500,
        #cpus_per_task=8,
        slurm_partition="shared",
    shell: "Rscript test_run_5.R {wildcards.effect_size} 1 {wildcards.cre} {wildcards.rep}"
