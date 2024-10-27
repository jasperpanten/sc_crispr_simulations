library(tidyverse)
library(SingleCellExperiment)
library(DropletUtils)


setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path)))

## create input so that read10xCounts can understand it
data_dir <- "~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/"

system(
  paste0("mkdir -p ", data_dir, "/morris_1/")
)

system(
  paste0("ln -sf ", data_dir, "GSE171452_RAW/GSM5225857_STINGseq-v1_cDNA.barcodes.tsv", " ", data_dir, "/morris_1/barcodes.tsv")
)
system(
  paste0("ln -sf ", data_dir, "GSE171452_RAW/GSM5225857_STINGseq-v1_cDNA.features.tsv", " ", data_dir, "/morris_1/features.tsv")
)
system(
  paste0("ln -sf ", data_dir, "GSE171452_RAW/GSM5225857_STINGseq-v1_cDNA.matrix.mtx", " ", data_dir, "/morris_1/matrix.mtx")
)

morris_dataset_1 <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_1/")


### dataset 2
for (x in c("A", "B", "C", "D")){
  
  system(
    paste0("mkdir -p ", data_dir, paste0("/morris_2-", x, "/"))
  )
  
  system(
    paste0("ln -sf ", data_dir, "GSE171452_RAW/GSM7108117_STINGseq-v2_cDNA-", x, "_barcodes.mtx", " ", data_dir, "/morris_2-", x, "/barcodes.tsv")
  )
  system(
    paste0("ln -sf ", data_dir, "GSE171452_RAW/GSM7108117_STINGseq-v2_cDNA-", x, "_features.mtx", " ", data_dir, "/morris_2-", x, "/features.tsv")
  )
  system(
    paste0("ln -sf ", data_dir, "GSE171452_RAW/GSM7108117_STINGseq-v2_cDNA-", x, "_matrix.mtx", " ", data_dir, "/morris_2-", x, "/matrix.mtx")
  )
}

morris_dataset_2 <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_2-A/")
