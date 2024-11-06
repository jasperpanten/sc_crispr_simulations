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

# GDOs: 
system(
  paste0("ln -sf ", data_dir, "GSE171452_RAW/GSM5225859_STINGseq-v1_GDO.barcodes.tsv", " ", data_dir, "/morris_1_gdo/barcodes.tsv")
)
system(
  paste0("ln -sf ", data_dir, "GSE171452_RAW/GSM5225859_STINGseq-v1_GDO.features.tsv", " ", data_dir, "/morris_1_gdo/features.tsv")
)
system(
  paste0("ln -sf ", data_dir, "GSE171452_RAW/GSM5225859_STINGseq-v1_GDO.matrix.mtx", " ", data_dir, "/morris_1_gdo/matrix.mtx")
)

morris_dataset_1_gdo <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_1_gdo/")
morris_dataset_1_gdo <- morris_dataset_1_gdo[grepl("SNP", rownames(morris_dataset_1_gdo)), ]

colnames(morris_dataset_1) <- paste0("SmallScreen_", morris_dataset_1$Barcode)
colnames(morris_dataset_1_gdo) <- paste0("SmallScreen_", morris_dataset_1_gdo$Barcode)

joint_cells <- intersect(colnames(morris_dataset_1), colnames(morris_dataset_1_gdo))

morris_dataset_1 <- morris_dataset_1[ , joint_cells]
morris_dataset_1_gdo <- morris_dataset_1_gdo[ , joint_cells]

altExps(morris_dataset_1)[["cre_pert"]] <- morris_dataset_1_gdo

saveRDS(morris_dataset_1, "../data/morris_raw_small_screen.rds")

### dataset 2
letter_to_gsm <- list("A" = "GSM7108117", "B" = "GSM7108118", "C" = "GSM7108119", "D" = "GSM7108120")

for (x in c("A", "B", "C", "D")){
  
  system(
    paste0("mkdir -p ", data_dir, paste0("/morris_2-", x, "/"))
  )
  
  system(
    paste0("ln -sf ", data_dir, "GSE171452_RAW/", letter_to_gsm[x], "_STINGseq-v2_cDNA-", x, "_barcodes.mtx", " ", data_dir, "/morris_2-", x, "/barcodes.tsv")
  )
  system(
    paste0("ln -sf ", data_dir, "GSE171452_RAW/", letter_to_gsm[x], "_STINGseq-v2_cDNA-", x, "_features.mtx", " ", data_dir, "/morris_2-", x, "/features.tsv")
  )
  system(
    paste0("ln -sf ", data_dir, "GSE171452_RAW/", letter_to_gsm[x], "_STINGseq-v2_cDNA-", x, "_matrix.mtx", " ", data_dir, "/morris_2-", x, "/matrix.mtx")
  )
}

morris_dataset_2_A <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_2-A/")
morris_dataset_2_B <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_2-B/")
morris_dataset_2_C <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_2-C/")
morris_dataset_2_D <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_2-D/")

morris_dataset_2 <- cbind(morris_dataset_2_A, morris_dataset_2_B, morris_dataset_2_C, morris_dataset_2_D)

# GDOs
letter_to_gsm <- list("A" = "GSM7108121", "B" = "GSM7108122", "C" = "GSM7108123", "D" = "GSM7108124")

for (x in c("A", "B", "C", "D")){
  
  system(
    paste0("mkdir -p ", data_dir, paste0("/morris_2_gdo-", x, "/"))
  )
  
  system(
    paste0("ln -sf ", data_dir, "GSE171452_RAW/", letter_to_gsm[x], "_STINGseq-v2_GDO-", x, "_barcodes.mtx", " ", data_dir, "/morris_2_gdo-", x, "/barcodes.tsv")
  )
  system(
    paste0("ln -sf ", data_dir, "GSE171452_RAW/", letter_to_gsm[x], "_STINGseq-v2_GDO-", x, "_features.mtx", " ", data_dir, "/morris_2_gdo-", x, "/features.tsv")
  )
  system(
    paste0("ln -sf ", data_dir, "GSE171452_RAW/", letter_to_gsm[x], "_STINGseq-v2_GDO-", x, "_matrix.mtx", " ", data_dir, "/morris_2_gdo-", x, "/matrix.mtx")
  )
}

morris_dataset_2_gdo_A <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_2_gdo-A/")
morris_dataset_2_gdo_B <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_2_gdo-B/")
morris_dataset_2_gdo_C <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_2_gdo-C/")
morris_dataset_2_gdo_D <- read10xCounts("~/Desktop/PostDoc_TL_Lab/Projects/Sam/data/morris_2_gdo-D/")

morris_dataset_2_gdo <- cbind(morris_dataset_2_gdo_A, morris_dataset_2_gdo_B, morris_dataset_2_gdo_C, morris_dataset_2_gdo_D)
morris_dataset_2_gdo <- morris_dataset_2_gdo[grepl("SNP", rownames(morris_dataset_2_gdo)), ]

colnames(morris_dataset_2) <- paste0(basename(morris_dataset_2$Sample), morris_dataset_2$Barcode)
colnames(morris_dataset_2_gdo) <- paste0(basename(morris_dataset_2_gdo$Sample), morris_dataset_2_gdo$Barcode)

joint_cells <- intersect(colnames(morris_dataset_2), colnames(morris_dataset_2_gdo))

morris_dataset_2 <- morris_dataset_2[ , joint_cells]
morris_dataset_2_gdo <- morris_dataset_2_gdo[ , joint_cells]

altExps(morris_dataset_2)[["cre_pert"]] <- morris_dataset_2_gdo

saveRDS(morris_dataset_2, "../data/morris_raw_large_screen.rds")
