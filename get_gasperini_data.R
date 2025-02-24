### download gasperini data

library(tidyverse)
library(SingleCellExperiment)

setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd("../")

# Gasperini data
raw_data_dir_gasp = "./data/"
# URL of data
remote <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120861&format=file&file="
# Gasperini et al results
all_deg_results_filename <- "GSE120861_all_deg_results.at_scale.txt"
# names of genes
genes_filename <- "GSE120861_at_scale_screen.genes.txt"
# names of cells
cells_filename <- "GSE120861_at_scale_screen.cells.txt"
# "reference cells" used by Gasperini et al for computational purposes
reference_cells_filename <- "GSE120861_50k_reference_cells.rds"
# all (gRNA, gene) pairs
gRNAgroup_pair_table_filename <- "GSE120861_gene_gRNAgroup_pair_table.at_scale.txt"
# list of gRNA groups used
gRNA_groups_filename <- "GSE120861_grna_groups.at_scale.txt"
# Monocle Cell Data Set object with all gRNA data
cds_filename <- "GSE120861_at_scale_screen.cds.rds"
# Expression data
expression_filename <- "GSE120861_at_scale_screen.exprs.mtx"
# list of files to download
filenames <- c(gRNAgroup_pair_table_filename,
               gRNA_groups_filename,
               all_deg_results_filename,
               genes_filename,
               cells_filename,
               reference_cells_filename,
               cds_filename
               # expression_filename, 
               )
options(timeout = 600)
# download files if not already present
for (filename in filenames) {
  if (!file.exists(paste0(raw_data_dir_gasp, "/", filename))) {
    cat(paste0("Downloading ", filename, "\n"))
    source <- paste0(remote, filename, ".gz")
    dest <- paste0(raw_data_dir_gasp, "/", filename, ".gz")
    download.file(source, dest)
    R.utils::gunzip(paste0(dest))
  }
}

### 
library(monocle)
cds_obj <- readRDS("./data/GSE120861_at_scale_screen.cds.rds")
gene_gRNAgroup_pair_table <- as.data.frame(readr::read_delim("./data/GSE120861_gene_gRNAgroup_pair_table.at_scale.txt", delim='\t'))
gene_gRNAgroup_pair_table_filtered <- gene_gRNAgroup_pair_table %>%
  dplyr::filter(gRNAgroup.chr != "NTC")

## get gRNA annotations:
grna_file <- read_tsv("./data/GSE120861_grna_groups.at_scale.txt", col_names = c(""))

# get all genes per perturbation: 
cre_target_genes <- gene_gRNAgroup_pair_table_filtered %>% 
  group_by(gRNAgroup) %>%
  summarize(target_genes = list(ENSG.targetgene), chr = unique(gRNAgroup.chr), start = unique(gRNAgroup.start), end = unique(gRNAgroup.stop)) %>%
  mutate(strand = "*") %>%
  mutate(cre_target = gRNAgroup) %>%
  column_to_rownames("gRNAgroup")

# here we parse annotations between target regions and guides, i think this is the level we want to do it on
# not sure where to get results per individual guide from
# 
# this takes a while... 
# guide_cell_annotations <- lapply( rownames(cre_target_genes), function(guide_here) {
#   which(grepl(paste0('(^|_)', guide_here, '(_|$)'), cds_obj@phenoData@data$gene ))
# })
# saveRDS(guide_cell_annotations, "~/Desktop/guide_cell_annotations.rds")

guide_cell_annotations <- readRDS("~/Desktop/dosage_paper/guide_cell_annotations.rds")

library(Matrix)

n_genes <- as.integer(length(grna_file$...1))
n_cells <- as.integer(ncol(cds_obj)[[1]])

# make sparse matrix to avoid massive matrix filled with 0s
is <- unlist(lapply(1:length(guide_cell_annotations), function(i){rep(i, length(guide_cell_annotations[[i]]))}))
js <- unlist(guide_cell_annotations)
guide_indicator_matrix <- sparseMatrix(is, js, x = rep(1, length(n_genes)))

# assemble SCE
# for now we leave out the guide-level information - not necessary i think? 
perturb_annotation <- SingleCellExperiment(assays = list("counts" = guide_indicator_matrix, "perts" = guide_indicator_matrix), 
                                           rowRanges = GenomicRanges::makeGRangesFromDataFrame(cre_target_genes, keep.extra.columns = T))
colnames(perturb_annotation) <- colnames(cds_obj)

sce_gasperini <- SingleCellExperiment(assays = list("counts" = cds_obj@assayData$exprs),
                                    colData = pData(cds_obj)[ , 1:18],
                                    rowData = fData(cds_obj),
                                    altExps = list("cre_pert" = perturb_annotation)
                                    )

perturbations_test <- c("ACTB_TSS", "ACTG1_TSS", "ACYP1_TSS")
cells_use <- colSums(counts(altExps(sce_gasperini)[["cre_pert"]][perturbations_test, ])) > 0
sce_gasperini_subset <- sce_gasperini[ , cells_use]
altExps(sce_gasperini_subset)[["cre_pert"]] <- altExps(sce_gasperini_subset)[["cre_pert"]][perturbations_test, ]

saveRDS(sce_gasperini_subset, "./sce_gasperini_subset.rds")

#perturbations_sam <- read.csv("~/Desktop/dosage_paper/gasperini_target_sites_Jasper_long.txt", header = F)[[1]]
perturbations_sam <- read.csv("~/Desktop/gasperini_target_sites_Jasper.txt", header = F)[[1]]
perturbations_sam_possible <- c(paste0(perturbations_sam, "_top_two"), paste0(perturbations_sam, "_second_two"))
perturbations_sam_existing <- perturbations_sam_possible[perturbations_sam_possible %in% rownames(altExps(sce_gasperini)[["cre_pert"]])]
cells_use <- colSums(counts(altExps(sce_gasperini)[["cre_pert"]][perturbations_sam_existing, ])) > 0
sce_gasperini_subset <- sce_gasperini[ , cells_use]
altExps(sce_gasperini_subset)[["cre_pert"]] <- altExps(sce_gasperini_subset)[["cre_pert"]][perturbations_sam_existing, ]

saveRDS(sce_gasperini_subset, "./sce_gasperini_sam.rds")

saveRDS(sce_gasperini, "./sce_gasperini.rds")
