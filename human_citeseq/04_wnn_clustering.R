################################################################################
# CITE-seq: weighted-nearest-neighbour (RNA + protein) clustering
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
#
# Purpose
#   Builds a WNN graph from Harmony-corrected RNA PCs and dsb protein values,
#   clusters at several resolutions and computes WNN, RNA-only and protein-only UMAPs.
#
# Inputs
#   normalized_harmony.rds  (from 03_harmony_integration.Rmd)
#
# Outputs
#   clustered2.rds
#
# Run order
#   Upstream  : 03_harmony_integration.Rmd
#   Downstream: global annotation (not included) -> 05a
################################################################################

library(Seurat)
library(dplyr)

## ---- Paths (EDIT) ----
in_rds <- "path/to/normalized_harmony.rds"   # EDIT: from 03_harmony_integration.Rmd

sample <- readRDS(in_rds)

# Use dsb protein values directly as the protein 'reduction' (excluding
# isotype controls and failed antibodies, as in 02)
pseudo = t(sample@assays$CITE@data)[,-c(150,181,220:222,242:245)]
pseudo_colnames = paste('pseudo', c(1:149,151:180,182:219,223:241,246:279), sep = "_")
colnames(pseudo) = pseudo_colnames
# add to object 
sample@reductions$pdsb@cell.embeddings = pseudo

###### GEX Dimensionality
# Determine the dimensionality of the dataset
pct <- sample[["harmony_rna"]]@stdev / sum(sample[["harmony_rna"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
n_pcs_to_use_RNA <- which(cumu > 90 & pct < 5)[1]

sample = FindMultiModalNeighbors(
  object = sample,
  reduction.list = list("harmony_rna", "pdsb"),
  weighted.nn.name = "sct.dsb_wnn", 
  knn.graph.name = "sct.dsb_knn",
  modality.weight.name = "sct.dsb_weight",
  snn.graph.name = "sct.dsb_snn",
  dims.list = list(1:n_pcs_to_use_RNA, 1:270), 
  verbose = TRUE
)

sample <- FindClusters(sample, graph.name = "sct.dsb_knn", algorithm = 3, 
            resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), verbose = TRUE)
sample <- RunUMAP(sample, nn.name = "sct.dsb_wnn", reduction.name = "sct.dsb_wnn_umap", 
            reduction.key = "sct.dsb_wnnUMAP_", verbose = TRUE, return.model = TRUE)

# Run just RNA or ADT
sample <- RunUMAP(sample, reduction = 'harmony_rna', dims = 1:n_pcs_to_use_RNA, assay = 'SCT', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

sample <- RunUMAP(sample, reduction = 'pdsb', dims = 1:270, assay = 'CITE', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')


saveRDS(sample, "./clustered2.rds")