################################################################################
# Visium: load, normalize, merge and Harmony-integrate all sections
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
#
# Purpose
#   Batch version of 01_spatial_merge_normalize.Rmd: loads the 28 sections in
#   metadata/spatial_samples.csv, SCTransform-normalizes and merges them, then
#   integrates sections with Harmony and clusters the spots.
#
# Inputs
#   <spatial_dir>/*.h5seurat  (Visium sections from Kuppe et al., Nature 2022)
#   metadata/spatial_samples.csv
#
# Outputs
#   spatial_merged_normalized.rds, kuppe_spatial_integrated.rds
#
# Run order
#   Upstream  : none
#   Downstream: 03_spatial_fibroblast_niches.Rmd
################################################################################

library(Seurat)
library(dplyr)
library(harmony)
library(SeuratDisk)

## ---- Paths (EDIT) ----
repo_dir       <- ".."                               # EDIT: root of this repository
spatial_dir    <- "path/to/Spatial/converted_for_seurat"  # EDIT: h5seurat files
integrated_rds <- "kuppe_spatial_integrated.rds"     # EDIT: output file
source(file.path(repo_dir, "R", "utils.R"))

## ---- Load, normalize and merge all sections ----
sheet <- read_sample_sheet(file.path(repo_dir, "metadata", "spatial_samples.csv"))
sections <- load_visium_samples(sheet, spatial_dir)

sample <- merge(sections[[1]], y = sections[-1])
saveRDS(sample, "./spatial_merged_normalized.rds")

# Union of the per-section variable features, in section order
DefaultAssay(sample) <- "SCT"
VariableFeatures(sample) <- unlist(lapply(sections, VariableFeatures), use.names = FALSE)

## ---- Integrate sections (Harmony) and cluster ----
sample <- RunPCA(sample, npcs = 100, verbose = TRUE)
sample <- RunHarmony(sample, c("sample"), reduction = "pca", reduction.save = "harmony", assay.use = "SCT")
sample <- RunUMAP(sample, reduction = "harmony", dims = 1:50)
sample <- FindNeighbors(sample, reduction = "harmony", dims = 1:50)
sample <- FindClusters(sample, graph.name = "SCT_snn", algorithm = 3, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5), verbose = TRUE)

saveRDS(sample, integrated_rds)
