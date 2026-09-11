################################################################################
# Multiome RNA: load 23 samples, QC, Harmony integration and clustering
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
# Part  : Single-nucleus multiome (paired RNA + ATAC)
#
# Purpose
#   Loads the gene-expression counts of every multiome sample (only nuclei that
#   passed ATAC QC), merges them, filters on mitochondrial fraction and UMIs,
#   normalizes with SCTransform, integrates samples with Harmony, clusters at
#   several resolutions and writes cluster markers.
#
# Inputs
#   metadata/multiome_samples.csv         sample -> Cell Ranger ARC batch folder
#   <counts_dir>/<batch>/<sample>/outs/filtered_feature_bc_matrix/
#   <cells_dir>/<sample>_cells            barcodes passing ATAC QC (RDS)
#
# Outputs (in out_dir)
#   RNA_merged_postQC.rds, integrated_preClean.rds
#   DE_RNA_SCT_snn_res.{0.5,0.7,1.0}.csv
#
# Run order
#   Upstream  : Cell Ranger ARC; per-sample ATAC QC (ArrowFiles)
#   Downstream: cell-type annotation (myocardium_multiome.rds) -> 02_archr_project.R
################################################################################

library(Seurat)
library(dplyr)
library(harmony)

## ---- Paths (EDIT) ----
repo_dir   <- ".."                              # EDIT: root of this repository
counts_dir <- "path/to/Multiome/counts"         # EDIT: Cell Ranger ARC runs (Results/, batch3/, ...)
cells_dir  <- "path/to/samplePostQC"            # EDIT: <sample>_cells barcode lists
out_dir    <- "path/to/RNA/globalObjectConstruction_batch1_2_3"  # EDIT: output folder
source(file.path(repo_dir, "R", "utils.R"))

## ---- Load all samples and merge ----
samples <- read_sample_sheet(file.path(repo_dir, "metadata", "multiome_samples.csv"))
rna <- load_multiome_rna(samples, counts_dir, cells_dir, project = "HF")

# Cell barcodes are prefixed s1 ... s23 (sample order in the sheet)
merged <- merge(rna[[1]], y = rna[-1], add.cell.ids = paste0("s", seq_along(rna)), project = "HF")
rm(rna)

## ---- QC ----
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged <- subset(merged, subset = percent.mt < 5 & nCount_RNA < 15000)
saveRDS(merged, file.path(out_dir, "RNA_merged_postQC.rds"))

## ---- Normalize, integrate across samples, cluster ----
DefaultAssay(merged) <- 'RNA'
merged <- SCTransform(merged, vars.to.regress = c("percent.mt", "nCount_RNA"))
merged <- RunPCA(merged, npcs = 100, verbose = TRUE)
merged <- RunHarmony(merged, c("sample"), reduction = "pca", reduction.save = "harmony", assay.use = "SCT")
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:50)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:50)
merged <- FindClusters(merged, graph.name = "SCT_snn", algorithm = 3,
                       resolution = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0), verbose = TRUE)
saveRDS(merged, file.path(out_dir, "integrated_preClean.rds"))

## ---- Cluster markers at three resolutions ----
# names = resolution in the metadata column, values = label in the file name
write_cluster_markers(merged, c("0.5" = "0.5", "0.7" = "0.7", "1" = "1.0"), out_dir)
