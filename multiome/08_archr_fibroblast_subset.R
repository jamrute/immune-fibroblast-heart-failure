################################################################################
# Multiome: fibroblast sub-project, joint RNA + ATAC clustering
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
# Part  : Single-nucleus multiome (paired RNA + ATAC)
#
# Purpose
#   (optional) Creates a fibroblast-only ArchR project, transfers fibroblast
#   states mapped from the CITE-seq reference, and recomputes the joint RNA +
#   ATAC embedding and clusters. Plots UMAPs and confusion matrices.
#
# Inputs  : ArchR project (Save-proj1); multiome fibroblasts mapped to the
#           CITE-seq fibroblast reference (multiome_fibroblast_mappedToCITEseq.rds)
# Outputs : ArchR project Save-projFibro; UMAP and confusion-matrix PDFs
#
# Run order
#   Upstream  : 05_archr_multiomic_clustering.R; human_citeseq/08_fibroblast_states.Rmd
#   Downstream: none
################################################################################

library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)
library(pheatmap)

## ---- Paths and options (EDIT) ----
analysis_dir <- "path/to/analysis_batch_1_2_3"   # EDIT: folder containing Save-proj1/
repo_dir     <- ".."                             # EDIT: root of this repository
fibro_rna_rds <- "path/to/RNA/fibroblast/multiome_fibroblast_mappedToCITEseq.rds"  # EDIT
build_fibro_project <- FALSE  # TRUE to (re)create Save-projFibro and its embeddings
setwd(analysis_dir)
source(file.path(repo_dir, "R", "utils.R"))

## ---- Fibroblast sub-project ----
if (build_fibro_project) {
  proj1 <- loadArchRProject(path = "Save-proj1")
  projFibro <- saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-projFibro", load = TRUE, dropCells = TRUE)
  idxSample <- BiocGenerics::which(projFibro$mappedCellType %in% "Fibroblast")
  projFibro <- projFibro[projFibro$cellNames[idxSample], ]

  # Fibroblast states mapped from the CITE-seq reference
  # NOTE: assumes nuclei are in the same order in scRNA and projFibro
  scRNA <- readRDS(fibro_rna_rds)
  projFibro$predictedcelltype <- as.character(scRNA$predicted.celltype)

  # Joint embedding, as in 05_archr_multiomic_clustering.R
  lsi_cluster_params <- list(resolution = 0.2, sampleCells = 10000, n.start = 10)
  projFibro <- addIterativeLSI(ArchRProj = projFibro, clusterParams = lsi_cluster_params, saveIterations = FALSE,
                               useMatrix = "GeneExpressionMatrix", depthCol = "Gex_nUMI", varFeatures = 2500,
                               firstSelection = "variable", binarize = FALSE, name = "LSI_RNA", force = TRUE)
  projFibro <- addIterativeLSI(ArchRProj = projFibro, clusterParams = lsi_cluster_params, saveIterations = FALSE,
                               useMatrix = "TileMatrix", depthCol = "nFrags", name = "LSI_ATAC", force = TRUE)
  projFibro <- addCombinedDims(projFibro, reducedDims = c("LSI_RNA", "LSI_ATAC"), name = "LSI_Combined")
  projFibro <- addHarmony(ArchRProj = projFibro, reducedDims = "LSI_Combined", name = "HAR_Combined",
                          groupBy = "Sample", force = TRUE)

  projFibro <- addUMAP(projFibro, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
  projFibro <- addUMAP(projFibro, reducedDims = "LSI_RNA", name = "UMAP_LSI_RNA", minDist = 0.8, force = TRUE)
  projFibro <- addUMAP(projFibro, reducedDims = "HAR_Combined", name = "UMAP_HAR_Combined", minDist = 0.8, force = TRUE)

  projFibro <- addClusters(projFibro, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.1, force = TRUE)
  projFibro <- addClusters(projFibro, reducedDims = "LSI_RNA", name = "Clusters_LSI_RNA", resolution = 0.1, force = TRUE)
  projFibro <- addClusters(projFibro, reducedDims = "HAR_Combined", name = "Clusters_HAR_Combined", resolution = 0.1, force = TRUE)

  saveArchRProject(ArchRProj = projFibro, outputDirectory = "Save-projFibro", load = TRUE)
}

projFibro <- loadArchRProject(path = "Save-projFibro")

## ---- UMAPs of each embedding ----
p1 <- plotEmbedding(projFibro, name = "Clusters_ATAC", embedding = "UMAP_ATAC", size = 1, labelAsFactors = F, labelMeans = F)
p2 <- plotEmbedding(projFibro, name = "Clusters_LSI_RNA", embedding = "UMAP_LSI_RNA", size = 1, labelAsFactors = F, labelMeans = F)
p3 <- plotEmbedding(projFibro, name = "Clusters_HAR_Combined", embedding = "UMAP_HAR_Combined", size = 1, labelAsFactors = F, labelMeans = F)
p4 <- plotEmbedding(projFibro, name = "predictedcelltype", embedding = "UMAP_HAR_Combined", size = 1, labelAsFactors = F, labelMeans = F)

p <- lapply(list(p1, p2, p3), function(x) {
  x + guides(color = "none", fill = "none") +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
})
do.call(cowplot::plot_grid, c(list(ncol = 3), p))
plotPDF(p1, p2, p3, p4, name = "UMAP-scATAC-scRNA-Combined.pdf", addDOC = FALSE, ArchRProj = projFibro)

## ---- Agreement between clusterings ----
plotPDF(plot_confusion_heatmap(projFibro$Clusters_ATAC, projFibro$Clusters_LSI_RNA),
        name = "Clusters_ATAC_Clusters_LSI_RNA_confusionMatrix.pdf", addDOC = FALSE, ArchRProj = projFibro)
plotPDF(plot_confusion_heatmap(projFibro$predictedcelltype, projFibro$Clusters_HAR_Combined),
        name = "predictedcelltype_HAR_Combined_confusionMatrix.pdf", addDOC = FALSE, ArchRProj = projFibro)

saveArchRProject(ArchRProj = projFibro, outputDirectory = "Save-projFibro", load = TRUE)
