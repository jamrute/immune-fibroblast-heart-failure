################################################################################
# Multiome: joint RNA + ATAC embedding and clustering
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
# Part  : Single-nucleus multiome (paired RNA + ATAC)
#
# Purpose
#   (optional) Computes iterative LSI on gene expression and on the tile matrix,
#   combines them, corrects with Harmony and clusters each embedding; plots the
#   UMAPs and confusion matrices between ATAC clusters, joint clusters and RNA
#   cell types.
#
# Inputs  : ArchR project (Save-proj1)
# Outputs : UMAP-scATAC-scRNA-Combined.pdf and confusion-matrix PDFs
#
# Run order
#   Upstream  : 02_archr_project.R
#   Downstream: 06_archr_peak2gene.R (uses the HAR_Combined embedding)
################################################################################

library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)
library(pheatmap)

## ---- Paths and options (EDIT) ----
analysis_dir <- "path/to/analysis_batch_1_2_3"   # EDIT: folder containing Save-proj1/
repo_dir     <- ".."                             # EDIT: root of this repository
run_multiomic_clustering <- FALSE  # TRUE to (re)compute the embeddings and clusters below
setwd(analysis_dir)
source(file.path(repo_dir, "R", "utils.R"))

proj1 <- loadArchRProject(path = "Save-proj1")

## ---- Joint embedding: LSI on RNA and ATAC, combined and Harmony-corrected ----
if (run_multiomic_clustering) {
  lsi_cluster_params <- list(resolution = 0.2, sampleCells = 10000, n.start = 10)
  proj1 <- addIterativeLSI(ArchRProj = proj1, clusterParams = lsi_cluster_params, saveIterations = FALSE,
                           useMatrix = "GeneExpressionMatrix", depthCol = "Gex_nUMI", varFeatures = 2500,
                           firstSelection = "variable", binarize = FALSE, name = "LSI_RNA", force = TRUE)
  proj1 <- addIterativeLSI(ArchRProj = proj1, clusterParams = lsi_cluster_params, saveIterations = FALSE,
                           useMatrix = "TileMatrix", depthCol = "nFrags", name = "LSI_ATAC", force = TRUE)
  proj1 <- addCombinedDims(proj1, reducedDims = c("LSI_RNA", "LSI_ATAC"), name = "LSI_Combined")
  proj1 <- addHarmony(ArchRProj = proj1, reducedDims = "LSI_Combined", name = "HAR_Combined",
                      groupBy = "Sample", force = TRUE)

  proj1 <- addUMAP(proj1, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
  proj1 <- addUMAP(proj1, reducedDims = "LSI_RNA", name = "UMAP_LSI_RNA", minDist = 0.8, force = TRUE)
  proj1 <- addUMAP(proj1, reducedDims = "HAR_Combined", name = "UMAP_HAR_Combined", minDist = 0.8, force = TRUE)

  proj1 <- addClusters(proj1, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.1, force = TRUE)
  proj1 <- addClusters(proj1, reducedDims = "LSI_RNA", name = "Clusters_LSI_RNA", resolution = 0.1, force = TRUE)
  proj1 <- addClusters(proj1, reducedDims = "HAR_Combined", name = "Clusters_HAR_Combined", resolution = 0.1, force = TRUE)

  saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)
}

## ---- UMAPs of each embedding ----
p1 <- plotEmbedding(proj1, name = "Clusters_ATAC", embedding = "UMAP_ATAC", size = 1, labelAsFactors = F, labelMeans = F)
p2 <- plotEmbedding(proj1, name = "Clusters_LSI_RNA", embedding = "UMAP_LSI_RNA", size = 1, labelAsFactors = F, labelMeans = F)
p3 <- plotEmbedding(proj1, name = "Clusters_HAR_Combined", embedding = "UMAP_HAR_Combined", size = 1, labelAsFactors = F, labelMeans = F)
p4 <- plotEmbedding(proj1, name = "mappedCellType", embedding = "UMAP_HAR_Combined", size = 1, labelAsFactors = F, labelMeans = F)

p <- lapply(list(p1, p2, p3), function(x) {
  x + guides(color = "none", fill = "none") +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
})
do.call(cowplot::plot_grid, c(list(ncol = 3), p))
plotPDF(p1, p2, p3, p4, name = "UMAP-scATAC-scRNA-Combined.pdf", addDOC = FALSE, ArchRProj = proj1)

## ---- Agreement between clusterings (row-normalized confusion matrices) ----
plotPDF(plot_confusion_heatmap(proj1$Clusters_ATAC, proj1$mappedCellType),
        name = "Clusters_ATAC_Clusters_mappedCellType_confusionMatrix.pdf", addDOC = FALSE, ArchRProj = proj1)
plotPDF(plot_confusion_heatmap(proj1$mappedCellType, proj1$Clusters_HAR_Combined),
        name = "mappedCellType_HAR_Combined_confusionMatrix.pdf", addDOC = FALSE, ArchRProj = proj1)
