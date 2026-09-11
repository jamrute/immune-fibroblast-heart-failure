################################################################################
# Multiome ATAC: peak calling, cell-type marker peaks and ATAC-only embeddings
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
# Part  : Single-nucleus multiome (paired RNA + ATAC)
#
# Purpose
#   (optional) Calls reproducible peaks per cell type with MACS2 and adds the
#   PeakMatrix; finds marker peaks per cell type (heatmap); (optional) runs
#   ATAC-only iterative LSI, clustering, Harmony and UMAP with confusion
#   matrices against the RNA cell types.
#
# Inputs  : ArchR project (Save-proj1) from 02_archr_project.R
# Outputs : markerPeak_mappedCellType_final.rds, marker-peak heatmap and
#           UMAP / confusion-matrix PDFs (ArchR Plots/ folder)
#
# Run order
#   Upstream  : 02_archr_project.R
#   Downstream: 04-08
################################################################################

library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

## ---- Paths and options (EDIT) ----
analysis_dir <- "path/to/analysis_batch_1_2_3"   # EDIT: folder containing Save-proj1/
repo_dir     <- ".."                             # EDIT: root of this repository
run_peak_calling    <- FALSE  # TRUE to (re)call peaks with MACS2 and rebuild the PeakMatrix
run_atac_clustering <- FALSE  # TRUE to (re)compute ATAC-only LSI, clusters and UMAPs
setwd(analysis_dir)
source(file.path(repo_dir, "R", "utils.R"))

proj1 <- loadArchRProject(path = "Save-proj1")

## ---- 1. Peak calling (MACS2, pseudobulk per cell type) ----
if (run_peak_calling) {
  proj1 <- addGroupCoverages(ArchRProj = proj1, groupBy = "mappedCellType")
  pathToMacs2 <- findMacs2()
  proj1 <- addReproduciblePeakSet(
    ArchRProj = proj1,
    groupBy = "mappedCellType",
    pathToMacs2 = pathToMacs2,
    maxPeaks = 10000000
  )
  proj1 <- addPeakMatrix(proj1)
  saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)
}

## ---- 2. Marker peaks per cell type ----
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj1,
  useMatrix = "PeakMatrix",
  useGroups = c("Cardiomyocyte", "Myeloid", "Pericyte", "Fibroblast", "Endothelium"),
  groupBy = "mappedCellType",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 5000
)
saveRDS(markersPeaks, file = "markerPeak_mappedCellType_final.rds")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker_noMax-Heatmap_mappedCellType", width = 8, height = 6,
        ArchRProj = proj1, addDOC = FALSE)

## ---- 3. ATAC-only dimensionality reduction and clustering ----
if (run_atac_clustering) {
  proj1 <- addIterativeLSI(
    ArchRProj = proj1,
    useMatrix = "PeakMatrix",
    name = "LSI_ATAC",
    iterations = 2,
    clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10),  # see Seurat::FindClusters
    varFeatures = 25000,
    dimsToUse = 1:50,
    seed = 1, force = TRUE
  )
  proj1 <- addClusters(input = proj1, reducedDims = "LSI_ATAC", method = "Seurat",
                       name = "Clusters_LSI", resolution = 0.5, force = TRUE, seed = 1)
  proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "LSI_ATAC", name = "UMAP_ATAC",
                   nNeighbors = 50, minDist = 0.5, metric = "cosine", force = TRUE)

  # Same after Harmony batch correction across samples
  proj1 <- addHarmony(ArchRProj = proj1, reducedDims = "LSI_ATAC", name = "Harmony",
                      groupBy = "Sample", force = TRUE)
  proj1 <- addClusters(input = proj1, reducedDims = "Harmony", method = "Seurat",
                       name = "Clusters_HAR", resolution = 0.5, force = TRUE, seed = 1)
  proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "Harmony", name = "UMAP_ATAC_HARMONY",
                   nNeighbors = 50, minDist = 0.5, metric = "cosine", force = TRUE)

  # UMAPs colored by sample, ATAC cluster and RNA cell type
  embeddings <- list(LSI = c(umap = "UMAP_ATAC", clusters = "Clusters_LSI"),
                     HAR = c(umap = "UMAP_ATAC_HARMONY", clusters = "Clusters_HAR"))
  for (tag in names(embeddings)) {
    emb <- embeddings[[tag]]
    p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = emb[["umap"]])
    p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = emb[["clusters"]], embedding = emb[["umap"]])
    p3 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "mappedCellType", embedding = emb[["umap"]])
    plotPDF(p1, p2, p3, name = paste0("Plot-", emb[["umap"]], "-Sample-", emb[["clusters"]], "-mappedCellType-", tag, ".pdf"),
            ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)
  }

  # Confusion matrices: ATAC clusters vs sample and vs RNA cell type
  for (clusters in c("Clusters_LSI", "Clusters_HAR")) {
    plotPDF(plot_confusion_heatmap(getCellColData(proj1, select = clusters, drop = TRUE), proj1$Sample),
            name = paste0("confusionMap_heatmap_", clusters, "_Sample.pdf"), ArchRProj = proj1, addDOC = FALSE)
    plotPDF(plot_confusion_heatmap(getCellColData(proj1, select = clusters, drop = TRUE), proj1$mappedCellType),
            name = paste0("confusionMap_heatmap_", clusters, "_celltypeRNA.pdf"), ArchRProj = proj1, addDOC = FALSE)
  }

  saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = FALSE)
}
