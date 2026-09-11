library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/", sep=""))

scRNA <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/integrated/myocardium_multiome.rds")
proj1 <- loadArchRProject(path = "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-proj1")


markersPeaks <- getMarkerFeatures(
    ArchRProj = proj1, 
    useMatrix = "PeakMatrix", 
    useGroups = c("Cardiomyocyte","Myeloid","Pericyte","Fibroblast","Endothelium"),
    groupBy = "mappedCellType",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    maxCells = 5000
)

saveRDS(markersPeaks, file="markerPeak_mappedCellType_final.rds")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

#Plot heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker_noMax-Heatmap_mappedCellType", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)


################### Peak Calling with MACS2

# #Pseudobulk ATAC and Call Peaks
# proj1 <- addGroupCoverages(ArchRProj = proj1, groupBy = "mappedCellType")
# pathToMacs2 <- findMacs2()
# proj1 <- addReproduciblePeakSet(
#     ArchRProj = proj1, 
#     groupBy = "mappedCellType", 
#     pathToMacs2 = pathToMacs2,
#     maxPeaks = 10000000
# )

# proj1 <- addPeakMatrix(proj1)

# saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)

# getAvailableMatrices(proj1)

# markersPeaks <- getMarkerFeatures(
#     ArchRProj = proj1, 
#     useMatrix = "PeakMatrix", 
#     useGroups = c("Fibroblast","Cardiomyocyte","Pericyte","Endothelium","Myeloid"),
#     groupBy = "mappedCellType",
#     bias = c("TSSEnrichment", "log10(nFrags)"),
#     testMethod = "wilcoxon",
#     maxCells = 5000
# )

# saveRDS(markersPeaks, file="markerPeak_mappedCellType_final.rds")

# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

# #Plot heatmap
# heatmapPeaks <- plotMarkerHeatmap(
#   seMarker = markersPeaks, 
#   cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
#   transpose = TRUE
# )
# draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# plotPDF(heatmapPeaks, name = "Peak-Marker_noMax-Heatmap_mappedCellType", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)

# ########## Dimensional Reduction and clustering

# proj1 <- addIterativeLSI(
#     ArchRProj = proj1,
#     useMatrix = "PeakMatrix", 
#     name = "LSI_ATAC", 
#     iterations = 2, 
#     clusterParams = list( #See Seurat::FindClusters
#         resolution = c(0.2), 
#         sampleCells = 10000, 
#         n.start = 10
#     ), 
#     varFeatures = 25000, 
#     dimsToUse = 1:50,
#     seed=1,force=T
# )

# # basic clustering 
# proj1 <- addClusters(
#     input = proj1,
#     reducedDims = "LSI_ATAC",
#     method = "Seurat",
#     name = "Clusters_LSI",
#     resolution = 0.5,
#     force=T,seed=1
# )

# # UMAP embedding
# proj1 <- addUMAP(
#     ArchRProj = proj1, 
#     reducedDims = "LSI_ATAC", 
#     name = "UMAP_ATAC", 
#     nNeighbors = 50, 
#     minDist = 0.5, 
#     metric = "cosine",force=TRUE
# )

# p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAP_ATAC")
# p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters_LSI", embedding = "UMAP_ATAC")
# p3 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "mappedCellType", embedding = "UMAP_ATAC")
# plotPDF(p1,p2,p3, name = "Plot-UMAP_ATAC-Sample-Clusters_LSI-mappedCellType-LSI.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

# ############ Repeat clustering and dim reduction with harmony batch correction
# proj1 <- addHarmony(
#     ArchRProj = proj1,
#     reducedDims = "LSI_ATAC",
#     name = "Harmony",
#     groupBy = "Sample",force=TRUE
# )
# proj1 <- addClusters(
#     input = proj1,
#     reducedDims = "Harmony",
#     method = "Seurat",
#     name = "Clusters_HAR",
#     resolution = 0.5,
#     force=T,seed=1
# )

# proj1 <- addUMAP(
#     ArchRProj = proj1, 
#     reducedDims = "Harmony", 
#     name = "UMAP_ATAC_HARMONY", 
#     nNeighbors = 50, 
#     minDist = 0.5, 
#     metric = "cosine",force=TRUE
# )

# p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAP_ATAC_HARMONY")
# p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters_HAR", embedding = "UMAP_ATAC_HARMONY")
# p3 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "mappedCellType", embedding = "UMAP_ATAC_HARMONY")
# plotPDF(p1,p2,p3, name = "Plot-UMAP_ATAC_HARMONY-Sample-Clusters_HAR-mappedCellType-HAR.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

# # Create confusion matrices: LSI
# cM <- confusionMatrix(paste0(proj1$Clusters_LSI), paste0(proj1$Sample))
# cM <- cM / Matrix::rowSums(cM)
# p1 <- pheatmap::pheatmap(
#     mat = as.matrix(cM), 
#     color = paletteContinuous("whiteBlue"), 
#     border_color = "black"
# )

# cM <- confusionMatrix(paste0(proj1$Clusters_LSI), paste0(proj1$mappedCellType))
# cM <- cM / Matrix::rowSums(cM)
# p2 <- pheatmap::pheatmap(
#     mat = as.matrix(cM), 
#     color = paletteContinuous("whiteBlue"), 
#     border_color = "black"
# )
# plotPDF(p1, name = "confusionMap_heatmap_Clusters_LSI_Sample.pdf", ArchRProj = proj1, addDOC = FALSE)
# plotPDF(p2, name = "confusionMap_heatmap_Clusters_LSI_celltypeRNA.pdf", ArchRProj = proj1, addDOC = FALSE)


# # Create confusion matrices: HAR
# cM <- confusionMatrix(paste0(proj1$Clusters_HAR), paste0(proj1$Sample))
# cM <- cM / Matrix::rowSums(cM)
# p1 <- pheatmap::pheatmap(
#     mat = as.matrix(cM), 
#     color = paletteContinuous("whiteBlue"), 
#     border_color = "black"
# )

# cM <- confusionMatrix(paste0(proj1$Clusters_HAR), paste0(proj1$mappedCellType))
# cM <- cM / Matrix::rowSums(cM)
# p2 <- pheatmap::pheatmap(
#     mat = as.matrix(cM), 
#     color = paletteContinuous("whiteBlue"), 
#     border_color = "black"
# )
# plotPDF(p1, name = "confusionMap_heatmap_Clusters_HAR_Sample.pdf", ArchRProj = proj1, addDOC = FALSE)
# plotPDF(p2, name = "confusionMap_heatmap_Clusters_HAR_celltypeRNA.pdf", ArchRProj = proj1, addDOC = FALSE)

# saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = FALSE)

