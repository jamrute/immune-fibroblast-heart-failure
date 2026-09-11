library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/", sep=""))

scRNA <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/integrated/myocardium_multiome.rds")
proj1 <- loadArchRProject(path = "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-proj1")

################### Peak Calling with MACS2

###### Multiomic Clustering

# # RNA
# proj1 <- addIterativeLSI(
#   ArchRProj = proj1, 
#   clusterParams = list(
#     resolution = 0.2, 
#     sampleCells = 10000,
#     n.start = 10
#   ),
#   saveIterations = FALSE,
#   useMatrix = "GeneExpressionMatrix", 
#   depthCol = "Gex_nUMI",
#   varFeatures = 2500,
#   firstSelection = "variable",
#   binarize = FALSE,
#   name = "LSI_RNA",
#   force = TRUE
# )

# # ATAC
# proj1 <- addIterativeLSI(
#   ArchRProj = proj1, 
#   clusterParams = list(
#     resolution = 0.2, 
#     sampleCells = 10000,
#     n.start = 10
#   ),
#   saveIterations = FALSE,
#   useMatrix = "TileMatrix", 
#   depthCol = "nFrags",
#   name = "LSI_ATAC",
#   force = TRUE
# )

# # Combined embedding
# proj1 <- addCombinedDims(proj1, reducedDims = c("LSI_RNA", "LSI_ATAC"), name =  "LSI_Combined")

# # Harmony
# proj1 <- addHarmony(
#     ArchRProj = proj1,
#     reducedDims = "LSI_Combined",
#     name = "HAR_Combined",
#     groupBy = "Sample",force=TRUE
# )

# proj1 <- addUMAP(proj1, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
# proj1 <- addUMAP(proj1, reducedDims = "LSI_RNA", name = "UMAP_LSI_RNA", minDist = 0.8, force = TRUE)
# proj1 <- addUMAP(proj1, reducedDims = "HAR_Combined", name = "UMAP_HAR_Combined", minDist = 0.8, force = TRUE)

# proj1 <- addClusters(proj1, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.1, force = TRUE)
# proj1 <- addClusters(proj1, reducedDims = "LSI_RNA", name = "Clusters_LSI_RNA", resolution = 0.1, force = TRUE)
# proj1 <- addClusters(proj1, reducedDims = "HAR_Combined", name = "Clusters_HAR_Combined", resolution = 0.1, force = TRUE)

# saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)

# We can plot how each of these dimensionality reductions look with respect to the clusters called in "LSI_Combined".
p1 <- plotEmbedding(proj1, name = "Clusters_ATAC", embedding = "UMAP_ATAC", size = 1, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj1, name = "Clusters_LSI_RNA", embedding = "UMAP_LSI_RNA", size = 1, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj1, name = "Clusters_HAR_Combined", embedding = "UMAP_HAR_Combined", size = 1, labelAsFactors=F, labelMeans=F)
p4 <- plotEmbedding(proj1, name = "mappedCellType", embedding = "UMAP_HAR_Combined", size = 1, labelAsFactors=F, labelMeans=F)

p <- lapply(list(p1,p2,p3), function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p))

plotPDF(p1, p2, p3, p4, name = "UMAP-scATAC-scRNA-Combined.pdf", addDOC = FALSE, ArchRProj = proj1)

# Clusters_ATAC vs mappedCellType
cM_atac_rna <- confusionMatrix(paste0(proj1$Clusters_ATAC), paste0(proj1$mappedCellType))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)
library(pheatmap)
p_atac_rna <- pheatmap::pheatmap(
  mat = as.matrix(cM_atac_rna), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

plotPDF(p_atac_rna, name = "Clusters_ATAC_Clusters_mappedCellType_confusionMatrix.pdf", addDOC = FALSE, ArchRProj = proj1)


# CellType_Seurat_RNA vs Clusters_HAR_Combined
cM_atac_rna <- confusionMatrix(paste0(proj1$mappedCellType), paste0(proj1$Clusters_HAR_Combined))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)
library(pheatmap)
p_atac_rna <- pheatmap::pheatmap(
  mat = as.matrix(cM_atac_rna), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

plotPDF(p_atac_rna, name = "mappedCellType_HAR_Combined_confusionMatrix.pdf", addDOC = FALSE, ArchRProj = proj1)

#saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)



