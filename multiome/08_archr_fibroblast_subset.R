library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/", sep=""))

# proj1 <- loadArchRProject(path = "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-proj1")
# projFibro <- saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-projFibro", load = TRUE, dropCells = TRUE)

# idxSample <- BiocGenerics::which(projFibro$mappedCellType %in% "Fibroblast")
# cellsSample <- projFibro$cellNames[idxSample]
# projFibro <- projFibro[cellsSample, ]

# scRNA <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/fibroblast/multiome_fibroblast_mappedToCITEseq.rds")

# # Transfer cluster and identities from scRNA
# predictedcelltype <- as.character(scRNA$predicted.celltype)
# projFibro$predictedcelltype <- predictedcelltype

# ###### Multiomic Clustering

# # RNA
# projFibro <- addIterativeLSI(
#   ArchRProj = projFibro, 
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
# projFibro <- addIterativeLSI(
#   ArchRProj = projFibro, 
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
# projFibro <- addCombinedDims(projFibro, reducedDims = c("LSI_RNA", "LSI_ATAC"), name =  "LSI_Combined")

# # Harmony
# projFibro <- addHarmony(
#     ArchRProj = projFibro,
#     reducedDims = "LSI_Combined",
#     name = "HAR_Combined",
#     groupBy = "Sample",force=TRUE
# )

# projFibro <- addUMAP(projFibro, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
# projFibro <- addUMAP(projFibro, reducedDims = "LSI_RNA", name = "UMAP_LSI_RNA", minDist = 0.8, force = TRUE)
# projFibro <- addUMAP(projFibro, reducedDims = "HAR_Combined", name = "UMAP_HAR_Combined", minDist = 0.8, force = TRUE)

# projFibro <- addClusters(projFibro, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.1, force = TRUE)
# projFibro <- addClusters(projFibro, reducedDims = "LSI_RNA", name = "Clusters_LSI_RNA", resolution = 0.1, force = TRUE)
# projFibro <- addClusters(projFibro, reducedDims = "HAR_Combined", name = "Clusters_HAR_Combined", resolution = 0.1, force = TRUE)

# saveArchRProject(ArchRProj = projFibro, outputDirectory = "Save-projFibro", load = TRUE)

projFibro <- loadArchRProject(path = "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-projFibro")

# We can plot how each of these dimensionality reductions look with respect to the clusters called in "LSI_Combined".
p1 <- plotEmbedding(projFibro, name = "Clusters_ATAC", embedding = "UMAP_ATAC", size = 1, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(projFibro, name = "Clusters_LSI_RNA", embedding = "UMAP_LSI_RNA", size = 1, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(projFibro, name = "Clusters_HAR_Combined", embedding = "UMAP_HAR_Combined", size = 1, labelAsFactors=F, labelMeans=F)
p4 <- plotEmbedding(projFibro, name = "predictedcelltype", embedding = "UMAP_HAR_Combined", size = 1, labelAsFactors=F, labelMeans=F)

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

plotPDF(p1, p2, p3, p4, name = "UMAP-scATAC-scRNA-Combined.pdf", addDOC = FALSE, ArchRProj = projFibro)

# Clusters_ATAC vs Clusters_LSI_RNA
cM_atac_rna <- confusionMatrix(paste0(projFibro$Clusters_ATAC), paste0(projFibro$Clusters_LSI_RNA))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)
library(pheatmap)
p_atac_rna <- pheatmap::pheatmap(
  mat = as.matrix(cM_atac_rna), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

plotPDF(p_atac_rna, name = "Clusters_ATAC_Clusters_LSI_RNA_confusionMatrix.pdf", addDOC = FALSE, ArchRProj = projFibro)


# CellType_Seurat_RNA vs Clusters_HAR_Combined
cM_atac_rna <- confusionMatrix(paste0(projFibro$predictedcelltype), paste0(projFibro$Clusters_HAR_Combined))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)
library(pheatmap)
p_atac_rna <- pheatmap::pheatmap(
  mat = as.matrix(cM_atac_rna), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

plotPDF(p_atac_rna, name = "predictedcelltype_HAR_Combined_confusionMatrix.pdf", addDOC = FALSE, ArchRProj = projFibro)

saveArchRProject(ArchRProj = projFibro, outputDirectory = "Save-projFibro", load = TRUE)




