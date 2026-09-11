library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/", sep=""))

scRNA <- readRDS("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/RNA/globalObjectConstruction_batch1_2_3/integrated/myocardium_multiome.rds")
proj1 <- loadArchRProject(path = "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-proj1")

############################################### Peak 2 Gene ###########################################################

proj1 <- addPeak2GeneLinks(ArchRProj = proj1, reducedDims = "HAR_Combined", useMatrix = "GeneExpressionMatrix")

p <- plotPeak2GeneHeatmap(ArchRProj = proj1, groupBy = "mappedCellType", corCutOff = 0.3)
plotPDF(p, name = "plotPeak2GeneHeatmap.pdf", width = 4, height = 8, ArchRProj = proj1, addDOC = FALSE)

#Co-accessibility analysis
test <- addCoAccessibility(
    ArchRProj = proj1,
    reducedDims = "HAR_Combined"
)
saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)