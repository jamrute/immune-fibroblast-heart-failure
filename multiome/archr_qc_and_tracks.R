library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/", sep=""))

proj1 <- loadArchRProject(path = "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-proj1")

p1 <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

p2 <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

plotPDF(p1,p2, name = "QC-Sample-Statistics.pdf", ArchRProj = proj1, addDOC = FALSE, width = 60, height = 10)
