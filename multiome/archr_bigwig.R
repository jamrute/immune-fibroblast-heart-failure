library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/", sep=""))

proj1 <- loadArchRProject(path = "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-proj1")

getGroupBW(
  ArchRProj = proj1,
  groupBy = "mappedCellType",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 30000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)
