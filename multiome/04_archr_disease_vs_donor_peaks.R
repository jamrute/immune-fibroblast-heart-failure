library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

setwd(paste("/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/", sep=""))

proj1 <- loadArchRProject(path = "/data/Junedh/Amgen_Epigenomics/analysis_batch_1_2_3/Save-proj1")

peakSet <- proj1@peakSet
saveRDS(peakSet, "./Save-proj1/ArchR_output/peakSet.rds")

# Create a meta-column for ancestry + cell type
proj1$cellType_condition <- paste(proj1$mappedCellType, proj1$condition,sep="_")

# Fibroblast
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_condition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Fibroblast_Donor",
  bgdGroups = "Fibroblast_HF",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Fibroblast_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Fibroblast_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Fibroblast_Donor-vs-Fibroblast_HF-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Fibroblast_Donor-vs-Fibroblast_HF.rds")

# Myeloid
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_condition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Myeloid_Donor",
  bgdGroups = "Myeloid_HF",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Myeloid_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Myeloid_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Myeloid_Donor-vs-Myeloid_HF-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Myeloid_Donor-vs-Myeloid_HF.rds")



# Endothelium
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_condition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Endothelium_Donor",
  bgdGroups = "Endothelium_HF",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Endothelium_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Endothelium_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Endothelium_Donor-vs-Endothelium_HF-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Endothelium_Donor-vs-Endothelium_HF.rds")


# Cardiomyocyte
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_condition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiomyocyte_Donor",
  bgdGroups = "Cardiomyocyte_HF",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Cardiomyocyte_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Cardiomyocyte_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Cardiomyocyte_Donor-vs-Cardiomyocyte_HF-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Cardiomyocyte_Donor-vs-Cardiomyocyte_HF.rds")

# Pericyte
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_condition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Pericyte_Donor",
  bgdGroups = "Pericyte_HF",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Pericyte_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Pericyte_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Pericyte_Donor-vs-Pericyte_HF-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Pericyte_Donor-vs-Pericyte_HF.rds")


# SMC
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_condition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "SMC_Donor",
  bgdGroups = "SMC_HF",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "SMC_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "SMC_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "SMC_Donor-vs-SMC_HF-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/SMC_Donor-vs-SMC_HF.rds")



# Endocardium
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_condition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Endocardium_Donor",
  bgdGroups = "Endocardium_HF",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Endocardium_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Endocardium_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Endocardium_Donor-vs-Endocardium_HF-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Endocardium_Donor-vs-Endocardium_HF.rds")



# TNKCells
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_condition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "TNKCells_Donor",
  bgdGroups = "TNKCells_HF",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "TNKCells_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "TNKCells_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "TNKCells_Donor-vs-TNKCells_HF-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/TNKCells_Donor-vs-TNKCells_HF.rds")

################### AMI
# Create a meta-column for ancestry + cell type
proj1$cellType_HFetiology <- paste(proj1$mappedCellType, proj1$HFetiology,sep="_")

# Fibroblast
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_HFetiology",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Fibroblast_Donor",
  bgdGroups = "Fibroblast_AMI",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Fibroblast_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Fibroblast_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Fibroblast_Donor-vs-Fibroblast_AMI-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Fibroblast_Donor-vs-Fibroblast_AMI.rds")

# Myeloid
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_HFetiology",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Myeloid_Donor",
  bgdGroups = "Myeloid_AMI",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Myeloid_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Myeloid_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Myeloid_Donor-vs-Myeloid_AMI-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Myeloid_Donor-vs-Myeloid_AMI.rds")





# Endothelium
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_HFetiology",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Endothelium_Donor",
  bgdGroups = "Endothelium_AMI",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Endothelium_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Endothelium_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Endothelium_Donor-vs-Endothelium_AMI-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Endothelium_Donor-vs-Endothelium_AMI.rds")



# Cardiomyocyte
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_HFetiology",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiomyocyte_Donor",
  bgdGroups = "Cardiomyocyte_AMI",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Cardiomyocyte_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Cardiomyocyte_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Cardiomyocyte_Donor-vs-Cardiomyocyte_AMI-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Cardiomyocyte_Donor-vs-Cardiomyocyte_AMI.rds")



# Pericyte
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_HFetiology",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Pericyte_Donor",
  bgdGroups = "Pericyte_AMI",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Pericyte_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Pericyte_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Pericyte_Donor-vs-Pericyte_AMI-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Pericyte_Donor-vs-Pericyte_AMI.rds")


# SMC
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_HFetiology",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "SMC_Donor",
  bgdGroups = "SMC_AMI",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "SMC_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "SMC_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "SMC_Donor-vs-SMC_AMI-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/SMC_Donor-vs-SMC_AMI.rds")


# Endocardium
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_HFetiology",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Endocardium_Donor",
  bgdGroups = "Endocardium_AMI",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "Endocardium_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "Endocardium_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "Endocardium_Donor-vs-Endocardium_AMI-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/Endocardium_Donor-vs-Endocardium_AMI.rds")




# TNKCells
markerTest <- getMarkerFeatures(
  ArchRProj = proj1, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_HFetiology",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "TNKCells_Donor",
  bgdGroups = "TNKCells_AMI",
  maxCells = 5000
)

pma <- markerPlot(seMarker = markerTest, name = "TNKCells_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "TNKCells_Donor", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")

plotPDF(pma, pv, name = "TNKCells_Donor-vs-TNKCells_AMI-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj1, addDOC = FALSE)
saveRDS(markerTest, "./Save-proj1/ArchR_output/TNKCells_Donor-vs-TNKCells_AMI.rds")





