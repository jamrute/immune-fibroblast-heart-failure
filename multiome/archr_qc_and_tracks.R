################################################################################
# Multiome ATAC: per-sample QC plots and cell-type bigWig tracks
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
# Part  : Single-nucleus multiome (paired RNA + ATAC)
#
# Purpose
#   Plots TSS enrichment and fragment counts per sample, and exports
#   pseudobulk accessibility tracks (bigWig, ReadsInTSS-normalized) per cell type
#   for genome browsers.
#
# Inputs  : ArchR project (Save-proj1)
# Outputs : QC-Sample-Statistics.pdf; Save-proj1/GroupBigWigs/mappedCellType/*.bw
#
# Run order
#   Upstream  : 02_archr_project.R
#   Downstream: none
################################################################################

library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

## ---- Paths (EDIT) ----
analysis_dir <- "path/to/analysis_batch_1_2_3"   # EDIT: folder containing Save-proj1/
setwd(analysis_dir)

proj1 <- loadArchRProject(path = "Save-proj1")

## ---- QC per sample ----
qc_plots <- lapply(c("TSSEnrichment", "log10(nFrags)"), function(metric) {
  plotGroups(ArchRProj = proj1, groupBy = "Sample", colorBy = "cellColData", name = metric,
             plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
})
plotPDF(qc_plots[[1]], qc_plots[[2]], name = "QC-Sample-Statistics.pdf", ArchRProj = proj1,
        addDOC = FALSE, width = 60, height = 10)

## ---- Cell-type bigWig tracks ----
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
