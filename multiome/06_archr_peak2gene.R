################################################################################
# Multiome: peak-to-gene links and co-accessibility
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
# Part  : Single-nucleus multiome (paired RNA + ATAC)
#
# Purpose
#   Links peaks to genes by correlating accessibility with paired expression
#   (HAR_Combined embedding), plots the peak-to-gene heatmap by cell type and
#   computes peak co-accessibility.
#
# Inputs  : ArchR project (Save-proj1) with the HAR_Combined embedding
# Outputs : plotPeak2GeneHeatmap.pdf; peak2gene links and co-accessibility
#           stored in the project
#
# Run order
#   Upstream  : 05_archr_multiomic_clustering.R
#   Downstream: 07_archr_peak2gene_export.R
################################################################################

library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

## ---- Paths (EDIT) ----
analysis_dir <- "path/to/analysis_batch_1_2_3"   # EDIT: folder containing Save-proj1/
setwd(analysis_dir)

proj1 <- loadArchRProject(path = "Save-proj1")

## ---- Peak-to-gene links ----
proj1 <- addPeak2GeneLinks(ArchRProj = proj1, reducedDims = "HAR_Combined", useMatrix = "GeneExpressionMatrix")

p <- plotPeak2GeneHeatmap(ArchRProj = proj1, groupBy = "mappedCellType", corCutOff = 0.3)
plotPDF(p, name = "plotPeak2GeneHeatmap.pdf", width = 4, height = 8, ArchRProj = proj1, addDOC = FALSE)

## ---- Co-accessibility ----
proj1 <- addCoAccessibility(ArchRProj = proj1, reducedDims = "HAR_Combined")

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = TRUE)
