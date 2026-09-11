################################################################################
# Multiome ATAC: differential accessibility, donor vs HF and donor vs AMI
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
# Part  : Single-nucleus multiome (paired RNA + ATAC)
#
# Purpose
#   For each major cell type, tests differential peak accessibility between
#   donor and heart failure (condition) and between donor and acute MI
#   (HF etiology), saving MA/volcano plots and the test results.
#
# Inputs  : ArchR project (Save-proj1) with a PeakMatrix
# Outputs : Save-proj1/ArchR_output/peakSet.rds
#           Save-proj1/ArchR_output/<CellType>_Donor-vs-<CellType>_{HF,AMI}.rds
#           <...>-Markers-MA-Volcano.pdf (ArchR Plots/ folder)
#
# Run order
#   Upstream  : 03_archr_peak_calling_markers.R
#   Downstream: none
################################################################################

library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

## ---- Paths (EDIT) ----
analysis_dir <- "path/to/analysis_batch_1_2_3"   # EDIT: folder containing Save-proj1/
repo_dir     <- ".."                             # EDIT: root of this repository
setwd(analysis_dir)
source(file.path(repo_dir, "R", "utils.R"))

proj1 <- loadArchRProject(path = "Save-proj1")
saveRDS(proj1@peakSet, "./Save-proj1/ArchR_output/peakSet.rds")

cell_types <- c("Fibroblast", "Myeloid", "Endothelium", "Cardiomyocyte",
                "Pericyte", "SMC", "Endocardium", "TNKCells")

## ---- Donor vs HF (condition) ----
proj1$cellType_condition <- paste(proj1$mappedCellType, proj1$condition, sep = "_")
for (ct in cell_types) {
  archr_differential_peaks(proj1, group_by = "cellType_condition",
                           use_group = paste0(ct, "_Donor"), bgd_group = paste0(ct, "_HF"))
}

## ---- Donor vs acute MI (HF etiology) ----
proj1$cellType_HFetiology <- paste(proj1$mappedCellType, proj1$HFetiology, sep = "_")
for (ct in cell_types) {
  archr_differential_peaks(proj1, group_by = "cellType_HFetiology",
                           use_group = paste0(ct, "_Donor"), bgd_group = paste0(ct, "_AMI"))
}
