################################################################################
# Multiome ATAC: build the ArchR project and add paired gene expression
#
# Paper : Amrute, Luo et al., Nature 635, 423-433 (2024) | doi:10.1038/s41586-024-08008-5
# Part  : Single-nucleus multiome (paired RNA + ATAC)
#
# Purpose
#   Creates the ArchR project from the per-sample ArrowFiles, keeps the nuclei
#   retained in the annotated RNA object, adds the paired gene-expression matrix
#   and transfers cell type, HF etiology and condition from the RNA analysis.
#
# Inputs
#   metadata/multiome_samples.csv
#   <arrow_dir>/<sample>/<sample>.arrow
#   <counts_dir>/<batch>/<sample>/outs/filtered_feature_bc_matrix.h5
#   Annotated multiome RNA object and its cell list (see Paths)
#
# Outputs
#   ArchR project in <proj_dir> (Save-proj1)
#
# Run order
#   Upstream  : 01_rna_merge_qc_integration.R (+ annotation)
#   Downstream: 03-08 ArchR scripts
################################################################################

library(ArchR)
addArchRGenome("hg38")
library(Seurat)
library(dplyr)

## ---- Paths (EDIT) ----
repo_dir      <- ".."                                         # EDIT
arrow_dir     <- "path/to/samplePostQC"                       # EDIT: <sample>/<sample>.arrow
counts_dir    <- "path/to/Multiome/counts"                    # EDIT: Cell Ranger ARC runs
rna_rds       <- "path/to/integrated/myocardium_multiome.rds" # EDIT: annotated multiome RNA
cell_list_rds <- "path/to/integrated/myocardium_multiome_cellList"  # EDIT
proj_dir      <- "path/to/analysis_batch_1_2_3/Save-proj1"    # EDIT: ArchR project output
source(file.path(repo_dir, "R", "utils.R"))

samples <- read_sample_sheet(file.path(repo_dir, "metadata", "multiome_samples.csv"))

## ---- ArchR project from ArrowFiles (QC already applied when creating arrows) ----
ArrowFiles <- file.path(arrow_dir, samples$sample, paste0(samples$sample, ".arrow"))
proj1 <- ArchRProject(ArrowFiles, copyArrows = TRUE)

# Keep only nuclei present in the annotated RNA object
scRNA <- readRDS(rna_rds)
col <- readRDS(cell_list_rds)
proj1 <- subsetCells(ArchRProj = proj1, cellNames = col)

## ---- Paired gene-expression matrix ----
seRNA <- import10xFeatureMatrix(input = multiome_matrix_paths(samples, counts_dir, h5 = TRUE),
                                names = samples$sample)
# One count matrix per sample, combined in sample order
seRNAcombined <- do.call(cbind, lapply(seRNA, assay))
seRNA2 <- SummarizedExperiment(assays = list(counts = seRNAcombined), rowRanges = rowRanges(seRNA[[1]]))

proj1 <- addGeneExpressionMatrix(
  input = proj1,
  seRNA = seRNA2,
  chromSizes = getChromSizes(proj1),
  excludeChr = c("chrM", "chrY"),
  scaleTo = 10000,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addGeneExpressionMatrix")
)

## ---- Transfer annotations from the RNA object ----
# NOTE: assumes nuclei are in the same order in scRNA and proj1
proj1$mappedCellType <- as.character(scRNA$predicted.celltype)
proj1$HFetiology     <- as.character(scRNA$HF.etiology)
proj1$condition      <- as.character(scRNA$condition)

proj1 <- saveArchRProject(ArchRProj = proj1, outputDirectory = proj_dir, load = TRUE)
